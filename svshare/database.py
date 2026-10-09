# Write analyze results to a SQLite database.
#
# A call is one caller's report of an SV in one sample (one record in the VCF an
# SV caller wrote). An SV is the variant those calls describe: Jasmine decides
# which calls are the same SV, and keeps only one of them in its merged record.
# The database keeps both, so nothing is lost:
#
#   svs              one row per SV Jasmine identified
#   calls            one row per original call, linked to its SV by sv_id
#   commands         every third-party command run, in order, with the tool
#                    and the sample it was run on, if applicable
#
# How closely an SV's calls agree is kept as each call's own position and length,
# plus the SV's start_spread and length_spread (max - min across its calls, in bp),
# so the matching can later be made stricter without rerunning Jasmine.
#
# An SV's pos, end_pos, and svlen are those of the call Jasmine kept, and
# gnomad_id and gnomad_af are those of the matching gnomAD record, if any.
import shlex
import sqlite3
from pathlib import Path

SCHEMA = """
CREATE TABLE svs (
    sv_id INTEGER PRIMARY KEY,
    chrom TEXT,
    pos INTEGER,
    end_pos INTEGER,
    svtype TEXT,
    svlen INTEGER,
    start_spread INTEGER,
    length_spread INTEGER,
    gnomad_id TEXT,
    gnomad_af REAL
);
CREATE TABLE calls (
    call_id INTEGER PRIMARY KEY,
    sv_id INTEGER NOT NULL REFERENCES svs(sv_id),
    sample TEXT NOT NULL,
    caller TEXT NOT NULL,
    vcf_id TEXT NOT NULL,
    chrom TEXT,
    pos INTEGER,
    end_pos INTEGER,
    svtype TEXT,
    svlen INTEGER,
    genotype TEXT,
    read_support INTEGER
);
CREATE INDEX calls_sv_id ON calls(sv_id);
CREATE INDEX calls_sample ON calls(sample);
CREATE INDEX calls_read_support ON calls(read_support);
CREATE TABLE commands (
    command_id INTEGER PRIMARY KEY,
    sample TEXT,
    tool TEXT NOT NULL,
    command TEXT NOT NULL
);
"""


# Yield (id, chrom, pos, info, format_values) for each record in a VCF, where
# info and format_values are dicts of the INFO and first sample's FORMAT fields.
def read_vcf_records(vcf):
    with open(vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info = dict(
                item.split("=", 1) if "=" in item else (item, True)
                for item in fields[7].split(";")
            )
            fmt = {}
            if len(fields) > 9:
                fmt = dict(zip(fields[8].split(":"), fields[9].split(":")))
            yield fields[2], fields[0], int(fields[1]), info, fmt


def _int_or_none(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


# For each record in Jasmine's merged VCF, yield (chrom, pos, end_pos, svtype,
# svlen, members), where members lists the (input index, call ID) pairs merged
# into it. SUPP_VEC has one digit per input VCF, in file_list order, and IDLIST
# names the merged calls in that same order.
def read_jasmine_groups(merged_vcf):
    for _id, chrom, pos, info, _fmt in read_vcf_records(merged_vcf):
        inputs = [i for i, bit in enumerate(info["SUPP_VEC"]) if bit == "1"]
        members = list(zip(inputs, info["IDLIST"].split(",")))
        yield (chrom, pos, _int_or_none(info.get("END")), info.get("SVTYPE"),
               _int_or_none(info.get("SVLEN")), members)


# Write the database to db_path, replacing any earlier one. inputs lists
# (sample, caller, vcf) in the order the VCFs were given to Jasmine; commands
# lists (sample, tool, command) for each command run, in the order they ran.
def write_database(merged_vcf, inputs, db_path, commands):
    calls_by_input = [
        {call_id: (chrom, pos, info, fmt)
         for call_id, chrom, pos, info, fmt in read_vcf_records(vcf)}
        for _sample, _caller, vcf in inputs
    ]

    db_path = Path(db_path)
    db_path.unlink(missing_ok=True)
    with sqlite3.connect(db_path) as db:
        db.executescript(SCHEMA)
        db.executemany(
            "INSERT INTO commands (sample, tool, command) VALUES (?, ?, ?)",
            [(sample, tool, shlex.join(cmd)) for sample, tool, cmd in commands],
        )

        for sv_id, (chrom, pos, end_pos, svtype, svlen, members) in enumerate(
            read_jasmine_groups(merged_vcf), start=1
        ):
            calls = []
            for index, vcf_id in members:
                sample, caller, _vcf = inputs[index]
                call_chrom, call_pos, info, fmt = calls_by_input[index][vcf_id]
                calls.append((
                    sv_id, sample, caller, vcf_id, call_chrom, call_pos,
                    _int_or_none(info.get("END")), info.get("SVTYPE"),
                    _int_or_none(info.get("SVLEN")), fmt.get("GT"),
                    # Both sniffles2 and cuteSV report variant-supporting reads as DV.
                    _int_or_none(fmt.get("DV")),
                ))

            starts = [call[5] for call in calls]
            lengths = [abs(call[8]) for call in calls if call[8] is not None]
            db.execute(
                "INSERT INTO svs (sv_id, chrom, pos, end_pos, svtype, svlen, start_spread,"
                " length_spread) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                (sv_id, chrom, pos, end_pos, svtype, svlen, max(starts) - min(starts),
                 max(lengths) - min(lengths) if lengths else None),
            )
            db.executemany(
                "INSERT INTO calls (sv_id, sample, caller, vcf_id, chrom, pos, end_pos,"
                " svtype, svlen, genotype, read_support)"
                " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                calls,
            )
    db.close()
    return db_path
