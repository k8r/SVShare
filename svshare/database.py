# Read and write analyze results in a SQLite database.
#
# A call is one caller's report of an SV in one sample (one record in the VCF an
# SV caller wrote). An SV is the variant those calls describe: Jasmine decides
# which calls are the same SV, and keeps only one of them in its merged record.
# The database keeps both, so nothing is lost:
#
#   svs              one row per SV Jasmine identified
#   calls            one row per original call, linked to its SV by sv_id
#   transcripts      one row per transcript near any SV, with its gene, type,
#                    and whether it's its gene's MANE Select transcript
#   sv_transcript_parts
#                    one row per part of a transcript (CDS, start or stop codon,
#                    UTR, exon, intron, or promoter) near an SV, with its exon
#                    or intron number, its position, and whether the SV overlaps it
#   genes            one row per gene with a transcript near any SV, with its
#                    readable name, a summary of what it does, and its ClinGen
#                    dosage sensitivity scores, where known
#   commands         every third-party command run, in order, with the tool
#                    and the sample it was run on, if applicable
#
# How closely an SV's calls agree is kept as each call's own position and length,
# plus the SV's start_spread and length_spread (max - min across its calls, in bp),
# so the matching can later be made stricter without rerunning Jasmine.
#
# An SV's pos, vcf_end, and svlen are those of the call Jasmine kept. The
# gnomad_ and colorsdb_ columns hold the ID and allele frequency of the matching
# gnomAD and CoLoRSdb record, if any. near_reference_gap describes how close the
# SV is to a gap in the reference, if it's at or near one.
import sqlite3
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

SCHEMA = """
CREATE TABLE svs (
    sv_id INTEGER PRIMARY KEY,
    chrom TEXT,
    pos INTEGER,
    vcf_end INTEGER,
    svtype TEXT,
    svlen INTEGER,
    start_spread INTEGER,
    length_spread INTEGER,
    gnomad_id TEXT,
    gnomad_af REAL,
    colorsdb_id TEXT,
    colorsdb_af REAL,
    near_reference_gap TEXT
);
CREATE TABLE calls (
    call_id INTEGER PRIMARY KEY,
    sv_id INTEGER NOT NULL REFERENCES svs(sv_id),
    sample TEXT NOT NULL,
    caller TEXT NOT NULL,
    vcf_id TEXT NOT NULL,
    chrom TEXT,
    pos INTEGER,
    vcf_end INTEGER,
    svtype TEXT,
    svlen INTEGER,
    genotype TEXT,
    read_support INTEGER
);
CREATE INDEX calls_sv_id ON calls(sv_id);
CREATE INDEX calls_sample ON calls(sample);
CREATE INDEX calls_read_support ON calls(read_support);
CREATE TABLE transcripts (
    transcript_id TEXT PRIMARY KEY,
    gene_name TEXT NOT NULL,
    gene_id TEXT NOT NULL,
    transcript_type TEXT NOT NULL,
    mane_select INTEGER NOT NULL,
    chrom TEXT NOT NULL,
    strand TEXT NOT NULL
);
CREATE INDEX transcripts_gene_name ON transcripts(gene_name);
CREATE TABLE sv_transcript_parts (
    sv_id INTEGER NOT NULL REFERENCES svs(sv_id),
    transcript_id TEXT NOT NULL REFERENCES transcripts(transcript_id),
    region TEXT NOT NULL,
    number INTEGER,
    pos INTEGER NOT NULL,
    end_pos INTEGER NOT NULL,
    overlaps_sv INTEGER NOT NULL
);
CREATE INDEX sv_transcript_parts_sv_id ON sv_transcript_parts(sv_id);
CREATE TABLE genes (
    gene_id TEXT PRIMARY KEY,
    hgnc_id TEXT,
    full_name TEXT,
    summary TEXT,
    haploinsufficiency_score INTEGER,
    haploinsufficiency TEXT,
    triplosensitivity_score INTEGER,
    triplosensitivity TEXT
);
CREATE TABLE commands (
    command_id INTEGER PRIMARY KEY,
    sample TEXT,
    tool TEXT NOT NULL,
    command TEXT NOT NULL
);
"""


# Where an SV is, and its type and length.
@dataclass
class SVPosition:
    sv_id: int
    chrom: str
    pos: int
    vcf_end: int | None
    svtype: str
    svlen: int | None


# Which SVs a report includes. Each field left at its default keeps every SV.
@dataclass
class SVFilter:
    # A sample counts toward min_samples if every caller in callers called the
    # SV in it, or any caller if callers is empty.
    min_samples: int = 1
    callers: tuple[str, ...] = ()
    svtype: str | None = None
    # SVs not found in a population dataset still pass.
    population_datasets: tuple[str, ...] = ()
    max_allele_frequency: float | None = None
    genes: tuple[str, ...] = ()
    region: tuple[str, int, int] | None = None
    # True keeps only SVs at or near a reference gap, False only the rest.
    near_reference_gap: bool | None = None


# Open the database, saving any changes and closing it when done.
@contextmanager
def connect(db_path):
    db = sqlite3.connect(db_path)
    try:
        with db:
            yield db
    finally:
        db.close()


# Return every SV's position, type, and length.
def sv_positions(db_path):
    with connect(db_path) as db:
        return [SVPosition(*row) for row in db.execute(
            "SELECT sv_id, chrom, pos, vcf_end, svtype, svlen FROM svs"
        )]


# Return the ID of every gene with a saved transcript, sorted.
def gene_ids(db_path):
    with connect(db_path) as db:
        return [gene_id for (gene_id,) in db.execute(
            "SELECT DISTINCT gene_id FROM transcripts ORDER BY gene_id"
        )]


# Return the name of every gene with a saved transcript.
def gene_names(db_path):
    with connect(db_path) as db:
        return {gene_name for (gene_name,) in db.execute(
            "SELECT DISTINCT gene_name FROM transcripts"
        )}


def _sample_support_condition(callers, min_samples):
    samples = "SELECT sv_id, sample FROM calls"
    params = []
    if callers:
        samples += (f" WHERE caller IN ({', '.join('?' * len(callers))})"
                    " GROUP BY sv_id, sample HAVING COUNT(DISTINCT caller) = ?")
        params = list(callers) + [len(callers)]
    else:
        samples += " GROUP BY sv_id, sample"
    return (
        f"sv_id IN (SELECT sv_id FROM ({samples}) GROUP BY sv_id HAVING COUNT(*) >= ?)",
        params + [min_samples],
    )


def _svtype_condition(svtype):
    if svtype is None:
        return None
    return "svtype = ?", [svtype]


def _allele_frequency_condition(population_datasets, max_allele_frequency):
    if max_allele_frequency is None:
        return None
    return (
        " AND ".join(f"({dataset}_af IS NULL OR {dataset}_af <= ?)"
                     for dataset in population_datasets),
        [max_allele_frequency] * len(population_datasets),
    )


def _genes_condition(genes):
    if not genes:
        return None
    return (
        "sv_id IN (SELECT sv_id FROM sv_transcript_parts JOIN transcripts USING (transcript_id)"
        f" WHERE overlaps_sv = 1 AND gene_name COLLATE NOCASE IN ({', '.join('?' * len(genes))}))",
        list(genes),
    )


def _region_condition(region):
    if region is None:
        return None
    chrom, start, end = region
    # An insertion's or breakend's end position isn't where it ends, so only its
    # start position is checked.
    return (
        "chrom = ? AND pos <= ?"
        " AND (CASE WHEN svtype IN ('INS', 'BND') THEN pos ELSE MAX(vcf_end, pos) END) >= ?",
        [chrom, end, start],
    )


def _reference_gap_condition(near_reference_gap):
    if near_reference_gap is None:
        return None
    return f"near_reference_gap IS {'NOT ' if near_reference_gap else ''}NULL", []


# Return the IDs, in order, of the SVs in the database that pass sv_filter.
def select_sv_ids(db_path, sv_filter):
    conditions, params = [], []
    for condition in (
        _sample_support_condition(sv_filter.callers, sv_filter.min_samples),
        _svtype_condition(sv_filter.svtype),
        _allele_frequency_condition(sv_filter.population_datasets,
                                    sv_filter.max_allele_frequency),
        _genes_condition(sv_filter.genes),
        _region_condition(sv_filter.region),
        _reference_gap_condition(sv_filter.near_reference_gap),
    ):
        if condition:
            conditions.append(condition[0])
            params += condition[1]
    where = f" WHERE {' AND '.join(conditions)}" if conditions else ""
    with connect(db_path) as db:
        return [sv_id for (sv_id,) in db.execute(f"SELECT sv_id FROM svs{where} ORDER BY sv_id", params)]


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


# For each record in Jasmine's merged VCF, yield (chrom, pos, vcf_end, svtype,
# svlen, members), where members lists the (input index, call ID) pairs merged
# into it. SUPP_VEC has one digit per input VCF, in file_list order, and IDLIST
# names the merged calls in that same order.
def read_jasmine_groups(merged_vcf):
    for _id, chrom, pos, info, _fmt in read_vcf_records(merged_vcf):
        inputs = [i for i, bit in enumerate(info["SUPP_VEC"]) if bit == "1"]
        members = list(zip(inputs, info["IDLIST"].split(",")))
        yield (chrom, pos, _int_or_none(info.get("END")), info.get("SVTYPE"),
               _int_or_none(info.get("SVLEN")), members)


# Write the database to db_path, replacing any earlier one. prepped_vcfs lists
# the VCFs prepped for merging, in the order Jasmine was given them.
def write_database(merged_vcf, prepped_vcfs, db_path):
    calls_by_input = [
        {call_id: (chrom, pos, info, fmt)
         for call_id, chrom, pos, info, fmt in read_vcf_records(prepped_vcf.vcf)}
        for prepped_vcf in prepped_vcfs
    ]

    db_path = Path(db_path)
    db_path.unlink(missing_ok=True)
    with connect(db_path) as db:
        db.executescript(SCHEMA)

        for sv_id, (chrom, pos, vcf_end, svtype, svlen, members) in enumerate(
            read_jasmine_groups(merged_vcf), start=1
        ):
            calls = []
            for index, vcf_id in members:
                prepped_vcf = prepped_vcfs[index]
                call_chrom, call_pos, info, fmt = calls_by_input[index][vcf_id]
                calls.append((
                    sv_id, prepped_vcf.sample, prepped_vcf.caller, vcf_id, call_chrom, call_pos,
                    _int_or_none(info.get("END")), info.get("SVTYPE"),
                    _int_or_none(info.get("SVLEN")), fmt.get("GT"),
                    # Both sniffles2 and cuteSV report variant-supporting reads as DV.
                    _int_or_none(fmt.get("DV")),
                ))

            starts = [call[5] for call in calls]
            lengths = [abs(call[8]) for call in calls if call[8] is not None]
            db.execute(
                "INSERT INTO svs (sv_id, chrom, pos, vcf_end, svtype, svlen, start_spread,"
                " length_spread) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                (sv_id, chrom, pos, vcf_end, svtype, svlen, max(starts) - min(starts),
                 max(lengths) - min(lengths) if lengths else None),
            )
            db.executemany(
                "INSERT INTO calls (sv_id, sample, caller, vcf_id, chrom, pos, vcf_end,"
                " svtype, svlen, genotype, read_support)"
                " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                calls,
            )
    return db_path


# Save each gene's names, summary, and dosage sensitivity scores.
def save_genes(db_path, rows):
    with connect(db_path) as db:
        db.executemany(
            "INSERT INTO genes (gene_id, hgnc_id, full_name, summary,"
            " haploinsufficiency_score, haploinsufficiency, triplosensitivity_score,"
            " triplosensitivity) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            rows,
        )


# Save the transcripts near SVs, by ID, and the parts of them near each SV.
def save_transcripts(db_path, transcripts, parts):
    with connect(db_path) as db:
        db.executemany(
            "INSERT INTO transcripts (transcript_id, gene_name, gene_id, transcript_type,"
            " mane_select, chrom, strand) VALUES (?, ?, ?, ?, ?, ?, ?)",
            [(transcript_id, *values) for transcript_id, values in sorted(transcripts.items())],
        )
        db.executemany(
            "INSERT INTO sv_transcript_parts (sv_id, transcript_id, region, number, pos,"
            " end_pos, overlaps_sv) VALUES (?, ?, ?, ?, ?, ?, ?)",
            sorted(parts, key=lambda part: (part[0], part[1], part[4])),
        )


# Save the ID and allele frequency of the population dataset record each SV in
# matches, by ID, matched. dataset is the prefix of its columns in svs.
def save_population_matches(db_path, dataset, matches):
    with connect(db_path) as db:
        db.executemany(
            f"UPDATE svs SET {dataset}_id = ?, {dataset}_af = ? WHERE sv_id = ?",
            [(record_id, af, sv_id) for sv_id, (record_id, af) in matches.items()],
        )


# Save how close each SV in descriptions, by ID, is to a gap in the reference.
def save_reference_gaps(db_path, descriptions):
    with connect(db_path) as db:
        db.executemany(
            "UPDATE svs SET near_reference_gap = ? WHERE sv_id = ?",
            [(description, sv_id) for sv_id, description in descriptions.items()],
        )


# Add each outside tool run to the database, in the order they ran.
def save_tool_runs(db_path, tool_runs):
    with connect(db_path) as db:
        db.executemany(
            "INSERT INTO commands (sample, tool, command) VALUES (?, ?, ?)",
            [(tool_run.sample, tool_run.tool, str(tool_run)) for tool_run in tool_runs],
        )
