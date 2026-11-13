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
# An SV's pos, vcf_end, and svlen are those of the call Jasmine kept, and
# reference_end is where it ends on the reference, which isn't always its VCF
# END. The gnomad_ and colorsdb_ columns hold the ID and allele frequency of the matching
# gnomAD and CoLoRSdb record, if any. near_reference_gap describes how close the
# SV is to a gap in the reference, if it's at or near one.
import json
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
    reference_end INTEGER,
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


# Open the database, saving any changes and closing it when done. Rows can be
# read by column name.
@contextmanager
def connect(db_path):
    db = sqlite3.connect(db_path)
    db.row_factory = sqlite3.Row
    try:
        with db:
            yield db
    finally:
        db.close()


# Return the name of every gene with a saved transcript.
def gene_names(db_path):
    with connect(db_path) as db:
        return {row["gene_name"] for row in db.execute(
            "SELECT DISTINCT gene_name FROM transcripts"
        )}


# Matches rows whose sv_id is in a JSON list of IDs, given as one parameter, so
# there's no limit on how many IDs can be passed.
SV_ID_IN = "sv_id IN (SELECT value FROM json_each(?))"


# Return each SV in sv_ids's position, type, length, allele frequencies, and
# how close it is to a gap in the reference.
def sv_rows(db_path, sv_ids):
    with connect(db_path) as db:
        return db.execute(
            "SELECT sv_id, chrom, pos, vcf_end, reference_end, svtype, svlen, gnomad_af,"
            f" colorsdb_af, near_reference_gap FROM svs WHERE {SV_ID_IN}",
            [json.dumps(list(sv_ids))],
        ).fetchall()


# Return the sample, caller, genotype, and read support of each call of the SVs
# in sv_ids, by sample and then caller.
def sv_calls(db_path, sv_ids):
    with connect(db_path) as db:
        return db.execute(
            "SELECT sv_id, sample, caller, genotype, read_support FROM calls"
            f" WHERE {SV_ID_IN} ORDER BY sample, caller",
            [json.dumps(list(sv_ids))],
        ).fetchall()


# Return each transcript part near the SVs in sv_ids, with its transcript and its
# gene's details, by SV and then transcript and position.
def sv_transcript_parts(db_path, sv_ids):
    with connect(db_path) as db:
        return db.execute(
            "SELECT sv_id, region, number, pos, end_pos, overlaps_sv, transcript_id,"
            " gene_name, gene_id, transcript_type, mane_select, chrom, strand, hgnc_id,"
            " full_name, summary, haploinsufficiency_score, haploinsufficiency,"
            " triplosensitivity_score, triplosensitivity"
            " FROM sv_transcript_parts JOIN transcripts USING (transcript_id)"
            f" LEFT JOIN genes USING (gene_id) WHERE {SV_ID_IN}"
            " ORDER BY sv_id, transcript_id, pos",
            [json.dumps(list(sv_ids))],
        ).fetchall()


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
    return "chrom = ? AND pos <= ? AND reference_end >= ?", [chrom, end, start]


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
        return [row["sv_id"] for row in db.execute(
            f"SELECT sv_id FROM svs{where} ORDER BY sv_id", params
        )]


# Write the SVs, everything found about them, and the outside tool runs to a
# new database at db_path, replacing any earlier one.
def write_database(db_path, svs, tool_runs):
    db_path = Path(db_path)
    db_path.unlink(missing_ok=True)
    with connect(db_path) as db:
        db.executescript(SCHEMA)
        _insert_svs(db, svs)
        _insert_transcripts(db, svs)
        _insert_genes(db, svs)
        _insert_tool_runs(db, tool_runs)
    return db_path


# Insert each SV, with how far apart its calls' starts and lengths are, and its
# calls.
def _insert_svs(db, svs):
    for sv in svs:
        position = sv.position
        starts = [call.pos for call in sv.calls]
        lengths = [abs(call.svlen) for call in sv.calls if call.svlen is not None]
        db.execute(
            "INSERT INTO svs (sv_id, chrom, pos, vcf_end, reference_end, svtype, svlen,"
            " start_spread, length_spread, gnomad_id, gnomad_af, colorsdb_id, colorsdb_af,"
            " near_reference_gap) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (sv.sv_id, position.chrom, position.pos, position.vcf_end,
             position.reference_end, position.svtype, position.svlen,
             max(starts) - min(starts), max(lengths) - min(lengths) if lengths else None,
             sv.gnomad_id, sv.gnomad_af, sv.colorsdb_id, sv.colorsdb_af,
             sv.near_reference_gap),
        )
        db.executemany(
            "INSERT INTO calls (sv_id, sample, caller, vcf_id, chrom, pos, vcf_end,"
            " svtype, svlen, genotype, read_support)"
            " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            [(sv.sv_id, call.sample, call.caller, call.vcf_id, call.chrom, call.pos,
              call.vcf_end, call.svtype, call.svlen, call.genotype, call.read_support)
             for call in sv.calls],
        )


# Insert the transcripts near the SVs, once each, and the parts of them near
# each SV.
def _insert_transcripts(db, svs):
    transcripts = {part.transcript.transcript_id: part.transcript
                   for sv in svs for part in sv.transcript_parts}
    db.executemany(
        "INSERT INTO transcripts (transcript_id, gene_name, gene_id, transcript_type,"
        " mane_select, chrom, strand) VALUES (?, ?, ?, ?, ?, ?, ?)",
        [(transcript.transcript_id, transcript.gene_name, transcript.gene_id,
          transcript.transcript_type, transcript.mane_select, transcript.chrom,
          transcript.strand)
         for _, transcript in sorted(transcripts.items())],
    )
    db.executemany(
        "INSERT INTO sv_transcript_parts (sv_id, transcript_id, region, number, pos,"
        " end_pos, overlaps_sv) VALUES (?, ?, ?, ?, ?, ?, ?)",
        [(sv.sv_id, part.transcript.transcript_id, part.region, part.number, part.pos,
          part.end, part.overlaps_sv)
         for sv in svs for part in sv.transcript_parts],
    )


# Insert the details of each gene near the SVs, once each.
def _insert_genes(db, svs):
    genes = {part.transcript.gene_details.gene_id: part.transcript.gene_details
             for sv in svs for part in sv.transcript_parts if part.transcript.gene_details}
    db.executemany(
        "INSERT INTO genes (gene_id, hgnc_id, full_name, summary,"
        " haploinsufficiency_score, haploinsufficiency, triplosensitivity_score,"
        " triplosensitivity) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
        [(gene.gene_id, gene.hgnc_id, gene.full_name, gene.summary,
          gene.haploinsufficiency_score, gene.haploinsufficiency,
          gene.triplosensitivity_score, gene.triplosensitivity)
         for _, gene in sorted(genes.items())],
    )


# Insert each outside tool run, in the order they ran.
def _insert_tool_runs(db, tool_runs):
    db.executemany(
        "INSERT INTO commands (sample, tool, command) VALUES (?, ?, ?)",
        [(tool_run.sample, tool_run.tool, str(tool_run)) for tool_run in tool_runs],
    )
