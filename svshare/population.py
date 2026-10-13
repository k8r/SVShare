# Look up each SV in population SV datasets, gnomAD and CoLoRSdb, to see how
# common it is. A record matches an SV when it's a compatible type and close enough.
import math
import sqlite3

import pysam

# The population datasets SVs can be looked up in, by the prefix of their
# columns in the svs table (e.g. gnomad_af), with the name to show for each.
DATASETS = {"gnomad": "gnomAD", "colorsdb": "CoLoRSdb"}

# How far apart (in bp) a record and an SV can be to count as the same SV:
# MAX_DIST_LINEAR times the SV's length, but at least MIN_DIST.
MAX_DIST_LINEAR = 0.5
MIN_DIST = 100

# Which population SV types each SV type can match. Short reads often show a
# duplication as a DUP where long reads show it as an INS, so the two can match.
# Other types, such as BND, aren't looked up.
COMPATIBLE_TYPES = {
    "DEL": {"DEL"},
    "INS": {"INS", "DUP"},
    "DUP": {"DUP", "INS"},
    "INV": {"INV"},
}


# Return the described SV's start and end. An insertion's end is taken as
# start + length, so insertions are compared by length as well as position and
# can match a duplication of the same sequence.
def sv_point(svtype, pos, end, svlen):
    if svtype == "INS" or end is None:
        return pos, pos + abs(svlen or 0)
    return pos, end


def max_distance(svlen):
    return max(MAX_DIST_LINEAR * abs(svlen or 0), MIN_DIST)


# The value of an INFO field that some VCFs declare as a list (one per ALT).
def _first(value):
    return value[0] if isinstance(value, tuple) else value


# Return (id, af) of the VCF record matching the SV with the highest allele
# frequency, or None if no record matches.
def best_match(vcf, chrom, pos, end, svtype, svlen):
    if svtype not in COMPATIBLE_TYPES or chrom not in vcf.header.contigs:
        return None
    point = sv_point(svtype, pos, end, svlen)
    limit = max_distance(svlen)

    best = None
    # A matching record starts within limit bp of the SV; fetch takes 0-based starts.
    for record in vcf.fetch(chrom, max(0, pos - 1 - int(limit)), pos + int(limit)):
        record_type = record.info.get("SVTYPE")
        # FILTER is PASS in gnomAD and empty in CoLoRSdb, which doesn't filter.
        if set(record.filter.keys()) - {"PASS"} or record_type not in COMPATIBLE_TYPES[svtype]:
            continue
        af = _first(record.info.get("AF"))
        if af is None:
            continue
        distance = math.dist(
            point,
            sv_point(record_type, record.pos, record.stop, _first(record.info.get("SVLEN"))),
        )
        if distance <= limit and (best is None or (af, -distance) > (best[1], -best[2])):
            # CoLoRSdb records have no ID, so label those by position and type.
            record_id = record.id or f"{record.chrom}:{record.pos}:{record_type}"
            best = (record_id, af, distance)
    return best[:2] if best else None


# Look up each SV in the database in a population dataset and save the matching
# record's ID and allele frequency, e.g. in gnomad_id and gnomad_af. dataset is
# a key of DATASETS, and vcf_path needs its .tbi index next to it.
# Returns (SVs matched, total SVs).
def add_population_frequencies(db_path, dataset, vcf_path):
    if dataset not in DATASETS:
        raise ValueError(f"Unknown population dataset: {dataset}")
    with pysam.VariantFile(str(vcf_path)) as vcf, sqlite3.connect(db_path) as db:
        svs = db.execute("SELECT sv_id, chrom, pos, end_pos, svtype, svlen FROM svs").fetchall()
        matches = []
        for sv_id, *sv in svs:
            match = best_match(vcf, *sv)
            if match:
                matches.append((*match, sv_id))
        db.executemany(
            f"UPDATE svs SET {dataset}_id = ?, {dataset}_af = ? WHERE sv_id = ?", matches
        )
    db.close()
    return len(matches), len(svs)
