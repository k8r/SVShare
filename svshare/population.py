# Look up each SV in a population SV dataset, such as gnomAD, to see how common
# it is. A record matches an SV when it's a compatible type and close enough.
import math
import sqlite3

import pysam

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


# Return (id, af) of the PASS record in the VCF that matches the SV, taking the
# highest allele frequency if several do (the closest, if tied), or None if none
# do.
def best_match(vcf, chrom, pos, end, svtype, svlen):
    if svtype not in COMPATIBLE_TYPES or chrom not in vcf.header.contigs:
        return None
    point = sv_point(svtype, pos, end, svlen)
    limit = max_distance(svlen)

    best = None
    # A matching record starts within limit bp of the SV; fetch takes 0-based starts.
    for record in vcf.fetch(chrom, max(0, pos - 1 - int(limit)), pos + int(limit)):
        record_type = record.info.get("SVTYPE")
        if list(record.filter) != ["PASS"] or record_type not in COMPATIBLE_TYPES[svtype]:
            continue
        af = record.info.get("AF")
        if isinstance(af, tuple):
            af = af[0]
        if af is None:
            continue
        distance = math.dist(
            point, sv_point(record_type, record.pos, record.stop, record.info.get("SVLEN"))
        )
        if distance <= limit and (best is None or (af, -distance) > (best[1], -best[2])):
            best = (record.id, af, distance)
    return best[:2] if best else None


# Fill in gnomad_id and gnomad_af for every SV in the database from a gnomAD
# SV sites VCF, which needs its .tbi index next to it.
# Returns (SVs matched, total SVs).
def add_gnomad_frequencies(db_path, gnomad_vcf):
    with pysam.VariantFile(str(gnomad_vcf)) as vcf, sqlite3.connect(db_path) as db:
        svs = db.execute("SELECT sv_id, chrom, pos, end_pos, svtype, svlen FROM svs").fetchall()
        matches = []
        for sv_id, *sv in svs:
            match = best_match(vcf, *sv)
            if match:
                matches.append((*match, sv_id))
        db.executemany("UPDATE svs SET gnomad_id = ?, gnomad_af = ? WHERE sv_id = ?", matches)
    db.close()
    return len(matches), len(svs)
