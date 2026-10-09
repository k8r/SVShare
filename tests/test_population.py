import sqlite3

import pysam
import pytest

from svshare.database import SCHEMA
from svshare.population import add_gnomad_frequencies, best_match

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=UNRESOLVED,Description="Variant is unresolved">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


# A gnomAD-style sites VCF, compressed and indexed like the real one, with
# (id, pos, filter, svtype, end, svlen, af) records on chr1.
def write_gnomad(tmp_path, records):
    path = tmp_path / "gnomad.vcf"
    path.write_text(HEADER + "".join(
        f"chr1\t{pos}\t{id_}\tN\t<{svtype}>\t.\t{filter_}\t"
        f"SVTYPE={svtype};END={end};SVLEN={svlen};AF={af}\n"
        # Indexing needs the records in position order.
        for id_, pos, filter_, svtype, end, svlen, af in sorted(records, key=lambda r: r[1])
    ))
    return pysam.VariantFile(pysam.tabix_index(str(path), preset="vcf", force=True))


# A 150 bp deletion at 1000-1150 can be up to 100 bp (min_dist) from a match.
# Checks that of the nearby PASS deletions, the one with the highest frequency
# wins, and that records too far away, of another type, or not PASS are ignored.
def test_best_match_takes_highest_frequency_among_nearby_matches(tmp_path):
    gnomad = write_gnomad(tmp_path, [
        ("near_rare", 1010, "PASS", "DEL", 1160, 150, 0.01),
        ("near_common", 1040, "PASS", "DEL", 1190, 150, 0.2),
        ("not_pass", 1000, "UNRESOLVED", "DEL", 1150, 150, 0.5),
        ("other_type", 1000, "PASS", "INV", 1150, 150, 0.6),
        ("too_far", 1100, "PASS", "DEL", 1250, 150, 0.9),
    ])

    match_id, af = best_match(gnomad, "chr1", 1000, 1150, "DEL", -150)
    assert match_id == "near_common"
    assert af == pytest.approx(0.2)


# Checks that a long-read insertion matches a short-read duplication of the same
# sequence, which starts at the insertion and spans its length.
def test_best_match_matches_insertion_to_duplication(tmp_path):
    gnomad = write_gnomad(tmp_path, [("dup", 5000, "PASS", "DUP", 5300, 300, 0.05)])

    assert best_match(gnomad, "chr1", 5000, 5000, "INS", 300)[0] == "dup"


# Checks that matched SVs get gnomAD's ID and frequency, and that the rest,
# including a BND, which isn't looked up, are left NULL.
def test_add_gnomad_frequencies_fills_matched_svs_only(tmp_path):
    write_gnomad(tmp_path, [("del", 1000, "PASS", "DEL", 1150, 150, 0.25)])
    db_path = tmp_path / "svshare.db"
    with sqlite3.connect(db_path) as db:
        db.executescript(SCHEMA)
        db.executemany(
            "INSERT INTO svs (sv_id, chrom, pos, end_pos, svtype, svlen) VALUES (?, ?, ?, ?, ?, ?)",
            [(1, "chr1", 1000, 1150, "DEL", -150),
             (2, "chr1", 50000, 50150, "DEL", -150),
             (3, "chr1", 1000, 1000, "BND", None)],
        )
    db.close()

    assert add_gnomad_frequencies(db_path, tmp_path / "gnomad.vcf.gz") == (1, 3)
    rows = sqlite3.connect(db_path).execute(
        "SELECT sv_id, gnomad_id, gnomad_af FROM svs ORDER BY sv_id"
    ).fetchall()
    assert rows == [(1, "del", 0.25), (2, None, None), (3, None, None)]
