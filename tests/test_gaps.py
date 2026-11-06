import sqlite3

import pysam

from svshare.database import SCHEMA
from svshare.gaps import add_reference_gaps, describe_nearest_gap, find_gaps


# chr1 has a 20-base gap at 101-120 and a run of 5 N's too short to count; chr2
# has a gap at its start.
def write_reference(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(
        ">chr1\n" + "A" * 100 + "N" * 20 + "C" * 100 + "NNNNN" + "G" * 2000 + "\n"
        ">chr2\n" + "n" * 15 + "T" * 50 + "\n"
    )
    pysam.faidx(str(reference))
    return reference


def test_find_gaps(tmp_path):
    assert find_gaps(write_reference(tmp_path)) == {"chr1": [(101, 120)], "chr2": [(1, 15)]}


def test_describe_nearest_gap():
    gaps = [(101, 120), (5001, 5100)]
    assert describe_nearest_gap(gaps, 90, 105) == "overlaps 20 unknown bases"
    assert describe_nearest_gap(gaps, 50, 5050) == "overlaps 120 unknown bases"
    assert describe_nearest_gap(gaps, 470, 470) == "350 bp from 20 unknown bases"
    assert describe_nearest_gap(gaps, 4000, 4500) == "501 bp from 100 unknown bases"
    assert describe_nearest_gap(gaps, 2000, 3000) is None
    assert describe_nearest_gap([], 100, 200) is None


# Insertions are checked at their position only, since their end position
# isn't where they end.
def test_add_reference_gaps(tmp_path):
    db_path = tmp_path / "svshare.db"
    with sqlite3.connect(db_path) as db:
        db.executescript(SCHEMA)
        db.executemany(
            "INSERT INTO svs (sv_id, chrom, pos, end_pos, svtype) VALUES (?, ?, ?, ?, ?)",
            [(1, "chr1", 50, 150, "DEL"), (2, "chr1", 800, 5000, "INS"),
             (3, "chr1", 1200, 1200, "INS"), (4, "chr3", 10, 20, "DEL")],
        )
    db.close()

    assert add_reference_gaps(db_path, write_reference(tmp_path)) == (2, 4)
    with sqlite3.connect(db_path) as db:
        near = dict(db.execute("SELECT sv_id, near_reference_gap FROM svs"))
    db.close()
    assert near == {1: "overlaps 20 unknown bases", 2: "680 bp from 20 unknown bases",
                    3: None, 4: None}
