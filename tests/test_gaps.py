import pysam

from helpers import make_database, read_rows
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
    db_path = make_database(tmp_path, svs=[
        dict(sv_id=1, chrom="chr1", pos=50, vcf_end=150, svtype="DEL"),
        dict(sv_id=2, chrom="chr1", pos=800, vcf_end=5000, svtype="INS"),
        dict(sv_id=3, chrom="chr1", pos=1200, vcf_end=1200, svtype="INS"),
        dict(sv_id=4, chrom="chr3", pos=10, vcf_end=20, svtype="DEL"),
    ])

    assert add_reference_gaps(db_path, write_reference(tmp_path)) == (2, 4)
    assert read_rows(db_path, "svs", "sv_id", "near_reference_gap") == [
        (1, "overlaps 20 unknown bases"), (2, "680 bp from 20 unknown bases"), (3, None), (4, None),
    ]
