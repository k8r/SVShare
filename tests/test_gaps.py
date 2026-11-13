import pysam

from svshare.gaps import add_reference_gaps, describe_nearest_gap, find_gaps
from svshare.svs import SV, SVPosition


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


# Each SV is checked from its position to where it ends on the reference, not
# to its VCF END.
def test_add_reference_gaps(tmp_path):
    svs = [
        SV(1, SVPosition("chr1", 50, 150, 150, "DEL", -100), []),
        SV(2, SVPosition("chr1", 800, 5000, 800, "INS", 300), []),
        SV(3, SVPosition("chr1", 1200, 1200, 1200, "INS", 300), []),
        SV(4, SVPosition("chr3", 10, 20, 20, "DEL", -10), []),
    ]

    assert add_reference_gaps(svs, write_reference(tmp_path)) == 2
    assert [sv.near_reference_gap for sv in svs] == [
        "overlaps 20 unknown bases", "680 bp from 20 unknown bases", None, None,
    ]
