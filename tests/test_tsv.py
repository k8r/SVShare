from svshare.svs import SV, Call, GeneDetails, SVPosition, Transcript, TranscriptPart
from svshare.tsv import COLUMNS, _row, write_tsv

GENE1 = GeneDetails("G1", full_name="gene one", haploinsufficiency="Sufficient evidence")
T1 = Transcript("T1", "GENE1", "G1", "protein_coding", True, "chr2", "+", GENE1)
T1B = Transcript("T1b", "GENE1", "G1", "lncRNA", False, "chr2", "+", GENE1)
T2 = Transcript("T2", "GENE2", "G2", "protein_coding", True, "chr2", "+", GeneDetails("G2"))


def make_sv(calls=(), transcript_parts=(), near_reference_gap=None):
    return SV(1, SVPosition("chr2", 300, 300, 300, "INS", 50), list(calls),
              colorsdb_af=0.005, near_reference_gap=near_reference_gap,
              transcript_parts=list(transcript_parts))


# Callers' genotypes are listed separately only where they differ, and unknown
# genotypes are left out.
def test_samples_list_each_sample_with_its_genotype():
    row = _row(make_sv(calls=[
        Call("HG002", "cutesv", "1/1"), Call("HG002", "sniffles2", "0/1"),
        Call("HG003", "cutesv", "./."),
        Call("HG004", "sniffles2", "0|1"), Call("HG004", "cutesv", "0|1"),
    ]))
    assert row["samples"] == "HG002 (0/1 sniffles2, 1/1 cutesv), HG003, HG004 (0|1)"
    assert row["callers"] == "sniffles2,cutesv"


def test_read_support_lists_each_call():
    row = _row(make_sv(calls=[
        Call("HG003", "sniffles2", "0/1", 6),
        Call("HG002", "cutesv", "1/1", 4), Call("HG002", "sniffles2", "0/1", 3),
        Call("HG004", "cutesv", "./.", None),
    ]))
    assert row["read_support"] == (
        "3 from sniffles2 for HG002; 4 from cutesv for HG002; 6 from sniffles2 for HG003")


# Each sample's best-supported call counts, and the lowest of those is reported.
def test_weakest_sample_support():
    row = _row(make_sv(calls=[
        Call("HG002", "cutesv", "1/1", 4), Call("HG002", "sniffles2", "0/1", 3),
        Call("HG003", "sniffles2", "0/1", 6), Call("HG004", "cutesv", "./.", None),
    ]))
    assert row["weakest_sample_support"] == 4
    assert _row(make_sv())["weakest_sample_support"] is None


# Each gene is listed once, with the kinds of part the SV overlaps across its
# transcripts, upstream first; parts it doesn't overlap are left out.
def test_genes_list_overlapped_parts_and_gene_info():
    row = _row(make_sv(transcript_parts=[
        TranscriptPart(T2, "promoter", None, 1, 400, True),
        TranscriptPart(T1, "intron", 1, 250, 400, True),
        TranscriptPart(T1B, "UTR", 1, 280, 320, True),
        TranscriptPart(T1B, "intron", 1, 321, 500, True),
        TranscriptPart(T1, "CDS", 2, 401, 600, False),
    ]))
    assert row["genes"] == "GENE1 (UTR, intron), GENE2 (promoter)"
    assert row["gene_names"] == "GENE1: gene one"
    assert row["haploinsufficiency"] == "GENE1: Sufficient evidence"
    assert row["triplosensitivity"] is None


# Kinds of part REGION_ORDER doesn't know are listed last instead of failing.
def test_unknown_regions_are_listed_last():
    row = _row(make_sv(transcript_parts=[
        TranscriptPart(T1, "exon", 1, 250, 300, True),
        TranscriptPart(T1, "intron", 1, 301, 400, True),
        TranscriptPart(T1, "CDS", 1, 401, 500, True),
    ]))
    assert row["genes"] == "GENE1 (CDS, intron, exon)"


def test_no_genes():
    row = _row(make_sv())
    assert row["genes"] is None
    assert row["gene_names"] is None


def test_write_tsv(tmp_path):
    output = tmp_path / "report.tsv"
    write_tsv([
        make_sv(calls=[Call("HG002", "sniffles2", "0|1", 7)],
                transcript_parts=[TranscriptPart(T1, "intron", 1, 250, 400, True)],
                near_reference_gap="350 bp from 20 unknown bases"),
        make_sv(),
    ], output)
    header, first, second = output.read_text().splitlines()
    assert header.split("\t") == COLUMNS
    assert first.split("\t") == [
        "1", "chr2", "300", "300", "INS", "50", "HG002 (0|1)", "7 from sniffles2 for HG002",
        "7", "350 bp from 20 unknown bases", "sniffles2", "", "0.005", "GENE1 (intron)",
        "GENE1: gene one", "GENE1: Sufficient evidence", "",
    ]
    assert second.split("\t")[0] == "2"
