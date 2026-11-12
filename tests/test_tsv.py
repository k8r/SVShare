from svshare.svs import SV, Call, Gene, SVPosition
from svshare.tsv import COLUMNS, _row, write_tsv


def make_sv(calls=(), genes=(), near_reference_gap=None):
    return SV(SVPosition(1, "chr2", 300, 300, 300, "INS", 50), None, 0.005, list(calls),
              list(genes), near_reference_gap)


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


def test_genes_list_overlapped_parts_and_gene_info():
    row = _row(make_sv(genes=[
        Gene("GENE1", "gene one", "Sufficient evidence", None, ["UTR", "intron"]),
        Gene("GENE2", None, None, None, ["promoter"]),
    ]))
    assert row["genes"] == "GENE1 (UTR, intron), GENE2 (promoter)"
    assert row["gene_names"] == "GENE1: gene one"
    assert row["haploinsufficiency"] == "GENE1: Sufficient evidence"
    assert row["triplosensitivity"] is None


def test_no_genes():
    row = _row(make_sv())
    assert row["genes"] is None
    assert row["gene_names"] is None


def test_write_tsv(tmp_path):
    output = tmp_path / "report.tsv"
    write_tsv([
        make_sv(calls=[Call("HG002", "sniffles2", "0|1", 7)],
                genes=[Gene("GENE1", "gene one", None, None, ["intron"])],
                near_reference_gap="350 bp from 20 unknown bases"),
        make_sv(),
    ], output)
    header, first, second = output.read_text().splitlines()
    assert header.split("\t") == COLUMNS
    assert first.split("\t") == [
        "1", "chr2", "300", "300", "INS", "50", "HG002 (0|1)", "7 from sniffles2 for HG002",
        "7", "350 bp from 20 unknown bases", "sniffles2", "", "0.005", "GENE1 (intron)",
        "GENE1: gene one", "", "",
    ]
    assert second.split("\t")[0] == "2"
