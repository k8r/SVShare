from svshare.svs import SV, Call, Gene
from svshare.tsv import COLUMNS, _row, write_tsv


def make_sv(calls=(), genes=()):
    return SV(1, "chr2", 300, None, "INS", 50, None, 0.005, list(calls), list(genes))


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
    write_tsv([make_sv(
        calls=[Call("HG002", "sniffles2", "0|1")],
        genes=[Gene("GENE1", "gene one", None, None, ["intron"])],
    )], output)
    header, row = output.read_text().splitlines()
    assert header.split("\t") == COLUMNS
    assert row.split("\t") == [
        "chr2", "300", "", "INS", "50", "HG002 (0|1)", "sniffles2", "", "0.005",
        "GENE1 (intron)", "GENE1: gene one", "", "",
    ]
