import sqlite3

from helpers import make_database
from svshare.genes import add_genes

HEADER = [
    "ensembl_gene_id", "hgnc_id", "symbol", "name", "summary",
    "haploinsufficiency_score", "haploinsufficiency",
    "triplosensitivity_score", "triplosensitivity",
]

# A gene info file with two genes: GENE1, with a summary and ClinGen scores, and
# GENE2, with neither, as for most genes.
GENE_INFO_ROWS = [
    ["ENSG01", "HGNC:1", "GENE1", "gene one", "Gene one does things.",
     "3", "Sufficient evidence for dosage pathogenicity", "Not yet evaluated",
     "Not yet evaluated"],
    ["ENSG02", "HGNC:2", "GENE2", "gene two", "", "", "", "", ""],
]


def write_gene_info(tmp_path):
    path = tmp_path / "gene_info.tsv"
    with open(path, "w") as f:
        f.write("# Built for a test.\n")
        for row in [HEADER, *GENE_INFO_ROWS]:
            f.write("\t".join(row) + "\n")
    return path


# Checks that add_genes stores one row per gene in transcripts (GENE1 has two
# transcripts), and what each row gets:
#   GENE1  everything in the gene info file, matched even though transcripts
#          gives its ID with a version (ENSG01.5); its "Not yet evaluated"
#          triplosensitivity score is stored as no score
#   GENE2  its name only, since the file has no summary or scores for it
#   GENE3  only its ID, since it isn't in the file
def test_add_genes_stores_info_for_each_gene_near_an_sv(tmp_path):
    db_path = make_database(tmp_path, transcripts=[
        dict(transcript_id="T1", gene_name="GENE1", gene_id="ENSG01.5",
             transcript_type="protein_coding", mane_select=1, chrom="chr1", strand="+"),
        dict(transcript_id="T1b", gene_name="GENE1", gene_id="ENSG01.5",
             transcript_type="protein_coding", mane_select=1, chrom="chr1", strand="+"),
        dict(transcript_id="T2", gene_name="GENE2", gene_id="ENSG02.1",
             transcript_type="protein_coding", mane_select=1, chrom="chr1", strand="+"),
        dict(transcript_id="T3", gene_name="GENE3", gene_id="ENSG03.2",
             transcript_type="protein_coding", mane_select=1, chrom="chr1", strand="+"),
    ])

    assert add_genes(db_path, write_gene_info(tmp_path)) == (1, 1, 3)
    assert sqlite3.connect(db_path).execute("SELECT * FROM genes ORDER BY gene_id").fetchall() == [
        ("ENSG01.5", "HGNC:1", "gene one", "Gene one does things.",
         3, "Sufficient evidence for dosage pathogenicity", None, "Not yet evaluated"),
        ("ENSG02.1", "HGNC:2", "gene two", None, None, None, None, None),
        ("ENSG03.2", None, None, None, None, None, None, None),
    ]
