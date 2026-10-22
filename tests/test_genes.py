import sqlite3

from svshare.database import SCHEMA
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
    db_path = tmp_path / "svshare.db"
    with sqlite3.connect(db_path) as db:
        db.executescript(SCHEMA)
        db.executemany(
            "INSERT INTO transcripts (transcript_id, gene_name, gene_id, transcript_type,"
            " mane_select, chrom, strand) VALUES (?, ?, ?, 'protein_coding', 1, 'chr1', '+')",
            [("T1", "GENE1", "ENSG01.5"), ("T1b", "GENE1", "ENSG01.5"),
             ("T2", "GENE2", "ENSG02.1"), ("T3", "GENE3", "ENSG03.2")],
        )
    db.close()

    assert add_genes(db_path, write_gene_info(tmp_path)) == (1, 1, 3)
    assert sqlite3.connect(db_path).execute("SELECT * FROM genes ORDER BY gene_id").fetchall() == [
        ("ENSG01.5", "HGNC:1", "gene one", "Gene one does things.",
         3, "Sufficient evidence for dosage pathogenicity", None, "Not yet evaluated"),
        ("ENSG02.1", "HGNC:2", "gene two", None, None, None, None, None),
        ("ENSG03.2", None, None, None, None, None, None, None),
    ]
