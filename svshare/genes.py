# Functions for gene information: names, summaries, and ClinGen dosage
# sensitivity scores.
import csv
import sqlite3
from pathlib import Path

GENE_INFO = Path(__file__).parent / "data" / "gene_info.tsv"


# ClinGen scores are 0-3 for how strong the evidence is, or 30 (autosomal
# recessive) or 40 (dosage sensitivity unlikely); anything else, such as "Not
# yet evaluated", is stored as None.
def _score(value):
    return int(value) if value.isdigit() else None


# Return {Ensembl gene ID without its version: row} from a gene info file.
def read_gene_info(path):
    with open(path, newline="") as f:
        lines = (line for line in f if not line.startswith("#"))
        return {row["ensembl_gene_id"]: row for row in csv.DictReader(lines, delimiter="\t")}


# Store a row in genes for each gene in transcripts, with whatever the gene info
# file has for it. Returns (genes with a summary, genes with ClinGen scores,
# total genes).
def add_genes(db_path, gene_info=GENE_INFO):
    info = read_gene_info(gene_info)
    rows = []
    with sqlite3.connect(db_path) as db:
        for (gene_id,) in db.execute("SELECT DISTINCT gene_id FROM transcripts ORDER BY gene_id"):
            # GENCODE's gene IDs carry a version, e.g. ENSG00000159842.18.
            gene = info.get(gene_id.split(".")[0])
            if gene is None:
                rows.append((gene_id, None, None, None, None, None, None, None))
                continue
            rows.append((
                gene_id, gene["hgnc_id"], gene["name"], gene["summary"] or None,
                _score(gene["haploinsufficiency_score"]), gene["haploinsufficiency"] or None,
                _score(gene["triplosensitivity_score"]), gene["triplosensitivity"] or None,
            ))
        db.executemany(
            "INSERT INTO genes (gene_id, hgnc_id, full_name, summary,"
            " haploinsufficiency_score, haploinsufficiency, triplosensitivity_score,"
            " triplosensitivity) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            rows,
        )
    db.close()
    with_summary = sum(1 for row in rows if row[3])
    with_dosage = sum(1 for row in rows if row[5])
    return with_summary, with_dosage, len(rows)
