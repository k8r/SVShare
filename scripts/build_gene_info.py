# Build svshare/data/gene_info.tsv, the gene descriptions and dosage
# sensitivity scores svshare/genes.py adds to each gene near an SV. Rerun it
# to refresh them.
#
# It joins three sources by gene:
#   HGNC          each gene's approved name, e.g. "ABR activator of RhoGEF and
#                 GTPase", and its Ensembl and NCBI Gene / Entrez IDs
#   NCBI Gene     a RefSeq summary of what the gene does
#   ClinGen       haploinsufficiency and triplosensitivity scores

import csv
import io
import json
import urllib.request
from datetime import date
from pathlib import Path

HGNC_URL = (
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
)
CLINGEN_URL = "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv"
NCBI_GENE_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/gene"
NCBI_BATCH_SIZE = 1000

OUT_PATH = Path(__file__).resolve().parent.parent / "svshare" / "data" / "gene_info.tsv"

COLUMNS = [
    "ensembl_gene_id", "hgnc_id", "symbol", "name", "summary",
    "haploinsufficiency_score", "haploinsufficiency",
    "triplosensitivity_score", "triplosensitivity",
]


def fetch_text(url, data=None):
    headers = {"Content-Type": "application/json"} if data else {}
    request = urllib.request.Request(url, data=data, headers=headers)
    with urllib.request.urlopen(request, timeout=120) as response:
        return response.read().decode("utf-8")


# Return HGNC's genes that have an Ensembl gene ID, as dicts of its columns.
def read_hgnc():
    rows = csv.DictReader(io.StringIO(fetch_text(HGNC_URL)), delimiter="\t")
    return [row for row in rows if row["ensembl_gene_id"]]


def read_clingen():
    lines = fetch_text(CLINGEN_URL).splitlines()
    # The header is the last of the comment lines at the top.
    header = max(i for i, line in enumerate(lines) if line.startswith("#"))
    rows = csv.DictReader([lines[header].lstrip("#"), *lines[header + 1:]], delimiter="\t")
    return {row["Gene ID"]: row for row in rows}


# Return {NCBI Gene ID / Entrez ID: summary} for the given IDs, asking NCBI in
# batches.
def read_ncbi_summaries(entrez_ids):
    summaries = {}
    for i in range(0, len(entrez_ids), NCBI_BATCH_SIZE):
        batch = entrez_ids[i:i + NCBI_BATCH_SIZE]
        body = json.dumps({"gene_ids": [int(x) for x in batch], "page_size": NCBI_BATCH_SIZE})
        reports = json.loads(fetch_text(NCBI_GENE_URL, body.encode()))["reports"]
        for report in reports:
            gene = report["gene"]
            text = " ".join(s.get("description", "") for s in gene.get("summary", []))
            if text:
                summaries[gene["gene_id"]] = text
        print(f"  NCBI summaries: {min(i + NCBI_BATCH_SIZE, len(entrez_ids))} of "
              f"{len(entrez_ids)} genes")
    return summaries


# Tabs and newlines would break the TSV.
def clean(value):
    return " ".join((value or "").split())


def main():
    print("Downloading HGNC")
    hgnc = read_hgnc()
    print("Downloading ClinGen")
    clingen = read_clingen()
    print("Downloading NCBI summaries")
    summaries = read_ncbi_summaries(sorted({row["entrez_id"] for row in hgnc if row["entrez_id"]}))

    OUT_PATH.parent.mkdir(exist_ok=True)
    with open(OUT_PATH, "w", newline="") as f:
        f.write(f"# Built {date.today()} by scripts/build_gene_info.py from HGNC, NCBI Gene "
                "and ClinGen.\n")
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        for row in sorted(hgnc, key=lambda row: row["ensembl_gene_id"]):
            # NCBI Gene and ClinGen identify genes by NCBI Gene ID / Entrez ID,
            # which HGNC gives alongside each gene's Ensembl ID.
            dosage = clingen.get(row["entrez_id"], {})
            writer.writerow([clean(value) for value in (
                row["ensembl_gene_id"], row["hgnc_id"], row["symbol"], row["name"],
                summaries.get(row["entrez_id"]),
                dosage.get("Haploinsufficiency Score"),
                dosage.get("Haploinsufficiency Description"),
                dosage.get("Triplosensitivity Score"),
                dosage.get("Triplosensitivity Description"),
            )])
    print(f"Wrote {len(hgnc)} genes ({len(summaries)} with a summary, "
          f"{sum(1 for row in hgnc if row['entrez_id'] in clingen)} with ClinGen scores) "
          f"to {OUT_PATH}")


if __name__ == "__main__":
    main()
