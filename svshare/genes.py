# Functions for gene information: names, summaries, and ClinGen dosage
# sensitivity scores.
import csv
from pathlib import Path

from .svs import GeneDetails

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


# Return the gene symbols in a gene info file.
def read_gene_symbols(path=GENE_INFO):
    return {row["symbol"] for row in read_gene_info(path).values()}


# Set gene_details on each transcript near the SVs, with whatever the gene info
# file has for its gene. Transcripts of the same gene share one GeneDetails.
# Returns (genes with a summary, genes with ClinGen scores, total genes).
def add_genes(svs, gene_info=GENE_INFO):
    info = read_gene_info(gene_info)
    genes = {}
    for sv in svs:
        for part in sv.transcript_parts:
            transcript = part.transcript
            if transcript.gene_id not in genes:
                genes[transcript.gene_id] = _gene_details(transcript.gene_id, info)
            transcript.gene_details = genes[transcript.gene_id]
    with_summary = sum(1 for details in genes.values() if details.summary)
    with_dosage = sum(1 for details in genes.values() if details.haploinsufficiency)
    return with_summary, with_dosage, len(genes)


# Return what info, from read_gene_info, says about the gene with gene_id.
def _gene_details(gene_id, info):
    # GENCODE's gene IDs carry a version, e.g. ENSG00000159842.18.
    gene = info.get(gene_id.split(".")[0])
    if gene is None:
        return GeneDetails(gene_id)
    return GeneDetails(
        gene_id, gene["hgnc_id"], gene["name"], gene["summary"] or None,
        _score(gene["haploinsufficiency_score"]), gene["haploinsufficiency"] or None,
        _score(gene["triplosensitivity_score"]), gene["triplosensitivity"] or None,
    )
