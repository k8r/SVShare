from svshare.genes import add_genes
from svshare.svs import SV, GeneDetails, SVPosition, Transcript, TranscriptPart

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


def make_transcript(transcript_id, gene_name, gene_id):
    return Transcript(transcript_id, gene_name, gene_id, "protein_coding", True, "chr1", "+")


# Checks that add_genes gives each transcript near an SV its gene's details,
# shared by transcripts of the same gene (GENE1 has two), and what each gets:
#   GENE1  everything in the gene info file, matched even though the transcript
#          gives its ID with a version (ENSG01.5); its "Not yet evaluated"
#          triplosensitivity score is kept as no score
#   GENE2  its name only, since the file has no summary or scores for it
#   GENE3  only its ID, since it isn't in the file
def test_add_genes_sets_details_for_each_gene_near_an_sv(tmp_path):
    t1 = make_transcript("T1", "GENE1", "ENSG01.5")
    t1b = make_transcript("T1b", "GENE1", "ENSG01.5")
    t2 = make_transcript("T2", "GENE2", "ENSG02.1")
    t3 = make_transcript("T3", "GENE3", "ENSG03.2")
    position = SVPosition("chr1", 1000, 2000, 2000, "DEL", -1000)
    svs = [
        SV(1, position, [], transcript_parts=[
            TranscriptPart(t1, "CDS", 1, 1000, 1100, True),
            TranscriptPart(t1b, "CDS", 1, 1000, 1100, True),
            TranscriptPart(t2, "intron", 1, 1500, 1600, True),
        ]),
        SV(2, position, [], transcript_parts=[TranscriptPart(t3, "CDS", 1, 1000, 1100, True)]),
    ]

    assert add_genes(svs, write_gene_info(tmp_path)) == (1, 1, 3)
    assert t1.gene_details == GeneDetails(
        "ENSG01.5", "HGNC:1", "gene one", "Gene one does things.",
        3, "Sufficient evidence for dosage pathogenicity", None, "Not yet evaluated",
    )
    assert t1b.gene_details is t1.gene_details
    assert t2.gene_details == GeneDetails("ENSG02.1", "HGNC:2", "gene two")
    assert t3.gene_details == GeneDetails("ENSG03.2")
