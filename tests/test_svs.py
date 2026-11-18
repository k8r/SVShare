from svshare.svs import SV, SVPosition, Transcript, TranscriptPart, reference_end


# Insertions, breakends, and SVs without an end sit at their position.
def test_reference_end():
    assert reference_end("DEL", 100, 250) == 250
    assert reference_end("INS", 100, 5000) == 100
    assert reference_end("BND", 100, 9000) == 100
    assert reference_end("DUP", 100, None) == 100
    assert reference_end("INV", 100, 90) == 100


# An SV with GENE1's intron overlapped and GENE2's promoter only nearby.
def make_sv_near_genes():
    gene1 = Transcript("T1", "GENE1", "ENSG01.5", "protein_coding", True, "chr1", "+")
    gene2 = Transcript("T2", "GENE2", "G2", "protein_coding", True, "chr1", "+")
    return SV(1, SVPosition("chr1", 1000, 2000, 2000, "DEL", -1000), [],
              overlapped_parts=[TranscriptPart(gene1, "intron", 1, 1500, 2500)],
              nearby_parts=[TranscriptPart(gene2, "promoter", None, 2500, 4499)])


def test_all_parts_has_overlapped_and_nearby_parts():
    sv = make_sv_near_genes()
    assert sv.all_parts() == [sv.overlapped_parts[0], sv.nearby_parts[0]]


# Genes match by name or Ensembl ID, in any case and with or without the ID's
# version, and a gene only nearby doesn't count.
def test_overlaps_gene():
    sv = make_sv_near_genes()
    assert sv.overlaps_gene("GENE1")
    assert sv.overlaps_gene("gene1")
    assert sv.overlaps_gene("ENSG01")
    assert sv.overlaps_gene("ensg01.5")
    assert sv.overlaps_gene("ENSG01.9")
    assert not sv.overlaps_gene("GENE2")
    assert not sv.overlaps_gene("GENE3")
