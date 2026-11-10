from helpers import make_database
from svshare.svs import SV, Call, Gene, _chrom_order, load_svs, sv_end


# A database with three SVs:
#   1  a DEL at chr10:500-600, called in HG002 and HG003
#   2  an INS at chr2:300, in the UTR and an intron of GENE1, at 0.5% in CoLoRSdb,
#      near a gap in the reference
#   3  a BND at chr2:100 with its mate at 9000, overlapping GENE2's promoter and
#      near GENE1, which it doesn't overlap
def write_test_database(tmp_path):
    return make_database(
        tmp_path,
        svs=[
            dict(sv_id=1, chrom="chr10", pos=500, vcf_end=600, svtype="DEL", svlen=-100),
            dict(sv_id=2, chrom="chr2", pos=300, vcf_end=300, svtype="INS", svlen=50,
                 colorsdb_af=0.005, near_reference_gap="350 bp from 20 unknown bases"),
            dict(sv_id=3, chrom="chr2", pos=100, vcf_end=9000, svtype="BND", svlen=0),
        ],
        calls=[
            dict(sv_id=1, sample="HG002", caller="cutesv", vcf_id="a", genotype="1/1",
                 read_support=4),
            dict(sv_id=1, sample="HG002", caller="sniffles2", vcf_id="b", genotype="0/1",
                 read_support=3),
            dict(sv_id=1, sample="HG003", caller="cutesv", vcf_id="c", genotype="./."),
            dict(sv_id=2, sample="HG002", caller="sniffles2", vcf_id="d", genotype="0|1",
                 read_support=6),
            dict(sv_id=2, sample="HG002", caller="cutesv", vcf_id="e", genotype="0|1",
                 read_support=6),
            dict(sv_id=3, sample="HG003", caller="sniffles2", vcf_id="f", genotype="0/1",
                 read_support=5),
        ],
        transcripts=[
            dict(transcript_id="T1", gene_name="GENE1", gene_id="G1",
                 transcript_type="protein_coding", mane_select=1, chrom="chr2", strand="+"),
            dict(transcript_id="T1b", gene_name="GENE1", gene_id="G1",
                 transcript_type="protein_coding", mane_select=1, chrom="chr2", strand="+"),
            dict(transcript_id="T2", gene_name="GENE2", gene_id="G2",
                 transcript_type="protein_coding", mane_select=1, chrom="chr2", strand="+"),
        ],
        sv_transcript_parts=[
            dict(sv_id=2, transcript_id="T1", region="intron", number=1, pos=250, end_pos=400,
                 overlaps_sv=1),
            dict(sv_id=2, transcript_id="T1b", region="UTR", number=1, pos=280, end_pos=320,
                 overlaps_sv=1),
            dict(sv_id=3, transcript_id="T1", region="promoter", pos=150, end_pos=249,
                 overlaps_sv=0),
            dict(sv_id=3, transcript_id="T2", region="promoter", pos=50, end_pos=150,
                 overlaps_sv=1),
        ],
        genes=[
            dict(gene_id="G1", full_name="gene one", haploinsufficiency="Sufficient evidence"),
            dict(gene_id="G2"),
        ],
    )


def test_svs_are_in_genome_order_and_only_the_given_ones(tmp_path):
    db_path = write_test_database(tmp_path)
    assert [sv.sv_id for sv in load_svs(db_path, [1, 2, 3])] == [3, 2, 1]
    assert [sv.sv_id for sv in load_svs(db_path, [1, 2])] == [2, 1]


# Only genes the SV overlaps are included, with the parts it overlaps across
# their transcripts.
def test_load_svs(tmp_path):
    svs = load_svs(write_test_database(tmp_path), [1, 2, 3])
    assert svs == [
        SV(3, "chr2", 100, None, "BND", 0, None, None,
           [Call("HG003", "sniffles2", "0/1", 5)],
           [Gene("GENE2", None, None, None, ["promoter"])]),
        SV(2, "chr2", 300, None, "INS", 50, None, 0.005,
           [Call("HG002", "cutesv", "0|1", 6), Call("HG002", "sniffles2", "0|1", 6)],
           [Gene("GENE1", "gene one", "Sufficient evidence", None, ["UTR", "intron"])],
           "350 bp from 20 unknown bases"),
        SV(1, "chr10", 500, 600, "DEL", -100, None, None,
           [Call("HG002", "cutesv", "1/1", 4), Call("HG002", "sniffles2", "0/1", 3),
            Call("HG003", "cutesv", "./.", None)],
           []),
    ]


def test_chromosomes_sort_by_number_then_x_y_m_then_name():
    chroms = ["chrUn_KI270302v1", "chrM", "chr10", "chrY", "chr1_KI270706v1_random", "chrX", "chr2"]
    assert sorted(chroms, key=_chrom_order) == [
        "chr2", "chr10", "chrX", "chrY", "chrM", "chr1_KI270706v1_random", "chrUn_KI270302v1"]


# Insertions, breakends, and SVs without an end sit at their position.
def test_sv_end():
    assert sv_end("DEL", 100, 250) == 250
    assert sv_end("INS", 100, 5000) == 100
    assert sv_end("BND", 100, 9000) == 100
    assert sv_end("DUP", 100, None) == 100
    assert sv_end("INV", 100, 90) == 100
