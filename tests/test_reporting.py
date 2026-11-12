from helpers import make_database
from svshare.callers import CALLERS
from svshare.database import SVFilter, select_sv_ids
from svshare.population import DATASETS
from svshare.reporting import (_chrom_order, find_unknown_genes, load_svs, make_sv_filter,
                              parse_gene_list)
from svshare.svs import SV, Call, Gene, SVPosition


def test_genes_only():
    assert parse_gene_list(["TP53", "BRCA1"], None) == ["BRCA1", "TP53"]


def test_gene_file_only(tmp_path):
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text("BRCA2\nPALB2\n")

    assert parse_gene_list(None, gene_file) == ["BRCA2", "PALB2"]


def test_genes_and_gene_file_are_combined_and_deduplicated(tmp_path):
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text("BRCA1\nPALB2\n")

    assert parse_gene_list(["TP53", "BRCA1"], gene_file) == ["BRCA1", "PALB2", "TP53"]


def test_gene_file_skips_blank_lines_and_strips_whitespace(tmp_path):
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text("  BRCA1  \n\nTP53\n\n")

    assert parse_gene_list(None, gene_file) == ["BRCA1", "TP53"]


def test_no_genes_or_gene_file_returns_empty_list():
    assert parse_gene_list(None, None) == []


# A database with four SVs:
#   1  a DEL at chr1:1000-2000 called only by sniffles2, in neither population
#      database, overlapping an intron of GENE1
#   2  an INS at chr1:5000 called only by cuteSV (in two samples), at 1% in
#      gnomAD, near GENE1's promoter but not overlapping it, and near a gap in
#      the reference
#   3  a DEL at chr2:1000-3000 called by both callers in the same sample, at 20%
#      in CoLoRSdb, overlapping GENE2's promoter
#   4  a BND at chr1:8000 called by both callers in different samples, at 2% in
#      gnomAD and 0.1% in CoLoRSdb
def write_test_database(tmp_path):
    return make_database(
        tmp_path,
        svs=[
            dict(sv_id=1, chrom="chr1", pos=1000, vcf_end=2000, reference_end=2000, svtype="DEL"),
            dict(sv_id=2, chrom="chr1", pos=5000, vcf_end=5001, reference_end=5000, svtype="INS",
                 gnomad_af=0.01, near_reference_gap="overlaps 20 unknown bases"),
            dict(sv_id=3, chrom="chr2", pos=1000, vcf_end=3000, reference_end=3000, svtype="DEL",
                 colorsdb_af=0.2),
            # A breakend's END is its mate's position, so it's before pos here.
            dict(sv_id=4, chrom="chr1", pos=8000, vcf_end=100, reference_end=8000, svtype="BND",
                 gnomad_af=0.02, colorsdb_af=0.001),
        ],
        transcripts=[
            dict(transcript_id="T1", gene_name="GENE1", gene_id="G1",
                 transcript_type="protein_coding", mane_select=1, chrom="chr1", strand="+"),
            dict(transcript_id="T2", gene_name="GENE2", gene_id="G2",
                 transcript_type="protein_coding", mane_select=1, chrom="chr2", strand="+"),
        ],
        sv_transcript_parts=[
            dict(sv_id=1, transcript_id="T1", region="intron", number=1, pos=1500, end_pos=2500,
                 overlaps_sv=1),
            dict(sv_id=2, transcript_id="T1", region="promoter", pos=5500, end_pos=7499,
                 overlaps_sv=0),
            dict(sv_id=3, transcript_id="T2", region="promoter", pos=2500, end_pos=4499,
                 overlaps_sv=1),
        ],
        calls=[
            dict(sv_id=1, sample="HG002", caller="sniffles2", vcf_id="a"),
            dict(sv_id=2, sample="HG002", caller="cutesv", vcf_id="b"),
            dict(sv_id=2, sample="HG003", caller="cutesv", vcf_id="c"),
            dict(sv_id=3, sample="HG002", caller="sniffles2", vcf_id="d"),
            dict(sv_id=3, sample="HG002", caller="cutesv", vcf_id="e"),
            dict(sv_id=4, sample="HG002", caller="sniffles2", vcf_id="f"),
            dict(sv_id=4, sample="HG003", caller="cutesv", vcf_id="g"),
        ],
    )


def test_no_callers_keeps_every_sv(tmp_path):
    assert select_sv_ids(write_test_database(tmp_path), SVFilter()) == [1, 2, 3, 4]


# Naming one caller keeps SVs the other caller also called.
def test_one_caller_keeps_svs_it_called(tmp_path):
    db_path = write_test_database(tmp_path)
    assert select_sv_ids(db_path, SVFilter(callers=("sniffles2",))) == [1, 3, 4]
    assert select_sv_ids(db_path, SVFilter(callers=("cutesv",))) == [2, 3, 4]


# Both callers must call the SV in the same sample, so 4 (called by each in a
# different sample) doesn't pass.
def test_every_caller_must_call_the_sv_in_one_sample(tmp_path):
    sv_filter = SVFilter(callers=("sniffles2", "cutesv"))
    assert select_sv_ids(write_test_database(tmp_path), sv_filter) == [3]


def test_min_samples_with_any_caller(tmp_path):
    db_path = write_test_database(tmp_path)
    assert select_sv_ids(db_path, SVFilter(min_samples=1)) == [1, 2, 3, 4]
    assert select_sv_ids(db_path, SVFilter(min_samples=2)) == [2, 4]
    assert select_sv_ids(db_path, SVFilter(min_samples=3)) == []


# Only samples where every listed caller called the SV count toward min_samples.
def test_min_samples_counts_samples_where_the_callers_called_it(tmp_path):
    db_path = write_test_database(tmp_path)
    assert select_sv_ids(db_path, SVFilter(min_samples=2, callers=("cutesv",))) == [2]
    assert select_sv_ids(db_path, SVFilter(min_samples=2, callers=("sniffles2",))) == []
    assert select_sv_ids(db_path, SVFilter(min_samples=2, callers=("sniffles2", "cutesv"))) == []


def test_svtype_keeps_svs_of_that_type(tmp_path):
    db_path = write_test_database(tmp_path)
    assert select_sv_ids(db_path, SVFilter(svtype="DEL")) == [1, 3]
    assert select_sv_ids(db_path, SVFilter(svtype="DUP")) == []


def test_filters_combine(tmp_path):
    sv_filter = SVFilter(callers=("sniffles2", "cutesv"), svtype="DEL")
    assert select_sv_ids(write_test_database(tmp_path), sv_filter) == [3]


# An SV passes if it's at most the frequency in every database that has a match
# for it, so 1 (no match) and 2 (exactly 1%) pass, and 3 and 4 (common in one
# database) don't.
def test_max_allele_frequency_keeps_svs_rare_in_every_database(tmp_path):
    sv_filter = SVFilter(population_datasets=("gnomad", "colorsdb"), max_allele_frequency=0.01)
    assert select_sv_ids(write_test_database(tmp_path), sv_filter) == [1, 2]


# Only parts the SV overlaps count, so 2 (near GENE1's promoter) doesn't pass.
def test_genes_keeps_svs_overlapping_a_listed_gene(tmp_path):
    db_path = write_test_database(tmp_path)
    assert select_sv_ids(db_path, SVFilter(genes=("GENE1",))) == [1]
    assert select_sv_ids(db_path, SVFilter(genes=("gene2",))) == [3]
    assert select_sv_ids(db_path, SVFilter(genes=("GENE1", "GENE2"))) == [1, 3]


# Regions are compared with where each SV ends on the reference, so 4's END is
# ignored.
def test_region_keeps_svs_overlapping_it(tmp_path):
    db_path = write_test_database(tmp_path)
    assert select_sv_ids(db_path, SVFilter(region=("chr1", 1500, 5000))) == [1, 2]
    assert select_sv_ids(db_path, SVFilter(region=("chr1", 2001, 4999))) == []
    assert select_sv_ids(db_path, SVFilter(region=("chr1", 1, 999))) == []
    assert select_sv_ids(db_path, SVFilter(region=("chr1", 8000, 8000))) == [4]
    assert select_sv_ids(db_path, SVFilter(region=("chr2", 1, 1000))) == [3]


def write_test_gene_info(tmp_path):
    gene_info = tmp_path / "gene_info.tsv"
    gene_info.write_text(
        "# test\n"
        "ensembl_gene_id\thgnc_id\tsymbol\tname\tsummary\thaploinsufficiency_score"
        "\thaploinsufficiency\ttriplosensitivity_score\ttriplosensitivity\n"
        "ENSG1\tHGNC:1\tBRCA1\t\t\t\t\t\t\n"
        "ENSG2\tHGNC:2\tTP53\t\t\t\t\t\t\n"
    )
    return gene_info


# Symbols from the gene info file and gene names in the database are both
# known, in any case.
def test_find_unknown_genes_accepts_known_symbols(tmp_path):
    db_path, gene_info = write_test_database(tmp_path), write_test_gene_info(tmp_path)
    assert find_unknown_genes(["brca1", "TP53", "GENE1"], db_path, gene_info) == []


def test_find_unknown_genes_returns_unknown_symbols(tmp_path):
    db_path, gene_info = write_test_database(tmp_path), write_test_gene_info(tmp_path)
    assert find_unknown_genes(["BRCA1", "brac1", "TP35"], db_path, gene_info) == ["brac1", "TP35"]


def test_near_reference_gap_keeps_only_those_svs_or_only_the_rest(tmp_path):
    db_path = write_test_database(tmp_path)
    assert select_sv_ids(db_path, SVFilter(near_reference_gap=True)) == [2]
    assert select_sv_ids(db_path, SVFilter(near_reference_gap=False)) == [1, 3, 4]


# Call make_sv_filter with the given options, and options that keep every SV for
# the rest.
def sv_filter_from(caller_support="any", min_samples=1, sv_type="any", max_frequency=None,
                   genes=(), region=None, near_reference_gap="include"):
    return make_sv_filter(caller_support=caller_support, min_samples=min_samples,
                          sv_type=sv_type, max_frequency=max_frequency, genes=genes,
                          region=region, near_reference_gap=near_reference_gap)


def test_make_sv_filter_turns_caller_support_into_callers():
    assert sv_filter_from(caller_support="any").callers == ()
    assert sv_filter_from(caller_support="cutesv").callers == ("cutesv",)
    assert sv_filter_from(caller_support="both").callers == tuple(CALLERS)


def test_make_sv_filter_turns_options_into_plain_values():
    sv_filter = sv_filter_from(sv_type="DEL", max_frequency=0.1, genes=["TP53"])
    assert (sv_filter.svtype, sv_filter.max_allele_frequency, sv_filter.genes) == (
        "DEL", 0.1, ("TP53",))
    assert sv_filter.population_datasets == tuple(DATASETS)
    assert sv_filter_from(sv_type="any").svtype is None
    assert sv_filter_from(near_reference_gap="include").near_reference_gap is None
    assert sv_filter_from(near_reference_gap="only").near_reference_gap is True
    assert sv_filter_from(near_reference_gap="exclude").near_reference_gap is False


# A database with three SVs:
#   1  a DEL at chr10:500-600, called in HG002 and HG003
#   2  an INS at chr2:300, in the UTR and an intron of GENE1, at 0.5% in CoLoRSdb,
#      near a gap in the reference
#   3  a BND at chr2:100 with its mate at 9000, overlapping GENE2's promoter and
#      near GENE1, which it doesn't overlap
def write_load_test_database(tmp_path):
    return make_database(
        tmp_path,
        svs=[
            dict(sv_id=1, chrom="chr10", pos=500, vcf_end=600, reference_end=600, svtype="DEL",
                 svlen=-100),
            dict(sv_id=2, chrom="chr2", pos=300, vcf_end=300, reference_end=300, svtype="INS",
                 svlen=50, colorsdb_af=0.005, near_reference_gap="350 bp from 20 unknown bases"),
            dict(sv_id=3, chrom="chr2", pos=100, vcf_end=9000, reference_end=100, svtype="BND",
                 svlen=0),
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
    db_path = write_load_test_database(tmp_path)
    assert [sv.position.sv_id for sv in load_svs(db_path, [1, 2, 3])] == [3, 2, 1]
    assert [sv.position.sv_id for sv in load_svs(db_path, [1, 2])] == [2, 1]


# Only genes the SV overlaps are included, with the parts it overlaps across
# their transcripts.
def test_load_svs(tmp_path):
    svs = load_svs(write_load_test_database(tmp_path), [1, 2, 3])
    assert svs == [
        SV(SVPosition(3, "chr2", 100, 9000, 100, "BND", 0), None, None,
           [Call("HG003", "sniffles2", "0/1", 5)],
           [Gene("GENE2", None, None, None, ["promoter"])]),
        SV(SVPosition(2, "chr2", 300, 300, 300, "INS", 50), None, 0.005,
           [Call("HG002", "cutesv", "0|1", 6), Call("HG002", "sniffles2", "0|1", 6)],
           [Gene("GENE1", "gene one", "Sufficient evidence", None, ["UTR", "intron"])],
           "350 bp from 20 unknown bases"),
        SV(SVPosition(1, "chr10", 500, 600, 600, "DEL", -100), None, None,
           [Call("HG002", "cutesv", "1/1", 4), Call("HG002", "sniffles2", "0/1", 3),
            Call("HG003", "cutesv", "./.", None)],
           []),
    ]


def test_chromosomes_sort_by_number_then_x_y_m_then_name():
    chroms = ["chrUn_KI270302v1", "chrM", "chr10", "chrY", "chr1_KI270706v1_random", "chrX", "chr2"]
    assert sorted(chroms, key=_chrom_order) == [
        "chr2", "chr10", "chrX", "chrY", "chrM", "chr1_KI270706v1_random", "chrUn_KI270302v1"]
