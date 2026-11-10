import sqlite3

from helpers import make_database
from svshare.callers import CALLERS
from svshare.database import SVFilter, select_sv_ids
from svshare.population import DATASETS
from svshare.reporting import find_unknown_genes, make_sv_filter, parse_gene_list


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
#      gnomAD, near GENE1's promoter but not overlapping it
#   3  a DEL at chr2:1000-3000 called by both callers in the same sample, at 20%
#      in CoLoRSdb, overlapping GENE2's promoter
#   4  a BND at chr1:8000 called by both callers in different samples, at 2% in
#      gnomAD and 0.1% in CoLoRSdb
def write_test_database(tmp_path):
    return make_database(
        tmp_path,
        svs=[
            dict(sv_id=1, chrom="chr1", pos=1000, end_pos=2000, svtype="DEL"),
            dict(sv_id=2, chrom="chr1", pos=5000, end_pos=5001, svtype="INS", gnomad_af=0.01),
            dict(sv_id=3, chrom="chr2", pos=1000, end_pos=3000, svtype="DEL", colorsdb_af=0.2),
            # A breakend's END is its mate's position, so it's before pos here.
            dict(sv_id=4, chrom="chr1", pos=8000, end_pos=100, svtype="BND", gnomad_af=0.02,
                 colorsdb_af=0.001),
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


# Insertions and breakends cover just their position, so 4's END is ignored.
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
    with sqlite3.connect(db_path) as db:
        db.execute("UPDATE svs SET near_reference_gap = 'overlaps 20 unknown bases'"
                   " WHERE sv_id = 2")
    db.close()
    assert select_sv_ids(db_path, SVFilter(near_reference_gap=True)) == [2]
    assert select_sv_ids(db_path, SVFilter(near_reference_gap=False)) == [1, 3, 4]


# The report's defaults keep every SV, except that they check every population dataset.
def test_make_sv_filter_defaults():
    assert make_sv_filter() == SVFilter(population_datasets=tuple(DATASETS))


def test_make_sv_filter_turns_caller_support_into_callers():
    assert make_sv_filter(caller_support="any").callers == ()
    assert make_sv_filter(caller_support="cutesv").callers == ("cutesv",)
    assert make_sv_filter(caller_support="both").callers == tuple(CALLERS)


def test_make_sv_filter_turns_options_into_plain_values():
    sv_filter = make_sv_filter(sv_type="DEL", max_frequency=0.1, genes=["TP53"])
    assert (sv_filter.svtype, sv_filter.max_allele_frequency, sv_filter.genes) == (
        "DEL", 0.1, ("TP53",))
    assert make_sv_filter(sv_type="any").svtype is None
    assert make_sv_filter(near_reference_gap="include").near_reference_gap is None
    assert make_sv_filter(near_reference_gap="only").near_reference_gap is True
    assert make_sv_filter(near_reference_gap="exclude").near_reference_gap is False
