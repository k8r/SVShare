import sqlite3

from svshare.database import SCHEMA
from svshare.filtering import filter_sv_ids, parse_gene_list


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


# A database with four SVs: 1, a DEL, called only by sniffles2; 2, an INS, only
# by cuteSV (in two samples); 3, a DEL, by both in the same sample; and 4, a BND,
# by both in different samples.
def write_test_database(tmp_path):
    db_path = tmp_path / "svshare.db"
    with sqlite3.connect(db_path) as db:
        db.executescript(SCHEMA)
        db.executemany(
            "INSERT INTO svs (sv_id, svtype) VALUES (?, ?)",
            [(1, "DEL"), (2, "INS"), (3, "DEL"), (4, "BND")],
        )
        db.executemany(
            "INSERT INTO calls (sv_id, sample, caller, vcf_id) VALUES (?, ?, ?, ?)",
            [(1, "HG002", "sniffles2", "a"),
             (2, "HG002", "cutesv", "b"), (2, "HG003", "cutesv", "c"),
             (3, "HG002", "sniffles2", "d"), (3, "HG002", "cutesv", "e"),
             (4, "HG002", "sniffles2", "f"), (4, "HG003", "cutesv", "g")],
        )
    db.close()
    return db_path


def test_caller_support_any_keeps_every_sv(tmp_path):
    assert filter_sv_ids(write_test_database(tmp_path), caller_support="any") == [1, 2, 3, 4]


# Naming one caller keeps SVs the other caller also called.
def test_caller_support_one_caller_keeps_svs_it_called(tmp_path):
    db_path = write_test_database(tmp_path)
    assert filter_sv_ids(db_path, caller_support="sniffles2") == [1, 3, 4]
    assert filter_sv_ids(db_path, caller_support="cutesv") == [2, 3, 4]


# Both callers count whether they called the SV in the same sample or not.
def test_caller_support_both_keeps_svs_every_caller_called(tmp_path):
    assert filter_sv_ids(write_test_database(tmp_path), caller_support="both") == [3, 4]


def test_sv_type_keeps_svs_of_that_type(tmp_path):
    db_path = write_test_database(tmp_path)
    assert filter_sv_ids(db_path, sv_type="any") == [1, 2, 3, 4]
    assert filter_sv_ids(db_path, sv_type="DEL") == [1, 3]
    assert filter_sv_ids(db_path, sv_type="DUP") == []


def test_filters_combine(tmp_path):
    assert filter_sv_ids(write_test_database(tmp_path), caller_support="both", sv_type="DEL") == [3]
