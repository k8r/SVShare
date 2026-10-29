import sqlite3

from svshare.database import SCHEMA
from svshare.svs import SV, Call, Gene, _chrom_order, load_svs


# A database with three SVs:
#   1  a DEL at chr10:500-600, called in HG002 and HG003
#   2  an INS at chr2:300, in the UTR and an intron of GENE1, at 0.5% in CoLoRSdb
#   3  a BND at chr2:100 with its mate at 9000, overlapping GENE2's promoter and
#      near GENE1, which it doesn't overlap
def write_test_database(tmp_path):
    db_path = tmp_path / "svshare.db"
    with sqlite3.connect(db_path) as db:
        db.executescript(SCHEMA)
        db.executemany(
            "INSERT INTO svs (sv_id, chrom, pos, end_pos, svtype, svlen, gnomad_af, colorsdb_af)"
            " VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            [(1, "chr10", 500, 600, "DEL", -100, None, None),
             (2, "chr2", 300, 300, "INS", 50, None, 0.005),
             (3, "chr2", 100, 9000, "BND", 0, None, None)],
        )
        db.executemany(
            "INSERT INTO calls (sv_id, sample, caller, vcf_id, genotype) VALUES (?, ?, ?, ?, ?)",
            [(1, "HG002", "cutesv", "a", "1/1"), (1, "HG002", "sniffles2", "b", "0/1"),
             (1, "HG003", "cutesv", "c", "./."),
             (2, "HG002", "sniffles2", "d", "0|1"), (2, "HG002", "cutesv", "e", "0|1"),
             (3, "HG003", "sniffles2", "f", "0/1")],
        )
        db.executemany(
            "INSERT INTO transcripts (transcript_id, gene_name, gene_id, transcript_type,"
            " mane_select, chrom, strand) VALUES (?, ?, ?, 'protein_coding', 1, 'chr2', '+')",
            [("T1", "GENE1", "G1"), ("T1b", "GENE1", "G1"), ("T2", "GENE2", "G2")],
        )
        db.executemany(
            "INSERT INTO sv_transcript_parts (sv_id, transcript_id, region, number, pos,"
            " end_pos, overlaps_sv) VALUES (?, ?, ?, ?, ?, ?, ?)",
            [(2, "T1", "intron", 1, 250, 400, 1), (2, "T1b", "UTR", 1, 280, 320, 1),
             (3, "T1", "promoter", None, 150, 249, 0),
             (3, "T2", "promoter", None, 50, 150, 1)],
        )
        db.executemany(
            "INSERT INTO genes (gene_id, full_name, haploinsufficiency, triplosensitivity)"
            " VALUES (?, ?, ?, ?)",
            [("G1", "gene one", "Sufficient evidence", None), ("G2", None, None, None)],
        )
    db.close()
    return db_path


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
           [Call("HG003", "sniffles2", "0/1")],
           [Gene("GENE2", None, None, None, ["promoter"])]),
        SV(2, "chr2", 300, None, "INS", 50, None, 0.005,
           [Call("HG002", "cutesv", "0|1"), Call("HG002", "sniffles2", "0|1")],
           [Gene("GENE1", "gene one", "Sufficient evidence", None, ["UTR", "intron"])]),
        SV(1, "chr10", 500, 600, "DEL", -100, None, None,
           [Call("HG002", "cutesv", "1/1"), Call("HG002", "sniffles2", "0/1"),
            Call("HG003", "cutesv", "./.")],
           []),
    ]


def test_chromosomes_sort_by_number_then_x_y_m_then_name():
    chroms = ["chrUn_KI270302v1", "chrM", "chr10", "chrY", "chr1_KI270706v1_random", "chrX", "chr2"]
    assert sorted(chroms, key=_chrom_order) == [
        "chr2", "chr10", "chrX", "chrY", "chrM", "chr1_KI270706v1_random", "chrUn_KI270302v1"]
