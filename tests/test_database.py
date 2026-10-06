import sqlite3

from svshare.database import write_database

HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"


def write_vcf(path, records):
    path.write_text(HEADER + "".join("\t".join(map(str, r)) + "\n" for r in records))
    return path


# A Sniffles2 VCF, a cuteSV VCF, and a Jasmine merge of them: a deletion called
# in both HG002 (sniffles2) and HG003 (cuteSV), 10 bp apart, and an insertion
# only in HG003.
def write_inputs(tmp_path):
    sniffles = write_vcf(tmp_path / "HG002.sniffles2.vcf", [
        ["chr1", 1000, "Sniffles2.DEL.1", "N", "<DEL>", ".", "PASS",
         "SVTYPE=DEL;SVLEN=-150;END=1150", "GT:GQ:DR:DV", "0/1:60:40:22"],
    ])
    cutesv = write_vcf(tmp_path / "HG003.cutesv.vcf", [
        ["chr1", 1010, "cuteSV.DEL.0", "N", "<DEL>", ".", "PASS",
         "SVTYPE=DEL;SVLEN=-140;END=1150", "GT:DR:DV:PL:GQ", "./.:.:18:.,.,.:."],
        ["chr2", 5000, "cuteSV.INS.0", "N", "<INS>", ".", "PASS",
         "SVTYPE=INS;SVLEN=300;END=5000", "GT:DR:DV:PL:GQ", "0/1:20:9:.,.,.:5"],
    ])
    merged = write_vcf(tmp_path / "merged.vcf", [
        ["chr1", 1000, "0_Sniffles2.DEL.1", "N", "<DEL>", ".", "PASS",
         "SVTYPE=DEL;SUPP_VEC=11;SUPP=2;IDLIST=Sniffles2.DEL.1,cuteSV.DEL.0",
         "GT", "0/1"],
        ["chr2", 5000, "1_cuteSV.INS.0", "N", "<INS>", ".", "PASS",
         "SVTYPE=INS;SUPP_VEC=01;SUPP=1;IDLIST=cuteSV.INS.0", "GT", "0/1"],
    ])
    inputs = [("HG002", "sniffles2", sniffles), ("HG003", "cutesv", cutesv)]
    return merged, inputs


def write_test_database(tmp_path):
    merged, inputs = write_inputs(tmp_path)
    return write_database(merged, inputs, tmp_path / "svshare.db", [])


# Checks that each Jasmine record becomes one SV, with how far apart its calls'
# starts and lengths are.
def test_write_database_has_one_sv_per_jasmine_record(tmp_path):
    db = sqlite3.connect(write_test_database(tmp_path))

    assert db.execute("SELECT * FROM svs ORDER BY sv_id").fetchall() == [
        (1, "chr1", 1000, "DEL", 10, 10),
        (2, "chr2", 5000, "INS", 0, 0),
    ]


# Checks that every original call gets its own row, tagged with the sample and
# caller of the input it came from and linked to the SV it was merged into.
def test_write_database_keeps_one_row_per_call(tmp_path):
    db = sqlite3.connect(write_test_database(tmp_path))

    rows = db.execute(
        "SELECT sv_id, sample, caller, vcf_id, pos, svlen, genotype, read_support"
        " FROM calls ORDER BY call_id"
    ).fetchall()
    assert rows == [
        (1, "HG002", "sniffles2", "Sniffles2.DEL.1", 1000, -150, "0/1", 22),
        (1, "HG003", "cutesv", "cuteSV.DEL.0", 1010, -140, "./.", 18),
        (2, "HG003", "cutesv", "cuteSV.INS.0", 5000, 300, "0/1", 9),
    ]
