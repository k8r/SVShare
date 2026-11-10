import sqlite3
from pathlib import Path

import pytest

from helpers import make_database, read_rows
from svshare.commands import ToolRun
from svshare.database import connect, save_tool_runs, write_database
from svshare.merging import VCFPreppedForMerging

HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"


def write_vcf(path, records):
    path.write_text(HEADER + "".join("\t".join(map(str, r)) + "\n" for r in records))
    return path


# A Sniffles2 VCF, a cuteSV VCF, and a Jasmine merge of them: a deletion called
# in both HG002 (sniffles2) and HG003 (cuteSV), 10 bp apart, and an insertion
# only in HG003.
def write_merge_vcfs(tmp_path):
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
         "SVTYPE=DEL;SVLEN=-150;END=1150;SUPP_VEC=11;SUPP=2;IDLIST=Sniffles2.DEL.1,cuteSV.DEL.0",
         "GT", "0/1"],
        ["chr2", 5000, "1_cuteSV.INS.0", "N", "<INS>", ".", "PASS",
         "SVTYPE=INS;SVLEN=300;END=5000;SUPP_VEC=01;SUPP=1;IDLIST=cuteSV.INS.0",
         "GT", "0/1"],
    ])
    prepped_vcfs = [VCFPreppedForMerging("HG002", "sniffles2", sniffles),
                    VCFPreppedForMerging("HG003", "cutesv", cutesv)]
    return merged, prepped_vcfs


# Checks that each Jasmine record becomes one SV, with the kept call's end and
# length and how far apart its calls' starts and lengths are.
def test_write_database_has_one_sv_per_jasmine_record(tmp_path):
    merged_vcf, prepped_vcfs = write_merge_vcfs(tmp_path)
    db_path = write_database(merged_vcf, prepped_vcfs, tmp_path / "svshare.db")

    columns = ("sv_id", "chrom", "pos", "end_pos", "svtype", "svlen",
               "start_spread", "length_spread")
    assert read_rows(db_path, "svs", *columns) == [
        (1, "chr1", 1000, 1150, "DEL", -150, 10, 10),
        (2, "chr2", 5000, 5000, "INS", 300, 0, 0),
    ]


# Checks that every original call gets its own row, tagged with the sample and
# caller of the input it came from and linked to the SV it was merged into.
def test_write_database_keeps_one_row_per_call(tmp_path):
    merged_vcf, prepped_vcfs = write_merge_vcfs(tmp_path)
    db_path = write_database(merged_vcf, prepped_vcfs, tmp_path / "svshare.db")

    columns = ("sv_id", "sample", "caller", "vcf_id", "pos", "svlen", "genotype", "read_support")
    assert read_rows(db_path, "calls", *columns) == [
        (1, "HG002", "sniffles2", "Sniffles2.DEL.1", 1000, -150, "0/1", 22),
        (1, "HG003", "cutesv", "cuteSV.DEL.0", 1010, -140, "./.", 18),
        (2, "HG003", "cutesv", "cuteSV.INS.0", 5000, 300, "0/1", 9),
    ]


def test_save_tool_runs(tmp_path):
    db_path = make_database(tmp_path)
    save_tool_runs(db_path, [
        ToolRun("HG002", "mosdepth", ["mosdepth", "out prefix"]),
        ToolRun(None, "bedtools", ["bedtools", "intersect"], Path("overlaps.bed")),
    ])
    assert read_rows(db_path, "commands", "sample", "tool", "command") == [
        ("HG002", "mosdepth", "mosdepth 'out prefix'"),
        (None, "bedtools", "bedtools intersect > overlaps.bed"),
    ]


# Changes are saved, and the database is closed even when something fails.
def test_connect_saves_and_closes(tmp_path):
    db_path = tmp_path / "test.db"
    with connect(db_path) as db:
        db.execute("CREATE TABLE t (x INTEGER)")
        db.execute("INSERT INTO t VALUES (1)")
    with pytest.raises(sqlite3.ProgrammingError):
        db.execute("SELECT 1")

    with pytest.raises(ValueError):
        with connect(db_path) as db:
            db.execute("INSERT INTO t VALUES (2)")
            raise ValueError
    with pytest.raises(sqlite3.ProgrammingError):
        db.execute("SELECT 1")
    with connect(db_path) as db:
        assert db.execute("SELECT x FROM t").fetchall() == [(1,)]
