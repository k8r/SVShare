from pathlib import Path
from types import SimpleNamespace

import pysam
import pytest

from svshare import analyze
from svshare.reference import check_bam_reference_compatibility, ensure_reference_index

TEST_BAM = Path(__file__).parent.parent / "DevData" / "HG002.test.bam"
TEST_REFERENCE = Path(__file__).parent / "fixtures" / "test_reference.fasta"


def read_fai(fai_path):
    contigs = {}
    with open(fai_path) as f:
        for line in f:
            name, length = line.split("\t")[:2]
            contigs[name] = int(length)
    return contigs


def write_fai(fai_path, contigs):
    with open(fai_path, "w") as f:
        for name, length in contigs.items():
            f.write(f"{name}\t{length}\t0\t60\t61\n")


def test_matching_reference_reports_no_issues():
    assert check_bam_reference_compatibility(TEST_BAM, TEST_REFERENCE) == []


def test_missing_contig_is_reported(tmp_path):
    contigs = read_fai(f"{TEST_REFERENCE}.fai")
    missing_name = "chr1"
    del contigs[missing_name]

    reference_path = tmp_path / "reference.fasta"
    write_fai(f"{reference_path}.fai", contigs)

    issues = check_bam_reference_compatibility(TEST_BAM, reference_path)

    assert len(issues) == 1
    assert missing_name in issues[0]
    assert "not found" in issues[0]


def test_length_mismatch_is_reported(tmp_path):
    contigs = read_fai(f"{TEST_REFERENCE}.fai")
    mismatched_name = "chr1"
    contigs[mismatched_name] += 1

    reference_path = tmp_path / "reference.fasta"
    write_fai(f"{reference_path}.fai", contigs)

    issues = check_bam_reference_compatibility(TEST_BAM, reference_path)

    assert len(issues) == 1
    assert mismatched_name in issues[0]
    assert "length mismatch" in issues[0]


def test_ensure_reference_index_noop_when_fai_exists(monkeypatch):
    calls = []
    monkeypatch.setattr(pysam, "faidx", lambda path: calls.append(path))

    ensure_reference_index(TEST_REFERENCE)

    assert calls == []


def test_ensure_reference_index_builds_when_missing_and_writable(tmp_path, monkeypatch):
    reference_path = tmp_path / "reference.fasta"
    reference_path.write_text("placeholder")

    calls = []
    monkeypatch.setattr(pysam, "faidx", lambda path: calls.append(path))

    ensure_reference_index(reference_path)

    assert calls == [str(reference_path)]


def test_ensure_reference_index_raises_when_directory_not_writable(tmp_path, monkeypatch):
    reference_path = tmp_path / "reference.fasta"
    reference_path.write_text("placeholder")

    monkeypatch.setattr("svshare.reference.os.access", lambda path, mode: False)
    calls = []
    monkeypatch.setattr(pysam, "faidx", lambda path: calls.append(path))

    with pytest.raises(PermissionError):
        ensure_reference_index(reference_path)

    assert calls == []


def test_run_aborts_when_reference_incompatible(tmp_path):
    contigs = read_fai(f"{TEST_REFERENCE}.fai")
    del contigs["chr1"]

    reference_path = tmp_path / "reference.fasta"
    write_fai(f"{reference_path}.fai", contigs)

    args = SimpleNamespace(
        samples=[TEST_BAM], reference=reference_path, output=tmp_path / "results", vcf_dir=None
    )

    with pytest.raises(SystemExit):
        analyze.run(args)
