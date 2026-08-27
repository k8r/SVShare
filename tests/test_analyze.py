# Tests for the analyze stage: currently just BAM/reference compatibility
# checking (svshare/reference.py's check_bam_reference_compatibility).
from pathlib import Path

from svshare.reference import check_bam_reference_compatibility

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
