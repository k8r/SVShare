# Helpers for working with the reference genome and its compatibility with sample data.
import os
from pathlib import Path

import pysam

# Ensure the reference FASTA has a .fai index, building it once if missing.
def ensure_reference_index(reference_path):
    fai_path = Path(f"{reference_path}.fai")
    if fai_path.exists():
        return
    if not os.access(reference_path.parent, os.W_OK):
        raise PermissionError(
            f"{reference_path.parent} is read-only, so {fai_path} cannot be created here. "
            f"Build the index first with `samtools faidx {reference_path}` somewhere writable."
        )
    print(f"Indexing reference {reference_path} ...")
    pysam.faidx(str(reference_path))

# Compare a BAM's header contigs (name + length) against the reference .fai,
# returning a list of mismatch descriptions (empty means compatible).
def check_bam_reference_compatibility(bam_path, reference_path):
    fai_path = Path(f"{reference_path}.fai")
    if not fai_path.exists():
        raise FileNotFoundError(
            f"Reference index {fai_path} not found; run ensure_reference_index() first."
        )

    ref_contigs = {}
    with open(fai_path) as f:
        for line in f:
            name, length = line.split("\t")[:2]
            ref_contigs[name] = int(length)

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        bam_contigs = dict(zip(bam.references, bam.lengths))

    issues = []
    for name, length in bam_contigs.items():
        if name not in ref_contigs:
            issues.append(f"Contig '{name}' in {bam_path} not found in reference {reference_path}.")
        elif ref_contigs[name] != length:
            issues.append(
                f"Contig '{name}' length mismatch: BAM={length}, reference={ref_contigs[name]} ({bam_path})."
            )

    return issues
