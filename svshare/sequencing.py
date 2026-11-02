# Facts about how a sample was sequenced: its platform and its read depth.
import re
from pathlib import Path

import pysam

from .commands import run_command

PLATFORMS = ["hifi", "ont", "clr"]

# Aligner presets that name the platform, as minimap2 (-x) and pbmm2 (--preset)
# record them in the BAM header.
PRESET_PLATFORMS = {
    "map-hifi": "hifi",
    "map-ont": "ont",
    "lr:hq": "ont",
    "lr:hqae": "ont",
    "map-pb": "clr",
    "hifi": "hifi",
    "ccs": "hifi",
    "subread": "clr",
}


# Return the platform the header's aligner commands name, or None if they name
# none or more than one.
def _platform_from_header(header):
    platforms = set()
    for program in header.get("PG", []):
        for preset in re.findall(r"(?:-\w*x|--preset)[ =](\S+)", program.get("CL", "")):
            if preset.lower() in PRESET_PLATFORMS:
                platforms.add(PRESET_PLATFORMS[preset.lower()])
    return platforms.pop() if len(platforms) == 1 else None


def detect_platform(bam):
    with pysam.AlignmentFile(str(bam)) as f:
        return _platform_from_header(f.header.to_dict())


# Return the mean depth across the autosomes from a mosdepth summary, or across
# everything if there are no autosomes. The sex chromosomes and other contigs
# are left out because their depth differs from the rest of the genome's.
def _autosome_depth(summary):
    length = bases = 0
    with open(summary) as f:
        next(f)
        for line in f:
            chrom, chrom_length, chrom_bases, chrom_mean = line.split("\t")[:4]
            if chrom == "total":
                total_mean = float(chrom_mean)
            elif re.fullmatch(r"(chr)?([1-9]|1\d|2[0-2])", chrom):
                length += int(chrom_length)
                bases += int(chrom_bases)
    return bases / length if length else total_mean


# Measure a BAM's mean depth with mosdepth, writing its summary to out_dir.
def mean_depth(bam, out_dir, threads=4):
    prefix = Path(out_dir) / Path(bam).stem
    cmd = ["mosdepth", "--no-per-base", "--fast-mode", "--threads", str(threads),
           str(prefix), str(bam)]
    run_command(cmd)
    return _autosome_depth(f"{prefix}.mosdepth.summary.txt"), cmd
