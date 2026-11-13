# Find SVs at or near gaps in the reference, where unknown sequence is written
# as N's. Reads can't align to a gap, so SVs found there are often artifacts.
import bisect
import re

import pysam

# The fewest N's in a row that count as a gap; shorter runs are usually just
# uncertain bases.
MIN_GAP_LENGTH = 10
# How far from a gap (in bp) an SV can be and still count as near it.
MAX_DISTANCE = 1000


# Return {chrom: [(start, end), ...]} of the reference's gaps, 1-based and
# inclusive, in order.
def find_gaps(reference, min_length=MIN_GAP_LENGTH):
    pattern = re.compile(f"[Nn]{{{min_length},}}")
    gaps = {}
    with pysam.FastaFile(str(reference)) as fasta:
        for chrom in fasta.references:
            found = [(match.start() + 1, match.end())
                     for match in pattern.finditer(fasta.fetch(chrom))]
            if found:
                gaps[chrom] = found
    return gaps


# Describe how close the SV spanning start to end is to the nearest gap in
# chrom_gaps, or return None if none is within max_distance.
def describe_nearest_gap(chrom_gaps, start, end, max_distance=MAX_DISTANCE):
    # Gaps don't overlap, so their ends are in order too.
    ends = [gap_end for _gap_start, gap_end in chrom_gaps]
    i = bisect.bisect_left(ends, start - max_distance)
    overlapped, nearest = 0, None
    for gap_start, gap_end in chrom_gaps[i:]:
        if gap_start > end + max_distance:
            break
        if gap_start <= end and gap_end >= start:
            overlapped += gap_end - gap_start + 1
            continue
        distance = gap_start - end if gap_start > end else start - gap_end
        if nearest is None or distance < nearest[0]:
            nearest = (distance, gap_end - gap_start + 1)
    if overlapped:
        return f"overlaps {overlapped:,} unknown bases"
    if nearest:
        return f"{nearest[0]:,} bp from {nearest[1]:,} unknown bases"
    return None


# Set near_reference_gap on each SV at or near a gap in the reference. Returns
# how many are.
def add_reference_gaps(svs, reference):
    gaps = find_gaps(reference)
    near = 0
    for sv in svs:
        position = sv.position
        sv.near_reference_gap = describe_nearest_gap(
            gaps.get(position.chrom, []), position.pos, position.reference_end
        )
        if sv.near_reference_gap:
            near += 1
    return near
