# Find SVs at or near gaps in the reference, where unknown sequence is written
# as N's. Reads can't align to a gap, so SVs found there are often artifacts.
import bisect
import re

import pysam

from .database import connect, sv_positions
from .svs import sv_end

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


# Save, for each SV in the database at or near a gap in the reference, how
# close it is. Returns (SVs near a gap, total SVs).
def add_reference_gaps(db_path, reference):
    gaps = find_gaps(reference)
    svs = sv_positions(db_path)
    near = []
    for sv in svs:
        description = describe_nearest_gap(gaps.get(sv.chrom, []), sv.pos,
                                           sv_end(sv.svtype, sv.pos, sv.end_pos))
        if description:
            near.append((description, sv.sv_id))
    with connect(db_path) as db:
        db.executemany("UPDATE svs SET near_reference_gap = ? WHERE sv_id = ?", near)
    return len(near), len(svs)
