# What an SV is: where it is, the calls of it, and the genes it overlaps.
from dataclasses import dataclass

# SV types that sit at one position, whatever end their caller gave.
POINT_TYPES = ("INS", "BND")


@dataclass
class Call:
    sample: str
    caller: str
    genotype: str | None
    read_support: int | None = None


# A gene an SV overlaps, with the parts of it the SV overlaps in any transcript.
@dataclass
class Gene:
    name: str
    full_name: str | None
    haploinsufficiency: str | None
    triplosensitivity: str | None
    regions: list[str]


@dataclass
class SV:
    sv_id: int
    chrom: str
    pos: int
    end: int | None
    svtype: str
    svlen: int | None
    gnomad_af: float | None
    colorsdb_af: float | None
    calls: list[Call]
    genes: list[Gene]
    near_reference_gap: str | None = None


# Where an SV ends: at its position if it's an insertion or breakend or has no
# end, otherwise at its end.
def sv_end(svtype, pos, vcf_end):
    if svtype in POINT_TYPES or vcf_end is None:
        return pos
    return max(vcf_end, pos)
