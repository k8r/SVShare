# What an SV is: where it is, the calls of it, and the genes and transcripts near it.
from dataclasses import dataclass, field

# SV types that sit at one position, whatever end their caller gave.
POINT_TYPES = ("INS", "BND")


# One caller's call of an SV in one sample, as its VCF recorded it.
@dataclass
class Call:
    sample: str
    caller: str
    genotype: str | None
    read_support: int | None = None
    vcf_id: str | None = None
    chrom: str | None = None
    pos: int | None = None
    vcf_end: int | None = None
    svtype: str | None = None
    svlen: int | None = None


# A gene an SV overlaps, with the parts of it the SV overlaps in any transcript.
@dataclass
class Gene:
    name: str
    full_name: str | None
    haploinsufficiency: str | None
    triplosensitivity: str | None
    regions: list[str]


# A transcript near an SV.
@dataclass
class Transcript:
    transcript_id: str
    gene_name: str
    gene_id: str
    transcript_type: str
    mane_select: bool
    chrom: str
    strand: str


# A part of a transcript near an SV, and whether the SV overlaps it.
@dataclass
class TranscriptPart:
    transcript: Transcript
    region: str
    number: int | None
    pos: int
    end: int
    overlaps_sv: bool


# Where an SV is, and its type and length.
@dataclass
class SVPosition:
    chrom: str
    pos: int
    vcf_end: int | None
    reference_end: int | None
    svtype: str
    svlen: int | None


@dataclass
class SV:
    sv_id: int
    position: SVPosition
    calls: list[Call]
    genes: list[Gene] = field(default_factory=list)
    transcript_parts: list[TranscriptPart] = field(default_factory=list)
    gnomad_id: str | None = None
    gnomad_af: float | None = None
    colorsdb_id: str | None = None
    colorsdb_af: float | None = None
    near_reference_gap: str | None = None


# Where an SV ends on the reference: at its position if it's an insertion or
# breakend or its VCF has no END, otherwise at its VCF END.
def reference_end(svtype, pos, vcf_end):
    if svtype in POINT_TYPES or vcf_end is None:
        return pos
    return max(vcf_end, pos)
