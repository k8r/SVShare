# What an SV is: where it is, the calls of it, and the transcripts and genes near it.
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


# What the gene info file says about a gene near an SV.
@dataclass
class GeneDetails:
    gene_id: str
    hgnc_id: str | None = None
    full_name: str | None = None
    summary: str | None = None
    # The ClinGen dosage sensitivity scores below are 0-3 for how strong the
    # evidence is, or 30 (autosomal recessive) or 40 (dosage sensitivity
    # unlikely), each with ClinGen's description.
    haploinsufficiency_score: int | None = None
    haploinsufficiency: str | None = None
    triplosensitivity_score: int | None = None
    triplosensitivity: str | None = None


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
    gene_details: GeneDetails | None = None


# A part of a transcript near an SV.
@dataclass
class TranscriptPart:
    transcript: Transcript
    region: str
    number: int | None
    pos: int
    end: int


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
    # The transcript parts this SV overlaps.
    overlapped_parts: list[TranscriptPart] = field(default_factory=list)
    # The transcript parts near this SV that it doesn't overlap, for showing it in
    # context.
    nearby_parts: list[TranscriptPart] = field(default_factory=list)
    gnomad_id: str | None = None
    gnomad_af: float | None = None
    colorsdb_id: str | None = None
    colorsdb_af: float | None = None
    near_reference_gap: str | None = None

    # Return every transcript part near this SV, overlapped or not.
    def all_parts(self):
        return self.overlapped_parts + self.nearby_parts

    # Return whether this SV overlaps a part of the given gene, by its name or its
    # Ensembl ID, ignoring case and any ID version.
    def overlaps_gene(self, gene):
        gene_id = gene.split(".")[0].upper()
        for part in self.overlapped_parts:
            transcript = part.transcript
            if transcript.gene_name.upper() == gene.upper():
                return True
            if transcript.gene_id.split(".")[0].upper() == gene_id:
                return True
        return False


# Where an SV ends on the reference: at its position if it's an insertion or
# breakend or its VCF has no END, otherwise at its VCF END.
def reference_end(svtype, pos, vcf_end):
    if svtype in POINT_TYPES or vcf_end is None:
        return pos
    return max(vcf_end, pos)
