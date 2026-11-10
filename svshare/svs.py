# The SVs analyze found, loaded from the database for reports.
from dataclasses import dataclass

from .database import gene_details, overlapped_gene_regions, sv_calls, sv_rows

# The order a gene's parts are listed in, roughly from upstream to downstream.
REGION_ORDER = ["promoter", "UTR", "start_codon", "CDS", "stop_codon", "noncoding_exon", "intron"]

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


# Sort key for chromosomes: numbered ones in number order, then X, Y, and M,
# then the rest by name.
def _chrom_order(chrom):
    name = chrom.removeprefix("chr")
    if name.isdigit():
        return (0, int(name), "")
    if name in ("X", "Y", "M"):
        return (1, "XYM".index(name), "")
    return (2, 0, name)


# Return the SVs in sv_ids, in genome order.
def load_svs(db_path, sv_ids):
    calls, regions = {}, {}
    for row in sv_calls(db_path, sv_ids):
        calls.setdefault(row["sv_id"], []).append(
            Call(row["sample"], row["caller"], row["genotype"], row["read_support"])
        )
    for row in overlapped_gene_regions(db_path, sv_ids):
        gene = (row["gene_name"], row["gene_id"])
        regions.setdefault(row["sv_id"], {}).setdefault(gene, set()).add(row["region"])
    gene_info = gene_details(db_path)

    svs = []
    for row in sv_rows(db_path, sv_ids):
        genes = []
        for (name, gene_id), gene_regions in sorted(regions.get(row["sv_id"], {}).items()):
            info = gene_info.get(gene_id, {})
            genes.append(Gene(
                name, info.get("full_name"), info.get("haploinsufficiency"),
                info.get("triplosensitivity"), sorted(gene_regions, key=REGION_ORDER.index),
            ))
        svs.append(SV(
            row["sv_id"], row["chrom"], row["pos"],
            None if row["svtype"] in POINT_TYPES else row["vcf_end"],
            row["svtype"], row["svlen"], row["gnomad_af"], row["colorsdb_af"],
            calls.get(row["sv_id"], []), genes, row["near_reference_gap"],
        ))
    svs.sort(key=lambda sv: (_chrom_order(sv.chrom), sv.pos, sv.sv_id))
    return svs
