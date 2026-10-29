# The SVs analyze found, loaded from the database for reports.
import sqlite3
from dataclasses import dataclass

# The order a gene's parts are listed in, roughly from upstream to downstream.
REGION_ORDER = ["promoter", "UTR", "start_codon", "CDS", "stop_codon", "noncoding_exon", "intron"]


@dataclass
class Call:
    sample: str
    caller: str
    genotype: str | None


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
    wanted = set(sv_ids)
    calls, regions = {}, {}
    with sqlite3.connect(db_path) as db:
        for sv_id, sample, caller, genotype in db.execute(
            "SELECT sv_id, sample, caller, genotype FROM calls ORDER BY sample, caller"
        ):
            if sv_id in wanted:
                calls.setdefault(sv_id, []).append(Call(sample, caller, genotype))
        for sv_id, gene_name, gene_id, region in db.execute(
            "SELECT DISTINCT sv_id, gene_name, gene_id, region FROM sv_transcript_parts"
            " JOIN transcripts USING (transcript_id) WHERE overlaps_sv = 1"
        ):
            if sv_id in wanted:
                regions.setdefault(sv_id, {}).setdefault((gene_name, gene_id), set()).add(region)
        gene_info = {gene_id: info for gene_id, *info in db.execute(
            "SELECT gene_id, full_name, haploinsufficiency, triplosensitivity FROM genes"
        )}
        rows = [row for row in db.execute(
            "SELECT sv_id, chrom, pos, end_pos, svtype, svlen, gnomad_af, colorsdb_af FROM svs"
        ) if row[0] in wanted]
    db.close()

    svs = []
    for sv_id, chrom, pos, end_pos, svtype, svlen, gnomad_af, colorsdb_af in rows:
        genes = [
            Gene(name, *gene_info.get(gene_id, (None, None, None)),
                 sorted(gene_regions, key=REGION_ORDER.index))
            for (name, gene_id), gene_regions in sorted(regions.get(sv_id, {}).items())
        ]
        svs.append(SV(
            sv_id, chrom, pos,
            # An insertion's or breakend's end position isn't where it ends.
            None if svtype in ("INS", "BND") else end_pos,
            svtype, svlen, gnomad_af, colorsdb_af, calls.get(sv_id, []), genes,
        ))
    svs.sort(key=lambda sv: (_chrom_order(sv.chrom), sv.pos, sv.sv_id))
    return svs
