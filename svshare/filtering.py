# Functions for choosing which SVs a report includes.
from dataclasses import dataclass

from .callers import CALLERS
from .database import connect, gene_names
from .genes import GENE_INFO, read_gene_symbols
from .population import DATASETS


# Combine --genes values and --gene-file lines into a sorted, deduplicated list.
def parse_gene_list(genes, gene_file):
    gene_set = set(genes or [])
    if gene_file:
        with open(gene_file) as f:
            gene_set.update(line.strip() for line in f if line.strip())
    return sorted(gene_set)


# Return the genes that aren't known symbols, ignoring case. Known symbols are
# HGNC's plus the database's gene names, which include genes HGNC doesn't name.
def find_unknown_genes(genes, db_path, gene_info=GENE_INFO):
    known = {symbol.upper() for symbol in read_gene_symbols(gene_info)}
    known.update(name.upper() for name in gene_names(db_path))
    return [gene for gene in genes if gene.upper() not in known]


# Which SVs a report includes. Each field left at its default keeps every SV.
@dataclass
class SVFilter:
    # A sample counts toward min_samples if every caller in callers called the
    # SV in it, or any caller if callers is empty.
    min_samples: int = 1
    callers: tuple[str, ...] = ()
    svtype: str | None = None
    # SVs not found in a population dataset still pass.
    population_datasets: tuple[str, ...] = ()
    max_allele_frequency: float | None = None
    genes: tuple[str, ...] = ()
    region: tuple[str, int, int] | None = None
    # True keeps only SVs at or near a reference gap, False only the rest.
    near_reference_gap: bool | None = None


# Turn the report's filter options into an SVFilter.
def make_sv_filter(caller_support="any", min_samples=1, sv_type="any", max_frequency=None,
                   genes=None, region=None, near_reference_gap="include"):
    if caller_support == "any":
        callers = ()
    elif caller_support == "both":
        callers = tuple(CALLERS)
    else:
        callers = (caller_support,)
    return SVFilter(
        callers=callers,
        min_samples=min_samples,
        svtype=None if sv_type == "any" else sv_type,
        population_datasets=tuple(DATASETS),
        max_allele_frequency=max_frequency,
        genes=tuple(genes or ()),
        region=region,
        near_reference_gap={"include": None, "only": True, "exclude": False}[near_reference_gap],
    )


def _sample_support_condition(callers, min_samples):
    samples = "SELECT sv_id, sample FROM calls"
    params = []
    if callers:
        samples += (f" WHERE caller IN ({', '.join('?' * len(callers))})"
                    " GROUP BY sv_id, sample HAVING COUNT(DISTINCT caller) = ?")
        params = list(callers) + [len(callers)]
    else:
        samples += " GROUP BY sv_id, sample"
    return (
        f"sv_id IN (SELECT sv_id FROM ({samples}) GROUP BY sv_id HAVING COUNT(*) >= ?)",
        params + [min_samples],
    )


def _svtype_condition(svtype):
    if svtype is None:
        return None
    return "svtype = ?", [svtype]


def _allele_frequency_condition(population_datasets, max_allele_frequency):
    if max_allele_frequency is None:
        return None
    return (
        " AND ".join(f"({dataset}_af IS NULL OR {dataset}_af <= ?)"
                     for dataset in population_datasets),
        [max_allele_frequency] * len(population_datasets),
    )


def _genes_condition(genes):
    if not genes:
        return None
    return (
        "sv_id IN (SELECT sv_id FROM sv_transcript_parts JOIN transcripts USING (transcript_id)"
        f" WHERE overlaps_sv = 1 AND gene_name COLLATE NOCASE IN ({', '.join('?' * len(genes))}))",
        list(genes),
    )


def _region_condition(region):
    if region is None:
        return None
    chrom, start, end = region
    # An insertion's or breakend's end position isn't where it ends, so only its
    # start position is checked.
    return (
        "chrom = ? AND pos <= ?"
        " AND (CASE WHEN svtype IN ('INS', 'BND') THEN pos ELSE MAX(end_pos, pos) END) >= ?",
        [chrom, end, start],
    )


def _reference_gap_condition(near_reference_gap):
    if near_reference_gap is None:
        return None
    return f"near_reference_gap IS {'NOT ' if near_reference_gap else ''}NULL", []


# Return the IDs, in order, of the SVs in the database that pass sv_filter.
def select_sv_ids(db_path, sv_filter):
    conditions, params = [], []
    for condition in (
        _sample_support_condition(sv_filter.callers, sv_filter.min_samples),
        _svtype_condition(sv_filter.svtype),
        _allele_frequency_condition(sv_filter.population_datasets,
                                    sv_filter.max_allele_frequency),
        _genes_condition(sv_filter.genes),
        _region_condition(sv_filter.region),
        _reference_gap_condition(sv_filter.near_reference_gap),
    ):
        if condition:
            conditions.append(condition[0])
            params += condition[1]
    where = f" WHERE {' AND '.join(conditions)}" if conditions else ""
    with connect(db_path) as db:
        return [sv_id for (sv_id,) in db.execute(f"SELECT sv_id FROM svs{where} ORDER BY sv_id", params)]


# Return the IDs, in order, of the SVs in the database that pass the report's
# filter options.
def filter_sv_ids(db_path, **options):
    return select_sv_ids(db_path, make_sv_filter(**options))
