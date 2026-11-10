# Functions for choosing which SVs a report includes.
from .callers import CALLERS
from .database import SVFilter, gene_names, select_sv_ids
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


# Return the IDs, in order, of the SVs in the database that pass the report's
# filter options.
def filter_sv_ids(db_path, **options):
    return select_sv_ids(db_path, make_sv_filter(**options))
