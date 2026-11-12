# `report` subcommand: filter and report analysis results.
import sys

from .callers import CALLERS
from .database import (SVFilter, gene_details, gene_names, overlapped_gene_regions,
                       select_sv_ids, sv_calls, sv_rows)
from .genes import GENE_INFO, read_gene_symbols
from .population import DATASETS
from .svs import SV, Call, Gene
from .tsv import write_tsv

# The order a gene's parts are listed in, roughly from upstream to downstream.
REGION_ORDER = ["promoter", "UTR", "start_codon", "CDS", "stop_codon", "noncoding_exon", "intron"]


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
def make_sv_filter(caller_support, min_samples, sv_type, max_frequency, genes, region,
                   near_reference_gap):
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
            row["sv_id"], row["chrom"], row["pos"], row["reference_end"], row["svtype"],
            row["svlen"], row["gnomad_af"], row["colorsdb_af"],
            calls.get(row["sv_id"], []), genes, row["near_reference_gap"],
        ))
    svs.sort(key=lambda sv: (_chrom_order(sv.chrom), sv.pos, sv.sv_id))
    return svs


# Ask a yes/no question, defaulting to no.
def _confirm(question):
    try:
        return input(f"{question} [y/N] ").strip().lower() in ("y", "yes")
    # There's no one to answer, so treat it as no.
    except EOFError:
        print()
        return False


def run(args):
    # SQLite would otherwise create an empty database at a mistyped path.
    if not args.database.is_file():
        sys.exit(f"svshare report: {args.database} isn't a file; pass the database "
                 "analyze wrote, e.g. results/svshare.db")
    genes = parse_gene_list(args.genes, args.gene_file)
    unknown = find_unknown_genes(genes, args.database)
    if unknown:
        print(f"These aren't known gene symbols, so no SVs will match them: {', '.join(unknown)}")
        if not _confirm("Continue anyway?"):
            sys.exit("svshare report: stopped; check the gene symbols and try again")
    output = args.output or (args.database.parent / f"report.{args.format}")
    print(f"Reporting on {args.database}")
    print(f"  caller-support: {args.caller_support}")
    print(f"  sv-type: {args.sv_type}")
    print(f"  max-frequency: {args.max_frequency}")
    print(f"  min-samples: {args.min_samples}")
    print(f"  near-reference-gap: {args.near_reference_gap}")
    if genes:
        print(f"  genes: {', '.join(genes)}")
    if args.region:
        chrom, start, end = args.region
        print(f"  region: {chrom}:{start}-{end}")

    sv_filter = make_sv_filter(caller_support=args.caller_support,
                               min_samples=args.min_samples, sv_type=args.sv_type,
                               max_frequency=args.max_frequency, genes=genes,
                               region=args.region,
                               near_reference_gap=args.near_reference_gap)
    sv_ids = select_sv_ids(args.database, sv_filter)
    print(f"{len(sv_ids)} SVs pass the filters")
    if args.format == "tsv":
        write_tsv(load_svs(args.database, sv_ids), output)
        print(f"Wrote {output}")
    else:
        print("HTML reports aren't written yet; pass --format tsv for a table")
