# `report` subcommand: filter and report analysis results.
import sys

from .callers import CALLERS
from .database import (SVFilter, gene_names, select_sv_ids, sv_calls, sv_rows,
                       sv_transcript_parts)
from .genes import GENE_INFO, read_gene_symbols
from .population import DATASETS
from .svs import SV, Call, GeneDetails, SVPosition, Transcript, TranscriptPart
from .tsv import write_tsv

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
    calls = {}
    for row in sv_calls(db_path, sv_ids):
        calls.setdefault(row["sv_id"], []).append(
            Call(row["sample"], row["caller"], row["genotype"], row["read_support"])
        )
    parts = _load_transcript_parts(db_path, sv_ids)

    svs = []
    for row in sv_rows(db_path, sv_ids):
        position = SVPosition(
            chrom=row["chrom"], pos=row["pos"], vcf_end=row["vcf_end"],
            reference_end=row["reference_end"], svtype=row["svtype"], svlen=row["svlen"],
        )
        svs.append(SV(
            row["sv_id"], position, calls.get(row["sv_id"], []),
            gnomad_af=row["gnomad_af"], colorsdb_af=row["colorsdb_af"],
            near_reference_gap=row["near_reference_gap"],
            transcript_parts=parts.get(row["sv_id"], []),
        ))
    svs.sort(key=lambda sv: (_chrom_order(sv.position.chrom), sv.position.pos, sv.sv_id))
    return svs


# Return the transcript parts near each SV in sv_ids, by SV ID. Parts of the same
# transcript share one Transcript, and transcripts of the same gene share one
# GeneDetails.
def _load_transcript_parts(db_path, sv_ids):
    genes, transcripts, parts = {}, {}, {}
    for row in sv_transcript_parts(db_path, sv_ids):
        gene_id, transcript_id = row["gene_id"], row["transcript_id"]
        if gene_id not in genes:
            genes[gene_id] = GeneDetails(
                gene_id, row["hgnc_id"], row["full_name"], row["summary"],
                row["haploinsufficiency_score"], row["haploinsufficiency"],
                row["triplosensitivity_score"], row["triplosensitivity"],
            )
        if transcript_id not in transcripts:
            transcripts[transcript_id] = Transcript(
                transcript_id, row["gene_name"], gene_id, row["transcript_type"],
                bool(row["mane_select"]), row["chrom"], row["strand"], genes[gene_id],
            )
        parts.setdefault(row["sv_id"], []).append(TranscriptPart(
            transcripts[transcript_id], row["region"], row["number"], row["pos"],
            row["end_pos"], bool(row["overlaps_sv"]),
        ))
    return parts


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
