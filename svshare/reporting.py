# `report` subcommand: filter analysis results and write them as HTML or TSV.
import sys

from .filtering import filter_sv_ids, find_unknown_genes, parse_gene_list


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
    if genes:
        print(f"  genes: {', '.join(genes)}")
    if args.region:
        chrom, start, end = args.region
        print(f"  region: {chrom}:{start}-{end}")

    sv_ids = filter_sv_ids(args.database, caller_support=args.caller_support,
                           min_samples=args.min_samples, sv_type=args.sv_type,
                           max_frequency=args.max_frequency, genes=genes,
                           region=args.region)
    print(f"{len(sv_ids)} SVs pass the filters")
    print(f"Output: {output}")
