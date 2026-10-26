# `report` subcommand: filter analysis results and write them as HTML or TSV.
import sys

from .filtering import filter_sv_ids, parse_gene_list


def run(args):
    # SQLite would otherwise create an empty database at a mistyped path.
    if not args.database.is_file():
        sys.exit(f"svshare report: {args.database} isn't a file; pass the database "
                 "analyze wrote, e.g. results/svshare.db")
    genes = parse_gene_list(args.genes, args.gene_file)
    output = args.output or (args.database.parent / f"report.{args.format}")
    print(f"Reporting on {args.database}")
    print(f"  caller-support: {args.caller_support}")
    print(f"  sv-type: {args.sv_type}")
    print(f"  max-gnomad-frequency: {args.max_gnomad_frequency}")
    print(f"  min-samples: {args.min_samples}")
    if genes:
        print(f"  genes: {', '.join(genes)}")
    if args.region:
        print(f"  region: {args.region}")

    sv_ids = filter_sv_ids(args.database, caller_support=args.caller_support,
                           sv_type=args.sv_type)
    print(f"{len(sv_ids)} SVs pass the filters")
    print(f"Output: {output}")
