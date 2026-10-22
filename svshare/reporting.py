# `report` subcommand: filter analysis results and write them as HTML or TSV.
from .analyze import DATABASE
from .filtering import filter_sv_ids, parse_gene_list


def run(args):
    genes = parse_gene_list(args.genes, args.gene_file)
    output = args.output or (args.results / f"report.{args.format}")
    print(f"Reporting results in {args.results}")
    print(f"  caller-support: {args.caller_support}")
    print(f"  sv-type: {args.sv_type}")
    print(f"  max-gnomad-frequency: {args.max_gnomad_frequency}")
    print(f"  min-samples: {args.min_samples}")
    if genes:
        print(f"  genes: {', '.join(genes)}")
    if args.region:
        print(f"  region: {args.region}")

    sv_ids = filter_sv_ids(args.results / DATABASE, caller_support=args.caller_support)
    print(f"{len(sv_ids)} SVs pass the filters")
    print(f"Output: {output}")
