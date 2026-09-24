# `filter` subcommand: narrow analysis results by caller, type, frequency, genes, or region.


# Combine --genes values and --gene-file lines into a sorted, deduplicated list.
def parse_gene_list(genes, gene_file):
    gene_set = set(genes or [])
    if gene_file:
        with open(gene_file) as f:
            gene_set.update(line.strip() for line in f if line.strip())
    return sorted(gene_set)


def run(args):
    genes = parse_gene_list(args.genes, args.gene_file)
    print(f"Filtering results in {args.results}")
    print(f"  caller-support: {args.caller_support}")
    print(f"  sv-type: {args.sv_type}")
    print(f"  max-gnomad-frequency: {args.max_gnomad_frequency}")
    print(f"  min-samples: {args.min_samples}")
    if genes:
        print(f"  genes: {', '.join(genes)}")
    if args.region:
        print(f"  region: {args.region}")
