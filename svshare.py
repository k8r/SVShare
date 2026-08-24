#!/usr/bin/env python3
import argparse
from pathlib import Path

def parse_gene_list(genes, gene_file):
    gene_set = set(genes or [])
    if gene_file:
        with open(gene_file) as f:
            gene_set.update(line.strip() for line in f if line.strip())
    return sorted(gene_set)

def cmd_analyze(args):
    print(f"Analyzing {len(args.samples)} sample(s):")
    for sample in args.samples:
        print(f"  - {sample}")
    print(f"Reference: {args.reference}")
    print(f"Output directory: {args.output}")

def cmd_filter(args):
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
    
def cmd_report(args):
    print(f"Generating report from {args.results}")
    print(f"Output: {args.output}")
    
def build_parser():
    parser = argparse.ArgumentParser(
        prog="svshare",
        description="Compare, filter, and annotate structural variants across samples.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    analyze_parser = subparsers.add_parser(
        "analyze", help="Call and merge structural variants across samples."
    )
    analyze_parser.add_argument(
        "--samples", nargs="+", required=True, type=Path,
        help="BAM or CRAM files, one per sample.",
    )
    analyze_parser.add_argument(
        "--reference", required=True, type=Path, help="Reference genome FASTA."
    )
    analyze_parser.add_argument(
        "--output", required=True, type=Path, help="Directory to write analysis results to."
    )
    analyze_parser.set_defaults(func=cmd_analyze)

    filter_parser = subparsers.add_parser("filter", help="Filter analysis results.")
    filter_parser.add_argument(
        "--results", required=True, type=Path, help="Directory containing analyze results."
    )
    filter_parser.add_argument(
        "--caller-support", choices=["both", "sniffles2", "cutesv", "any"], default="any"
    )
    filter_parser.add_argument(
        "--sv-type", choices=["DEL", "DUP", "INS", "INV", "BND", "any"], default="any"
    )
    filter_parser.add_argument("--max-gnomad-frequency", type=float, default=0.01)
    filter_parser.add_argument("--min-samples", type=int, default=2)
    filter_parser.add_argument(
        "--genes", nargs="+", default=None, help="Gene symbols to filter to."
    )
    filter_parser.add_argument(
        "--gene-file", type=Path, default=None,
        help="Text file of gene symbols, one per line, no header.",
    )
    filter_parser.add_argument(
        "--region", default=None, help="Genomic region to filter to, e.g. chr1:1000-2000."
    )
    filter_parser.set_defaults(func=cmd_filter)

    report_parser = subparsers.add_parser(
        "report", help="Generate a formatted report from results."
    )
    report_parser.add_argument(
        "--results", required=True, type=Path, help="Directory containing analyze or filter results."
    )
    report_parser.add_argument(
        "--output", required=True, type=Path, help="Report file to write, e.g. report.html."
    )
    report_parser.set_defaults(func=cmd_report)

    return parser

def main():
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
