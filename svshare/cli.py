# Command-line interface: parse arguments and run the chosen subcommand.
import argparse
from pathlib import Path

from . import analyze, filtering, reporting
from .callers import CALLERS


# Build the parser with the analyze, filter, and report subcommands.
def build_parser():
    parser = argparse.ArgumentParser(
        prog="svshare",
        description="Compare, filter, and annotate structural variants across samples.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    analyze_parser = subparsers.add_parser(
        "analyze", help="Call, compare, and annotate structural variants across samples."
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
    # Development only (hidden from --help): reuse caller VCFs from an earlier run,
    # named <sample>.sniffles2.vcf and <sample>.cutesv.vcf, instead of rerunning the callers.
    analyze_parser.add_argument(
        "--vcf-dir", type=Path, default=None, help=argparse.SUPPRESS
    )
    analyze_parser.set_defaults(func=analyze.run)

    filter_parser = subparsers.add_parser("filter", help="Filter analysis results.")
    filter_parser.add_argument(
        "--results", required=True, type=Path, help="Directory containing analyze results."
    )
    filter_parser.add_argument(
        "--caller-support", choices=[*CALLERS, "both", "any"], default="any"
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
    filter_parser.set_defaults(func=filtering.run)

    report_parser = subparsers.add_parser(
        "report", help="Generate a formatted report from results."
    )
    report_parser.add_argument(
        "--results", required=True, type=Path, help="Directory containing analyze or filter results."
    )
    report_parser.add_argument(
        "--output", type=Path, default=None,
        help="Report file to write (default: <results>/report.html).",
    )
    report_parser.set_defaults(func=reporting.run)

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
