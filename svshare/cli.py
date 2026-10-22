# Command-line interface: parse arguments and run the chosen subcommand.
import argparse
from pathlib import Path

import svshare

from . import analyze, reporting
from .callers import CALLERS


# Build the parser with the analyze and report subcommands.
def build_parser():
    parser = argparse.ArgumentParser(
        prog="svshare",
        description=svshare.__doc__,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    analyze_parser = subparsers.add_parser(
        "analyze",
        help="Call, compare, and annotate structural variants across samples and callers.",
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
    analyze_parser.add_argument(
        "--gnomad", type=Path, default=None,
        help="gnomAD SV sites VCF (.vcf.gz with its .tbi index) to look up how common "
        "each SV is in the population. Without it, the lookup is skipped.",
    )
    analyze_parser.add_argument(
        "--colorsdb", type=Path, default=None,
        help="CoLoRSdb SV VCF (.vcf.gz with its .tbi index), a long-read population "
        "dataset, to look up how common each SV is. Without it, the lookup is skipped.",
    )
    analyze_parser.add_argument(
        "--gencode", type=Path, default=None,
        help="GENCODE gene annotation GTF (.gtf or .gtf.gz) to find the genes each SV "
        "affects. Without it, SVs aren't matched to genes.",
    )
    # Development only (hidden from --help): reuse caller VCFs from an earlier run,
    # named <sample>.sniffles2.vcf and <sample>.cutesv.vcf, instead of rerunning the callers.
    analyze_parser.add_argument(
        "--vcf-dir", type=Path, default=None, help=argparse.SUPPRESS
    )
    # Development only (hidden from --help): cache the transcript regions BED that
    # --gencode builds at this path, reusing it if it exists. Delete it when the GTF
    # or the region rules in transcripts.py change.
    analyze_parser.add_argument(
        "--regions-bed", type=Path, default=None, help=argparse.SUPPRESS
    )
    analyze_parser.set_defaults(func=analyze.run)

    report_parser = subparsers.add_parser(
        "report", help="Filter analysis results and write them as a report."
    )
    report_parser.add_argument(
        "--results", required=True, type=Path, help="Directory containing analyze results."
    )
    report_parser.add_argument(
        "--format", choices=["html", "tsv"], default="html",
        help="html for a report, tsv for a table (default: html).",
    )
    report_parser.add_argument(
        "--output", type=Path, default=None,
        help="File to write (default: <results>/report.html or <results>/report.tsv).",
    )
    report_parser.add_argument(
        "--caller-support", choices=[*CALLERS, "both", "any"], default="any"
    )
    report_parser.add_argument(
        "--sv-type", choices=["DEL", "DUP", "INS", "INV", "BND", "any"], default="any"
    )
    report_parser.add_argument("--max-gnomad-frequency", type=float, default=0.01)
    report_parser.add_argument("--min-samples", type=int, default=2)
    report_parser.add_argument(
        "--genes", nargs="+", default=None, help="Gene symbols to filter to."
    )
    report_parser.add_argument(
        "--gene-file", type=Path, default=None,
        help="Text file of gene symbols, one per line, no header.",
    )
    report_parser.add_argument(
        "--region", default=None, help="Genomic region to filter to, e.g. chr1:1000-2000."
    )
    report_parser.set_defaults(func=reporting.run)

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
