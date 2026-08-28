# `analyze` subcommand: call, compare, and annotate SVs across samples.
import sys

from .callers import CALLERS
from .reference import check_bam_reference_compatibility, ensure_reference_index

def run(args):
    ensure_reference_index(args.reference)

    all_issues = []
    print(f"Analyzing {len(args.samples)} sample(s):")
    for sample in args.samples:
        print(f"  - {sample}")
        issues = check_bam_reference_compatibility(sample, args.reference)
        for issue in issues:
            print(f"    ! {issue}")
        all_issues.extend(issues)

    if all_issues:
        sys.exit(f"{len(all_issues)} contig mismatch(es) found; aborting.")

    print(f"Reference: {args.reference}")
    print(f"Output directory: {args.output}")

    args.output.mkdir(parents=True, exist_ok=True)
    for sample in args.samples:
        for name, caller in CALLERS.items():
            print(f"  calling SVs with {name}: {sample}")
            caller(sample, args.reference, args.output)
