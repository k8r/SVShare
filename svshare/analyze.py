# `analyze` subcommand: call, compare, and annotate SVs across samples.
import shlex
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

    args.output.mkdir(parents=True, exist_ok=True)
    commands = {sample: [] for sample in args.samples}
    for i, sample in enumerate(args.samples, start=1):
        print(f"\nSample {i}/{len(args.samples)}: {sample}")
        for name, caller in CALLERS.items():
            print(f"  calling SVs with {name}")
            _vcf, cmd = caller(sample, args.reference, args.output)
            commands[sample].append(cmd)

    print_summary(args, commands)


# Print the samples, reference, and output location, then the numbered caller
# commands run for each sample.
def print_summary(args, commands):
    print()
    print(f"Analyzed {len(args.samples)} sample(s):")
    for sample in args.samples:
        print(f"  - {sample}")
    print(f"Reference: {args.reference}")
    print(f"Output directory: {args.output}")

    print()
    print("Commands run:")
    for sample, cmds in commands.items():
        print()
        print(f"Sample: {sample}")
        for n, cmd in enumerate(cmds, start=1):
            print(f"  {n}) {shlex.join(cmd)}")
