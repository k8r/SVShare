# `analyze` subcommand: call, compare, and annotate SVs across samples.
import shlex
import sys
from pathlib import Path

from .callers import CALLERS, caller_vcf_path
from .merging import run_jasmine
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
        caller_vcfs = []
        for name, caller in CALLERS.items():
            if args.vcf_dir:
                vcf = caller_vcf_path(sample, args.vcf_dir, name)
                print(f"  using existing {name} VCF: {vcf}")
            else:
                print(f"  calling SVs with {name}")
                vcf, cmd = caller(sample, args.reference, args.output)
                commands[sample].append(cmd)
            caller_vcfs.append(vcf)

        print("  merging caller results with Jasmine")
        merged_vcf = args.output / f"{Path(sample).stem}.merged.vcf"
        _vcf, cmd = run_jasmine(caller_vcfs, merged_vcf, args.output)
        commands[sample].append(cmd)

    print_summary(args, commands)


# Print the samples, reference, and output location, then the numbered
# commands (callers and Jasmine) run for each sample.
def print_summary(args, commands):
    print()
    print(f"Analyzed {len(args.samples)} sample(s):")
    for sample in args.samples:
        print(f"  - {sample}")
    print(f"Reference: {args.reference}")
    print(f"Output directory: {args.output}")
    if args.vcf_dir:
        print(f"Existing caller VCFs from: {args.vcf_dir}")

    print()
    print("Commands run:")
    for sample, cmds in commands.items():
        print()
        print(f"Sample: {sample}")
        for n, cmd in enumerate(cmds, start=1):
            print(f"  {n}) {shlex.join(cmd)}")
