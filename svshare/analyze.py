# `analyze` subcommand: call, compare, and annotate SVs across samples and callers.
import shlex
import sys
from pathlib import Path

from .callers import CALLERS, caller_vcf_path
from .database import write_database
from .merging import run_jasmine
from .reference import check_bam_reference_compatibility, ensure_reference_index

# Jasmine's merge of every VCF from every sample across callers.
MERGED_VCF = "all_calls.merged.vcf"
DATABASE = "svshare.db"


def run(args):
    # Use absolute paths so the commands run, shown in the summary and stored in
    # the database, work from any directory. absolute() rather than resolve(), so
    # a symlinked reference still finds the .fai next to the link.
    args.samples = [sample.absolute() for sample in args.samples]
    args.reference = args.reference.absolute()
    args.output = args.output.absolute()
    if args.vcf_dir:
        args.vcf_dir = args.vcf_dir.absolute()

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
    # (sample, tool, command) for every command run, in order, for the database.
    run_commands = []
    # (sample, caller, vcf) for the VCF from every SV caller for every sample, in
    # the order Jasmine reads them.
    inputs = []
    for i, sample in enumerate(args.samples, start=1):
        print(f"\nSample {i}/{len(args.samples)}: {sample}")
        for name, caller in CALLERS.items():
            if args.vcf_dir:
                vcf = caller_vcf_path(sample, args.vcf_dir, name)
                print(f"  using existing {name} VCF: {vcf}")
            else:
                print(f"  calling SVs with {name}")
                vcf, cmd = caller(sample, args.reference, args.output)
                commands[sample].append(cmd)
                run_commands.append((Path(sample).stem, name, cmd))
            inputs.append((Path(sample).stem, name, vcf))

    # Group calls of the same SV across all samples and callers.
    print("\nMerging all samples and callers with Jasmine")
    merged_vcf = args.output / MERGED_VCF
    _vcf, merge_cmd = run_jasmine([vcf for _s, _c, vcf in inputs], merged_vcf, args.output)
    run_commands.append((None, "jasmine", merge_cmd))
    write_database(merged_vcf, inputs, args.output / DATABASE, run_commands)

    print_summary(args, commands, merge_cmd)


# Print the samples, reference, and output location, then the caller commands
# run for each sample, and the Jasmine merge, one per line so each can be copied
# and rerun.
def print_summary(args, commands, merge_cmd):
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
        if not cmds:
            continue
        print()
        print(f"Sample: {sample}")
        for cmd in cmds:
            print(f"  {shlex.join(cmd)}")

    print()
    print("All samples and callers:")
    print(f"  {shlex.join(merge_cmd)}")
