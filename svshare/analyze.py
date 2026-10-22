# `analyze` subcommand: call, compare, and annotate SVs across samples and callers.
import shlex
import sys
from pathlib import Path

from .callers import CALLERS, caller_vcf_path
from .database import write_database
from .genes import add_genes
from .transcripts import add_sv_transcripts
from .merging import run_jasmine
from .population import DATASETS, add_population_frequencies
from .reference import check_bam_reference_compatibility, ensure_reference_index

# Jasmine's merge of every VCF from every sample across callers.
MERGED_VCF = "all_calls.merged.vcf"
DATABASE = "svshare.db"


def run(args):
    # Make paths absolute so the paths and commands in the summary and database
    # work from any directory.
    args.samples = [sample.absolute() for sample in args.samples]
    args.reference = args.reference.absolute()
    args.output = args.output.absolute()
    if args.vcf_dir:
        args.vcf_dir = args.vcf_dir.absolute()
    if args.gnomad:
        args.gnomad = args.gnomad.absolute()
    if args.colorsdb:
        args.colorsdb = args.colorsdb.absolute()
    if args.gencode:
        args.gencode = args.gencode.absolute()
    if args.regions_bed:
        args.regions_bed = args.regions_bed.absolute()

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
    db_path = write_database(merged_vcf, inputs, args.output / DATABASE, run_commands)

    population_vcfs = {"gnomad": args.gnomad, "colorsdb": args.colorsdb}
    for dataset, vcf in population_vcfs.items():
        name = DATASETS[dataset]
        if vcf:
            print(f"\nLooking up SVs in {name}")
            matched, total = add_population_frequencies(db_path, dataset, vcf)
            print(f"  {matched} of {total} SVs found in {name}")
        else:
            print(f"\nNo --{dataset} file given, so SVs weren't looked up in {name}")

    gene_command = None
    if args.gencode:
        print("\nFinding the genes each SV affects")
        if args.regions_bed and args.regions_bed.exists():
            print(f"  using existing transcript regions: {args.regions_bed}")
        gene_command, affected, total = add_sv_transcripts(
            db_path, args.gencode, args.output, args.regions_bed
        )
        print(f"  {affected} of {total} SVs overlap a gene")
        with_summary, with_dosage, genes = add_genes(db_path)
        print(f"  {genes} genes nearby: {with_summary} with a summary, "
              f"{with_dosage} with ClinGen dosage scores")
    else:
        print("\nNo --gencode file given, so SVs weren't matched to genes")

    print_summary(args, commands, merge_cmd, gene_command)


# Print paths and commands, one per line for easy rerunning.
def print_summary(args, commands, merge_cmd, gene_command):
    print()
    print(f"Analyzed {len(args.samples)} sample(s):")
    for sample in args.samples:
        print(f"  - {sample}")
    print(f"Reference: {args.reference}")
    print(f"Output directory: {args.output}")
    if args.vcf_dir:
        print(f"Existing caller VCFs from: {args.vcf_dir}")
    if args.gnomad:
        print(f"gnomAD: {args.gnomad}")
    if args.colorsdb:
        print(f"CoLoRSdb: {args.colorsdb}")
    if args.gencode:
        print(f"GENCODE: {args.gencode}")

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

    if gene_command:
        print()
        print("Gene lookup:")
        print(f"  {gene_command}")
