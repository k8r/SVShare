# `analyze` subcommand: call, compare, and annotate SVs across samples and callers.
import sys
from pathlib import Path

from .callers import CALLERS, caller_vcf_path, cutesv_min_support
from .database import save_tool_runs, write_database
from .gaps import add_reference_gaps
from .genes import add_genes
from .transcripts import add_sv_transcripts
from .merging import VCFPreppedForMerging, drop_0_0_calls, run_jasmine
from .population import DATASETS, add_population_frequencies
from .reference import check_bam_reference_compatibility, ensure_reference_index
from .sequencing import detect_platform, mean_depth

# Jasmine's merge of every VCF from every sample across callers.
MERGED_VCF = "all_calls.merged.vcf"
DATABASE = "svshare.db"


def run(args):
    _make_paths_absolute(args)
    ensure_reference_index(args.reference)
    platforms = _check_samples(args)

    args.output.mkdir(parents=True, exist_ok=True)
    tool_runs, prepped_vcfs = _call_svs(args, platforms)
    tool_runs.append(_merge_calls(args, prepped_vcfs))
    db_path = write_database(args.output / MERGED_VCF, prepped_vcfs, args.output / DATABASE)

    _find_reference_gaps(db_path, args.reference)
    _look_up_population_frequencies(db_path, args)
    gene_run = _find_genes(db_path, args)
    if gene_run:
        tool_runs.append(gene_run)

    save_tool_runs(db_path, tool_runs)
    print_summary(args, tool_runs)


# Make paths absolute so the paths and commands in the summary and database
# work from any directory.
def _make_paths_absolute(args):
    args.samples = [bam.absolute() for bam in args.samples]
    for option in ("reference", "output", "vcf_dir", "gnomad", "colorsdb", "gencode",
                   "regions_bed"):
        path = getattr(args, option)
        if path:
            setattr(args, option, path.absolute())


# Check each sample matches the reference and find its sequencing platform,
# exiting if either fails. Returns each sample's platform.
def _check_samples(args):
    all_issues = []
    platforms = {}
    print(f"Analyzing {len(args.samples)} sample(s):")
    for bam in args.samples:
        print(f"  - {bam}")
        issues = check_bam_reference_compatibility(bam, args.reference)
        for issue in issues:
            print(f"    ! {issue}")
        all_issues.extend(issues)
        # The platform only matters when the callers run.
        if not args.vcf_dir:
            platforms[bam] = args.platform or detect_platform(bam)
            print(f"    platform: {platforms[bam] or 'unknown'}")

    if all_issues:
        sys.exit(f"{len(all_issues)} contig mismatch(es) found; aborting.")
    unknown = [str(bam) for bam, platform in platforms.items() if not platform]
    if unknown:
        sys.exit(f"Couldn't tell which platform sequenced {', '.join(unknown)}; "
                 "pass --platform hifi, ont, or clr.")
    return platforms


# Call SVs in every sample with every caller. Returns every outside tool run, in
# order, and the VCFs prepped for merging.
def _call_svs(args, platforms):
    tool_runs, prepped_vcfs = [], []
    for i, bam in enumerate(args.samples, start=1):
        print(f"\nSample {i}/{len(args.samples)}: {bam}")
        sample_runs, sample_vcfs = _call_sample_svs(args, bam, platforms.get(bam))
        tool_runs.extend(sample_runs)
        prepped_vcfs.extend(sample_vcfs)
    return tool_runs, prepped_vcfs


# Call SVs in one sample with each caller, or use existing VCFs from --vcf-dir.
# Returns its tool runs and VCFs.
def _call_sample_svs(args, bam, platform):
    sample = Path(bam).stem
    tool_runs = []
    if not args.vcf_dir:
        print("  measuring read depth with mosdepth")
        depth, tool_run = mean_depth(bam, sample, args.output)
        tool_runs.append(tool_run)
        min_support = cutesv_min_support(depth, args.cutesv_min_support)
        print(f"  mean depth {depth:.1f}x, so cuteSV needs {min_support} reads to call an SV")
        # Settings only some callers take.
        settings = {"cutesv": {"platform": platform, "min_support": min_support}}

    prepped_vcfs = []
    for name, caller in CALLERS.items():
        if args.vcf_dir:
            vcf = caller_vcf_path(bam, args.vcf_dir, name)
            print(f"  using existing {name} VCF: {vcf}")
        else:
            print(f"  calling SVs with {name}")
            vcf, tool_run = caller(bam, sample, args.reference, args.output,
                                   **settings.get(name, {}))
            tool_runs.append(tool_run)
        vcf, tool_run = drop_0_0_calls(vcf, sample, args.output)
        tool_runs.append(tool_run)
        prepped_vcfs.append(VCFPreppedForMerging(sample, name, vcf))
    return tool_runs, prepped_vcfs


# Group calls of the same SV across all samples and callers. Returns the tool run.
def _merge_calls(args, prepped_vcfs):
    print("\nMerging all samples and callers with Jasmine")
    vcfs = [prepped_vcf.vcf for prepped_vcf in prepped_vcfs]
    _vcf, tool_run = run_jasmine(vcfs, args.output / MERGED_VCF, args.output)
    return tool_run


def _find_reference_gaps(db_path, reference):
    print("\nFinding SVs near gaps in the reference")
    near, total = add_reference_gaps(db_path, reference)
    print(f"  {near} of {total} SVs are at or near a gap")


def _look_up_population_frequencies(db_path, args):
    population_vcfs = {"gnomad": args.gnomad, "colorsdb": args.colorsdb}
    for dataset, vcf in population_vcfs.items():
        name = DATASETS[dataset]
        if vcf:
            print(f"\nLooking up SVs in {name}")
            matched, total = add_population_frequencies(db_path, dataset, vcf)
            print(f"  {matched} of {total} SVs found in {name}")
        else:
            print(f"\nNo --{dataset} file given, so SVs weren't looked up in {name}")


# Find the genes each SV overlaps, if --gencode was given. Returns the tool run,
# or None.
def _find_genes(db_path, args):
    if not args.gencode:
        print("\nNo --gencode file given, so SVs weren't matched to genes")
        return None
    print("\nFinding the genes each SV affects")
    if args.regions_bed and args.regions_bed.exists():
        print(f"  using existing transcript regions: {args.regions_bed}")
    gene_run, affected, total = add_sv_transcripts(
        db_path, args.gencode, args.output, args.regions_bed
    )
    print(f"  {affected} of {total} SVs overlap a gene")
    with_summary, with_dosage, genes = add_genes(db_path)
    print(f"  {genes} genes nearby: {with_summary} with a summary, "
          f"{with_dosage} with ClinGen dosage scores")
    return gene_run


# Print paths and commands, one per line for easy rerunning.
def print_summary(args, tool_runs):
    print()
    print(f"Analyzed {len(args.samples)} sample(s):")
    for bam in args.samples:
        print(f"  - {bam}")
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
    samples = list(dict.fromkeys(tool_run.sample for tool_run in tool_runs))
    for sample in samples:
        print()
        print(f"Sample: {sample}" if sample else "All samples and callers:")
        for tool_run in tool_runs:
            if tool_run.sample == sample:
                print(f"  {tool_run}")
