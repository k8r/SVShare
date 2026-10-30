# Wrappers around bcftools, which prepares each caller's VCF for merging, and
# Jasmine, which merges SV calls from several VCFs into one VCF.
import shutil
from pathlib import Path

from .commands import run_command

# How close two calls must be for Jasmine to count them as the same SV. Set
# here so they're recorded in the command.
JASMINE_SETTINGS = {
    "max_dist_linear": "0.5",
    "min_dist": "100",
    "max_dist": "1000000000",  # no limit
    "kd_tree_norm": "2",
    "min_seq_id": "0",
    "min_overlap": "0",
}


# Copy a caller's VCF into out_dir without its 0/0 calls, so Jasmine doesn't
# merge them. A 0/0 genotype means the sample doesn't have the SV.
def drop_0_0_calls(vcf, out_dir):
    out_vcf = Path(out_dir) / f"{Path(vcf).stem}.called.vcf"
    cmd = ["bcftools", "view", "--exclude", 'GT="ref"', str(vcf), "-o", str(out_vcf)]
    run_command(cmd)
    return out_vcf, cmd


# Merge the given VCFs into out_vcf.
def run_jasmine(vcfs, out_vcf, out_dir, threads=4):
    out_vcf = Path(out_vcf)
    # Jasmine reads its input VCFs from a file listing one per line.
    input_list = out_vcf.with_suffix(".jasmine_inputs.txt")
    input_list.write_text("".join(f"{vcf}\n" for vcf in vcfs))
    cmd = [
        "jasmine",
        f"file_list={input_list}",
        f"out_file={out_vcf}",
        f"out_dir={out_dir}",
        f"threads={threads}",
        *(f"{name}={value}" for name, value in JASMINE_SETTINGS.items()),
        # Without this, Jasmine recalculates each call's end and length from its
        # sequence, assuming the sequence starts with the reference base just
        # before the SV.
        # cuteSV's insertions don't, so their lengths would come out one base
        # short.
        "--leave_breakpoints",
    ]
    # Bioconda's jasmine launcher is a bash script with no #! line, so it can't be
    # executed directly from Python; run it with bash. The returned cmd stays the
    # plain `jasmine ...` form, which works as-is from a terminal.
    run_command(["bash", shutil.which("jasmine") or "jasmine", *cmd[1:]])
    return out_vcf, cmd
