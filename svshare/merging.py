# Wrapper around Jasmine, which merges SV calls from several VCFs into one VCF.
# Like the callers, run_jasmine returns (vcf_path, command) so analyze can
# report exactly what was run.
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


# Merge the given VCFs into out_vcf. Jasmine reads its inputs from a text file
# listing one VCF path per line; it's saved in the same directory as out_vcf and
# named after it (e.g. the input list for HG002.merged.vcf is
# HG002.merged.jasmine_inputs.txt).
def run_jasmine(vcfs, out_vcf, out_dir, threads=4):
    out_vcf = Path(out_vcf)
    input_list = out_vcf.with_suffix(".jasmine_inputs.txt")
    input_list.write_text("".join(f"{vcf}\n" for vcf in vcfs))
    cmd = [
        "jasmine",
        f"file_list={input_list}",
        f"out_file={out_vcf}",
        f"out_dir={out_dir}",
        f"threads={threads}",
        *(f"{name}={value}" for name, value in JASMINE_SETTINGS.items()),
    ]
    # Bioconda's jasmine launcher is a bash script with no #! line, so it can't be
    # executed directly from Python; run it with bash. The returned cmd stays the
    # plain `jasmine ...` form, which works as-is from a terminal.
    run_command(["bash", shutil.which("jasmine") or "jasmine", *cmd[1:]])
    return out_vcf, cmd
