# Wrappers around the external structural-variant callers (sniffles2, cuteSV).
# Each caller function has the same signature -- (bam, reference, out_dir, threads)
# -- and returns (vcf_path, command), so analyze can treat them uniformly and
# report exactly what was run.
from pathlib import Path

from .commands import run_command


# Path of the VCF a caller writes for a BAM, e.g. results/HG002.sniffles2.vcf.
# analyze --vcf-dir expects existing VCFs to follow the same naming.
def caller_vcf_path(bam, out_dir, caller_name):
    return Path(out_dir) / f"{Path(bam).stem}.{caller_name}.vcf"


def run_sniffles2(bam, reference, out_dir, threads=4):
    out_vcf = caller_vcf_path(bam, out_dir, "sniffles2")
    cmd = [
        "sniffles",
        "--input", str(bam),
        "--reference", str(reference),
        "--vcf", str(out_vcf),
        "--threads", str(threads),
        "--allow-overwrite",
    ]
    run_command(cmd)
    return out_vcf, cmd


def run_cutesv(bam, reference, out_dir, threads=4):
    out_vcf = caller_vcf_path(bam, out_dir, "cutesv")
    cmd = [
        "cuteSV",
        str(bam),
        str(reference),
        str(out_vcf),
        str(out_dir),
        "--threads", str(threads),
    ]
    run_command(cmd)
    return out_vcf, cmd


CALLERS = {
    "sniffles2": run_sniffles2,
    "cutesv": run_cutesv,
}
