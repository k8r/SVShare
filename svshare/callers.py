# Wrappers around the external structural-variant callers (sniffles2, cuteSV).
# Each caller function has the same signature -- (bam, reference, out_dir, threads)
# -- and returns (vcf_path, command), so analyze can treat them uniformly and
# report exactly what was run.
import subprocess
from pathlib import Path


def _run(cmd):
    result = subprocess.run(cmd, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"{cmd[0]} failed (exit {result.returncode})\n{result.stderr.strip()}"
        )


def run_sniffles2(bam, reference, out_dir, threads=4):
    out_vcf = Path(out_dir) / f"{Path(bam).stem}.sniffles2.vcf"
    cmd = [
        "sniffles",
        "--input", str(bam),
        "--reference", str(reference),
        "--vcf", str(out_vcf),
        "--threads", str(threads),
        "--allow-overwrite",
    ]
    _run(cmd)
    return out_vcf, cmd


def run_cutesv(bam, reference, out_dir, threads=4):
    out_dir = Path(out_dir)
    out_vcf = out_dir / f"{Path(bam).stem}.cutesv.vcf"
    cmd = [
        "cuteSV",
        str(bam),
        str(reference),
        str(out_vcf),
        str(out_dir),
        "--threads", str(threads),
    ]
    _run(cmd)
    return out_vcf, cmd


CALLERS = {
    "sniffles2": run_sniffles2,
    "cutesv": run_cutesv,
}
