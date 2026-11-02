# Wrappers around the external structural-variant callers (sniffles2, cuteSV).
import shutil
from pathlib import Path

from .commands import run_command

# cuteSV's recommended settings for each platform.
CUTESV_PLATFORM_SETTINGS = {
    "hifi": {"max_cluster_bias_INS": 1000, "diff_ratio_merging_INS": 0.9,
             "max_cluster_bias_DEL": 1000, "diff_ratio_merging_DEL": 0.5},
    "ont": {"max_cluster_bias_INS": 100, "diff_ratio_merging_INS": 0.3,
            "max_cluster_bias_DEL": 100, "diff_ratio_merging_DEL": 0.3},
    "clr": {"max_cluster_bias_INS": 100, "diff_ratio_merging_INS": 0.3,
            "max_cluster_bias_DEL": 200, "diff_ratio_merging_DEL": 0.5},
}


# Path of the VCF a caller writes for a BAM, e.g. results/HG002.sniffles2.vcf.
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


# The number of reads cuteSV needs to call an SV: a fraction of the sample's
# mean depth, but at least 2.
def cutesv_min_support(depth, fraction):
    return max(2, round(depth * fraction))


# Delete a directory and everything in it, if it exists.
def _remove_dir(path):
    # On some drives, macOS's hidden ._ files disappear with the file they belong
    # to, so a file can be gone by the time it's deleted.
    def ignore_missing(_function, _path, exc_info):
        if not issubclass(exc_info[0], FileNotFoundError):
            raise exc_info[1]
    shutil.rmtree(path, onerror=ignore_missing)


def run_cutesv(bam, reference, out_dir, platform, min_support, threads=4):
    out_vcf = caller_vcf_path(bam, out_dir, "cutesv")
    # cuteSV won't start if an interrupted run left its working files behind, so
    # it gets a folder of its own that's cleared first.
    work_dir = Path(out_dir) / f"{Path(bam).stem}.cutesv_work"
    _remove_dir(work_dir)
    work_dir.mkdir()
    cmd = [
        "cuteSV",
        str(bam),
        str(reference),
        str(out_vcf),
        str(work_dir),
        "--threads", str(threads),
        "--genotype",
        "--min_support", str(min_support),
        *(arg for name, value in CUTESV_PLATFORM_SETTINGS[platform].items()
          for arg in (f"--{name}", str(value))),
    ]
    run_command(cmd)
    _remove_dir(work_dir)
    return out_vcf, cmd


CALLERS = {
    "sniffles2": run_sniffles2,
    "cutesv": run_cutesv,
}
