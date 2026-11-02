from svshare.sequencing import _autosome_depth, _platform_from_header


def header(*command_lines):
    return {"PG": [{"ID": f"p{i}", "CL": cl} for i, cl in enumerate(command_lines)]}


def test_platform_from_minimap2_and_pbmm2_presets():
    assert _platform_from_header(header("minimap2 -ax map-hifi -t 8 ref.fa -")) == "hifi"
    assert _platform_from_header(header("minimap2 -x lr:hq ref.fa reads.fq")) == "ont"
    assert _platform_from_header(header("pbmm2 align --preset SUBREAD ref.fa in.bam")) == "clr"


# Repeating the same preset is fine; presets that disagree, or none, aren't.
def test_platform_is_none_unless_the_header_names_exactly_one():
    assert _platform_from_header(header("minimap2 -ax map-hifi a", "minimap2 -ax map-hifi b")) == "hifi"
    assert _platform_from_header(header("minimap2 -ax map-hifi a", "minimap2 -ax map-ont b")) is None
    assert _platform_from_header(header("samtools merge out.bam a.bam")) is None
    assert _platform_from_header({}) is None


def test_autosome_depth_leaves_out_other_chromosomes(tmp_path):
    summary = tmp_path / "HG002.mosdepth.summary.txt"
    summary.write_text(
        "chrom\tlength\tbases\tmean\tmin\tmax\n"
        "chr1\t100\t1000\t10.00\t0\t50\n"
        "chr2\t300\t6000\t20.00\t0\t50\n"
        "chrX\t100\t500\t5.00\t0\t50\n"
        "chrUn_KI270302v1\t10\t1000\t100.00\t0\t500\n"
        "total\t510\t8500\t16.67\t0\t500\n"
    )
    assert _autosome_depth(summary) == 17.5
