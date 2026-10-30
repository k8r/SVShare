from svshare import merging


# Checks that the input list keeps the VCFs in the order they were passed in.
# Jasmine records which VCF each SV came from by its position in this list.
def test_run_jasmine_lists_input_vcfs_in_order(tmp_path, monkeypatch):
    monkeypatch.setattr(merging, "run_command", lambda cmd: None)
    vcfs = [tmp_path / "HG002.sniffles2.vcf", tmp_path / "HG002.cutesv.vcf"]

    merging.run_jasmine(vcfs, tmp_path / "HG002.merged.vcf", tmp_path)

    # run_jasmine saves its input list in the same directory as the output VCF.
    input_list = tmp_path / "HG002.merged.jasmine_inputs.txt"
    assert input_list.read_text().splitlines() == [str(vcf) for vcf in vcfs]


# Only 0/0 calls are dropped; calls with an unknown genotype are kept.
def test_drop_0_0_calls(tmp_path):
    vcf = tmp_path / "vcfs" / "HG002.cutesv.vcf"
    vcf.parent.mkdir()
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=100000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        + "".join(f"chr1\t{pos}\t{vcf_id}\tN\t<DEL>\t.\tPASS\t.\tGT\t{gt}\n"
                  for pos, vcf_id, gt in [(100, "a", "0/1"), (200, "b", "0/0"),
                                          (300, "c", "./."), (400, "d", "1/1"),
                                          (500, "e", "0|0")])
    )

    out_vcf, _cmd = merging.drop_0_0_calls(vcf, tmp_path)

    assert out_vcf == tmp_path / "HG002.cutesv.called.vcf"
    ids = [line.split("\t")[2] for line in out_vcf.read_text().splitlines()
           if not line.startswith("#")]
    assert ids == ["a", "c", "d"]
