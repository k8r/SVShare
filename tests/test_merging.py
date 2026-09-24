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
