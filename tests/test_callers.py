from svshare import callers


# Working files left by an interrupted run are cleared first, and the folder is
# removed once cuteSV finishes.
def test_cutesv_gets_a_cleared_work_dir(tmp_path, monkeypatch):
    work_dir = tmp_path / "HG002.cutesv_work"
    (work_dir / "signatures").mkdir(parents=True)
    (work_dir / "signatures" / "_chr1_0_10000.bed").write_text("left over")

    seen = []
    def fake_run_command(cmd):
        seen.append(sorted(path.name for path in work_dir.iterdir()))
    monkeypatch.setattr(callers, "run_command", fake_run_command)

    _vcf, cmd = callers.run_cutesv("HG002.bam", "reference.fa", tmp_path)

    assert cmd[4] == str(work_dir)
    assert seen == [[]]
    assert not work_dir.exists()
