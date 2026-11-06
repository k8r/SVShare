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

    _vcf, run = callers.run_cutesv("HG002.bam", "HG002", "reference.fa", tmp_path, "hifi", 3)

    assert run.command[4] == str(work_dir)
    assert seen == [[]]
    assert not work_dir.exists()


# cuteSV gets the minimum support and its settings for the platform.
def test_cutesv_settings(tmp_path, monkeypatch):
    monkeypatch.setattr(callers, "run_command", lambda cmd: None)

    _vcf, run = callers.run_cutesv("HG002.bam", "HG002", "reference.fa", tmp_path, "ont", 4)

    assert (run.sample, run.tool) == ("HG002", "cutesv")
    cmd = run.command
    assert cmd[cmd.index("--min_support") + 1] == "4"
    assert cmd[cmd.index("--max_cluster_bias_DEL") + 1] == "100"


def test_cutesv_min_support_is_a_fraction_of_depth_but_at_least_2():
    assert callers.cutesv_min_support(10.2, 0.3) == 3
    assert callers.cutesv_min_support(31.0, 0.3) == 9
    assert callers.cutesv_min_support(4.0, 0.3) == 2
