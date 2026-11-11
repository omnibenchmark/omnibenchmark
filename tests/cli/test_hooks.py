"""Run hooks (design 016)."""

from unittest.mock import MagicMock, patch

import pytest

import omnibenchmark.config
from omnibenchmark.cli.run import _run_benchmark
from omnibenchmark.config import ConfigAccessor
from omnibenchmark.core import hooks
from omnibenchmark.model import SoftwareBackendEnum


@pytest.fixture
def cfg(tmp_path, monkeypatch):
    c = ConfigAccessor(tmp_path / "ob.cfg")
    monkeypatch.setattr(omnibenchmark.config, "config", c)
    monkeypatch.delenv("OB_HOOKS_PATH", raising=False)
    return c


@pytest.fixture
def plan(tmp_path):
    p = tmp_path / "plan.yaml"
    p.write_text("id: x\n")
    return p


def script(path, body):
    path.write_text("#!/bin/sh\n" + body + "\n")
    path.chmod(0o755)
    return path


# Records its arguments and hook env into OB_HOOKS_OUT, then exits with $RC.
RECORD = 'mkdir -p "$OB_HOOKS_OUT"; echo "$OB_HOOK $* $OB_PLAN_ID" >> "$OB_HOOKS_OUT/calls"; echo "$OB_SELECTION" > "$OB_HOOKS_OUT/selection"; exit ${RC:-0}'


@pytest.mark.short
def test_nothing_configured_runs_nothing(cfg, plan, tmp_path):
    assert hooks.find("pre-run") is None
    with patch.object(hooks, "logger") as log:
        hooks.log_active()
        assert hooks.run("pre-run", plan, tmp_path / "out") == 0
    log.info.assert_not_called()
    assert not (tmp_path / "out").exists()


@pytest.mark.short
def test_config_command_wins_over_hooks_path(cfg, plan, tmp_path):
    d = tmp_path / "hooks"
    d.mkdir()
    script(d / "post-run", "exit 9")
    cfg.set("hooks", "path", str(d))
    assert hooks.find("post-run")[0] == [str(d / "post-run")]
    cfg.set("hooks", "post-run", f"sh -c '{RECORD}' hook")
    out = tmp_path / "out"
    with patch.object(hooks, "logger") as log:
        hooks.log_active()
        rc = hooks.run("post-run", plan, out, 2, env={"OB_PLAN_ID": "x"})
    assert rc == 0
    assert "hooks.post-run" in log.info.call_args[0][0]
    assert (out / ".hooks" / "calls").read_text().split() == [
        "post-run",
        str(plan),
        str(out.resolve()),
        "2",
        "x",
    ]
    log.warning.assert_not_called()  # writing to OB_HOOKS_OUT is fine


@pytest.mark.short
def test_timeout_and_plan_edit_count_as_failure(cfg, plan, tmp_path):
    cfg.set("hooks", "timeout", "0.2")
    cfg.set("hooks", "pre-run", "sleep 5")
    assert hooks.run("pre-run", plan, tmp_path / "out") != 0
    cfg.set("hooks", "pre-run", f"sh -c 'echo edit >> {plan}'")
    with patch.object(hooks, "logger") as log:
        assert hooks.run("pre-run", plan, tmp_path / "out") != 0
    assert "plan modified" in log.error.call_args[0][0]


@pytest.mark.short
def test_post_run_touching_outputs_warns(cfg, plan, tmp_path):
    out = tmp_path / "out"
    (out / "DATA").mkdir(parents=True)
    (out / "DATA" / "x.tsv").write_text("1\n")
    cfg.set("hooks", "post-run", f"sh -c 'echo 2 >> {out}/DATA/x.tsv'")
    with patch.object(hooks, "logger") as log:
        assert hooks.run("post-run", plan, out, 0) == 0
    assert "DATA/x.tsv" in log.warning.call_args[0][0]


def _run(tmp_path, dry=False, snakemake=None, **kw):
    with (
        patch("omnibenchmark.cli.run.BenchmarkExecution") as be,
        patch("omnibenchmark.cli.run.populate_git_cache"),
        patch("omnibenchmark.cli.run._generate_explicit_snakefile"),
        patch("omnibenchmark.cli.run.write_run_manifest") as manifest,
        patch("omnibenchmark.cli.run._run_snakemake", side_effect=snakemake),
    ):
        b = MagicMock()
        b.model.id, b.model.provenance = "x", None
        b.get_benchmark_software_backend.return_value = SoftwareBackendEnum.host
        be.return_value = b
        with pytest.raises((SystemExit, KeyboardInterrupt)) as exc:
            _run_benchmark(
                benchmark_path=str(tmp_path / "plan.yaml"),
                cores=1,
                dry=dry,
                continue_on_error=False,
                out_dir=str(tmp_path / "out"),
                debug=False,
                dirty=False,
                **kw,
            )
    return getattr(exc.value, "code", None), manifest


@pytest.mark.short
def test_pre_run_veto_exits_3_before_writing(cfg, plan, tmp_path):
    cfg.set("hooks", "pre-run", "false")
    code, manifest = _run(tmp_path)
    assert code == hooks.VETO
    manifest.assert_not_called()
    assert not (tmp_path / "out").exists()  # nothing written but the hook's own


@pytest.mark.short
def test_dry_runs_no_hooks(cfg, plan, tmp_path, monkeypatch):
    monkeypatch.setenv("RC", "1")
    cfg.set("hooks", "pre-run", f"sh -c '{RECORD}' hook")
    code, _ = _run(tmp_path, dry=True)
    assert code == 0
    assert not (tmp_path / "out" / ".hooks").exists()


@pytest.mark.short
def test_post_run_gets_snakemake_exit_code(cfg, plan, tmp_path):
    cfg.set("hooks", "post-run", f"sh -c '{RECORD}' hook")
    code, _ = _run(tmp_path, snakemake=SystemExit(2))
    assert code == 2
    calls = (tmp_path / "out" / ".hooks" / "calls").read_text().split()
    assert calls[0] == "post-run" and calls[3] == "2"


@pytest.mark.short
def test_post_run_sees_ctrl_c_as_130(cfg, plan, tmp_path):
    cfg.set("hooks", "post-run", f"sh -c '{RECORD}' hook")
    _run(tmp_path, snakemake=KeyboardInterrupt())  # re-raised, caught by _run
    calls = (tmp_path / "out" / ".hooks" / "calls").read_text().split()
    assert calls[3] == "130"


@pytest.mark.short
def test_snapshot_skips_envs_and_module_checkouts(tmp_path):
    for d in ("DATA", ".snakemake/conda/x", ".modules/m"):
        (tmp_path / d).mkdir(parents=True)
        (tmp_path / d / "f").write_text("1")
    assert set(hooks._snapshot(tmp_path)) == {tmp_path / "DATA" / "f"}


@pytest.mark.short
def test_hooks_get_the_selection(cfg, plan, tmp_path):
    cfg.set("hooks", "post-run", f"sh -c '{RECORD}' hook")
    _run(tmp_path, snakemake=SystemExit(0), selection={"until": "DATA"})
    sel = (tmp_path / "out" / ".hooks" / "selection").read_text()
    assert sel.strip() == '{"until": "DATA"}'
