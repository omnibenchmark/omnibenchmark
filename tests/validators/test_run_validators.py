"""In-process tests for run_validators and the `ob validate` commands.

Snakemake and the status machinery are stubbed: these cover the exit paths and
the verdict gate, not the generated workflow (test_validate_outputs does that).
"""

import json
import logging
from pathlib import Path
from types import SimpleNamespace

import pytest
from click.testing import CliRunner

from omnibenchmark.backend import validators as bv
from omnibenchmark.cli.describe import status as describe_status
from omnibenchmark.cli.validate import validate
from omnibenchmark.core.validators import build_snakefile

PLAN = """\
id: b
version: "1.0"
benchmarker: t
api_version: "0.8.0"
software_backend: host
software_environments:
  host:
    description: host
stages:
  - id: data
    modules:
      - id: D1
        software_environment: host
        repository: {url: https://example.org/x.git, commit: abc1234}
    outputs:
      - id: raw
        path: "{module.id}_data.json"
__VALIDATORS__"""

TWO_VALIDATORS = """\
        validators:
          - script: checks/a.py
          - script: checks/b.sh
"""


@pytest.fixture
def project(tmp_path, monkeypatch, caplog):
    """A plan with two validators on `data/raw` and one produced file."""
    caplog.set_level(logging.INFO, logger="omnibenchmark")
    (tmp_path / "checks").mkdir()
    (tmp_path / "checks" / "a.py").write_text("")
    (tmp_path / "checks" / "b.sh").write_text("")
    out = tmp_path / "out"
    produced = out / "data.D1" / "D1_data.json"
    produced.parent.mkdir(parents=True)
    produced.write_text("{}")
    _write_plan(tmp_path, TWO_VALIDATORS)
    _stub_status(monkeypatch, [produced])
    return SimpleNamespace(dir=tmp_path, out=out, rel="data.D1/D1_data.json")


def _write_plan(d: Path, validators: str, extra: str = "") -> None:
    (d / "benchmark.yaml").write_text(
        PLAN.replace("__VALIDATORS__", validators) + extra
    )


def _stub_status(monkeypatch, files):
    step = SimpleNamespace(
        stage_id="data",
        node=SimpleNamespace(get_outputs=lambda: ["raw"]),
        output_stage_files=[
            SimpleNamespace(output_file=f, output_file_path=str(f)) for f in files
        ],
    )
    exec_paths = {"p": SimpleNamespace(exec_path={"n": step})}
    monkeypatch.setattr(
        "omnibenchmark.core.status.status.prepare_status",
        lambda *a, **k: (None, None, exec_paths),
    )


def _stub_snakemake(monkeypatch, exit_codes=None, rc=0, log=""):
    """Fake snakemake: write a verdict per key in *exit_codes*, print *log*."""
    calls = []

    def run(cmd, cwd, stdout, stderr):
        calls.append(cmd)
        stdout.write(log)
        for key, code in (exit_codes or {}).items():
            result = Path(cwd) / ".validation" / f"{key}.json"
            result.parent.mkdir(parents=True, exist_ok=True)
            result.write_text(json.dumps({"exit_code": code}))
        return SimpleNamespace(returncode=rc)

    monkeypatch.setattr(bv.subprocess, "run", run)
    return calls


def _run(p, **kw):
    return bv.run_validators(str(p.dir / "benchmark.yaml"), str(p.out), **kw)


@pytest.mark.short
def test_all_pass(project, monkeypatch):
    r = project.rel
    calls = _stub_snakemake(monkeypatch, {f"{r}.0": 0, f"{r}.1": 0})
    assert _run(project) == 0
    assert "--forceall" not in calls[0]
    assert _run(project, force=True) == 0
    assert "--forceall" in calls[1]


@pytest.mark.short
def test_any_fail_or_error_fails_the_gate(project, monkeypatch, caplog):
    r = project.rel
    _stub_snakemake(monkeypatch, {f"{r}.0": 0, f"{r}.1": 1})
    assert _run(project) == 1
    assert "Validation failed: data.D1/D1_data.json by b.sh" in caplog.text

    _stub_snakemake(monkeypatch, {f"{r}.0": 0, f"{r}.1": 127})
    assert _run(project) == 1
    assert "Validator could not run" in caplog.text


@pytest.mark.short
def test_missing_verdict_is_an_error(project, monkeypatch, caplog):
    _stub_snakemake(monkeypatch, {f"{project.rel}.0": 0})
    assert _run(project) == 1
    assert "Validator could not run" in caplog.text


@pytest.mark.short
def test_orphan_results_are_dropped(project, monkeypatch):
    orphan = project.out / ".validation" / f"{project.rel}.2.json"
    orphan.parent.mkdir(parents=True)
    orphan.write_text("{}")
    r = project.rel
    _stub_snakemake(monkeypatch, {f"{r}.0": 0, f"{r}.1": 0})
    assert _run(project) == 0
    assert not orphan.exists()


@pytest.mark.short
def test_snakemake_failure_shows_log_tail(project, monkeypatch, caplog):
    _stub_snakemake(monkeypatch, rc=2, log="boom\n")
    assert _run(project) == 2
    assert "boom" in caplog.text
    assert "--unlock" not in caplog.text


@pytest.mark.short
def test_lock_error_gets_a_hint(project, monkeypatch, caplog):
    _stub_snakemake(monkeypatch, rc=1, log="Error: Directory cannot be locked.\n")
    assert _run(project) == 1
    assert "snakemake --unlock" in caplog.text


@pytest.mark.short
def test_env_directive_and_flag_reach_snakemake(project, monkeypatch):
    _write_plan(project.dir, TWO_VALIDATORS, "validators:\n  env: host\n")
    monkeypatch.setattr(
        bv, "_resolve_env_directive", lambda *a: ('conda: "/e.yaml"', "--use-conda")
    )
    r = project.rel
    calls = _stub_snakemake(monkeypatch, {f"{r}.0": 0, f"{r}.1": 0})
    assert _run(project) == 0
    assert "--use-conda" in calls[0]
    snakefile = (project.out / ".validation" / "Snakefile").read_text()
    assert '    conda: "/e.yaml"\n' in snakefile


@pytest.mark.short
def test_env_resolution_failure(project, monkeypatch, caplog):
    _write_plan(project.dir, TWO_VALIDATORS, "validators:\n  env: host\n")

    def fail(*a):
        raise ValueError("could not resolve validators env 'host'")

    monkeypatch.setattr(bv, "_resolve_env_directive", fail)
    calls = _stub_snakemake(monkeypatch)
    assert _run(project) == 1
    assert "could not resolve" in caplog.text
    assert not calls


@pytest.mark.short
def test_no_validators(project, caplog):
    _write_plan(project.dir, "")
    assert _run(project) == 1
    assert "declares no validators" in caplog.text


@pytest.mark.short
def test_missing_script(project, caplog):
    (project.dir / "checks" / "b.sh").unlink()
    assert _run(project) == 1
    assert "b.sh" in caplog.text


@pytest.mark.short
def test_missing_out_dir(project, caplog):
    assert bv.run_validators(str(project.dir / "benchmark.yaml"), "nope") == 1
    assert "not found" in caplog.text


@pytest.mark.short
def test_nothing_produced(project, monkeypatch, caplog):
    _stub_status(monkeypatch, [project.out / "absent.json"])
    assert _run(project) == 1
    assert "No produced outputs" in caplog.text


@pytest.mark.short
def test_build_snakefile_env_line():
    targets = {"a.txt.0": ["a.txt", "python3", "/v.py", ""]}
    with_env = build_snakefile(targets, 'conda: "/e.yaml"')
    assert '        ".validation/{key}.log"\n    conda: "/e.yaml"\n' in with_env
    assert '"a.txt.0"' in with_env
    assert "conda:" not in build_snakefile(targets, "")


@pytest.mark.short
def test_cli_validate_plan_exit_codes(project):
    plan = str(project.dir / "benchmark.yaml")
    ok = CliRunner().invoke(validate, ["plan", plan])
    assert ok.exit_code == 0, ok.output

    (project.dir / "checks" / "b.sh").unlink()
    bad = CliRunner().invoke(validate, ["plan", plan])
    assert bad.exit_code == 1
    assert "unexpected error" not in bad.output


@pytest.mark.short
def test_cli_validate_outputs_reports_a_broken_plan(project):
    (project.dir / "benchmark.yaml").write_text("id: b\n")
    res = CliRunner().invoke(validate, ["outputs", str(project.dir / "benchmark.yaml")])
    assert res.exit_code == 1
    assert res.exception is None or isinstance(res.exception, SystemExit)


@pytest.mark.short
def test_cli_validate_outputs_passes_exit_code(project, monkeypatch):
    seen = {}

    def fake(**kw):
        seen.update(kw)
        return 3

    monkeypatch.setattr(bv, "run_validators", fake)
    plan = str(project.dir / "benchmark.yaml")
    res = CliRunner().invoke(validate, ["outputs", plan, "--out-dir", "o", "-c", "2"])
    assert res.exit_code == 3
    assert seen["out_dir"] == "o" and seen["cores"] == 2 and not seen["force"]


@pytest.mark.short
def test_refuses_an_out_dir_of_another_version(project, monkeypatch, caplog):
    md = project.out / ".metadata"
    md.mkdir()
    bench = {"id": "b", "version": "0.9", "summary_hash": "x"}
    (md / "runs.jsonl").write_text(json.dumps({"benchmark": bench}) + "\n")
    calls = _stub_snakemake(monkeypatch)
    assert _run(project) == 1
    assert "version '0.9'" in caplog.text
    assert not calls

    res = CliRunner().invoke(
        describe_status,
        [str(project.dir / "benchmark.yaml"), "--out-dir", str(project.out)],
    )
    assert res.exit_code == 1
