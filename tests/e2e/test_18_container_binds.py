"""Container jobs see the plan directory and the host TMPDIR.

Snakemake runs a container job with ``--home <out dir>``, so a file the plan
refers to next to itself is missing inside the job unless ob binds the plan
directory. Snakemake also exports the host TMPDIR into the job, where it does
not exist unless bound (#319). The module calls mktemp and copies a
plan-relative asset to its output.
"""

import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest

from tests.cli.cli_setup import OmniCLISetup

RUN_SH = """\
#!/usr/bin/env bash
set -euo pipefail
while [[ $# -gt 0 ]]; do
  case $1 in
    --asset) asset=$2; shift 2;;
    --output) out=${2#*=}; shift 2;;
    *) shift 2;;
  esac
done
mktemp > /dev/null
cp "$asset" "$out"
"""

PLAN = """\
id: bind_repro
version: "1.0"
benchmarker: "test"
api_version: "0.8.0"
software_backend: apptainer
software_environments:
  debian: {{description: "debian", apptainer: debian.sandbox}}
stages:
  - id: data
    modules:
      - id: copy
        software_environment: debian
        repository: {{url: {mod}, commit: {commit}}}
        parameters:
          - asset: {asset}
    outputs:
      - id: copied
        path: greeting.txt
"""


def _git(cwd: Path, *args: str) -> str:
    return subprocess.run(
        ["git", "-c", "user.name=t", "-c", "user.email=t@t", *args],
        cwd=cwd,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


@pytest.mark.e2e
@pytest.mark.skipif(not shutil.which("apptainer"), reason="needs apptainer")
def test_container_job_sees_plan_dir_and_tmpdir(monkeypatch):
    # Not pytest's tmp_path: apptainer mounts /tmp by default, which would hide
    # both bugs.
    with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as d:
        (Path(d) / "hosttmp").mkdir()
        monkeypatch.setenv("TMPDIR", str(Path(d) / "hosttmp"))
        _run(Path(d))


def _run(root: Path):
    mod = root / "mod"
    mod.mkdir()
    (mod / "run.sh").write_text(RUN_SH)
    (mod / "run.sh").chmod(0o755)
    (mod / "omnibenchmark.yaml").write_text("entrypoints:\n  default: run.sh\n")
    _git(mod, "init", "-q")
    _git(mod, "add", ".")
    _git(mod, "commit", "-qm", "init")
    commit = _git(mod, "rev-parse", "--short", "HEAD")

    plan = root / "plan"
    (plan / "assets").mkdir(parents=True)
    asset = plan / "assets" / "greeting.txt"
    asset.write_text("hello\n")
    # A sandbox runs without /dev/fuse, unlike a .sif under apptainer >= 1.5.
    subprocess.run(
        [
            "apptainer",
            "build",
            "--sandbox",
            "debian.sandbox",
            "docker://debian:stable-slim",
        ],
        cwd=plan,
        check=True,
        capture_output=True,
    )
    (plan / "plan.yml").write_text(PLAN.format(mod=mod, commit=commit, asset=asset))

    with OmniCLISetup() as omni:
        run = omni.call(["run", "plan.yml", "--dirty", "--cores", "1"], cwd=plan)
    logs = "".join(p.read_text() for p in (plan / "out" / ".logs").glob("data_*.log"))
    assert run.returncode == 0, run.stdout + run.stderr + logs
    (copied,) = (plan / "out").glob("data.copy.*/greeting.txt")
    assert copied.read_text() == "hello\n"
