"""Output validators under the apptainer backend.

The container mounts only out/, so a validator script that lives next to the
plan is invisible inside it unless it is staged under out/. Runs the bundled
dummy module and one validator in a python image and expects a pass.
Before the copy, this verdict was exit 127.
"""

import json
import shutil
import tempfile
from pathlib import Path

import pytest

from tests.cli.cli_setup import OmniCLISetup

BENCHMARK = """\
id: ValidatorsApptainer
version: "1.0"
benchmarker: "test"
api_version: "0.8.0"
software_backend: apptainer
software_environments:
  py:
    apptainer: "docker://python:3.12-slim"
stages:
  - id: data
    modules:
      - id: D1
        software_environment: py
        repository:
          url: bundles/dummymodule_4ff8427.bundle
          commit: 4ff8427
        parameters:
          - evaluate: "1+1"
    outputs:
      - id: raw
        path: "{module.id}_data.json"
        validators:
          - script: checks/nonempty.py
validators:
  env: py
"""

NONEMPTY = """\
import os, sys
sys.exit(0 if os.path.getsize(sys.argv[1]) > 0 else 1)
"""


@pytest.mark.e2e
@pytest.mark.skipif(not shutil.which("apptainer"), reason="needs apptainer")
def test_validators_run_inside_apptainer():
    # Not pytest's tmp_path: apptainer mounts /tmp by default, which would hide
    # the bug. The repo checkout is outside both /tmp and out/.
    with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as d:
        _run(Path(d))


def _run(tmp_path: Path):
    bundles = Path(__file__).parent.parent / "data" / "bundles"
    (tmp_path / "bundles").symlink_to(bundles, target_is_directory=True)
    (tmp_path / "benchmark.yaml").write_text(BENCHMARK)
    (tmp_path / "checks").mkdir()
    (tmp_path / "checks" / "nonempty.py").write_text(NONEMPTY)
    out = tmp_path / "out"

    with OmniCLISetup() as omni:
        run = omni.call(["run", "benchmark.yaml", "--out-dir", "out"], cwd=tmp_path)
        rule_logs = "".join(p.read_text() for p in (out / ".logs").glob("data_*.log"))
        assert run.returncode == 0, run.stdout + run.stderr + rule_logs
        val = omni.call(
            ["validate", "outputs", "benchmark.yaml", "--out-dir", "out"],
            cwd=tmp_path,
        )
    assert val.returncode == 0, val.stdout + val.stderr

    (verdict,) = (out / ".validation").rglob("D1_data.json.0.json")
    data = json.loads(verdict.read_text())
    assert data["exit_code"] == 0, (out / data["log"]).read_text()
    assert (out / ".validation" / "scripts" / data["script_sha256"]).is_file()
