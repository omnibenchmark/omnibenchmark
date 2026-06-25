"""Integration test for `ob validate outputs`.

Runs a real host-backend benchmark whose bundled dummy modules each produce a
JSON file, then runs the output validators (a non-empty check) and asserts the
recorded pass/fail verdicts. No network or conda — the module is a local bundle.
"""

import json
from pathlib import Path

import pytest

from tests.cli.cli_setup import OmniCLISetup

# Host-backend plan: two dummy modules under stage `data`, each writing
# {dataset}_data.json. Validators run in the `host` env.
BENCHMARK = """\
id: ValidatorsIntegration
description: validators integration test
version: "1.0"
benchmarker: "test"
benchmark_yaml_spec: "0.3"
software_backend: host
software_environments:
  host:
    description: "Host environment for direct execution"
stages:
  - id: data
    modules:
      - id: D1
        name: "Dataset 1"
        software_environment: "host"
        repository:
          url: bundles/dummymodule_4ff8427.bundle
          commit: 4ff8427
        parameters:
          - evaluate: "1+1"
      - id: D2
        name: "Dataset 2"
        software_environment: "host"
        repository:
          url: bundles/dummymodule_4ff8427.bundle
          commit: 4ff8427
        parameters:
          - evaluate: "2+2"
    outputs:
      - id: data.raw
        path: "{dataset}_data.json"
validators:
  env: host
"""

# The validator a stage `data` output named `data.json` gets: pass iff non-empty.
NONEMPTY_VALIDATOR = """\
#!/usr/bin/env python3
import os, sys
sys.exit(0 if os.path.getsize(sys.argv[1]) > 0 else 1)
"""


@pytest.fixture
def project(tmp_path):
    """A benchmark dir with the bundled module, plan, and one validator."""
    bundles = Path(__file__).parent.parent / "data" / "bundles"
    (tmp_path / "bundles").symlink_to(bundles, target_is_directory=True)
    (tmp_path / "benchmark.yaml").write_text(BENCHMARK)
    # keyed by (stage_id, output_id): stage "data", output id "data.raw"
    vdir = tmp_path / "validators" / "data" / "data.raw"
    vdir.mkdir(parents=True)
    (vdir / "validate.py").write_text(NONEMPTY_VALIDATOR)
    return tmp_path


def _verdicts(out_dir: Path) -> dict:
    """Map result-file name -> exit_code from out/.validation/**.json."""
    return {
        p.name: json.loads(p.read_text())["exit_code"]
        for p in (out_dir / ".validation").rglob("*_data.json.json")
    }


@pytest.mark.integration
def test_validate_outputs_pass_and_fail(project):
    out = project / "out"

    with OmniCLISetup() as omni:
        run = omni.call(["run", "benchmark.yaml", "--out-dir", "out"], cwd=str(project))
    assert run.returncode == 0, run.stdout + run.stderr
    assert list(out.rglob("*_data.json")), "benchmark produced no outputs"

    # 1. Both outputs are non-empty -> both pass.
    with OmniCLISetup() as omni:
        val = omni.call(
            ["validate", "outputs", "benchmark.yaml", "--out-dir", "out"],
            cwd=str(project),
        )
    assert val.returncode == 0, val.stdout + val.stderr
    verdicts = _verdicts(out)
    assert len(verdicts) == 2, verdicts
    assert all(code == 0 for code in verdicts.values()), verdicts

    # 2. Empty the D2 output -> it must now fail, D1 still passes.
    for p in out.rglob("D2_data.json"):
        p.write_text("")
    with OmniCLISetup() as omni:
        val = omni.call(
            ["validate", "outputs", "benchmark.yaml", "--out-dir", "out", "--force"],
            cwd=str(project),
        )
    assert val.returncode == 0, val.stdout + val.stderr
    verdicts = _verdicts(out)
    assert verdicts, verdicts
    assert all(code == 1 for name, code in verdicts.items() if "D2" in name)
    assert all(code == 0 for name, code in verdicts.items() if "D1" in name)
