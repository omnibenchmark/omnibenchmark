"""`[conda]` config under the conda backend.

Two runs of the same plan into different out dirs share `conda.cache`: the
first builds the environment there, the second reuses it. `conda.channel_alias`
decides where named channels are fetched from, so the installed packages carry
its URL.
"""

import json
import shutil
import sys
from pathlib import Path

import pytest

from tests.cli.cli_setup import OmniCLISetup

ALIAS = "https://conda.anaconda.org"

BENCHMARK = """\
id: CondaCache
version: "1.0"
benchmarker: "test"
api_version: "0.8.0"
software_backend: conda
software_environments:
  py:
    conda: envs/py.yaml
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
"""

ENV = """\
channels:
  - conda-forge
dependencies:
  - python=3.12
"""


@pytest.mark.e2e
@pytest.mark.skipif(not shutil.which("conda"), reason="needs conda")
def test_conda_cache_shared_across_runs(tmp_path: Path, monkeypatch):
    cache = tmp_path / "cache"
    cache.mkdir()
    # The config dir is XDG on Linux, ~/Library/Application Support on macOS.
    home = tmp_path / "home"
    if sys.platform == "darwin":
        cfg_dir = home / "Library" / "Application Support" / "omnibenchmark"
    else:
        cfg_dir = home / ".config" / "omnibenchmark"
    cfg_dir.mkdir(parents=True)
    (cfg_dir / "omnibenchmark.cfg").write_text(
        f"[conda]\ncache = {cache}\nchannel_alias = {ALIAS}\n"
    )
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.delenv("XDG_CONFIG_HOME", raising=False)
    for var in ("SNAKEMAKE_CONDA_PREFIX", "CONDA_CHANNEL_ALIAS"):
        monkeypatch.delenv(var, raising=False)

    plan = tmp_path / "plan"
    (plan / "envs").mkdir(parents=True)
    (plan / "envs" / "py.yaml").write_text(ENV)
    (plan / "benchmark.yaml").write_text(BENCHMARK)
    bundles = Path(__file__).parent.parent / "data" / "bundles"
    (plan / "bundles").symlink_to(bundles, target_is_directory=True)

    def run(out: str) -> str:
        with OmniCLISetup() as omni:
            r = omni.call(["run", "benchmark.yaml", "--out-dir", out], cwd=str(plan))
        assert r.returncode == 0, r.stdout + r.stderr
        assert not (plan / out / ".snakemake" / "conda").exists()
        return "".join(
            p.read_text() for p in (plan / out / ".logs").glob("snakemake_*.log")
        )

    assert "Creating conda environment" in run("out1")
    assert "Creating conda environment" not in run("out2")

    (env_dir,) = [p for p in cache.iterdir() if p.is_dir()]
    (python,) = (env_dir / "conda-meta").glob("python-3*.json")
    assert json.loads(python.read_text())["url"].startswith(f"{ALIAS}/conda-forge/")
