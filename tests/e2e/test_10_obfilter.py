import re
import pytest
from pathlib import Path

from omnibenchmark import filter as obfilter
from omnibenchmark.model.benchmark import Benchmark
from tests.e2e.common import E2ETestRunner

CONFIG = Path(__file__).parent / "configs" / "09_transitive_exclude.yaml"
# Topology: data[D1,D2] -> preprocessing[P1] -> methods[M1,M2].


def _blob(picks):
    """Pack a v3 blob anchored on the config's real summary_hash (so no drift)."""
    parent = {"sha256": Benchmark.from_yaml(CONFIG.read_text()).summary_hash()}
    return obfilter.pack_blob(picks, parent)


@pytest.fixture
def config():
    return CONFIG


@pytest.mark.e2e
def test_filter_prunes_unpicked_module(config, tmp_path, bundled_repos, keep_files):
    """--filter keeps M1 (picked) and drops M2 (unpicked); '*' keeps upstream stages."""
    runner = E2ETestRunner(tmp_path, keep_files)
    config_file = runner.setup_test_environment(config, "09_transitive_exclude.yaml")

    blob = _blob(
        {"data": {"*": "all"}, "preprocessing": {"*": "all"}, "methods": {"M1": "all"}}
    )
    runner.execute_cli_command(config_file, ["--dry", "--filter", blob])

    snakefile = (runner.out_dir / "Snakefile").read_text()
    assert re.search(r"methods/M1", snakefile), "picked M1 should be present"
    assert not re.search(r"methods/M2", snakefile), "unpicked M2 should be pruned"
    # upstream stars keep both datasets
    assert re.search(r"data/D1", snakefile) and re.search(r"data/D2", snakefile)


@pytest.mark.e2e
def test_filter_orphaned_pick_fails(config, tmp_path, bundled_repos, keep_files):
    """A pick naming a non-existent module hard-fails (no --allow-drift)."""
    from tests.e2e.common import OmniCLISetup

    runner = E2ETestRunner(tmp_path, keep_files)
    config_file = runner.setup_test_environment(config, "09_transitive_exclude.yaml")
    blob = _blob({"methods": {"M_ghost": "all"}})

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(config_file),
                "--out-dir",
                str(runner.out_dir),
                "--dry",
                "--filter",
                blob,
            ],
            cwd=str(runner.tmp_path),
        )
    assert result.returncode != 0
    assert "do not resolve" in (result.stdout + result.stderr)
