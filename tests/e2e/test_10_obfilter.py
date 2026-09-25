import re
import pytest
from pathlib import Path

from omnibenchmark import filter as obfilter
from omnibenchmark.model.benchmark import Benchmark
from tests.e2e.common import E2ETestRunner

CONFIG = Path(__file__).parent / "configs" / "09_transitive_exclude.yaml"
# Topology: data[D1,D2] -> preprocessing[P1] -> methods[M1,M2].

# The example in docs/src/howto.md ("Write a filter by hand").
PICKS_YAML = """\
picks:
  data: {'*': all}             # every dataset
  preprocessing: {'*': all}
  methods: {M1: all}           # M1 only; M2 is pruned
"""


def _filter_arg(tmp_path, text, form):
    """The --filter argument for `text` (YAML picks), as a YAML file or a packed blob."""
    if form == "blob":
        # Anchored on the config's real summary_hash, so no drift.
        parent = {"sha256": Benchmark.from_yaml(CONFIG.read_text()).summary_hash()}
        return obfilter.pack_blob(obfilter.load_filter(text)["picks"], parent)
    path = tmp_path / "picks.obfilter.yaml"
    path.write_text(text)
    return str(path)


@pytest.mark.short
def test_howto_example_is_the_tested_one():
    howto = (Path(__file__).parents[2] / "docs" / "src" / "howto.md").read_text()
    block = re.search(r"```yaml\n# picks.obfilter.yaml\n(.*?)```", howto, re.S)
    assert block and block.group(1) == PICKS_YAML


@pytest.fixture
def config():
    return CONFIG


@pytest.mark.e2e
@pytest.mark.parametrize("form", ["yaml", "blob"])
def test_filter_prunes_unpicked_module(
    config, tmp_path, bundled_repos, keep_files, form
):
    """--filter runs M1 (picked), not M2 (unpicked); '*' keeps upstream stages."""
    runner = E2ETestRunner(tmp_path, keep_files)
    config_file = runner.setup_test_environment(config, "09_transitive_exclude.yaml")

    arg = _filter_arg(tmp_path, PICKS_YAML, form)
    runner.execute_cli_command(config_file, ["--filter", arg])

    out = runner.out_dir
    # upstream stars keep both datasets; only the picked M1 runs
    for d in ("D1", "D2"):
        assert list(
            out.glob(f"data/{d}/.*/preprocessing/P1/.*/methods/M1/.*/M1_method.json")
        )
    assert not list(out.glob("**/methods/M2")), "unpicked M2 should be pruned"
    assert (out / "metrics" / "metrics.json").exists()


@pytest.mark.e2e
def test_filter_orphaned_pick_fails(config, tmp_path, bundled_repos, keep_files):
    """A pick naming a non-existent module hard-fails (no --allow-drift)."""
    from tests.e2e.common import OmniCLISetup

    runner = E2ETestRunner(tmp_path, keep_files)
    config_file = runner.setup_test_environment(config, "09_transitive_exclude.yaml")
    arg = _filter_arg(tmp_path, "picks:\n  methods: {M_ghost: all}\n", "yaml")

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(config_file),
                "--out-dir",
                str(runner.out_dir),
                "--dry",
                "--filter",
                arg,
            ],
            cwd=str(runner.tmp_path),
        )
    assert result.returncode != 0
    assert "do not resolve" in (result.stdout + result.stderr)
