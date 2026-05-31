import re
import pytest
from pathlib import Path

from tests.e2e.common import E2ETestRunner


@pytest.fixture
def transitive_exclude_config():
    """Path to the transitive (cross-stage) exclude config."""
    return Path(__file__).parent / "configs" / "09_transitive_exclude.yaml"


@pytest.mark.e2e
def test_transitive_exclude_prunes_across_stages(
    transitive_exclude_config, tmp_path, bundled_repos, keep_files
):
    """A dataset excluding a method must prune the combination even when the
    two modules are separated by an intervening stage.

    Topology: data -> preprocessing -> methods. ``D2`` declares ``exclude: [M2]``.
    Because M2 is not D2's immediate downstream neighbour, the pre-fix run loop
    (which only compared the direct predecessor) failed to prune it and emitted
    the degenerate D2 -> preprocessing -> M2 rule into the Snakefile.

    We assert on the generated Snakefile (``--dry``), mirroring the original
    bug repro, so the check is independent of module execution.
    """
    runner = E2ETestRunner(tmp_path, keep_files)
    config_file = runner.setup_test_environment(
        transitive_exclude_config, "09_transitive_exclude.yaml"
    )

    runner.execute_cli_command(config_file, ["--dry"])

    snakefile = (runner.out_dir / "Snakefile").read_text()

    # Every D2 path that reaches a method must route through M1, never M2.
    d2_method_paths = re.findall(r"\bdata/D2\S*methods/M\d", snakefile)
    assert d2_method_paths, "expected at least one D2->methods path in the Snakefile"
    assert all(
        p.endswith("methods/M1") for p in d2_method_paths
    ), f"D2->M2 should be excluded across the preprocessing stage, got: {d2_method_paths}"

    # Sanity: the legitimate D1->M2 combination is still generated.
    assert re.search(
        r"\bdata/D1\S*methods/M2", snakefile
    ), "D1->M2 should still be present"
