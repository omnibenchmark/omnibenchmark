"""E2E tests for `Stage.provides` → `Module.requires` lineage gating.

Uses --dry to assert on the generated Snakefile rather than running
modules; this avoids a quirk where the legacy dummy bundle's filename
convention (--name = module_id) does not match path templates under
api_version >= 0.5 for non-data stages. The lineage-gate behavior we
care about (which rules are created vs pruned) is fully observable in
the generated Snakefile.
"""

from pathlib import Path

import pytest


@pytest.fixture
def lineage_provides_config():
    """Path to the lineage-gating benchmark config."""
    return Path(__file__).parent / "configs" / "11_lineage_provides.yaml"


def _read_snakefile(out_dir: Path) -> str:
    snakefile = out_dir / "Snakefile"
    assert snakefile.exists(), f"Snakefile not generated at {snakefile}"
    return snakefile.read_text()


@pytest.mark.e2e
def test_requires_label_match_keeps_only_matching_lineage(
    lineage_provides_config, tmp_path, bundled_repos, keep_files
):
    """`requires: {dataset_size: lg}` keeps only the `huge` execution path."""
    from tests.e2e.common import E2ETestRunner

    runner = E2ETestRunner(tmp_path, keep_files)
    config_in_tmp = runner.setup_test_environment(
        lineage_provides_config, "11_lineage_provides.yaml"
    )

    # --dry: generate Snakefile, do not execute.
    runner.execute_cli_command(config_in_tmp, ["--dry"])

    snakefile = _read_snakefile(runner.out_dir)

    # Both data nodes must be present (the gate is downstream).
    assert "data/small/" in snakefile, "small data rule missing"
    assert "data/huge/" in snakefile, "huge data rule missing"

    # The methods rule for the matching lineage must be present.
    assert "methods/M_lg_only" in snakefile, "M_lg_only rule missing entirely"

    # Critical: the methods rule must run only off the `huge` branch.
    # Rule output paths nest the lineage; the small branch must not appear
    # under any methods rule.
    methods_lines = [
        line for line in snakefile.splitlines() if "methods/M_lg_only" in line
    ]
    assert methods_lines, "expected at least one methods rule line"
    huge_branch = any("huge" in line for line in methods_lines)
    small_branch = any("/small/" in line for line in methods_lines)
    assert huge_branch, f"M_lg_only should run on huge lineage; saw: {methods_lines}"
    assert (
        not small_branch
    ), f"M_lg_only must NOT run on small lineage; saw: {methods_lines}"
