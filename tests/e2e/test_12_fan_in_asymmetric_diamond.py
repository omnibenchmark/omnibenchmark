"""Fan-in output paths must be unique per parent set.

A join's directory prefix (`base_path`) is the parent of its DEEPEST input —
one branch only. Two joins sharing that branch but differing in another parent
therefore resolved to the same output path, and Snakemake rejected the plan:

    AmbiguousRuleException: Rules join_J1_… and join_J1_… are ambiguous for
    the file …/deep/D1/.75d9b523/join/J1/.ac913dcd/R1_j.json

The fixture is an *asymmetric* diamond on purpose — see its description. The
`snakemake --dry-run` check below is the general guard: it fails on any
duplicate-output regression, not only on a revert of the segment change.
"""

import re
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest

from tests.e2e.common import E2ETestRunner


@pytest.fixture
def fan_in_diamond_config():
    """Path to the asymmetric fan-in diamond config."""
    return Path(__file__).parent / "configs" / "12_fan_in_asymmetric_diamond.yaml"


def _generate(config, tmp_path, keep_files):
    runner = E2ETestRunner(tmp_path, keep_files)
    config_file = runner.setup_test_environment(
        config, "12_fan_in_asymmetric_diamond.yaml"
    )
    runner.execute_cli_command(config_file, ["--dry"])
    return runner


def _declared_outputs(snakefile: str) -> list:
    """Every path listed in an `output:` block, excluding `rule all`."""
    outputs = []
    for block in re.split(r"^rule ", snakefile, flags=re.M)[1:]:
        if block.startswith("all:"):
            continue
        body = block.split("output:", 1)
        if len(body) < 2:
            continue
        # the output block runs until the next directive at rule-indent level
        section = re.split(r"\n    \w+:", body[1], maxsplit=1)[0]
        outputs.extend(re.findall(r'"([^"]+)"', section))
    return outputs


@pytest.mark.e2e
def test_no_two_rules_declare_the_same_output(
    fan_in_diamond_config, tmp_path, bundled_repos, keep_files
):
    """The direct assertion: no output path is claimed twice."""
    runner = _generate(fan_in_diamond_config, tmp_path, keep_files)
    snakefile = (runner.out_dir / "Snakefile").read_text()

    duplicates = [p for p, n in Counter(_declared_outputs(snakefile)).items() if n > 1]

    assert not duplicates, (
        "two rules declare the same output file; a fan-in node's path must "
        f"carry its parent-set digest. Duplicated: {duplicates}"
    )


@pytest.mark.e2e
def test_both_joins_are_emitted_with_distinct_paths(
    fan_in_diamond_config, tmp_path, bundled_repos, keep_files
):
    """Both bundles survive — the fix must disambiguate, not deduplicate."""
    runner = _generate(fan_in_diamond_config, tmp_path, keep_files)
    snakefile = (runner.out_dir / "Snakefile").read_text()

    join_outputs = [p for p in _declared_outputs(snakefile) if "/join/J1/" in p]

    assert (
        len(join_outputs) == 2
    ), f"expected one join per shallow module: {join_outputs}"
    assert len(set(join_outputs)) == 2, "join outputs collided"
    # each records the branch its path omits
    for rule_block in re.split(r"^rule ", snakefile, flags=re.M)[1:]:
        if "/join/J1/" in rule_block and "lineage.json" in rule_block:
            assert "shallow" in rule_block, (
                "the join's lineage.json must name the shallow branch, which "
                "its path prefix drops"
            )


@pytest.mark.e2e
def test_snakemake_accepts_the_generated_workflow(
    fan_in_diamond_config, tmp_path, bundled_repos, keep_files
):
    """Snakemake itself builds the DAG — the authoritative check.

    Asserting on generated text can only catch regressions we predicted; this
    catches any way of producing two rules for one file.
    """
    snakemake = shutil.which(
        "snakemake", path=str(Path(sys.executable).parent)
    ) or shutil.which("snakemake")
    if snakemake is None:
        pytest.skip("snakemake not on PATH")

    runner = _generate(fan_in_diamond_config, tmp_path, keep_files)
    result = subprocess.run(
        [snakemake, "--cores", "1", "--dry-run"],
        cwd=str(runner.out_dir),
        capture_output=True,
        text=True,
        timeout=300,
    )
    combined = result.stdout + result.stderr

    assert "AmbiguousRuleException" not in combined, combined[-2000:]
    assert result.returncode == 0, combined[-2000:]
