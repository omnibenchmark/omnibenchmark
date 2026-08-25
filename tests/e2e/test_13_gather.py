"""Grouped gather, end to end (design 010 §6 "Testing").

Everything else that covers `gather:` builds `Stage`/`GatherSpec` objects in
Python, so nothing exercised the path a user actually takes: YAML -> parsed
model -> plan -> executed workflow. That path has its own failure modes — the
`from` alias (a Python keyword, so it only works via the Pydantic alias), the
api gate, the topology edges a gather stage needs despite declaring no
`inputs:`, and the list-valued flag the module receives.

The fixture is the shape §5 asks for: two method stages sharing the
`methods.result` output id, a gather over it grouped by `data`, and a
downstream stage chaining off the gather (scatter -> gather -> scatter).
"""

import json
import re
from pathlib import Path

import pytest

from tests.e2e.common import E2ETestRunner


@pytest.fixture
def gather_config():
    """Path to the grouped-gather config."""
    return Path(__file__).parent / "configs" / "13_gather.yaml"


def _run(config, tmp_path, keep_files, extra=None):
    runner = E2ETestRunner(tmp_path, keep_files)
    config_file = runner.setup_test_environment(config, "13_gather.yaml")
    runner.execute_cli_command(config_file, extra or [])
    return runner


def _rule_blocks(snakefile: str) -> dict:
    """Map rule name -> rule body, excluding `rule all`."""
    blocks = {}
    for block in re.split(r"^rule ", snakefile, flags=re.M)[1:]:
        name, _, body = block.partition(":")
        if name != "all":
            blocks[name] = body
    return blocks


@pytest.mark.e2e
def test_gather_plan_groups_producers_from_both_method_stages(
    gather_config, tmp_path, bundled_repos, keep_files
):
    """One gather node per data module, each collecting both method stages.

    Four producers of `methods.result` (2 data modules x 2 method stages)
    partition into two groups of two. Asserting the *input count* is the point:
    membership is §3.1 rule 2 — every producer of the id, from any stage — so a
    regression that kept only one producer per output id would still emit
    plausible-looking gather nodes, just with one member each.
    """
    runner = _run(gather_config, tmp_path, keep_files, ["--dry"])
    snakefile = (runner.out_dir / "Snakefile").read_text()

    # `metrics_MC_<group>_default` exactly — downstream report rules prefix-compose
    # off the gather's id, so a substring match would also catch those.
    gather_rules = {
        name: body
        for name, body in _rule_blocks(snakefile).items()
        if re.fullmatch(r"metrics_MC_[^_]+_default", name)
    }
    assert (
        len(gather_rules) == 2
    ), f"expected one gather node per data module, got {sorted(gather_rules)}"

    for name, body in gather_rules.items():
        inputs = re.findall(r'input_\d+="([^"]+)"', body)
        assert len(inputs) == 2, f"{name} gathered {inputs}, expected both methods"
        assert {"method_a/MA" in p for p in inputs} == {
            True,
            False,
        }, f"{name} must draw one member from each method stage: {inputs}"

    # `prefix:` + the group value start the cut chain; the group is a data
    # module id, never a node id (§3.2 "grouping merges parameter expansions").
    outputs = re.findall(r'"(aggregated/[^"]+)"', snakefile)
    groups = {p.split("/")[1] for p in outputs}
    assert groups == {"D1", "D2"}, f"unexpected group segments: {groups}"


@pytest.mark.e2e
def test_downstream_stage_chains_off_the_gather(
    gather_config, tmp_path, bundled_repos, keep_files
):
    """Scatter after gather: `report` consumes `metrics.summary` and extends
    the gather's directory, which is what makes a gather an ordinary node
    (§2 "Gather nodes are ordinary nodes") rather than a terminal collector."""
    runner = _run(gather_config, tmp_path, keep_files, ["--dry"])
    snakefile = (runner.out_dir / "Snakefile").read_text()

    # Deduplicated: every output also appears in `rule all`'s input list.
    report_outputs = set(
        re.findall(r'"(aggregated/\S*report/R1/\S*_r\.json)"', snakefile)
    )
    assert (
        len(report_outputs) == 2
    ), f"expected one report per gather node, got {sorted(report_outputs)}"
    for path in report_outputs:
        assert path.startswith(
            "aggregated/D"
        ), f"report must extend the gather's directory, got {path}"


@pytest.mark.e2e
def test_gather_runs_and_records_its_members(
    gather_config, tmp_path, bundled_repos, keep_files
):
    """The executed workflow, which is the part no unit test reaches.

    Proves the emitted shell actually works: the collector receives its members
    as one list-valued `--methods.result` flag and writes into `output_dir`.
    The aggregate it produces is checked against what the two method modules
    computed, so this fails if the gather collected the wrong members rather
    than merely the wrong number of them. Also checks the `lineage.json`
    sidecar, the only record of what a gather collected once the chain is cut
    (§3.3).
    """
    runner = _run(gather_config, tmp_path, keep_files)

    # data D1 evaluates 1+1=2, D2 evaluates 2+2=4; MA is input+1 and MB input+2.
    expected = {"D1": [3.0, 4.0], "D2": [5.0, 6.0]}

    for group in ("D1", "D2"):
        base = runner.out_dir / "aggregated" / group / "metrics" / "MC" / ".default"
        summary = base / "metrics.json"
        assert summary.is_file(), f"gather produced no output for {group}"

        aggregate = json.loads(summary.read_text())
        assert aggregate["files_processed"] == 2, aggregate
        assert (
            sorted(aggregate["values"]) == expected[group]
        ), f"{group} aggregated the wrong members: {aggregate['values']}"

        lineage = json.loads((base / "lineage.json").read_text())
        assert lineage["kind"] == "gather"
        modules = sorted(m["module"] for m in lineage["members"])
        assert modules == [
            "MA",
            "MB",
        ], f"{group}'s sidecar must name every contributing node, got {modules}"
        assert {m["stage"] for m in lineage["members"]} == {"method_a", "method_b"}

    reports = list(runner.out_dir.glob("aggregated/*/**/report/R1/**/R1_r.json"))
    assert len(reports) == 2, f"expected one report per group, got {reports}"
