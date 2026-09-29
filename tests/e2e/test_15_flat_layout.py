"""Flat output layout, end to end (007 §3.1).

One 0.8.0 plan with every segment kind: chain, join (parent-set digest),
grouped gather, and a scatter after the gather. Checks the tree the run
leaves behind, the `human/` view, `ob collect` over it, and that a second run
executes nothing.
"""

import csv
import json
import re
from pathlib import Path

import pytest
from click.testing import CliRunner

from tests.e2e.common import E2ETestRunner

_ID = r"[A-Za-z_][A-Za-z0-9_]*"
_SEGMENT = re.compile(
    rf"{_ID}\.{_ID}(\.{_ID})?\.(default|[0-9a-f]{{8}})(-[0-9a-f]{{8}})?"
)


@pytest.mark.e2e
def test_flat_layout(tmp_path, bundled_repos, keep_files):
    runner = E2ETestRunner(tmp_path, keep_files)
    config = Path(__file__).parent / "configs" / "15_flat_layout.yaml"
    config_file = runner.setup_test_environment(config, "15_flat_layout.yaml")
    runner.execute_cli_command(config_file)
    out = runner.out_dir

    # Every result directory is one node segment; nothing hidden, no links.
    for path in out.rglob("*"):
        rel = path.relative_to(out)
        if rel.parts[0].startswith(".") or rel.parts[0] in ("human", "Snakefile"):
            continue
        assert not path.is_symlink(), rel
        if path.is_dir():
            assert _SEGMENT.fullmatch(path.name), f"not a node segment: {rel}"

    for group in ("D1", "D2"):
        (gather,) = out.glob(f"metrics.MC.{group}.*")
        assert (gather / "MC_g.json").is_file()
        lineage = json.loads((gather / "lineage.json").read_text())
        assert [m["module"] for m in lineage["members"]] == ["J1"]
        (report,) = gather.glob("report.R1.*/R1_r.json")

        joins = list(out.glob(f"data.{group}.*/*/join.J1.*-*/J1_j.json"))
        assert len(joins) == 1, joins

        # The readable mirror resolves to the same files.
        (readable,) = (out / "human").glob(
            f"metrics.MC.{group}.evaluate-1+1_kind-g/report.R1.*/R1_r.json"
        )
        assert readable.is_symlink()
        assert readable.resolve() == report.resolve()

    # `ob collect` reads the flat tree and does not count `human/` twice.
    from omnibenchmark.cli.collect import collect

    result = CliRunner().invoke(collect, ["performance", "-o", str(out)])
    assert result.exit_code == 0, result.output
    with open(out / "performances.tsv", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    # 2 data + 2 shallow + 2 deep + 2 join + 2 report; gathers have no
    # benchmark file.
    assert len(rows) == 10, [r["path"] for r in rows]
    reports = [r for r in rows if r["stage"] == "report"]
    assert sorted(r["dataset"] for r in reports) == ["D1", "D2"]

    # Rebuilding `human/` happens outside Snakemake, so a rerun is a no-op.
    rerun = runner.execute_cli_command(config_file, debug_label="rerun")
    assert "Nothing to be done" in rerun.stdout + rerun.stderr or (
        "Completed 0 jobs" in rerun.stdout
    ), rerun.stdout
