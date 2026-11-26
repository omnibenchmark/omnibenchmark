"""A gather grouped by a list of stages, with a `path` (010 §3.2.2).

`compare` groups by `[methods, data]`: the data node is the anchor, `prep`
merges, so each job pairs the `as_is` and `denoise` result of one method
parameter set. Its `path` nests one `methods.M1.<param>` directory per job.
Checks the tree, lineage.json `key`, the `by-params/` view and `ob collect`
below the gather.
"""

import csv
import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from tests.e2e.common import E2ETestRunner


@pytest.mark.e2e
def test_gather_group_by_list(tmp_path, bundled_repos, keep_files):
    runner = E2ETestRunner(tmp_path, keep_files)
    config = Path(__file__).parent / "configs" / "18_gather_group_by_list.yaml"
    config_file = runner.setup_test_environment(config, "18_gather_group_by_list.yaml")
    runner.execute_cli_command(config_file)
    out = runner.out_dir

    for data in ("D1", "D2"):
        jobs = sorted(out.glob(f"data.{data}.*/compare/methods.M1.*"))
        assert len(jobs) == 2, jobs
        for job in jobs:
            assert (job / "C_c.json").is_file()
            lineage = json.loads((job / "lineage.json").read_text())
            (key,) = lineage["key"]
            assert (key["stage"], key["module"]) == ("methods", "M1")
            assert job.name == f"methods.M1.{key['params']}"
            members = lineage["members"]
            assert sorted(m["lineage"][1]["module"] for m in members) == [
                "as_is",
                "denoise",
            ]
            assert {m["lineage"][-1]["params"] for m in members} == {key["params"]}
            assert {m["lineage"][0]["module"] for m in members} == {data}
            (report,) = job.glob("report.R1.*/R1_r.json")

            # The readable view names the key directory by the method's
            # parameters, not by the gather module's.
            readable = [
                p
                for p in (out / "by-params").glob(
                    f"data.{data}.*/compare/methods.M1.*/report.R1.*/R1_r.json"
                )
                if p.resolve() == report.resolve()
            ]
            assert len(readable) == 1
            assert "evaluate-input" in readable[0].parts[-3]

    from omnibenchmark.cli.collect import collect

    result = CliRunner().invoke(collect, ["performance", "-o", str(out)])
    assert result.exit_code == 0, result.output
    with open(out / "performances.tsv", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    reports = [r for r in rows if r["stage"] == "report"]
    assert len(reports) == 4
    # The key directory is read as the gather, not as a `methods` node.
    assert {r["lineage"].split("/", 2)[-1] for r in reports} == {"compare/C/report/R1"}
    assert sorted(r["dataset"] for r in reports) == ["D1", "D1", "D2", "D2"]
