"""Snapshots, end to end (design 017).

A run is packed up to a stage with `ob archive --until`. A plan that differs
after that stage is then run with `--use-snapshots` into a fresh out dir: the
nodes of the snapshot are taken from it, the rest is computed.
"""

import os
import re
import shutil
from pathlib import Path

import pytest

from tests.cli.cli_setup import OmniCLISetup
from tests.e2e.common import E2ETestRunner

_FUSE = (
    os.path.exists("/dev/fuse")
    and shutil.which("mksquashfs")
    and (shutil.which("squashfuse_ll") or shutil.which("squashfuse"))
)


def _ob(tmp_path: Path, *args: str) -> str:
    with OmniCLISetup() as omni:
        run = omni.call(list(args), cwd=str(tmp_path))
    assert run.returncode == 0, run.stdout + run.stderr
    return run.stdout + run.stderr


@pytest.mark.e2e
@pytest.mark.parametrize(
    "pack_args",
    [
        ["--until", "join"],
        pytest.param(
            ["--until", "join", "--squashfs"],
            marks=pytest.mark.skipif(not _FUSE, reason="needs fuse and squashfs"),
        ),
    ],
)
def test_run_from_snapshot(tmp_path, bundled_repos, keep_files, pack_args):
    runner = E2ETestRunner(tmp_path, keep_files)
    config = Path(__file__).parent / "configs" / "15_flat_layout.yaml"
    plan = runner.setup_test_environment(config, "plan.yaml")
    runner.execute_cli_command(plan)

    # A base and what lies above it, in two files.
    _ob(tmp_path, "archive", "plan.yaml", "--until", "data", *pack_args[2:])
    _ob(tmp_path, "archive", "plan.yaml", *pack_args, "--base-until", "data")
    assert len(list((tmp_path / "snapshots").iterdir())) == 2

    # The metrics stage changes: everything up to `join` keeps its keys.
    changed = tmp_path / "changed.yaml"
    text = plan.read_text()
    old = '["--evaluate", "1+1", "--kind", "g"]'
    assert text.count(old) == 1
    changed.write_text(text.replace(old, '["--evaluate", "1+2", "--kind", "g"]'))

    out2 = tmp_path / "out2"
    try:
        log = _ob(
            tmp_path, "run", "changed.yaml", "--out-dir", "out2", "--use-snapshots"
        )
        assert "snapshots: none" not in log, log
        ran = "".join(p.read_text() for p in (out2 / ".logs").glob("snakemake_*.log"))
        # Only the metrics nodes and what follows them ran.
        rules = set(re.findall(r"^(?:local)?rule (\S+):", ran, re.M)) - {"all"}
        assert rules, ran
        assert all(
            "_metrics_MC_" in r or r.startswith("summary_") for r in rules
        ), rules

        for first in runner.out_dir.glob("data.*/shallow.*/join.*/performance.txt"):
            second = out2 / first.relative_to(runner.out_dir)
            assert second.read_text() == first.read_text()
            assert second.is_symlink() == ("--squashfs" in pack_args)
        assert len(list(out2.glob("data.*/metrics.MC.*/report.R1.*/R1_r.json"))) == 2

        # A lineage.json from a snapshot names its files relative to out/; one
        # written by this run names them by absolute path.
        (packed,) = out2.glob("data.D1.*/shallow.*/join.*/lineage.json")
        (fresh,) = out2.glob("data.D1.*/metrics.MC.*/lineage.json")
        assert str(tmp_path) not in packed.read_text()
        assert "data.D1." in packed.read_text()
        assert str(out2) in fresh.read_text()
    finally:
        _ob(tmp_path, "slice", "unmount", "--out-dir", "out2")


@pytest.mark.e2e
def test_archive_refuses_nodes_another_plan_produced(
    tmp_path, bundled_repos, keep_files
):
    """Each rule leaves its node key in its directory. Resolving a changed plan
    into the same out dir must not let its new keys label the old files."""
    runner = E2ETestRunner(tmp_path, keep_files)
    config = Path(__file__).parent / "configs" / "15_flat_layout.yaml"
    plan = runner.setup_test_environment(config, "plan.yaml")
    runner.execute_cli_command(plan)

    old = '"input+1"'  # the shallow stage
    assert plan.read_text().count(old) == 1
    changed = tmp_path / "changed.yaml"
    changed.write_text(plan.read_text().replace(old, '"input+3"'))
    _ob(tmp_path, "run", "changed.yaml", "--dry")

    with OmniCLISetup() as omni:
        stale = omni.call(["archive", "changed.yaml", "--snapshot"], cwd=str(tmp_path))
    assert stale.returncode != 0
    assert "not produced by the plan" in stale.stdout + stale.stderr
    # The data stage is upstream of the change: its keys still hold.
    _ob(tmp_path, "archive", "changed.yaml", "--until", "data")
