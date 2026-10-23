"""foreach outputs (design 012): a producer picks its samples at run time, each
gets its own downstream node, and a gather collects them."""

import subprocess
from pathlib import Path

import pytest

from tests.cli.cli_setup import OmniCLISetup
from tests.e2e.common import E2ETestRunner

CONFIG_FILENAME = "15_foreach.yaml"

# One script, three roles, told apart by the flags it receives.
_MODULE = r"""#!/usr/bin/env python3
import os, sys

args, i = {}, 1
while i < len(sys.argv):
    key, val = sys.argv[i][2:], []
    i += 1
    while i < len(sys.argv) and not sys.argv[i].startswith("--"):
        val.append(sys.argv[i])
        i += 1
    args.setdefault(key, []).extend(val)

out = args.get("output", [None])[0]
if out and "{sample}" in out:                       # producer
    pattern = out.split("=", 1)[1]
    here = os.path.dirname(os.path.abspath(__file__))
    drop = os.path.join(here, "drop.txt")
    skip = open(drop).read().split() if os.path.exists(drop) else []
    for s in [s for s in args["samples"][0].split(",") if s not in skip]:
        path = pattern.replace("{sample}", s)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        open(path, "w").write(s + "\n")
elif "counts" in args:                              # per sample
    text = open(args["counts"][0]).read()
    open(out.split("=", 1)[1], "w").write(text.upper())
else:                                               # gather
    parts = [open(p).read() for p in args["processed"]]
    open(os.path.join(args["output_dir"][0], "summary.txt"), "w").write("".join(parts))
"""


def _commit(repo: Path, files: dict) -> str:
    for name, text in files.items():
        (repo / name).write_text(text)
    git = ["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@t"]
    if not (repo / ".git").exists():
        subprocess.run(git[:3] + ["init", "-q"], check=True)
    subprocess.run(git + ["add", "."], check=True)
    subprocess.run(git + ["commit", "-qm", "init"], check=True)
    out = subprocess.run(
        git[:3] + ["rev-parse", "HEAD"], check=True, capture_output=True, text=True
    )
    return out.stdout.strip()


def _setup(tmp_path, keep_files):
    repo = tmp_path / "foreach_module"
    repo.mkdir()
    commit = _commit(
        repo, {"run.py": _MODULE, "config.cfg": "[DEFAULT]\nSCRIPT=run.py\n"}
    )
    (repo / "run.py").chmod(0o755)
    config = Path(__file__).parent / "configs" / CONFIG_FILENAME
    runner = E2ETestRunner(tmp_path, keep_files)
    config_file = runner.setup_test_environment(config, CONFIG_FILENAME)
    config_file.write_text(config_file.read_text().replace("HEAD_COMMIT", commit))
    return runner, config_file, repo, commit


def _summary(runner) -> str:
    found = list(runner.out_dir.glob("per_dataset/**/summary.txt"))
    assert len(found) == 1, found
    return found[0].read_text()


@pytest.mark.e2e
def test_foreach_fans_out_and_gathers(tmp_path, keep_files):
    runner, config_file, _repo, _commit_id = _setup(tmp_path, keep_files)
    runner.execute_cli_command(config_file)

    assert _summary(runner) == "S1\nS2\nS3\n"
    # `.*`: the hash directory, not its readable symlink sibling.
    processed = sorted(runner.out_dir.glob("data/D1/.*/*/process/**/processed.txt"))
    assert [p.read_text() for p in processed] == ["S1\n", "S2\n", "S3\n"]

    snakefile = (runner.out_dir / "Snakefile").read_text()
    assert "checkpoint data_D1" in snakefile
    assert "_ob_foreach(checkpoints.data_D1" in snakefile


@pytest.mark.e2e
def test_foreach_rerun_drops_stale_values(tmp_path, keep_files):
    runner, config_file, repo, first = _setup(tmp_path, keep_files)
    runner.execute_cli_command(config_file)
    assert _summary(runner) == "S1\nS2\nS3\n"

    # Same parameters, hence the same node directory; the new commit makes the
    # producer re-run and choose fewer values. s3 must not be globbed again.
    second = _commit(repo, {"drop.txt": "s3\n"})
    config_file.write_text(config_file.read_text().replace(first, second))
    runner.execute_cli_command(config_file, debug_label="second run")
    assert _summary(runner) == "S1\nS2\n"
    assert not list(runner.out_dir.glob("data/D1/.*/s3/counts.txt"))


@pytest.mark.e2e
def test_nested_foreach_is_rejected(tmp_path, keep_files):
    runner, config_file, _repo, _commit_id = _setup(tmp_path, keep_files)
    text = config_file.read_text().replace(
        '        path: "processed.txt"',
        '        path: "{lane}/processed.txt"\n        foreach: lane',
    )
    config_file.write_text(text)
    with OmniCLISetup() as omni:
        result = omni.call(
            ["run", str(config_file), "--out-dir", str(runner.out_dir)],
            cwd=str(tmp_path),
        )
    assert result.returncode != 0
    assert "nesting is not supported" in result.stdout + result.stderr
