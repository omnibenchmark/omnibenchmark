"""Design 014, end to end: failed nodes are recorded, and with -k the rest of
the benchmark, gathers included, still runs.

One throwaway module repo, built here, whose `--mode` picks the outcome. Five
method nodes (ok, crash, sleep past the timeout, exit 0 without output, failed
import) each feed a `post` node, and a gather collects all five.
"""

import json
import subprocess
from pathlib import Path

import pytest

from omnibenchmark.model.failure import Failure
from tests.cli.cli_setup import OmniCLISetup

NODE_PY = """\
import json, sys, time
args = sys.argv[1:]
opt = lambda k: args[args.index(k) + 1] if k in args else None
mode = opt("--mode") or "ok"
outs = [a.split("=", 1)[1] for i, a in enumerate(args) if i and args[i - 1] == "--output"]
if mode == "crash":
    sys.exit("boom")
if mode == "sleep":
    time.sleep(60)
if mode == "import":
    import not_a_package
content = "ok"
if mode == "gather":
    m = json.load(open(opt("--from-files")))
    content = json.dumps({"members": len(m["members"]),
                          "failed": sorted(f["failure"]["kind"] for f in m.get("failed", []))})
if mode != "noout":
    for o in outs:
        open(o, "w").write(content)
"""

PLAN = """\
id: Failures
version: "1.0"
benchmarker: "E2E Test Suite"
api_version: "0.8.0"
software_backend: host
software_environments:
  host:
    description: "Host"
stages:
  - id: data
    modules:
      - {{id: D, software_environment: host, repository: {{url: {repo}, commit: {commit}}}}}
    outputs:
      - {{id: raw, path: "raw.txt"}}
  - id: methods
    inputs: [raw]
    modules:
      - id: M
        software_environment: host
        repository: {{url: {repo}, commit: {commit}}}
        parameters:
          - mode: [ok, crash, sleep, noout, import]
    outputs:
      - {{id: result, path: "result.txt"}}
  - id: post
    inputs: [result]
    modules:
      - {{id: P, software_environment: host, repository: {{url: {repo}, commit: {commit}}}}}
    outputs:
      - {{id: post, path: "post.txt"}}
  - id: summary
    gather: {{from: result}}
    modules:
      - id: S
        software_environment: host
        repository: {{url: {repo}, commit: {commit}}}
        parameters:
          - mode: gather
    outputs:
      - {{id: summary, path: "summary.json"}}
"""


@pytest.fixture
def plan(tmp_path):
    repo = tmp_path / "mod"
    repo.mkdir()
    (repo / "node.py").write_text(NODE_PY)
    (repo / "config.cfg").write_text("[DEFAULT]\nSCRIPT=node.py\n")
    git = ["git", "-c", "user.name=t", "-c", "user.email=t@t", "-C", str(repo)]
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    subprocess.run([*git, "add", "."], check=True)
    subprocess.run([*git, "commit", "-qm", "m"], check=True)
    commit = subprocess.run(
        [*git, "rev-parse", "--short", "HEAD"], capture_output=True, text=True
    ).stdout.strip()
    subprocess.run(
        [*git, "bundle", "create", "-q", "../mod.bundle", "HEAD"], check=True
    )
    path = tmp_path / "plan.yaml"
    path.write_text(PLAN.format(repo="mod.bundle", commit=commit))
    return path


def ob_run(plan, *args):
    with OmniCLISetup() as omni:
        out = plan.parent / "out"
        return omni.call(
            ["run", str(plan), "--out-dir", str(out), *args], cwd=str(plan.parent)
        )


def failures_by_mode(out: Path) -> dict:
    """{(stage, mode): kind} for every recorded failure."""
    found = {}
    for f in out.rglob("failure.json"):
        if "by-params" in f.parts:  # the readable view repeats every node
            continue
        stage = "post" if f.parent.name.startswith("post.") else "methods"
        method = f.parent.parent if stage == "post" else f.parent
        mode = json.loads((method / "parameters.json").read_text())["mode"]
        found[stage, mode] = Failure.load(f.parent).kind.value
    return found


@pytest.mark.e2e
def test_failures_are_recorded_and_the_gather_still_runs(plan):
    result = ob_run(plan, "-k", "--timeout", "3s", "--cores", "4")
    assert result.returncode == 1, result.stdout + result.stderr
    out = plan.parent / "out"

    kinds = failures_by_mode(out)
    assert {k: v for k, v in kinds.items() if k[0] == "methods"} == {
        ("methods", "crash"): "crashed",
        ("methods", "sleep"): "timed_out",
        ("methods", "noout"): "missing_output",
        ("methods", "import"): "failed_to_start",
    }
    # Each failed method's `post` child was not run.
    assert {k[1]: v for k, v in kinds.items() if k[0] == "post"} == dict.fromkeys(
        ["crash", "sleep", "noout", "import"], "upstream"
    )

    summary = json.loads(next(out.glob("summary.*/summary.json")).read_text())
    assert summary == {
        "members": 1,
        "failed": ["crashed", "failed_to_start", "missing_output", "timed_out"],
    }
    assert "8 of 12 nodes failed" in result.stdout + result.stderr


@pytest.mark.e2e
def test_without_k_the_run_stops_but_keeps_the_record(plan):
    result = ob_run(plan, "--timeout", "3s")
    assert result.returncode != 0
    assert not list((plan.parent / "out").rglob("summary.json"))
    assert list((plan.parent / "out").rglob("failure.json"))


@pytest.mark.e2e
def test_rerun_keeps_failures_and_retry_reruns_only_those(plan):
    ob_run(plan, "-k", "--timeout", "3s", "--cores", "4")
    out = plan.parent / "out"
    summary = next(out.glob("summary.*/summary.json"))
    before = summary.stat().st_mtime_ns

    # Another timeout and no -k change no rule: nothing reruns.
    again = ob_run(plan, "-k", "--timeout", "1h", "--cores", "4")
    assert again.returncode == 1 and summary.stat().st_mtime_ns == before

    ob_run(plan, "-k", "--timeout", "3s", "--retry-failed", "crashed")
    crash = [
        f
        for f in out.rglob("failure.json")
        if "crashed" in f.read_text() and "by-params" not in f.parts
    ]
    assert len(crash) == 1 and crash[0].stat().st_mtime_ns > before
    assert summary.stat().st_mtime_ns > before  # the gather above it reran
