"""Design 014: the node wrapper (ob_node.sh) and the Failure model."""

import json
import shutil
import subprocess
import time
from pathlib import Path
from types import SimpleNamespace

import pytest

from omnibenchmark.model import failure as failures
from omnibenchmark.model.failure import Failure, FailureKind

WRAPPER = Path(failures.__file__).parents[1] / "backend" / "ob_node.sh"


def run_node(tmp_path, cmd, upstream="", **env):
    """Run *cmd* under ob_run as a rule would; (exit code, failure or None)."""
    (tmp_path / "n").mkdir(exist_ok=True)
    out = tmp_path / "n" / "out.txt"
    script = (
        f"set -euo pipefail; source {WRAPPER}; "
        f'OB_NODE=n OB_UPSTREAM="{upstream}" OB_OUTPUTS="{out}" ob_run {cmd}'
    )
    rc = subprocess.run(
        ["bash", "-c", script], cwd=tmp_path, env={"PATH": "/usr/bin:/bin", **env}
    ).returncode
    return rc, Failure.load(Path("n"), tmp_path)


@pytest.mark.parametrize(
    "cmd, kind, reason",
    [
        ("false", "crashed", "exit code 1"),
        ("sh -c 'echo boom >&2; exit 3'", "crashed", "exit code 3: boom"),
        ("sh -c 'kill -9 $$'", "crashed", "killed by signal 9 (out of memory?)"),
        ("./missing.sh", "failed_to_start", "entrypoint or interpreter not found"),
        (
            "sh -c 'echo \"ModuleNotFoundError: No module named x\" >&2; exit 1'",
            "failed_to_start",
            "ModuleNotFoundError: No module named x",
        ),
        ("true", "missing_output", "declared output out.txt not written"),
    ],
)
def test_wrapper_classifies(tmp_path, cmd, kind, reason):
    rc, f = run_node(tmp_path, cmd)
    assert rc != 0 and f.kind == FailureKind(kind) and f.reason == reason
    assert (tmp_path / "n" / "out.txt").read_text() == ""  # emptied, not removed


def test_wrapper_success_leaves_no_record(tmp_path):
    rc, f = run_node(tmp_path, "sh -c 'echo ok > n/out.txt'")
    assert (rc, f) == (0, None)


def test_wrapper_times_out_and_continues(tmp_path):
    rc, f = run_node(tmp_path, "sleep 10", OB_TIMEOUT="1", OB_CONTINUE_ON_ERROR="1")
    assert rc == 0  # -k: Snakemake sees a finished job
    assert (f.kind, f.exit_code, f.timeout_s) == (FailureKind.TIMED_OUT, 124, 1)


@pytest.mark.skipif(shutil.which("setsid") is None, reason="needs setsid")
def test_wrapper_times_out_past_a_descendant_holding_stderr(tmp_path):
    # The grandchild leaves the process group, so the timeout's KILL misses
    # it, as it misses one stuck in the kernel; it keeps stderr open for 60 s.
    stray = "setsid sh -c 'echo $$ > stray; exec sleep 60'"
    t0 = time.monotonic()
    try:
        rc, f = run_node(tmp_path, f'sh -c "{stray} & sleep 30"', OB_TIMEOUT="1")
        assert time.monotonic() - t0 < 20
        assert f.kind == FailureKind.TIMED_OUT and rc == 124
    finally:
        pid = (tmp_path / "stray").read_text().strip()
        subprocess.run(["kill", pid])


def test_wrapper_skips_a_node_below_a_failure(tmp_path):
    (tmp_path / "up").mkdir()
    (tmp_path / "up" / "failure.json").write_text("{}")
    rc, f = run_node(tmp_path, "touch ran", upstream="up")
    assert f.kind == FailureKind.UPSTREAM and f.upstream == ("up",)
    assert f.exit_code is None and not (tmp_path / "ran").exists()


def test_wrapper_json_escapes_the_reason(tmp_path):
    _, f = run_node(tmp_path, """sh -c 'printf "a \\"q\\" \\\\ b\\n" >&2; exit 2'""")
    assert f.reason == 'exit code 2: a "q" \\ b'


def test_failure_round_trip_and_signal(tmp_path):
    f = Failure(FailureKind.CRASHED, "killed", exit_code=137, upstream=("a",))
    (tmp_path / "d").mkdir()
    (tmp_path / "d" / "failure.json").write_text(json.dumps(f.to_dict()))
    loaded = Failure.load(Path("d"), tmp_path)
    assert loaded == f and loaded.node_dir == "d" and loaded.signal == 9
    assert Failure.load(Path("none"), tmp_path) is None


@pytest.mark.parametrize(
    "text, seconds", [("90s", 90), ("30m", 1800), ("2h", 7200), ("1d", 86400)]
)
def test_parse_duration(text, seconds):
    assert failures.parse_duration(text) == seconds


@pytest.mark.parametrize("text", ["120", "0h", "2 hours", ""])
def test_parse_duration_needs_a_unit(text):
    with pytest.raises(ValueError):
        failures.parse_duration(text)


def test_collect_clear_and_summary(tmp_path):
    nodes = [
        SimpleNamespace(node_dir=d, outputs={"o": f"{d}/o.txt"})
        for d in ("a", "b", "c")
    ]
    for d, kind in (("a", "crashed"), ("b", "timed_out")):
        (tmp_path / d).mkdir()
        (tmp_path / d / "o.txt").write_text("")
        (tmp_path / d / "failure.json").write_text(
            json.dumps({"kind": kind, "reason": f"{d} failed"})
        )
    found = failures.collect(tmp_path, nodes)
    assert [f.node_dir for f in found] == ["a", "b"]
    text = failures.summary(found, 3)
    assert text.splitlines()[0] == "2 of 3 nodes failed" and "    a   a failed" in text

    assert failures.clear(tmp_path, nodes, {"timed_out"}) == 1
    assert not (tmp_path / "b" / "o.txt").exists()
    assert [f.node_dir for f in failures.collect(tmp_path, nodes)] == ["a"]


def test_generated_rules_run_modules_under_the_wrapper(tmp_path):
    from omnibenchmark.model.benchmark import APIVersion
    from tests.backend.test_snakemake_gen import _fanin_snakefile, _gen, _make_node

    text, gather = _fanin_snakefile(tmp_path, 3, APIVersion.V0_8_0)
    assert (tmp_path / ".ob_node.sh").exists()
    assert "source .ob_node.sh" in gather and 'OB_MANIFEST="$MANIFEST" ob_run' in gather
    # The manifest rule waits for its members and knows where their records go.
    rule = text[text.index("rule metrics_summ_default_manifest:") :]
    assert rule.index("input:") < rule.index("output:")
    assert '"_failure": {' in rule and 'record["failed"] = failed' in rule

    parent = _make_node(node_id="p", outputs={"r": "out/p/r.txt"})
    child = _make_node(
        node_id="c", inputs={"r": "out/p/r.txt"}, outputs={"o": "out/p/c/o.txt"}
    )
    path = tmp_path / "Snakefile"
    _gen(api_version=APIVersion.V0_8_0).generate_snakefile([parent, child], path)
    assert (
        'OB_NODE=out/p/c OB_UPSTREAM="out/p" OB_OUTPUTS="$OUTPUT_FILE_0" ob_run'
        in path.read_text()
    )

    _gen(api_version=APIVersion.V0_7_0).generate_snakefile([parent, child], path)
    assert (
        "ob_run" not in path.read_text()
        and "source .ob_node.sh" not in path.read_text()
    )
