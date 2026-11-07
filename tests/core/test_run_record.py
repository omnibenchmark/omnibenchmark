"""Short tests for omnibenchmark.core.run_record and its use by `ob run`."""

import json
import shutil
from pathlib import Path
from unittest.mock import patch

import pytest

from omnibenchmark.cli.run import _run_benchmark
from omnibenchmark.core import run_record as rr
from omnibenchmark.model import Benchmark

DATA = Path(__file__).parents[1] / "data"
PLAN = DATA / "Benchmark_001_trimmed.yaml"


@pytest.fixture
def model():
    return Benchmark.from_yaml(PLAN)


def _record(out: Path, model, selection=None, **bench):
    """A run line as `ob run` writes it, with *bench* overriding the identity."""
    rr.write_run_record(out, model, PLAN, ["n1", "n2"], selection)
    md = out / ".metadata"
    lines = (md / rr.RUNLOG).read_text().splitlines()
    last = json.loads(lines[-1])
    last["benchmark"].update(bench)
    lines[-1] = json.dumps(last)
    (md / rr.RUNLOG).write_text("\n".join(lines) + "\n")


@pytest.mark.short
def test_write_run_record(tmp_path, model):
    (tmp_path / ".metadata").mkdir()
    (tmp_path / ".metadata" / "manifest.json").write_text('{"run_id": "r1"}')
    blob = {"picks": {"s": ["m"]}, "parent": {"sha256": "abc"}}
    sel = rr.selection_of(None, None, blob)
    line = rr.write_run_record(tmp_path, model, PLAN, ["b", "a"], sel, blob, ["gpu"])

    md = tmp_path / ".metadata"
    assert line["run_id"] == "r1"
    assert line["benchmark"] == rr.benchmark_record(model)
    assert (md / line["plan"]).read_text() == PLAN.read_text()
    assert json.loads((md / line["nodes"]).read_text()) == ["a", "b"]
    assert sel["filter"]["parent"] == "abc"
    filt = md / f"filter-{sel['filter']['sha256'][:8]}.json"
    assert json.loads(filt.read_text()) == blob
    assert line["capabilities"] == ["gpu"]
    assert rr.read_runs(tmp_path) == [line]


@pytest.mark.short
def test_selection_of():
    assert rr.selection_of(None, None, None) is None
    assert rr.selection_of("M1", None, None) == {"module": "M1"}
    assert rr.selection_of(None, "PCA", None) == {"until": "PCA"}


@pytest.mark.short
def test_nothing_recorded_is_not_an_error(tmp_path, model):
    assert rr.recorded_benchmark(tmp_path) is None
    assert rr.identity_mismatch(model, tmp_path) is None
    assert rr.check_out_dir(model, tmp_path) == ([], [])


@pytest.mark.short
def test_same_plan_passes(tmp_path, model):
    _record(tmp_path, model)
    assert rr.check_out_dir(model, tmp_path) == ([], [])


@pytest.mark.short
@pytest.mark.parametrize("field", ["id", "version"])
def test_other_benchmark_is_refused(tmp_path, model, field):
    _record(tmp_path, model, **{field: "other"})
    msg = rr.identity_mismatch(model, tmp_path)
    assert f"{field} 'other'" in msg
    assert rr.check_out_dir(model, tmp_path) == ([msg], [])


@pytest.mark.short
def test_plan_edited_since_last_run(tmp_path, model):
    _record(tmp_path, model, summary_hash="f" * 64)
    errors, _ = rr.check_out_dir(model, tmp_path)
    assert "Run `ob run` first" in errors[0]


@pytest.mark.short
def test_earlier_revisions_and_partial_runs_are_noted(tmp_path, model):
    _record(tmp_path, model, summary_hash="f" * 64)
    _record(tmp_path, model, selection={"until": "data"})
    errors, notes = rr.check_out_dir(model, tmp_path)
    assert errors == []
    assert "earlier plan revisions (ffffffff)" in notes[0]
    assert "partial run of this plan: --until data" in notes[1]

    _record(tmp_path, model)  # a full run of this plan
    assert len(rr.check_out_dir(model, tmp_path)[1]) == 1


@pytest.mark.short
def test_out_dir_from_an_older_ob(tmp_path, model):
    (tmp_path / ".metadata").mkdir()
    shutil.copy(PLAN, tmp_path / ".metadata" / "benchmark.yaml")
    assert rr.check_out_dir(model, tmp_path) == ([], [])

    text = PLAN.read_text().replace('version: "1.0"', 'version: "1.1"')
    (tmp_path / ".metadata" / "benchmark.yaml").write_text(text)
    assert "version '1.1'" in rr.identity_mismatch(model, tmp_path)

    # A copy this ob cannot parse still gives its id and version.
    (tmp_path / ".metadata" / "benchmark.yaml").write_text("id: x\nversion: 2\n")
    assert rr.recorded_benchmark(tmp_path) == {
        "id": "x",
        "version": "2",
        "summary_hash": None,
    }


@pytest.mark.short
def test_describe_selection():
    assert rr._describe({"module": "M1"}) == "-m M1"
    assert rr._describe({"filter": {"sha256": "abcdef0123"}}) == "--filter abcdef01"


def _ob_run(out: Path, dry: bool, **kwargs):
    with (
        patch("omnibenchmark.cli.run.populate_git_cache"),
        patch("omnibenchmark.cli.run._generate_explicit_snakefile", return_value=[]),
        patch("omnibenchmark.cli.run._run_snakemake") as snakemake,
        pytest.raises(SystemExit) as exc,
    ):
        _run_benchmark(
            benchmark_path=str(PLAN),
            cores=1,
            dry=dry,
            continue_on_error=False,
            out_dir=str(out),
            debug=False,
            dirty=False,
            **kwargs,
        )
        raise SystemExit(0)  # a real run returns after snakemake
    return exc.value.code, snakemake


@pytest.mark.short
def test_dry_run_records_nothing(tmp_path):
    code, _ = _ob_run(tmp_path, dry=True)
    assert code == 0
    assert not (tmp_path / ".metadata").exists()


@pytest.mark.short
def test_ob_run_records_and_refuses_another_version(tmp_path, model):
    code, snakemake = _ob_run(tmp_path, dry=False)
    assert code == 0 and snakemake.called
    (line,) = rr.read_runs(tmp_path)
    manifest = json.loads((tmp_path / ".metadata" / "manifest.json").read_text())
    assert line["run_id"] == manifest["run_id"]
    assert manifest["benchmark"] == line["benchmark"] == rr.benchmark_record(model)

    _record(tmp_path, model, version="0.9")
    code, snakemake = _ob_run(tmp_path, dry=True)
    assert code == 1 and not snakemake.called
    assert len(rr.read_runs(tmp_path)) == 2


@pytest.mark.short
def test_run_id_is_the_trace_id(tmp_path):
    from omnibenchmark.telemetry import TelemetryEmitter

    made = []

    def spy(**kwargs):
        made.append(TelemetryEmitter(**kwargs))
        return made[-1]

    trace = tmp_path / "trace.jsonl"
    with patch("omnibenchmark.telemetry.TelemetryEmitter", side_effect=spy):
        code, _ = _ob_run(tmp_path / "out", dry=False, telemetry_output=str(trace))
    assert code == 0
    (line,) = rr.read_runs(tmp_path / "out")
    assert line["run_id"] == made[0].trace_id
    assert line["trace"] == str(trace.resolve())


@pytest.mark.short
def test_no_telemetry_no_trace(tmp_path):
    _ob_run(tmp_path, dry=False)
    (line,) = rr.read_runs(tmp_path)
    assert line["trace"] is None
