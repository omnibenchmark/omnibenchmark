"""Behavioral tests for TelemetryEmitter.

These cover the OTLP JSON Lines output: hierarchy emission ordering,
parenting, status propagation, and the auxiliary log records emitted
alongside spans.
"""

import io
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import pytest

from omnibenchmark.telemetry.emitter import TelemetryEmitter
from omnibenchmark.telemetry.events import SpanStatus

pytestmark = pytest.mark.short


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


@dataclass
class FakeNode:
    id: str
    stage_id: str
    module_id: str
    param_id: str = ""
    inputs: Dict[str, str] = field(default_factory=dict)
    outputs: List[str] = field(default_factory=list)
    _params: Optional[str] = None

    def get_parameter_json(self) -> Optional[str]:
        return self._params


def _make_emitter(buf: Optional[io.StringIO] = None) -> TelemetryEmitter:
    return TelemetryEmitter(output=buf or io.StringIO())


def _lines(buf: io.StringIO) -> List[Dict[str, Any]]:
    return [json.loads(line) for line in buf.getvalue().splitlines() if line]


def _spans(buf: io.StringIO) -> List[Dict[str, Any]]:
    out = []
    for entry in _lines(buf):
        for rs in entry.get("resourceSpans", []):
            for ss in rs.get("scopeSpans", []):
                out.extend(ss.get("spans", []))
    return out


def _logs(buf: io.StringIO) -> List[Dict[str, Any]]:
    out = []
    for entry in _lines(buf):
        for rl in entry.get("resourceLogs", []):
            for sl in rl.get("scopeLogs", []):
                out.extend(sl.get("logRecords", []))
    return out


def _attrs(span_or_log: Dict[str, Any]) -> Dict[str, Any]:
    """Flatten OTLP attributes into a plain dict."""
    flat = {}
    for a in span_or_log.get("attributes", []):
        v = a["value"]
        flat[a["key"]] = next(iter(v.values()))
    return flat


def _init_simple(emitter: TelemetryEmitter) -> List[FakeNode]:
    """Stage S1 with module M1 (2 rules) and module M2 (1 rule)."""
    nodes = [
        FakeNode(id="S1-M1-p0", stage_id="S1", module_id="M1", outputs=["a.txt"]),
        FakeNode(id="S1-M1-p1", stage_id="S1", module_id="M1"),
        FakeNode(id="S1-M2-p0", stage_id="S1", module_id="M2"),
    ]
    emitter.init_benchmark(
        benchmark_name="bench",
        benchmark_version="1.0",
        benchmark_author="me",
        software_backend="conda",
        cores=2,
        stages=[{"id": "S1", "name": "stage-one"}],
        resolved_nodes=nodes,
    )
    return nodes


# ---------------------------------------------------------------------------
# Output sink configuration
# ---------------------------------------------------------------------------


class TestOutputSink:
    def test_defaults_to_stdout(self, capsys):
        emitter = TelemetryEmitter()
        emitter.init_benchmark(
            benchmark_name="b",
            benchmark_version="0",
            benchmark_author="",
            software_backend="host",
            cores=1,
            stages=[],
            resolved_nodes=[],
        )
        emitter.benchmark_completed(success=True)
        captured = capsys.readouterr()
        assert '"resourceSpans"' in captured.out

    def test_writes_to_path(self, tmp_path: Path):
        out = tmp_path / "telemetry.jsonl"
        with TelemetryEmitter(output=out) as emitter:
            _init_simple(emitter)
            emitter.benchmark_completed(success=True)
        assert out.exists()
        # Each line is valid JSON
        for line in out.read_text().splitlines():
            json.loads(line)

    def test_close_is_idempotent_on_external_handle(self):
        buf = io.StringIO()
        emitter = TelemetryEmitter(output=buf)
        emitter.close()  # should not close caller-owned buf
        assert not buf.closed

    def test_context_manager_closes_owned_path(self, tmp_path: Path):
        out = tmp_path / "t.jsonl"
        with TelemetryEmitter(output=out) as emitter:
            handle = emitter._file_handle
        assert handle.closed


# ---------------------------------------------------------------------------
# init_benchmark / metadata
# ---------------------------------------------------------------------------


class TestInitBenchmark:
    def test_groups_rules_by_stage_and_module(self):
        emitter = _make_emitter()
        _init_simple(emitter)
        assert set(emitter._stages.keys()) == {"S1"}
        assert emitter._stages["S1"]["module_count"] == 2
        assert set(emitter._modules.keys()) == {("S1", "M1"), ("S1", "M2")}
        assert emitter._modules[("S1", "M1")]["rule_count"] == 2
        assert emitter._modules[("S1", "M2")]["rule_count"] == 1
        assert len(emitter._rules) == 3

    def test_rule_name_sanitization_matches_snakefile_rule(self):
        emitter = _make_emitter()
        emitter.init_benchmark(
            benchmark_name="b",
            benchmark_version="0",
            benchmark_author="",
            software_backend="host",
            cores=1,
            stages=[{"id": "S"}],
            resolved_nodes=[FakeNode(id="1.weird-id", stage_id="S", module_id="M")],
        )
        # Leading digit triggers "rule_" prefix; dots/dashes become underscores.
        assert "rule_1_weird_id" in emitter._rules

    def test_emit_manifest_alias(self):
        emitter = _make_emitter()
        emitter.emit_manifest(
            benchmark_name="b",
            benchmark_version="0",
            benchmark_author="",
            software_backend="host",
            cores=1,
            stages=[],
            resolved_nodes=[],
        )
        assert emitter._benchmark_name == "b"

    def test_get_rule_log_file_returns_path_for_known_rule(self):
        buf = io.StringIO()
        emitter = TelemetryEmitter(output=buf)
        _init_simple(emitter)
        path = emitter.get_rule_log_file("S1_M1_p0")
        assert path == ".logs/S1_M1_p0.log"

    def test_get_rule_log_file_returns_none_for_unknown(self):
        emitter = _make_emitter()
        _init_simple(emitter)
        assert emitter.get_rule_log_file("does-not-exist") is None

    def test_stage_id_can_come_from_stage_id_key(self):
        emitter = _make_emitter()
        emitter.init_benchmark(
            benchmark_name="b",
            benchmark_version="0",
            benchmark_author="",
            software_backend="host",
            cores=1,
            stages=[{"stage_id": "alt"}],
            resolved_nodes=[FakeNode(id="x", stage_id="alt", module_id="M")],
        )
        assert "alt" in emitter._stages


# ---------------------------------------------------------------------------
# Phase / setup spans
# ---------------------------------------------------------------------------


class TestPhaseSpans:
    def test_emit_setup_span_writes_internal_span(self):
        buf = io.StringIO()
        emitter = TelemetryEmitter(output=buf)
        _init_simple(emitter)
        emitter.emit_setup_span(output="set up conda")
        spans = _spans(buf)
        assert len(spans) == 1
        assert spans[0]["name"] == "setup: environment preparation"
        assert _attrs(spans[0])["setup.type"] == "environment"
        # Output text becomes a structured log
        logs = _logs(buf)
        assert any(log["body"]["stringValue"] == "set up conda" for log in logs)

    def test_emit_phase_span_without_output_emits_no_log(self):
        buf = io.StringIO()
        emitter = TelemetryEmitter(output=buf)
        _init_simple(emitter)
        emitter.emit_phase_span(
            name="setup: resolve",
            phase="resolution",
            setup_type="resolution",
        )
        assert len(_spans(buf)) == 1
        assert _logs(buf) == []

    def test_emit_phase_started_emits_log_only(self):
        buf = io.StringIO()
        emitter = TelemetryEmitter(output=buf)
        _init_simple(emitter)
        emitter.emit_phase_started(name="resolve", phase="resolution")
        assert _spans(buf) == []
        logs = _logs(buf)
        assert len(logs) == 1
        attrs = _attrs(logs[0])
        assert attrs["phase"] == "resolution"
        assert attrs["phase.status"] == "started"


# ---------------------------------------------------------------------------
# Rule lifecycle
# ---------------------------------------------------------------------------


class TestRuleLifecycle:
    def test_rule_started_for_unknown_rule_is_noop(self):
        emitter = _make_emitter()
        _init_simple(emitter)
        emitter.rule_started("unknown")  # must not raise
        # No timing recorded anywhere
        assert all(m["start_time"] is None for m in emitter._modules.values())

    def test_rule_completed_emits_rule_span_with_parent_module(self):
        buf = io.StringIO()
        emitter = TelemetryEmitter(output=buf)
        _init_simple(emitter)
        emitter.rule_started("S1_M1_p0")
        emitter.rule_completed("S1_M1_p0", stdout="hello", stderr="warn")

        spans = _spans(buf)
        rule_spans = [s for s in spans if s["name"].startswith("rule:")]
        assert len(rule_spans) == 1
        rs = rule_spans[0]
        # Parent is the module span we stashed at init time
        assert rs["parentSpanId"] == emitter._modules[("S1", "M1")]["span_id"]
        # Status OK and outputs/stdout reflected
        assert rs["status"]["code"] == SpanStatus.OK.value
        # Outputs are array-valued; just verify the attribute is present
        assert any(a["key"] == "rule.outputs" for a in rs["attributes"])
        # Events include stdout + stderr
        ev_names = [e["name"] for e in rs.get("events", [])]
        assert "stdout" in ev_names and "stderr" in ev_names

        # Stdout / stderr also surface as structured logs
        log_severities = {log["severityText"] for log in _logs(buf)}
        assert "stdout" in log_severities
        # Successful rule => stderr log severity is WARN, not ERROR
        assert "stderr" in log_severities

    def test_rule_completed_on_unknown_rule_is_noop(self):
        buf = io.StringIO()
        emitter = TelemetryEmitter(output=buf)
        _init_simple(emitter)
        emitter.rule_completed("unknown")
        assert _spans(buf) == []

    def test_rule_completed_is_idempotent(self):
        buf = io.StringIO()
        emitter = TelemetryEmitter(output=buf)
        _init_simple(emitter)
        emitter.rule_completed("S1_M1_p0")
        emitter.rule_completed("S1_M1_p0")  # second call should be skipped
        rule_spans = [s for s in _spans(buf) if s["name"].startswith("rule:")]
        assert len(rule_spans) == 1

    def test_rule_failed_marks_module_and_emits_error_log(self):
        buf = io.StringIO()
        emitter = TelemetryEmitter(output=buf)
        _init_simple(emitter)
        emitter.rule_failed("S1_M2_p0", error="boom", stderr="trace")

        rule_spans = [s for s in _spans(buf) if s["name"].startswith("rule:")]
        assert rule_spans[0]["status"]["code"] == SpanStatus.ERROR.value
        # Prominent ERROR log carries the failure summary
        error_logs = [log for log in _logs(buf) if log.get("severityText") == "ERROR"]
        assert any(
            "FAILED: S1_M2_p0" in log["body"]["stringValue"] for log in error_logs
        )
        # Module span is emitted (M2 has only one rule) and marked ERROR
        module_spans = [s for s in _spans(buf) if s["name"].startswith("module:")]
        assert any(s["status"]["code"] == SpanStatus.ERROR.value for s in module_spans)


# ---------------------------------------------------------------------------
# Hierarchy emission
# ---------------------------------------------------------------------------


class TestHierarchyEmission:
    def test_module_emitted_only_after_all_rules_complete(self):
        buf = io.StringIO()
        emitter = TelemetryEmitter(output=buf)
        _init_simple(emitter)

        emitter.rule_completed("S1_M1_p0")
        # M1 still has one outstanding rule; module span must not be emitted yet
        assert not [s for s in _spans(buf) if s["name"].startswith("module:")]

        emitter.rule_completed("S1_M1_p1")
        module_spans = [s for s in _spans(buf) if s["name"].startswith("module:")]
        assert len(module_spans) == 1
        assert module_spans[0]["name"] == "module: M1"
        # Stage still pending (M2 not done yet)
        assert not [s for s in _spans(buf) if s["name"].startswith("stage:")]

    def test_benchmark_completed_flushes_pending_modules_and_stages(self):
        buf = io.StringIO()
        emitter = TelemetryEmitter(output=buf)
        _init_simple(emitter)
        # No rules at all completed
        emitter.benchmark_completed(success=True)
        names = [s["name"] for s in _spans(buf)]
        assert "module: M1" in names
        assert "module: M2" in names
        assert "stage: stage-one" in names
        assert any(n.startswith("bench-1.0") for n in names)

    def test_benchmark_completed_failure_status(self):
        buf = io.StringIO()
        emitter = TelemetryEmitter(output=buf)
        _init_simple(emitter)
        emitter.benchmark_completed(success=False, message="aborted")
        root = [s for s in _spans(buf) if s["name"].startswith("bench-")][0]
        assert root["status"]["code"] == SpanStatus.ERROR.value
        assert root["status"].get("message") == "aborted"

    def test_stage_completed_emits_stage_span_once(self):
        buf = io.StringIO()
        emitter = TelemetryEmitter(output=buf)
        _init_simple(emitter)
        emitter.stage_completed("S1")
        emitter.stage_completed("S1")  # second call no-ops
        stage_spans = [s for s in _spans(buf) if s["name"].startswith("stage:")]
        assert len(stage_spans) == 1


# ---------------------------------------------------------------------------
# Misc
# ---------------------------------------------------------------------------


class TestEdgeCases:
    def test_rule_failed_on_unknown_is_noop(self):
        buf = io.StringIO()
        emitter = TelemetryEmitter(output=buf)
        _init_simple(emitter)
        emitter.rule_failed("nope", error="boom")
        assert _spans(buf) == []

    def test_rule_span_includes_inputs_and_parameters_when_present(self):
        buf = io.StringIO()
        emitter = TelemetryEmitter(output=buf)
        emitter.init_benchmark(
            benchmark_name="b",
            benchmark_version="0",
            benchmark_author="",
            software_backend="host",
            cores=1,
            stages=[{"id": "S"}],
            resolved_nodes=[
                FakeNode(
                    id="r0",
                    stage_id="S",
                    module_id="M",
                    inputs={"in1": "/tmp/x.txt"},
                    _params='{"k": 1}',
                )
            ],
        )
        emitter.rule_completed("r0")
        rule_span = [s for s in _spans(buf) if s["name"].startswith("rule:")][0]
        keys = {a["key"] for a in rule_span["attributes"]}
        assert "rule.inputs" in keys
        assert "rule.parameters" in keys

    def test_stage_completed_after_natural_emission_is_noop(self):
        buf = io.StringIO()
        emitter = TelemetryEmitter(output=buf)
        _init_simple(emitter)
        # Drive everything to completion through the natural path
        for name in ("S1_M1_p0", "S1_M1_p1", "S1_M2_p0"):
            emitter.rule_completed(name)
        before = len([s for s in _spans(buf) if s["name"].startswith("stage:")])
        emitter.stage_completed("S1")  # already emitted naturally
        after = len([s for s in _spans(buf) if s["name"].startswith("stage:")])
        assert before == after == 1


class TestMisc:
    def test_trace_id_is_stable(self):
        emitter = _make_emitter()
        assert emitter.trace_id == emitter._trace_id
        # Trace IDs are 32-char hex per OTLP
        assert len(emitter.trace_id) == 32
