"""Behavioral tests for SpanBuilder.

Focus on the observable outcomes — span hierarchy, attribute presence, and
status transitions — rather than internal storage details.
"""

import pytest

from omnibenchmark.telemetry.events import SpanStatus
from omnibenchmark.telemetry.spans import SpanBuilder

pytestmark = pytest.mark.short


def _attrs(span) -> dict:
    """Flatten a span's attributes to {key: value} for easy assertions."""
    return {a.key: a.value for a in span.attributes}


def _make_builder_with_root() -> SpanBuilder:
    b = SpanBuilder()
    b.create_benchmark_span(
        benchmark_name="bench",
        benchmark_version="1.0",
        benchmark_author="me",
        total_rules=3,
        software_backend="conda",
        cores=4,
    )
    return b


# ---------------------------------------------------------------------------
# Root benchmark span
# ---------------------------------------------------------------------------


class TestBenchmarkSpan:
    def test_root_has_no_parent_and_carries_metadata(self):
        b = SpanBuilder()
        root = b.create_benchmark_span(
            benchmark_name="bench",
            benchmark_version="1.0",
            benchmark_author="me",
            total_rules=3,
            software_backend="conda",
            cores=4,
        )

        assert root.parent_span_id is None
        assert root.trace_id == b.trace_id
        assert root.start_time_ns is not None
        assert root.status == SpanStatus.UNSET
        assert "bench" in root.name

        attrs = _attrs(root)
        assert attrs["benchmark.name"] == "bench"
        assert attrs["benchmark.version"] == "1.0"
        assert attrs["benchmark.author"] == "me"
        assert attrs["benchmark.total_rules"] == 3
        assert attrs["benchmark.software_backend"] == "conda"
        assert attrs["benchmark.cores"] == 4
        assert attrs["service.name"] == "omnibenchmark"

    def test_custom_service_name_is_used(self):
        b = SpanBuilder(service_name="my-svc")
        root = b.create_benchmark_span("b", "1", "a", 0, "host", 1)
        assert _attrs(root)["service.name"] == "my-svc"


# ---------------------------------------------------------------------------
# Stage spans
# ---------------------------------------------------------------------------


class TestStageSpan:
    def test_stage_requires_root(self):
        b = SpanBuilder()
        with pytest.raises(RuntimeError):
            b.create_stage_span("s1")

    def test_stage_is_child_of_root(self):
        b = _make_builder_with_root()
        stage = b.create_stage_span("s1", stage_name="Stage One")

        assert stage.parent_span_id == b.get_all_spans()[0].span_id
        assert stage.trace_id == b.trace_id
        assert "Stage One" in stage.name
        assert stage.status == SpanStatus.UNSET
        # Stage start is set lazily when first rule starts.
        assert stage.start_time_ns is None

    def test_stage_name_defaults_to_id(self):
        b = _make_builder_with_root()
        stage = b.create_stage_span("s1")
        assert "s1" in stage.name
        assert _attrs(stage)["stage.name"] == "s1"

    def test_stage_attributes_record_counts(self):
        b = _make_builder_with_root()
        stage = b.create_stage_span("s1", module_count=2, rule_count=5)
        attrs = _attrs(stage)
        assert attrs["stage.id"] == "s1"
        assert attrs["stage.module_count"] == 2
        assert attrs["stage.rule_count"] == 5

    def test_get_stage_span_round_trip(self):
        b = _make_builder_with_root()
        stage = b.create_stage_span("s1")
        assert b.get_stage_span("s1") is stage
        assert b.get_stage_span("missing") is None


# ---------------------------------------------------------------------------
# Rule spans
# ---------------------------------------------------------------------------


class TestRuleSpan:
    def test_rule_is_child_of_stage(self):
        b = _make_builder_with_root()
        stage = b.create_stage_span("s1")
        rule = b.create_rule_span(
            rule_name="r1",
            stage_id="s1",
            module_id="m1",
            node_id="n1",
            inputs={"in": "data.tsv"},
            outputs=["out.tsv"],
        )
        assert rule.parent_span_id == stage.span_id
        assert rule.trace_id == b.trace_id

    def test_rule_falls_back_to_root_when_stage_missing(self):
        b = _make_builder_with_root()
        root_id = b.get_all_spans()[0].span_id
        rule = b.create_rule_span(
            rule_name="r1",
            stage_id="unknown-stage",
            module_id="m1",
            node_id="n1",
            inputs={},
            outputs=["out.tsv"],
        )
        assert rule.parent_span_id == root_id

    def test_rule_without_root_raises(self):
        b = SpanBuilder()
        with pytest.raises(RuntimeError):
            b.create_rule_span("r1", "s1", "m1", "n1", {}, [])

    def test_rule_attributes_include_core_metadata(self):
        b = _make_builder_with_root()
        b.create_stage_span("s1")
        rule = b.create_rule_span(
            rule_name="r1",
            stage_id="s1",
            module_id="m1",
            node_id="n1",
            inputs={"in": "data.tsv"},
            outputs=["out.tsv"],
        )
        attrs = _attrs(rule)
        assert attrs["rule.name"] == "r1"
        assert attrs["rule.stage_id"] == "s1"
        assert attrs["rule.module_id"] == "m1"
        assert attrs["rule.node_id"] == "n1"
        assert attrs["rule.outputs"] == ["out.tsv"]

    def test_inputs_attribute_omitted_when_empty(self):
        b = _make_builder_with_root()
        b.create_stage_span("s1")
        rule = b.create_rule_span("r1", "s1", "m1", "n1", {}, ["out.tsv"])
        assert "rule.inputs" not in _attrs(rule)

    def test_inputs_attribute_records_values_only(self):
        b = _make_builder_with_root()
        b.create_stage_span("s1")
        rule = b.create_rule_span(
            "r1",
            "s1",
            "m1",
            "n1",
            inputs={"a": "x.tsv", "b": "y.tsv"},
            outputs=["out.tsv"],
        )
        # Behavior: inputs values are exposed; we don't care about order.
        assert set(_attrs(rule)["rule.inputs"]) == {"x.tsv", "y.tsv"}

    def test_parameters_attribute_only_when_supplied(self):
        b = _make_builder_with_root()
        b.create_stage_span("s1")
        without = b.create_rule_span("r1", "s1", "m1", "n1", {}, [])
        assert "rule.parameters" not in _attrs(without)

        b2 = _make_builder_with_root()
        b2.create_stage_span("s1")
        with_params = b2.create_rule_span(
            "r2", "s1", "m1", "n1", {}, [], parameters="k=1"
        )
        assert _attrs(with_params)["rule.parameters"] == "k=1"

    def test_get_rule_span_round_trip(self):
        b = _make_builder_with_root()
        b.create_stage_span("s1")
        rule = b.create_rule_span("r1", "s1", "m1", "n1", {}, [])
        assert b.get_rule_span("r1") is rule
        assert b.get_rule_span("missing") is None


# ---------------------------------------------------------------------------
# Hierarchy view
# ---------------------------------------------------------------------------


class TestGetAllSpans:
    def test_empty_builder_returns_no_spans(self):
        assert SpanBuilder().get_all_spans() == []

    def test_hierarchy_order_is_root_then_stages_then_rules(self):
        b = _make_builder_with_root()
        s1 = b.create_stage_span("s1")
        s2 = b.create_stage_span("s2")
        r1 = b.create_rule_span("r1", "s1", "m1", "n1", {}, [])
        r2 = b.create_rule_span("r2", "s2", "m1", "n2", {}, [])

        spans = b.get_all_spans()
        assert len(spans) == 5
        root = spans[0]
        assert root.parent_span_id is None
        # The two stages come before the two rules.
        stage_ids = {s.span_id for s in spans[1:3]}
        rule_ids = {s.span_id for s in spans[3:]}
        assert stage_ids == {s1.span_id, s2.span_id}
        assert rule_ids == {r1.span_id, r2.span_id}


# ---------------------------------------------------------------------------
# Lifecycle transitions
# ---------------------------------------------------------------------------


class TestRuleLifecycle:
    def test_started_sets_rule_start_and_propagates_to_stage(self):
        b = _make_builder_with_root()
        stage = b.create_stage_span("s1")
        rule = b.create_rule_span("r1", "s1", "m1", "n1", {}, [])
        assert rule.start_time_ns is None
        assert stage.start_time_ns is None

        b.mark_rule_started("r1")
        assert rule.start_time_ns is not None
        assert stage.start_time_ns is not None

    def test_started_does_not_overwrite_stage_start(self):
        b = _make_builder_with_root()
        stage = b.create_stage_span("s1")
        b.create_rule_span("r1", "s1", "m1", "n1", {}, [])
        b.create_rule_span("r2", "s1", "m1", "n2", {}, [])

        b.mark_rule_started("r1")
        first_start = stage.start_time_ns

        b.mark_rule_started("r2")
        assert stage.start_time_ns == first_start

    def test_started_on_unknown_rule_returns_none(self):
        b = _make_builder_with_root()
        assert b.mark_rule_started("nope") is None

    def test_completed_sets_ok_status_and_end_time(self):
        b = _make_builder_with_root()
        b.create_stage_span("s1")
        rule = b.create_rule_span("r1", "s1", "m1", "n1", {}, [])
        b.mark_rule_started("r1")

        b.mark_rule_completed("r1")
        assert rule.status == SpanStatus.OK
        assert rule.end_time_ns is not None
        assert rule.events == []

    def test_completed_records_optional_stdout_and_stderr(self):
        b = _make_builder_with_root()
        b.create_stage_span("s1")
        rule = b.create_rule_span("r1", "s1", "m1", "n1", {}, [])

        b.mark_rule_completed("r1", stdout="hi", stderr="warn")
        event_names = [e.name for e in rule.events]
        assert event_names == ["stdout", "stderr"]
        outputs = {e.name: e.attributes[0].value for e in rule.events}
        assert outputs["stdout"] == "hi"
        assert outputs["stderr"] == "warn"

    def test_completed_skips_empty_streams(self):
        b = _make_builder_with_root()
        b.create_stage_span("s1")
        rule = b.create_rule_span("r1", "s1", "m1", "n1", {}, [])
        b.mark_rule_completed("r1", stdout="", stderr=None)
        assert rule.events == []

    def test_completed_on_unknown_rule_returns_none(self):
        b = _make_builder_with_root()
        assert b.mark_rule_completed("nope") is None

    def test_failed_sets_error_status_and_exception_event(self):
        b = _make_builder_with_root()
        b.create_stage_span("s1")
        rule = b.create_rule_span("r1", "s1", "m1", "n1", {}, [])

        b.mark_rule_failed("r1", error="boom")
        assert rule.status == SpanStatus.ERROR
        assert rule.status_message == "boom"
        assert rule.end_time_ns is not None
        # An exception event must be recorded with the message.
        exc_events = [e for e in rule.events if e.name == "exception"]
        assert len(exc_events) == 1
        attrs = {a.key: a.value for a in exc_events[0].attributes}
        assert attrs["exception.message"] == "boom"

    def test_failed_attaches_stdout_and_stderr_when_present(self):
        b = _make_builder_with_root()
        b.create_stage_span("s1")
        rule = b.create_rule_span("r1", "s1", "m1", "n1", {}, [])

        b.mark_rule_failed("r1", error="boom", stdout="o", stderr="e")
        names = [e.name for e in rule.events]
        assert names == ["exception", "stdout", "stderr"]

    def test_failed_on_unknown_rule_returns_none(self):
        b = _make_builder_with_root()
        assert b.mark_rule_failed("nope", "err") is None


class TestStageLifecycle:
    def test_stage_ok_when_no_rules_failed(self):
        b = _make_builder_with_root()
        stage = b.create_stage_span("s1")
        b.create_rule_span("r1", "s1", "m1", "n1", {}, [])
        b.mark_rule_completed("r1")

        b.mark_stage_completed("s1")
        assert stage.status == SpanStatus.OK
        assert stage.end_time_ns is not None

    def test_stage_error_when_any_rule_failed(self):
        b = _make_builder_with_root()
        stage = b.create_stage_span("s1")
        b.create_rule_span("r1", "s1", "m1", "n1", {}, [])
        b.create_rule_span("r2", "s1", "m1", "n2", {}, [])
        b.mark_rule_completed("r1")
        b.mark_rule_failed("r2", error="x")

        b.mark_stage_completed("s1")
        assert stage.status == SpanStatus.ERROR

    def test_stage_failure_does_not_leak_across_stages(self):
        b = _make_builder_with_root()
        s1 = b.create_stage_span("s1")
        s2 = b.create_stage_span("s2")
        b.create_rule_span("r1", "s1", "m1", "n1", {}, [])
        b.create_rule_span("r2", "s2", "m1", "n2", {}, [])
        b.mark_rule_failed("r1", error="x")
        b.mark_rule_completed("r2")

        b.mark_stage_completed("s1")
        b.mark_stage_completed("s2")
        assert s1.status == SpanStatus.ERROR
        assert s2.status == SpanStatus.OK

    def test_stage_completed_on_unknown_stage_returns_none(self):
        b = _make_builder_with_root()
        assert b.mark_stage_completed("nope") is None


class TestBenchmarkLifecycle:
    def test_success_sets_ok(self):
        b = _make_builder_with_root()
        root = b.mark_benchmark_completed(success=True)
        assert root.status == SpanStatus.OK
        assert root.end_time_ns is not None
        assert root.status_message == ""

    def test_failure_sets_error_with_message(self):
        b = _make_builder_with_root()
        root = b.mark_benchmark_completed(success=False, message="aborted")
        assert root.status == SpanStatus.ERROR
        assert root.status_message == "aborted"

    def test_no_root_returns_none(self):
        b = SpanBuilder()
        assert b.mark_benchmark_completed(success=True) is None
