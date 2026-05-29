"""Behavioral tests for omnibenchmark.telemetry.events.

Tests focus on observable OTLP output and identifier guarantees, not the
specific dispatch ladder used inside to_otlp().
"""

import time

import pytest

from omnibenchmark.telemetry.events import (
    Attribute,
    LogRecord,
    SeverityNumber,
    Span,
    SpanEvent,
    SpanKind,
    SpanStatus,
    generate_span_id,
    generate_trace_id,
    now_ns,
)

pytestmark = pytest.mark.short


# ---------------------------------------------------------------------------
# Attribute serialization
# ---------------------------------------------------------------------------


class TestAttributeToOtlp:
    def test_string_value(self):
        out = Attribute("k", "hello").to_otlp()
        assert out == {"key": "k", "value": {"stringValue": "hello"}}

    def test_bool_is_not_serialized_as_int(self):
        # bool is an int subclass in Python; OTLP must distinguish them.
        assert Attribute("k", True).to_otlp() == {
            "key": "k",
            "value": {"boolValue": True},
        }
        assert Attribute("k", False).to_otlp() == {
            "key": "k",
            "value": {"boolValue": False},
        }

    def test_int_is_serialized_as_string(self):
        # OTLP requires intValue to be a string (int64 wire format).
        out = Attribute("count", 42).to_otlp()
        assert out == {"key": "count", "value": {"intValue": "42"}}
        assert isinstance(out["value"]["intValue"], str)

    def test_float_value(self):
        out = Attribute("ratio", 1.5).to_otlp()
        assert out == {"key": "ratio", "value": {"doubleValue": 1.5}}

    def test_list_serialized_as_string_array(self):
        out = Attribute("files", ["a.txt", "b.txt"]).to_otlp()
        assert out == {
            "key": "files",
            "value": {
                "arrayValue": {
                    "values": [
                        {"stringValue": "a.txt"},
                        {"stringValue": "b.txt"},
                    ]
                }
            },
        }

    def test_tuple_serialized_same_as_list(self):
        from_list = Attribute("x", ["a", "b"]).to_otlp()
        from_tuple = Attribute("x", ("a", "b")).to_otlp()
        assert from_list == from_tuple

    def test_list_of_ints_is_stringified(self):
        out = Attribute("ids", [1, 2, 3]).to_otlp()
        values = out["value"]["arrayValue"]["values"]
        assert values == [
            {"stringValue": "1"},
            {"stringValue": "2"},
            {"stringValue": "3"},
        ]

    def test_unknown_type_falls_back_to_string(self):
        class Custom:
            def __str__(self):
                return "custom!"

        out = Attribute("k", Custom()).to_otlp()
        assert out == {"key": "k", "value": {"stringValue": "custom!"}}

    def test_none_falls_back_to_string(self):
        out = Attribute("k", None).to_otlp()
        assert out == {"key": "k", "value": {"stringValue": "None"}}


# ---------------------------------------------------------------------------
# SpanEvent serialization
# ---------------------------------------------------------------------------


class TestSpanEventToOtlp:
    def test_round_trip_shape(self):
        ev = SpanEvent(
            name="stdout",
            timestamp_ns=123,
            attributes=[Attribute("output", "hi")],
        )
        out = ev.to_otlp()
        assert out["name"] == "stdout"
        assert out["timeUnixNano"] == "123"  # OTLP requires string
        assert out["attributes"] == [{"key": "output", "value": {"stringValue": "hi"}}]

    def test_no_attributes_yields_empty_list(self):
        out = SpanEvent(name="x", timestamp_ns=1).to_otlp()
        assert out["attributes"] == []


# ---------------------------------------------------------------------------
# Span serialization
# ---------------------------------------------------------------------------


def _make_span(**overrides) -> Span:
    base = dict(
        trace_id="t" * 32,
        span_id="s" * 16,
        parent_span_id=None,
        name="root",
    )
    base.update(overrides)
    return Span(**base)


class TestSpanToOtlp:
    def test_required_fields_always_present(self):
        out = _make_span().to_otlp()
        assert out["traceId"] == "t" * 32
        assert out["spanId"] == "s" * 16
        assert out["name"] == "root"
        assert out["kind"] == int(SpanKind.INTERNAL)
        assert out["status"] == {"code": int(SpanStatus.UNSET)}
        assert out["attributes"] == []

    def test_optional_fields_omitted_when_unset(self):
        out = _make_span().to_otlp()
        assert "parentSpanId" not in out
        assert "startTimeUnixNano" not in out
        assert "endTimeUnixNano" not in out
        assert "events" not in out
        assert "message" not in out["status"]

    def test_parent_included_when_set(self):
        out = _make_span(parent_span_id="p" * 16).to_otlp()
        assert out["parentSpanId"] == "p" * 16

    def test_times_serialized_as_strings(self):
        out = _make_span(start_time_ns=10, end_time_ns=20).to_otlp()
        assert out["startTimeUnixNano"] == "10"
        assert out["endTimeUnixNano"] == "20"

    def test_kind_uses_int_value(self):
        out = _make_span(kind=SpanKind.SERVER).to_otlp()
        assert out["kind"] == int(SpanKind.SERVER)

    def test_error_status_adds_semantic_attributes(self):
        out = _make_span(status=SpanStatus.ERROR, status_message="boom").to_otlp()
        attr_pairs = {a["key"]: a["value"] for a in out["attributes"]}
        assert attr_pairs["otel.status_code"] == {"stringValue": "ERROR"}
        assert attr_pairs["otel.status_description"] == {"stringValue": "boom"}
        assert out["status"]["code"] == int(SpanStatus.ERROR)
        assert out["status"]["message"] == "boom"

    def test_error_without_message_omits_description(self):
        out = _make_span(status=SpanStatus.ERROR).to_otlp()
        keys = {a["key"] for a in out["attributes"]}
        assert "otel.status_code" in keys
        assert "otel.status_description" not in keys
        assert "message" not in out["status"]

    def test_ok_status_adds_semantic_code(self):
        out = _make_span(status=SpanStatus.OK).to_otlp()
        keys = {a["key"]: a["value"] for a in out["attributes"]}
        assert keys["otel.status_code"] == {"stringValue": "OK"}

    def test_unset_status_adds_no_semantic_attributes(self):
        out = _make_span(status=SpanStatus.UNSET).to_otlp()
        assert out["attributes"] == []

    def test_user_attributes_preserved_alongside_semantic_ones(self):
        out = _make_span(
            status=SpanStatus.OK,
            attributes=[Attribute("custom.key", "v")],
        ).to_otlp()
        keys = [a["key"] for a in out["attributes"]]
        assert "custom.key" in keys
        assert "otel.status_code" in keys

    def test_events_serialized_when_present(self):
        ev = SpanEvent(name="stdout", timestamp_ns=5)
        out = _make_span(events=[ev]).to_otlp()
        assert out["events"] == [ev.to_otlp()]


# ---------------------------------------------------------------------------
# LogRecord serialization
# ---------------------------------------------------------------------------


class TestLogRecordToOtlp:
    def test_basic_shape(self):
        rec = LogRecord(
            trace_id="t" * 32,
            span_id="s" * 16,
            body="hello",
            severity=SeverityNumber.WARN,
            timestamp_ns=42,
        )
        out = rec.to_otlp()
        assert out["traceId"] == "t" * 32
        assert out["spanId"] == "s" * 16
        assert out["body"] == {"stringValue": "hello"}
        assert out["severityNumber"] == int(SeverityNumber.WARN)
        assert out["timeUnixNano"] == "42"

    def test_severity_text_only_emitted_when_set(self):
        out = LogRecord(
            trace_id="t" * 32, span_id="s" * 16, body="x", timestamp_ns=1
        ).to_otlp()
        assert "severityText" not in out

        out2 = LogRecord(
            trace_id="t" * 32,
            span_id="s" * 16,
            body="x",
            timestamp_ns=1,
            severity_text="WARN",
        ).to_otlp()
        assert out2["severityText"] == "WARN"

    def test_missing_timestamp_uses_current_time(self, monkeypatch):
        monkeypatch.setattr("omnibenchmark.telemetry.events.now_ns", lambda: 999)
        out = LogRecord(trace_id="t" * 32, span_id="s" * 16, body="x").to_otlp()
        assert out["timeUnixNano"] == "999"


# ---------------------------------------------------------------------------
# ID and time helpers
# ---------------------------------------------------------------------------


class TestIdGenerators:
    def test_trace_id_is_32_hex_chars(self):
        tid = generate_trace_id()
        assert len(tid) == 32
        int(tid, 16)  # parses as hex

    def test_span_id_is_16_hex_chars(self):
        sid = generate_span_id()
        assert len(sid) == 16
        int(sid, 16)

    def test_ids_are_unique_across_calls(self):
        traces = {generate_trace_id() for _ in range(50)}
        spans = {generate_span_id() for _ in range(50)}
        assert len(traces) == 50
        assert len(spans) == 50


class TestNowNs:
    def test_returns_nanoseconds_close_to_wall_clock(self):
        expected = time.time() * 1_000_000_000
        actual = now_ns()
        # Within 1 second of wall clock.
        assert abs(actual - expected) < 1_000_000_000
        assert isinstance(actual, int)
