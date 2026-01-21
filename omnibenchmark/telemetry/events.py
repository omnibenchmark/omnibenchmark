"""
Event types and schema for OTLP-compatible telemetry.

This module defines the data structures that map to OpenTelemetry's
trace model: spans, attributes, status, and events.
"""

from dataclasses import dataclass, field
from enum import IntEnum
from typing import Any, Optional
import time
import uuid


class SpanKind(IntEnum):
    """OpenTelemetry span kinds."""

    INTERNAL = 1
    SERVER = 2
    CLIENT = 3
    PRODUCER = 4
    CONSUMER = 5


class SpanStatus(IntEnum):
    """OpenTelemetry span status codes."""

    UNSET = 0  # Pending / not yet determined
    OK = 1  # Completed successfully
    ERROR = 2  # Failed


@dataclass
class Attribute:
    """A key-value attribute for spans."""

    key: str
    value: Any

    def to_otlp(self) -> dict:
        """Convert to OTLP attribute format."""
        v = self.value
        if isinstance(v, bool):
            return {"key": self.key, "value": {"boolValue": v}}
        elif isinstance(v, int):
            return {"key": self.key, "value": {"intValue": str(v)}}
        elif isinstance(v, float):
            return {"key": self.key, "value": {"doubleValue": v}}
        elif isinstance(v, (list, tuple)):
            # Array of strings
            return {
                "key": self.key,
                "value": {
                    "arrayValue": {"values": [{"stringValue": str(x)} for x in v]}
                },
            }
        else:
            return {"key": self.key, "value": {"stringValue": str(v)}}


@dataclass
class SpanEvent:
    """An event within a span (e.g., log line, error)."""

    name: str
    timestamp_ns: int
    attributes: list[Attribute] = field(default_factory=list)

    def to_otlp(self) -> dict:
        return {
            "timeUnixNano": str(self.timestamp_ns),
            "name": self.name,
            "attributes": [a.to_otlp() for a in self.attributes],
        }


@dataclass
class Span:
    """
    An OTLP-compatible span representing a unit of work.

    In omnibenchmark, spans represent:
    - The entire benchmark run (root)
    - Individual stages
    - Individual rules/jobs
    """

    trace_id: str
    span_id: str
    parent_span_id: Optional[str]
    name: str
    kind: SpanKind = SpanKind.INTERNAL
    start_time_ns: Optional[int] = None
    end_time_ns: Optional[int] = None
    status: SpanStatus = SpanStatus.UNSET
    status_message: str = ""
    attributes: list[Attribute] = field(default_factory=list)
    events: list[SpanEvent] = field(default_factory=list)

    def to_otlp(self) -> dict:
        """Convert to OTLP span format."""
        span = {
            "traceId": self.trace_id,
            "spanId": self.span_id,
            "name": self.name,
            "kind": int(self.kind),
            "status": {
                "code": int(self.status),
            },
            "attributes": [a.to_otlp() for a in self.attributes],
        }

        if self.parent_span_id:
            span["parentSpanId"] = self.parent_span_id

        if self.start_time_ns:
            span["startTimeUnixNano"] = str(self.start_time_ns)

        if self.end_time_ns:
            span["endTimeUnixNano"] = str(self.end_time_ns)

        if self.status_message:
            span["status"]["message"] = self.status_message

        if self.events:
            span["events"] = [e.to_otlp() for e in self.events]

        return span


def generate_trace_id() -> str:
    """Generate a 32-character hex trace ID (16 bytes)."""
    return uuid.uuid4().hex  # uuid4().hex is already 32 chars


def generate_span_id() -> str:
    """Generate a 16-character hex span ID (8 bytes)."""
    return uuid.uuid4().hex[:16]


def now_ns() -> int:
    """Current time in nanoseconds since Unix epoch."""
    return int(time.time() * 1_000_000_000)


class SeverityNumber(IntEnum):
    """OpenTelemetry log severity numbers."""

    UNSPECIFIED = 0
    TRACE = 1
    DEBUG = 5
    INFO = 9
    WARN = 13
    ERROR = 17
    FATAL = 21


@dataclass
class LogRecord:
    """
    An OTLP-compatible log record.

    Log records can be correlated with spans via trace_id and span_id.
    """

    trace_id: str
    span_id: str
    body: str
    severity: SeverityNumber = SeverityNumber.INFO
    severity_text: str = ""
    timestamp_ns: Optional[int] = None
    attributes: list[Attribute] = field(default_factory=list)

    def to_otlp(self) -> dict:
        """Convert to OTLP log record format."""
        record = {
            "timeUnixNano": str(self.timestamp_ns or now_ns()),
            "severityNumber": int(self.severity),
            "body": {"stringValue": self.body},
            "attributes": [a.to_otlp() for a in self.attributes],
            "traceId": self.trace_id,
            "spanId": self.span_id,
        }

        if self.severity_text:
            record["severityText"] = self.severity_text

        return record
