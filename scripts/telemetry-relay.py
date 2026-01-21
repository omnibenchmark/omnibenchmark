#!/usr/bin/env python3
"""
Relay OTLP JSON telemetry to protobuf-based backends (Aspire Dashboard, etc.)

This script reads JSON Lines from stdin (or a file) and forwards them as
protobuf to an OTLP gRPC endpoint, preserving span hierarchy (traceId, spanId, parentSpanId).

Usage:
    # Stream from file
    tail -f telemetry.jsonl | python scripts/telemetry-relay.py

    # With custom endpoint
    tail -f telemetry.jsonl | python scripts/telemetry-relay.py --endpoint localhost:4317

    # From file (batch)
    python scripts/telemetry-relay.py --input telemetry.jsonl

Requirements:
    pip install opentelemetry-proto grpcio

Aspire Dashboard setup:
    docker run --rm -d --name aspire \\
      -p 18888:18888 \\
      -p 4317:18889 \\
      mcr.microsoft.com/dotnet/aspire-dashboard:latest

    Then open: http://localhost:18888
    Get token from: docker logs aspire 2>&1 | grep "login?t="
"""

import argparse
import json
import sys

try:
    import grpc
    from opentelemetry.proto.collector.trace.v1 import trace_service_pb2
    from opentelemetry.proto.collector.trace.v1 import trace_service_pb2_grpc
    from opentelemetry.proto.collector.logs.v1 import logs_service_pb2
    from opentelemetry.proto.collector.logs.v1 import logs_service_pb2_grpc
    from opentelemetry.proto.trace.v1 import trace_pb2
    from opentelemetry.proto.logs.v1 import logs_pb2
    from opentelemetry.proto.common.v1 import common_pb2
    from opentelemetry.proto.resource.v1 import resource_pb2
except ImportError:
    print("Error: Required packages not installed.", file=sys.stderr)
    print("Install with: pip install opentelemetry-proto grpcio", file=sys.stderr)
    sys.exit(1)


def hex_to_bytes(hex_str: str) -> bytes:
    """Convert hex string to bytes."""
    return bytes.fromhex(hex_str)


def json_value_to_any_value(value_obj: dict) -> common_pb2.AnyValue:
    """Convert JSON attribute value to protobuf AnyValue."""
    if "stringValue" in value_obj:
        return common_pb2.AnyValue(string_value=value_obj["stringValue"])
    elif "intValue" in value_obj:
        return common_pb2.AnyValue(int_value=int(value_obj["intValue"]))
    elif "doubleValue" in value_obj:
        return common_pb2.AnyValue(double_value=value_obj["doubleValue"])
    elif "boolValue" in value_obj:
        return common_pb2.AnyValue(bool_value=value_obj["boolValue"])
    elif "arrayValue" in value_obj:
        values = value_obj["arrayValue"].get("values", [])
        array_values = [json_value_to_any_value(v) for v in values]
        return common_pb2.AnyValue(
            array_value=common_pb2.ArrayValue(values=array_values)
        )
    return common_pb2.AnyValue(string_value="")


def json_to_protobuf(json_data: dict) -> trace_service_pb2.ExportTraceServiceRequest:
    """Convert OTLP JSON to protobuf ExportTraceServiceRequest."""
    resource_spans_list = []

    for rs in json_data.get("resourceSpans", []):
        # Build resource
        resource_attrs = []
        resource_data = rs.get("resource", {})
        for attr in resource_data.get("attributes", []):
            key = attr.get("key", "")
            value_obj = attr.get("value", {})
            resource_attrs.append(
                common_pb2.KeyValue(key=key, value=json_value_to_any_value(value_obj))
            )
        resource = resource_pb2.Resource(attributes=resource_attrs)

        # Build scope spans
        scope_spans_list = []
        for ss in rs.get("scopeSpans", []):
            scope_data = ss.get("scope", {})
            scope = common_pb2.InstrumentationScope(
                name=scope_data.get("name", ""),
                version=scope_data.get("version", ""),
            )

            # Build spans
            spans = []
            for span_data in ss.get("spans", []):
                # Convert IDs from hex to bytes
                trace_id = hex_to_bytes(span_data.get("traceId", "0" * 32))
                span_id = hex_to_bytes(span_data.get("spanId", "0" * 16))
                parent_span_id = b""
                if span_data.get("parentSpanId"):
                    parent_span_id = hex_to_bytes(span_data["parentSpanId"])

                # Build attributes
                attributes = []
                for attr in span_data.get("attributes", []):
                    key = attr.get("key", "")
                    value_obj = attr.get("value", {})
                    attributes.append(
                        common_pb2.KeyValue(
                            key=key, value=json_value_to_any_value(value_obj)
                        )
                    )

                # Build events
                events = []
                for event_data in span_data.get("events", []):
                    event_attrs = []
                    for attr in event_data.get("attributes", []):
                        key = attr.get("key", "")
                        value_obj = attr.get("value", {})
                        event_attrs.append(
                            common_pb2.KeyValue(
                                key=key, value=json_value_to_any_value(value_obj)
                            )
                        )
                    events.append(
                        trace_pb2.Span.Event(
                            time_unix_nano=int(event_data.get("timeUnixNano", 0)),
                            name=event_data.get("name", ""),
                            attributes=event_attrs,
                        )
                    )

                # Build status
                status_data = span_data.get("status", {})
                status = trace_pb2.Status(
                    code=status_data.get("code", 0),
                    message=status_data.get("message", ""),
                )

                # Create span
                span = trace_pb2.Span(
                    trace_id=trace_id,
                    span_id=span_id,
                    parent_span_id=parent_span_id,
                    name=span_data.get("name", ""),
                    kind=span_data.get("kind", 1),
                    start_time_unix_nano=int(span_data.get("startTimeUnixNano", 0)),
                    end_time_unix_nano=int(span_data.get("endTimeUnixNano", 0)),
                    attributes=attributes,
                    events=events,
                    status=status,
                )
                spans.append(span)

            scope_spans_list.append(
                trace_pb2.ScopeSpans(
                    scope=scope,
                    spans=spans,
                )
            )

        resource_spans_list.append(
            trace_pb2.ResourceSpans(
                resource=resource,
                scope_spans=scope_spans_list,
            )
        )

    return trace_service_pb2.ExportTraceServiceRequest(
        resource_spans=resource_spans_list
    )


def json_to_logs_protobuf(json_data: dict) -> logs_service_pb2.ExportLogsServiceRequest:
    """Convert OTLP JSON logs to protobuf ExportLogsServiceRequest."""
    resource_logs_list = []

    for rl in json_data.get("resourceLogs", []):
        # Build resource
        resource_attrs = []
        resource_data = rl.get("resource", {})
        for attr in resource_data.get("attributes", []):
            key = attr.get("key", "")
            value_obj = attr.get("value", {})
            resource_attrs.append(
                common_pb2.KeyValue(key=key, value=json_value_to_any_value(value_obj))
            )
        resource = resource_pb2.Resource(attributes=resource_attrs)

        # Build scope logs
        scope_logs_list = []
        for sl in rl.get("scopeLogs", []):
            scope_data = sl.get("scope", {})
            scope = common_pb2.InstrumentationScope(
                name=scope_data.get("name", ""),
                version=scope_data.get("version", ""),
            )

            # Build log records
            log_records = []
            for log_data in sl.get("logRecords", []):
                # Convert IDs from hex to bytes
                trace_id = b""
                span_id = b""
                if log_data.get("traceId"):
                    trace_id = hex_to_bytes(log_data["traceId"])
                if log_data.get("spanId"):
                    span_id = hex_to_bytes(log_data["spanId"])

                # Build attributes
                attributes = []
                for attr in log_data.get("attributes", []):
                    key = attr.get("key", "")
                    value_obj = attr.get("value", {})
                    attributes.append(
                        common_pb2.KeyValue(
                            key=key, value=json_value_to_any_value(value_obj)
                        )
                    )

                # Build body
                body_data = log_data.get("body", {})
                body = json_value_to_any_value(body_data)

                # Create log record
                log_record = logs_pb2.LogRecord(
                    time_unix_nano=int(log_data.get("timeUnixNano", 0)),
                    severity_number=log_data.get("severityNumber", 0),
                    severity_text=log_data.get("severityText", ""),
                    body=body,
                    attributes=attributes,
                    trace_id=trace_id,
                    span_id=span_id,
                )
                log_records.append(log_record)

            scope_logs_list.append(
                logs_pb2.ScopeLogs(
                    scope=scope,
                    log_records=log_records,
                )
            )

        resource_logs_list.append(
            logs_pb2.ResourceLogs(
                resource=resource,
                scope_logs=scope_logs_list,
            )
        )

    return logs_service_pb2.ExportLogsServiceRequest(resource_logs=resource_logs_list)


class TelemetryRelay:
    """Relays OTLP JSON to protobuf gRPC endpoint."""

    def __init__(self, endpoint: str = "localhost:4317"):
        self.endpoint = endpoint
        self.channel = grpc.insecure_channel(endpoint)
        self.trace_stub = trace_service_pb2_grpc.TraceServiceStub(self.channel)
        self.logs_stub = logs_service_pb2_grpc.LogsServiceStub(self.channel)

    def send(self, json_data: dict) -> bool:
        """Send a single OTLP JSON object as protobuf (traces or logs)."""
        try:
            if "resourceSpans" in json_data:
                request = json_to_protobuf(json_data)
                self.trace_stub.Export(request)
            elif "resourceLogs" in json_data:
                request = json_to_logs_protobuf(json_data)
                self.logs_stub.Export(request)
            else:
                print(
                    f"Unknown telemetry format: {list(json_data.keys())}",
                    file=sys.stderr,
                )
                return False
            return True
        except grpc.RpcError as e:
            print(f"gRPC error: {e}", file=sys.stderr)
            return False
        except Exception as e:
            print(f"Error: {e}", file=sys.stderr)
            return False

    def close(self):
        """Close the gRPC channel."""
        self.channel.close()


def main():
    parser = argparse.ArgumentParser(
        description="Relay OTLP JSON telemetry to protobuf backends",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--endpoint",
        default="localhost:4317",
        help="OTLP gRPC endpoint (default: localhost:4317)",
    )
    parser.add_argument(
        "--input",
        "-i",
        type=str,
        default=None,
        help="Input file (default: stdin)",
    )
    args = parser.parse_args()

    relay = TelemetryRelay(endpoint=args.endpoint)

    # Open input
    if args.input:
        input_file = open(args.input, "r")
    else:
        input_file = sys.stdin

    try:
        print(f"Relaying telemetry to {args.endpoint}...", file=sys.stderr)
        spans = 0
        logs = 0
        errors = 0
        for line in input_file:
            line = line.strip()
            if not line:
                continue
            try:
                json_data = json.loads(line)
                if relay.send(json_data):
                    if "resourceSpans" in json_data:
                        spans += 1
                    elif "resourceLogs" in json_data:
                        logs += 1
                    total = spans + logs
                    if total % 10 == 0:
                        print(f"  Sent {spans} spans, {logs} logs...", file=sys.stderr)
                else:
                    errors += 1
            except json.JSONDecodeError as e:
                print(f"Invalid JSON: {e}", file=sys.stderr)
                errors += 1
        print(
            f"Done. Sent {spans} spans, {logs} logs, {errors} errors.", file=sys.stderr
        )
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
    finally:
        if args.input:
            input_file.close()
        relay.close()


if __name__ == "__main__":
    main()
