"""
OTLP JSON → protobuf relay for telemetry data.

Reads OTLP JSON Lines and forwards them as protobuf to an OTLP gRPC endpoint
(Aspire Dashboard, Jaeger, etc.), preserving span hierarchy.

Optional dependency: opentelemetry-proto, grpcio.
ImportError is raised at class instantiation time if these are missing.

CLI entry point: ``python -m omnibenchmark.telemetry.relay`` (or via the
thin shim at ``scripts/telemetry-relay.py``).

Usage (CLI):
    tail -f telemetry.jsonl | python -m omnibenchmark.telemetry.relay
    python -m omnibenchmark.telemetry.relay --input telemetry.jsonl --endpoint localhost:18889
"""

import argparse
import json
import sys


def _import_grpc():
    """Import gRPC and opentelemetry-proto; raise ImportError with a helpful message if missing."""
    try:
        import grpc
        from opentelemetry.proto.collector.trace.v1 import (
            trace_service_pb2,
            trace_service_pb2_grpc,
        )
        from opentelemetry.proto.collector.logs.v1 import (
            logs_service_pb2,
            logs_service_pb2_grpc,
        )
        from opentelemetry.proto.collector.metrics.v1 import (
            metrics_service_pb2,
            metrics_service_pb2_grpc,
        )
        from opentelemetry.proto.trace.v1 import trace_pb2
        from opentelemetry.proto.logs.v1 import logs_pb2
        from opentelemetry.proto.metrics.v1 import metrics_pb2
        from opentelemetry.proto.common.v1 import common_pb2
        from opentelemetry.proto.resource.v1 import resource_pb2
    except ImportError as exc:
        raise ImportError(
            "Telemetry relay requires opentelemetry-proto and grpcio. "
            "Install with: pip install opentelemetry-proto grpcio"
        ) from exc
    return (
        grpc,
        trace_service_pb2,
        trace_service_pb2_grpc,
        logs_service_pb2,
        logs_service_pb2_grpc,
        metrics_service_pb2,
        metrics_service_pb2_grpc,
        trace_pb2,
        logs_pb2,
        metrics_pb2,
        common_pb2,
        resource_pb2,
    )


def hex_to_bytes(hex_str: str) -> bytes:
    return bytes.fromhex(hex_str)


def json_value_to_any_value(value_obj: dict):
    """Convert a JSON attribute value object to a protobuf AnyValue."""
    from opentelemetry.proto.common.v1 import common_pb2

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


def json_to_protobuf(json_data: dict):
    """Convert OTLP JSON traces to protobuf ExportTraceServiceRequest."""
    from opentelemetry.proto.collector.trace.v1 import trace_service_pb2
    from opentelemetry.proto.trace.v1 import trace_pb2
    from opentelemetry.proto.common.v1 import common_pb2
    from opentelemetry.proto.resource.v1 import resource_pb2

    resource_spans_list = []
    for rs in json_data.get("resourceSpans", []):
        resource_attrs = [
            common_pb2.KeyValue(
                key=attr.get("key", ""),
                value=json_value_to_any_value(attr.get("value", {})),
            )
            for attr in rs.get("resource", {}).get("attributes", [])
        ]
        resource = resource_pb2.Resource(attributes=resource_attrs)

        scope_spans_list = []
        for ss in rs.get("scopeSpans", []):
            scope_data = ss.get("scope", {})
            scope = common_pb2.InstrumentationScope(
                name=scope_data.get("name", ""),
                version=scope_data.get("version", ""),
            )
            spans = []
            for span_data in ss.get("spans", []):
                trace_id = hex_to_bytes(span_data.get("traceId", "0" * 32))
                span_id = hex_to_bytes(span_data.get("spanId", "0" * 16))
                parent_span_id = b""
                if span_data.get("parentSpanId"):
                    parent_span_id = hex_to_bytes(span_data["parentSpanId"])

                attributes = [
                    common_pb2.KeyValue(
                        key=attr.get("key", ""),
                        value=json_value_to_any_value(attr.get("value", {})),
                    )
                    for attr in span_data.get("attributes", [])
                ]
                events = []
                for event_data in span_data.get("events", []):
                    event_attrs = [
                        common_pb2.KeyValue(
                            key=attr.get("key", ""),
                            value=json_value_to_any_value(attr.get("value", {})),
                        )
                        for attr in event_data.get("attributes", [])
                    ]
                    events.append(
                        trace_pb2.Span.Event(
                            time_unix_nano=int(event_data.get("timeUnixNano", 0)),
                            name=event_data.get("name", ""),
                            attributes=event_attrs,
                        )
                    )
                status_data = span_data.get("status", {})
                status = trace_pb2.Status(
                    code=status_data.get("code", 0),
                    message=status_data.get("message", ""),
                )
                spans.append(
                    trace_pb2.Span(
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
                )
            scope_spans_list.append(trace_pb2.ScopeSpans(scope=scope, spans=spans))

        resource_spans_list.append(
            trace_pb2.ResourceSpans(resource=resource, scope_spans=scope_spans_list)
        )

    return trace_service_pb2.ExportTraceServiceRequest(
        resource_spans=resource_spans_list
    )


def json_to_logs_protobuf(json_data: dict):
    """Convert OTLP JSON logs to protobuf ExportLogsServiceRequest."""
    from opentelemetry.proto.collector.logs.v1 import logs_service_pb2
    from opentelemetry.proto.logs.v1 import logs_pb2
    from opentelemetry.proto.common.v1 import common_pb2
    from opentelemetry.proto.resource.v1 import resource_pb2

    resource_logs_list = []
    for rl in json_data.get("resourceLogs", []):
        resource_attrs = [
            common_pb2.KeyValue(
                key=attr.get("key", ""),
                value=json_value_to_any_value(attr.get("value", {})),
            )
            for attr in rl.get("resource", {}).get("attributes", [])
        ]
        resource = resource_pb2.Resource(attributes=resource_attrs)

        scope_logs_list = []
        for sl in rl.get("scopeLogs", []):
            scope_data = sl.get("scope", {})
            scope = common_pb2.InstrumentationScope(
                name=scope_data.get("name", ""),
                version=scope_data.get("version", ""),
            )
            log_records = []
            for log_data in sl.get("logRecords", []):
                trace_id = (
                    hex_to_bytes(log_data["traceId"])
                    if log_data.get("traceId")
                    else b""
                )
                span_id = (
                    hex_to_bytes(log_data["spanId"]) if log_data.get("spanId") else b""
                )
                attributes = [
                    common_pb2.KeyValue(
                        key=attr.get("key", ""),
                        value=json_value_to_any_value(attr.get("value", {})),
                    )
                    for attr in log_data.get("attributes", [])
                ]
                body = json_value_to_any_value(log_data.get("body", {}))
                log_records.append(
                    logs_pb2.LogRecord(
                        time_unix_nano=int(log_data.get("timeUnixNano", 0)),
                        severity_number=log_data.get("severityNumber", 0),
                        severity_text=log_data.get("severityText", ""),
                        body=body,
                        attributes=attributes,
                        trace_id=trace_id,
                        span_id=span_id,
                    )
                )
            scope_logs_list.append(
                logs_pb2.ScopeLogs(scope=scope, log_records=log_records)
            )

        resource_logs_list.append(
            logs_pb2.ResourceLogs(resource=resource, scope_logs=scope_logs_list)
        )

    return logs_service_pb2.ExportLogsServiceRequest(resource_logs=resource_logs_list)


def _build_exemplars(exemplars_data: list) -> list:
    from opentelemetry.proto.metrics.v1 import metrics_pb2

    result = []
    for ex in exemplars_data:
        kwargs = {
            "time_unix_nano": int(ex.get("timeUnixNano", 0)),
            "as_double": float(ex.get("asDouble", 0.0)),
        }
        if ex.get("traceId"):
            kwargs["trace_id"] = hex_to_bytes(ex["traceId"])
        if ex.get("spanId"):
            kwargs["span_id"] = hex_to_bytes(ex["spanId"])
        result.append(metrics_pb2.Exemplar(**kwargs))
    return result


def _build_gauge_data_points(data_points: list) -> list:
    from opentelemetry.proto.metrics.v1 import metrics_pb2
    from opentelemetry.proto.common.v1 import common_pb2

    result = []
    for dp in data_points:
        attrs = [
            common_pb2.KeyValue(
                key=attr.get("key", ""),
                value=json_value_to_any_value(attr.get("value", {})),
            )
            for attr in dp.get("attributes", [])
        ]
        result.append(
            metrics_pb2.NumberDataPoint(
                attributes=attrs,
                time_unix_nano=int(dp.get("timeUnixNano", 0)),
                as_double=float(dp.get("asDouble", 0.0)),
                exemplars=_build_exemplars(dp.get("exemplars", [])),
            )
        )
    return result


def json_to_metrics_protobuf(json_data: dict):
    """Convert OTLP JSON metrics to protobuf ExportMetricsServiceRequest."""
    from opentelemetry.proto.collector.metrics.v1 import metrics_service_pb2
    from opentelemetry.proto.metrics.v1 import metrics_pb2
    from opentelemetry.proto.common.v1 import common_pb2
    from opentelemetry.proto.resource.v1 import resource_pb2

    resource_metrics_list = []
    for rm in json_data.get("resourceMetrics", []):
        resource_attrs = [
            common_pb2.KeyValue(
                key=attr.get("key", ""),
                value=json_value_to_any_value(attr.get("value", {})),
            )
            for attr in rm.get("resource", {}).get("attributes", [])
        ]
        resource = resource_pb2.Resource(attributes=resource_attrs)

        scope_metrics_list = []
        for sm in rm.get("scopeMetrics", []):
            scope_data = sm.get("scope", {})
            scope = common_pb2.InstrumentationScope(
                name=scope_data.get("name", ""),
                version=scope_data.get("version", ""),
            )
            metric_protos = []
            for metric in sm.get("metrics", []):
                name = metric.get("name", "")
                description = metric.get("description", "")
                unit = metric.get("unit", "")
                if "gauge" in metric:
                    data_points = _build_gauge_data_points(
                        metric["gauge"].get("dataPoints", [])
                    )
                    metric_protos.append(
                        metrics_pb2.Metric(
                            name=name,
                            description=description,
                            unit=unit,
                            gauge=metrics_pb2.Gauge(data_points=data_points),
                        )
                    )
                elif "sum" in metric:
                    data_points = _build_gauge_data_points(
                        metric["sum"].get("dataPoints", [])
                    )
                    metric_protos.append(
                        metrics_pb2.Metric(
                            name=name,
                            description=description,
                            unit=unit,
                            sum=metrics_pb2.Sum(
                                data_points=data_points,
                                aggregation_temporality=metric["sum"].get(
                                    "aggregationTemporality",
                                    metrics_pb2.AGGREGATION_TEMPORALITY_CUMULATIVE,
                                ),
                                is_monotonic=metric["sum"].get("isMonotonic", False),
                            ),
                        )
                    )
            scope_metrics_list.append(
                metrics_pb2.ScopeMetrics(scope=scope, metrics=metric_protos)
            )

        resource_metrics_list.append(
            metrics_pb2.ResourceMetrics(
                resource=resource, scope_metrics=scope_metrics_list
            )
        )

    return metrics_service_pb2.ExportMetricsServiceRequest(
        resource_metrics=resource_metrics_list
    )


class TelemetryRelay:
    """Relays OTLP JSON Lines to a protobuf gRPC endpoint."""

    def __init__(self, endpoint: str = "localhost:18889"):
        import grpc
        from opentelemetry.proto.collector.trace.v1 import trace_service_pb2_grpc
        from opentelemetry.proto.collector.logs.v1 import logs_service_pb2_grpc
        from opentelemetry.proto.collector.metrics.v1 import metrics_service_pb2_grpc

        self.endpoint = endpoint
        self.channel = grpc.insecure_channel(endpoint)
        self.trace_stub = trace_service_pb2_grpc.TraceServiceStub(self.channel)
        self.logs_stub = logs_service_pb2_grpc.LogsServiceStub(self.channel)
        self.metrics_stub = metrics_service_pb2_grpc.MetricsServiceStub(self.channel)

    def send(self, json_data: dict) -> bool:
        """Send a single OTLP JSON object as protobuf. Returns True on success."""
        import grpc

        try:
            if "resourceSpans" in json_data:
                self.trace_stub.Export(json_to_protobuf(json_data))
            elif "resourceLogs" in json_data:
                self.logs_stub.Export(json_to_logs_protobuf(json_data))
            elif "resourceMetrics" in json_data:
                self.metrics_stub.Export(json_to_metrics_protobuf(json_data))
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
        self.channel.close()


def main():
    parser = argparse.ArgumentParser(
        description="Relay OTLP JSON telemetry to protobuf backends",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--endpoint",
        default="localhost:18889",
        help="OTLP gRPC endpoint (default: localhost:18889)",
    )
    parser.add_argument(
        "--input", "-i", type=str, default=None, help="Input file (default: stdin)"
    )
    args = parser.parse_args()

    try:
        relay = TelemetryRelay(endpoint=args.endpoint)
    except ImportError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    input_file = open(args.input, "r") if args.input else sys.stdin
    try:
        print(f"Relaying telemetry to {args.endpoint}...", file=sys.stderr)
        spans = logs = metrics = errors = 0
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
                    elif "resourceMetrics" in json_data:
                        metrics += 1
                    total = spans + logs + metrics
                    if total % 10 == 0:
                        print(
                            f"  Sent {spans} spans, {logs} logs, {metrics} metrics...",
                            file=sys.stderr,
                        )
                else:
                    errors += 1
            except json.JSONDecodeError as e:
                print(f"Invalid JSON: {e}", file=sys.stderr)
                errors += 1
        print(
            f"Done. Sent {spans} spans, {logs} logs, {metrics} metrics, {errors} errors.",
            file=sys.stderr,
        )
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
    finally:
        if args.input:
            input_file.close()
        relay.close()


if __name__ == "__main__":
    main()
