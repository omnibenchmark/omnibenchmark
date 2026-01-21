"""
Telemetry module for omnibenchmark.

Provides OTLP-compatible structured logging that can be piped to
logfire, OpenTelemetry collectors, or other observability backends.

The telemetry output represents the benchmark DAG as a span hierarchy:

    benchmark_run (root span)
    ├── stage: datasets
    │   ├── rule: datasets_iris_abc123
    │   └── rule: datasets_wine_def456
    ├── stage: methods
    │   ├── rule: methods_kmeans_iris_...
    │   └── ...
    └── stage: metrics
        └── ...

Usage:
    from omnibenchmark.telemetry import TelemetryEmitter

    emitter = TelemetryEmitter(output=sys.stdout)
    emitter.emit_manifest(benchmark, resolved_nodes)

    # During execution:
    emitter.rule_started(rule_name, node_id)
    emitter.rule_completed(rule_name, node_id, stdout, stderr)
    emitter.rule_failed(rule_name, node_id, error, stdout, stderr)
"""

from omnibenchmark.telemetry.emitter import TelemetryEmitter
from omnibenchmark.telemetry.events import SpanStatus
from omnibenchmark.telemetry.spans import SpanBuilder

__all__ = ["TelemetryEmitter", "SpanBuilder", "SpanStatus"]
