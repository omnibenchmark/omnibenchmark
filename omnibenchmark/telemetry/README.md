# Omnibenchmark Telemetry

OTLP-compatible telemetry for observing benchmark execution in logfire or other OpenTelemetry backends.

## Overview

The telemetry module emits structured JSON Lines (NDJSON) following the OpenTelemetry Protocol (OTLP) format. This allows you to:

- View the benchmark DAG as a browseable hierarchy in logfire
- Monitor rule execution status (pending → running → completed/failed)
- Capture stdout/stderr from rule execution
- Pipe telemetry over SSH from remote execution (e.g., Slurm clusters)

## Span Hierarchy

The benchmark structure is represented as nested spans:

```
benchmark: my_benchmark (root span)
├── stage: datasets
│   ├── rule: datasets_iris_abc123
│   └── rule: datasets_wine_def456
├── stage: methods
│   ├── rule: methods_kmeans_iris_abc123
│   ├── rule: methods_kmeans_wine_def456
│   └── ...
└── stage: metrics
    └── ...
```

## Quick Start

### 1. Start the Aspire Dashboard

The .NET Aspire Dashboard provides excellent trace visualization with minimal setup.

```bash
docker run --rm -d --name aspire \
  -p 18888:18888 \
  -p 4317:18889 \
  mcr.microsoft.com/dotnet/aspire-dashboard:latest
```

- **UI**: http://localhost:18888
- **OTLP gRPC endpoint**: localhost:4317 (mapped from container port 18889)

**Get the login token:**
```bash
docker logs aspire 2>&1 | grep "login?t="
```

### 2. Run Benchmark with Telemetry

```bash
ob run benchmark.yaml --telemetry json --telemetry-output telemetry.jsonl
```

### 3. Stream Telemetry to Aspire

Aspire only accepts protobuf over gRPC, so use the relay script to convert JSON to protobuf:

```bash
# Install dependencies (if not already installed)
pip install opentelemetry-proto grpcio

# Stream telemetry to Aspire (run in separate terminal)
tail -f telemetry.jsonl | python scripts/telemetry-relay.py --endpoint localhost:4317
```

Or send a completed telemetry file:

```bash
python scripts/telemetry-relay.py --input telemetry.jsonl --endpoint localhost:4317
```

## Output Format

The telemetry output is JSON Lines (NDJSON) with one OTLP TracesData object per line:

```json
{"resourceSpans":[{"resource":{"attributes":[{"key":"service.name","value":{"stringValue":"omnibenchmark"}}]},"scopeSpans":[{"scope":{"name":"omnibenchmark.telemetry","version":"1.0.0"},"spans":[{"traceId":"abc123...","spanId":"def456...","name":"benchmark: my_benchmark","kind":2,"status":{"code":0},"attributes":[...]}]}]}]}
```

### Span Attributes

**Benchmark span:**
- `benchmark.name` - Benchmark identifier
- `benchmark.version` - Semantic version
- `benchmark.author` - Author name
- `benchmark.total_rules` - Total number of rules
- `benchmark.software_backend` - conda, apptainer, docker, etc.
- `benchmark.cores` - Number of CPU cores

**Stage span:**
- `stage.id` - Stage identifier
- `stage.name` - Human-readable name
- `stage.module_count` - Number of modules in stage
- `stage.rule_count` - Number of rules in stage

**Rule span:**
- `rule.name` - Snakemake rule name
- `rule.stage_id` - Parent stage
- `rule.module_id` - Module identifier
- `rule.node_id` - Full node identifier
- `rule.inputs` - Input file paths
- `rule.outputs` - Output file paths
- `rule.parameters` - JSON-encoded parameters

### Span Status

- `code: 0` (UNSET) - Pending, not yet started
- `code: 1` (OK) - Completed successfully
- `code: 2` (ERROR) - Failed with error

### Span Events

Rule spans include events for:
- `stdout` - Captured standard output
- `stderr` - Captured standard error
- `exception` - Error message on failure

## Programmatic Usage

```python
from omnibenchmark.telemetry import TelemetryEmitter

# Create emitter (defaults to stdout)
emitter = TelemetryEmitter()

# Or write to file
emitter = TelemetryEmitter(output=Path("telemetry.jsonl"))

# Emit the manifest (full DAG) upfront
emitter.emit_manifest(
    benchmark_name="my_benchmark",
    benchmark_version="1.0.0",
    benchmark_author="Jane Doe",
    software_backend="conda",
    cores=4,
    stages=[{"id": "datasets", "name": "Datasets"}, ...],
    resolved_nodes=resolved_nodes,
)

# Stream execution updates
emitter.rule_started("datasets_iris_abc123")
emitter.rule_completed("datasets_iris_abc123", stdout="...", stderr="...")
# or
emitter.rule_failed("datasets_wine_def456", error="FileNotFoundError", stderr="...")

# Mark completion
emitter.benchmark_completed(success=True)
```

## Troubleshooting

### No spans appearing in logfire

1. Check the collector is receiving data:
   ```bash
   # Use debug exporter to see incoming spans
   docker logs <collector_container>
   ```

2. Verify JSON format:
   ```bash
   ob run benchmark.yaml --telemetry json | head -1 | jq .
   ```

3. Check OTLP endpoint is reachable:
   ```bash
   curl -v http://localhost:4318/v1/traces
   ```

### Spans not showing hierarchy

Ensure parent-child relationships are preserved:
- All spans share the same `traceId`
- Child spans have `parentSpanId` pointing to their parent
- Spans are emitted in order (parent before children)

### Large stdout/stderr truncation

By default, stdout/stderr are captured in full. For very large outputs, consider:
- Configuring max capture size in the emitter
- Using log aggregation for detailed output
