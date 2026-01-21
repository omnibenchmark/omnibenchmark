"""
OTLP JSON emitter for telemetry data.

Outputs OTLP-compatible JSON Lines (NDJSON) that can be:
- Piped to logfire or other OTLP collectors
- Written to a .jsonl file
- Sent via HTTP to an OTLP endpoint

Spans are emitted only when they complete (with start and end times),
building up the trace hierarchy progressively as rules finish.

Hierarchy:
    benchmark
    └── stage
        └── module (groups rules by module_id within a stage)
            └── rule (one per parameter combination)
"""

import json
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO, Optional, Union

from omnibenchmark.telemetry.events import (
    Span,
    SpanStatus,
    SpanKind,
    Attribute,
    SpanEvent,
    LogRecord,
    SeverityNumber,
    generate_trace_id,
    generate_span_id,
    now_ns,
)


@dataclass
class TelemetryEmitter:
    """
    Emits OTLP-compatible telemetry as JSON Lines.

    Spans are emitted only when they complete, preserving hierarchy through
    parent span IDs. The trace builds up progressively:
    - Rules emit when they complete
    - Modules emit when all their rules complete
    - Stages emit when all their modules complete
    - Benchmark emits at the very end

    Usage:
        emitter = TelemetryEmitter(output=Path("telemetry.jsonl"))

        # Initialize with benchmark info (doesn't emit yet)
        emitter.init_benchmark(name, version, stages, resolved_nodes)

        # As rules execute:
        emitter.rule_started("rule_name")
        emitter.rule_completed("rule_name")  # Emits the rule span

        # At the end:
        emitter.benchmark_completed(success=True)  # Emits remaining hierarchy
    """

    output: Union[IO, Path, None] = None
    service_name: str = "omnibenchmark"
    service_version: str = "0.1.0"

    _file_handle: Optional[IO] = None
    _owns_handle: bool = False

    # Trace-wide ID
    _trace_id: str = field(default_factory=generate_trace_id)

    # Benchmark span info
    _benchmark_span_id: str = field(default_factory=generate_span_id)
    _benchmark_name: str = ""
    _benchmark_version: str = ""
    _benchmark_author: str = ""
    _benchmark_start_time: int = 0
    _total_rules: int = 0
    _software_backend: str = ""
    _cores: int = 1

    # Stage tracking: stage_id -> {span_id, name, start_time, modules: set, completed_modules: int, has_failures}
    _stages: dict = field(default_factory=dict)

    # Module tracking: (stage_id, module_id) -> {span_id, start_time, rule_count, completed_rules, has_failures}
    _modules: dict = field(default_factory=dict)

    # Rule tracking: rule_name -> {span_id, stage_id, module_id, node_id, inputs, outputs, params, start_time, ...}
    _rules: dict = field(default_factory=dict)

    # Track what's been emitted
    _emitted_modules: set = field(default_factory=set)
    _emitted_stages: set = field(default_factory=set)

    def __post_init__(self):
        if self.output is None:
            self._file_handle = sys.stdout
            self._owns_handle = False
        elif isinstance(self.output, Path):
            self._file_handle = open(self.output, "w")
            self._owns_handle = True
        else:
            self._file_handle = self.output
            self._owns_handle = False

    def close(self):
        """Close the output file if we own it."""
        if self._owns_handle and self._file_handle:
            self._file_handle.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def _emit_span(self, span: Span):
        """Emit a single span as an OTLP JSON line."""
        traces_data = {
            "resourceSpans": [
                {
                    "resource": {
                        "attributes": [
                            {
                                "key": "service.name",
                                "value": {"stringValue": self.service_name},
                            },
                            {
                                "key": "service.version",
                                "value": {"stringValue": self.service_version},
                            },
                        ]
                    },
                    "scopeSpans": [
                        {
                            "scope": {
                                "name": "omnibenchmark.telemetry",
                                "version": "1.0.0",
                            },
                            "spans": [span.to_otlp()],
                        }
                    ],
                }
            ]
        }

        line = json.dumps(traces_data, separators=(",", ":"))
        self._file_handle.write(line + "\n")
        self._file_handle.flush()

    def _emit_log(self, log: LogRecord):
        """Emit a single log record as an OTLP JSON line."""
        logs_data = {
            "resourceLogs": [
                {
                    "resource": {
                        "attributes": [
                            {
                                "key": "service.name",
                                "value": {"stringValue": self.service_name},
                            },
                            {
                                "key": "service.version",
                                "value": {"stringValue": self.service_version},
                            },
                        ]
                    },
                    "scopeLogs": [
                        {
                            "scope": {
                                "name": "omnibenchmark.telemetry",
                                "version": "1.0.0",
                            },
                            "logRecords": [log.to_otlp()],
                        }
                    ],
                }
            ]
        }

        line = json.dumps(logs_data, separators=(",", ":"))
        self._file_handle.write(line + "\n")
        self._file_handle.flush()

    def init_benchmark(
        self,
        benchmark_name: str,
        benchmark_version: str,
        benchmark_author: str,
        software_backend: str,
        cores: int,
        stages: list[dict],
        resolved_nodes: list,
    ):
        """
        Initialize benchmark metadata. Does NOT emit spans yet.

        Spans will be emitted as rules complete, building up the hierarchy.
        """
        self._benchmark_name = benchmark_name
        self._benchmark_version = benchmark_version
        self._benchmark_author = benchmark_author
        self._software_backend = software_backend
        self._cores = cores
        self._benchmark_start_time = now_ns()
        self._total_rules = len(resolved_nodes)

        # Group nodes by stage and module
        nodes_by_stage: dict[str, list] = {}
        nodes_by_module: dict[tuple, list] = {}  # (stage_id, module_id) -> nodes

        for node in resolved_nodes:
            stage_id = node.stage_id
            module_id = node.module_id

            if stage_id not in nodes_by_stage:
                nodes_by_stage[stage_id] = []
            nodes_by_stage[stage_id].append(node)

            module_key = (stage_id, module_id)
            if module_key not in nodes_by_module:
                nodes_by_module[module_key] = []
            nodes_by_module[module_key].append(node)

        # Initialize stage tracking
        for stage in stages:
            stage_id = stage.get("id", stage.get("stage_id"))
            stage_name = stage.get("name", stage_id)
            stage_nodes = nodes_by_stage.get(stage_id, [])
            module_ids = set(n.module_id for n in stage_nodes)

            self._stages[stage_id] = {
                "span_id": generate_span_id(),
                "name": stage_name,
                "start_time": None,
                "modules": module_ids,
                "module_count": len(module_ids),
                "completed_modules": 0,
                "has_failures": False,
            }

        # Initialize module tracking
        for (stage_id, module_id), nodes in nodes_by_module.items():
            self._modules[(stage_id, module_id)] = {
                "span_id": generate_span_id(),
                "stage_id": stage_id,
                "module_id": module_id,
                "start_time": None,
                "rule_count": len(nodes),
                "completed_rules": 0,
                "has_failures": False,
            }

        # Initialize rule tracking
        for node in resolved_nodes:
            # Rule name must match Snakefile's _sanitize_rule_name(node.id)
            rule_name = node.id.replace("-", "_").replace(".", "_")
            if rule_name and not rule_name[0].isalpha():
                rule_name = "rule_" + rule_name

            self._rules[rule_name] = {
                "span_id": generate_span_id(),
                "stage_id": node.stage_id,
                "module_id": node.module_id,
                "node_id": node.id,
                "param_id": node.param_id,
                "inputs": node.inputs if hasattr(node, "inputs") else {},
                "outputs": node.outputs if hasattr(node, "outputs") else [],
                "parameters": node.get_parameter_json()
                if hasattr(node, "get_parameter_json")
                else None,
                "log_file": f".logs/{rule_name}.log",  # Per-rule log file path
                "start_time": None,
                "end_time": None,
                "status": SpanStatus.UNSET,
                "error": None,
                "stdout": None,
                "stderr": None,
            }

    # Keep emit_manifest as alias for backward compatibility
    def emit_manifest(self, **kwargs):
        """Alias for init_benchmark (backward compatibility)."""
        self.init_benchmark(**kwargs)

    def get_rule_log_file(self, rule_name: str) -> Optional[str]:
        """Get the log file path for a rule (relative to out_dir)."""
        if rule_name in self._rules:
            return self._rules[rule_name]["log_file"]
        return None

    def emit_setup_span(
        self, output: str, start_time_ns: int = None, end_time_ns: int = None
    ):
        """Emit a setup span for environment preparation (conda, apptainer, etc.)."""
        span = Span(
            trace_id=self._trace_id,
            span_id=generate_span_id(),
            parent_span_id=self._benchmark_span_id,
            name="setup: environment preparation",
            kind=SpanKind.INTERNAL,
            start_time_ns=start_time_ns or self._benchmark_start_time,
            end_time_ns=end_time_ns or now_ns(),
            status=SpanStatus.OK,
            attributes=[
                Attribute("setup.type", "environment"),
            ],
        )
        self._emit_span(span)

        # Also emit as structured log for easy viewing
        if output:
            self._emit_log(
                LogRecord(
                    trace_id=self._trace_id,
                    span_id=span.span_id,
                    body=output,
                    severity=SeverityNumber.INFO,
                    severity_text="setup",
                    timestamp_ns=end_time_ns or now_ns(),
                    attributes=[
                        Attribute("phase", "environment_setup"),
                    ],
                )
            )

    def rule_started(self, rule_name: str):
        """Record when a rule starts (doesn't emit yet)."""
        if rule_name not in self._rules:
            return

        rule = self._rules[rule_name]
        rule["start_time"] = now_ns()

        stage_id = rule["stage_id"]
        module_id = rule["module_id"]
        module_key = (stage_id, module_id)

        # Mark module start time if not set
        if (
            module_key in self._modules
            and self._modules[module_key]["start_time"] is None
        ):
            self._modules[module_key]["start_time"] = now_ns()

        # Mark stage start time if not set
        if stage_id in self._stages and self._stages[stage_id]["start_time"] is None:
            self._stages[stage_id]["start_time"] = now_ns()

    def rule_completed(
        self,
        rule_name: str,
        stdout: Optional[str] = None,
        stderr: Optional[str] = None,
    ):
        """Mark rule as completed and emit its span."""
        if rule_name not in self._rules:
            return

        rule = self._rules[rule_name]
        rule["end_time"] = now_ns()
        rule["status"] = SpanStatus.OK
        rule["stdout"] = stdout
        rule["stderr"] = stderr

        if rule["start_time"] is None:
            rule["start_time"] = rule["end_time"]

        # Emit the rule span
        self._emit_rule_span(rule_name)

        # Update module tracking
        stage_id = rule["stage_id"]
        module_id = rule["module_id"]
        module_key = (stage_id, module_id)

        if module_key in self._modules:
            self._modules[module_key]["completed_rules"] += 1
            self._maybe_emit_module(module_key)

    def rule_failed(
        self,
        rule_name: str,
        error: str,
        stdout: Optional[str] = None,
        stderr: Optional[str] = None,
    ):
        """Mark rule as failed and emit its span."""
        if rule_name not in self._rules:
            return

        rule = self._rules[rule_name]
        rule["end_time"] = now_ns()
        rule["status"] = SpanStatus.ERROR
        rule["error"] = error
        rule["stdout"] = stdout
        rule["stderr"] = stderr

        if rule["start_time"] is None:
            rule["start_time"] = rule["end_time"]

        # Emit the rule span
        self._emit_rule_span(rule_name)

        # Update module tracking
        stage_id = rule["stage_id"]
        module_id = rule["module_id"]
        module_key = (stage_id, module_id)

        if module_key in self._modules:
            self._modules[module_key]["completed_rules"] += 1
            self._modules[module_key]["has_failures"] = True
            self._maybe_emit_module(module_key)

    def _emit_rule_span(self, rule_name: str):
        """Emit a completed rule span."""
        rule = self._rules[rule_name]
        stage_id = rule["stage_id"]
        module_id = rule["module_id"]
        module_key = (stage_id, module_id)

        # Parent is the module span
        parent_span_id = (
            self._modules[module_key]["span_id"]
            if module_key in self._modules
            else self._benchmark_span_id
        )

        attributes = [
            Attribute("rule.name", rule_name),
            Attribute("rule.stage_id", stage_id),
            Attribute("rule.module_id", module_id),
            Attribute("rule.node_id", rule["node_id"]),
            Attribute("rule.param_id", rule["param_id"]),
        ]

        if rule["outputs"]:
            attributes.append(Attribute("rule.outputs", rule["outputs"]))
        if rule["inputs"]:
            attributes.append(Attribute("rule.inputs", list(rule["inputs"].values())))
        if rule["parameters"]:
            attributes.append(Attribute("rule.parameters", rule["parameters"]))

        events = []
        if rule["stdout"]:
            events.append(
                SpanEvent(
                    name="stdout",
                    timestamp_ns=rule["end_time"],
                    attributes=[Attribute("output", rule["stdout"])],
                )
            )
        if rule["stderr"]:
            events.append(
                SpanEvent(
                    name="stderr",
                    timestamp_ns=rule["end_time"],
                    attributes=[Attribute("output", rule["stderr"])],
                )
            )
        if rule["error"]:
            events.append(
                SpanEvent(
                    name="exception",
                    timestamp_ns=rule["end_time"],
                    attributes=[Attribute("exception.message", rule["error"])],
                )
            )

        span = Span(
            trace_id=self._trace_id,
            span_id=rule["span_id"],
            parent_span_id=parent_span_id,
            name=f"rule: {rule_name}",
            kind=SpanKind.INTERNAL,
            start_time_ns=rule["start_time"],
            end_time_ns=rule["end_time"],
            status=rule["status"],
            status_message=rule["error"] or "",
            attributes=attributes,
            events=events,
        )

        self._emit_span(span)

        # Emit structured logs for stdout/stderr (correlated with the span)
        if rule["stdout"]:
            self._emit_log(
                LogRecord(
                    trace_id=self._trace_id,
                    span_id=rule["span_id"],
                    body=rule["stdout"],
                    severity=SeverityNumber.INFO,
                    severity_text="stdout",
                    timestamp_ns=rule["end_time"],
                    attributes=[
                        Attribute("rule.name", rule_name),
                        Attribute("stream", "stdout"),
                    ],
                )
            )

        if rule["stderr"]:
            # Use ERROR severity if the rule failed, WARN otherwise
            severity = (
                SeverityNumber.ERROR
                if rule["status"] == SpanStatus.ERROR
                else SeverityNumber.WARN
            )
            self._emit_log(
                LogRecord(
                    trace_id=self._trace_id,
                    span_id=rule["span_id"],
                    body=rule["stderr"],
                    severity=severity,
                    severity_text="stderr",
                    timestamp_ns=rule["end_time"],
                    attributes=[
                        Attribute("rule.name", rule_name),
                        Attribute("stream", "stderr"),
                    ],
                )
            )

    def _maybe_emit_module(self, module_key: tuple):
        """Emit module span if all its rules are complete."""
        if module_key in self._emitted_modules:
            return

        module = self._modules[module_key]
        if module["completed_rules"] >= module["rule_count"]:
            self._emit_module_span(module_key)
            self._emitted_modules.add(module_key)

            # Update stage tracking
            stage_id = module["stage_id"]
            if stage_id in self._stages:
                self._stages[stage_id]["completed_modules"] += 1
                if module["has_failures"]:
                    self._stages[stage_id]["has_failures"] = True
                self._maybe_emit_stage(stage_id)

    def _emit_module_span(self, module_key: tuple):
        """Emit a completed module span."""
        module = self._modules[module_key]
        stage_id = module["stage_id"]
        end_time = now_ns()

        # Parent is the stage span
        parent_span_id = (
            self._stages[stage_id]["span_id"]
            if stage_id in self._stages
            else self._benchmark_span_id
        )

        span = Span(
            trace_id=self._trace_id,
            span_id=module["span_id"],
            parent_span_id=parent_span_id,
            name=f"module: {module['module_id']}",
            kind=SpanKind.INTERNAL,
            start_time_ns=module["start_time"] or self._benchmark_start_time,
            end_time_ns=end_time,
            status=SpanStatus.ERROR if module["has_failures"] else SpanStatus.OK,
            attributes=[
                Attribute("module.id", module["module_id"]),
                Attribute("module.stage_id", stage_id),
                Attribute("module.rule_count", module["rule_count"]),
            ],
        )

        self._emit_span(span)

    def _maybe_emit_stage(self, stage_id: str):
        """Emit stage span if all its modules are complete."""
        if stage_id in self._emitted_stages:
            return

        stage = self._stages[stage_id]
        if stage["completed_modules"] >= stage["module_count"]:
            self._emit_stage_span(stage_id)
            self._emitted_stages.add(stage_id)

    def _emit_stage_span(self, stage_id: str):
        """Emit a completed stage span."""
        stage = self._stages[stage_id]
        end_time = now_ns()

        span = Span(
            trace_id=self._trace_id,
            span_id=stage["span_id"],
            parent_span_id=self._benchmark_span_id,
            name=f"stage: {stage['name']}",
            kind=SpanKind.INTERNAL,
            start_time_ns=stage["start_time"] or self._benchmark_start_time,
            end_time_ns=end_time,
            status=SpanStatus.ERROR if stage["has_failures"] else SpanStatus.OK,
            attributes=[
                Attribute("stage.id", stage_id),
                Attribute("stage.name", stage["name"]),
                Attribute("stage.module_count", stage["module_count"]),
            ],
        )

        self._emit_span(span)

    def stage_completed(self, stage_id: str):
        """Manually mark a stage as completed and emit it."""
        if stage_id not in self._emitted_stages:
            self._emit_stage_span(stage_id)
            self._emitted_stages.add(stage_id)

    def benchmark_completed(self, success: bool, message: str = ""):
        """Emit any remaining modules, stages, and the benchmark root span."""
        end_time = now_ns()

        # Emit any modules that haven't been emitted yet
        for module_key in self._modules:
            if module_key not in self._emitted_modules:
                self._emit_module_span(module_key)
                self._emitted_modules.add(module_key)

        # Emit any stages that haven't been emitted yet
        for stage_id in self._stages:
            if stage_id not in self._emitted_stages:
                self._emit_stage_span(stage_id)
                self._emitted_stages.add(stage_id)

        # Emit the benchmark root span
        span = Span(
            trace_id=self._trace_id,
            span_id=self._benchmark_span_id,
            parent_span_id=None,
            name=f"{self._benchmark_name}-{self._benchmark_version}",
            kind=SpanKind.SERVER,
            start_time_ns=self._benchmark_start_time,
            end_time_ns=end_time,
            status=SpanStatus.OK if success else SpanStatus.ERROR,
            status_message=message,
            attributes=[
                Attribute("benchmark.name", self._benchmark_name),
                Attribute("benchmark.version", self._benchmark_version),
                Attribute("benchmark.author", self._benchmark_author),
                Attribute("benchmark.total_rules", self._total_rules),
                Attribute("benchmark.software_backend", self._software_backend),
                Attribute("benchmark.cores", self._cores),
            ],
        )

        self._emit_span(span)

    @property
    def trace_id(self) -> str:
        """Get the trace ID for this run."""
        return self._trace_id
