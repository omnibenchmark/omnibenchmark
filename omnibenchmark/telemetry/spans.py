"""
Span hierarchy builder for benchmark DAG.

Converts the benchmark structure (stages, modules, rules) into
a hierarchy of OTLP spans that can be browsed in logfire.
"""

from dataclasses import dataclass, field
from typing import Optional

from omnibenchmark.telemetry.events import (
    Span,
    SpanKind,
    SpanStatus,
    Attribute,
    SpanEvent,
    generate_trace_id,
    generate_span_id,
    now_ns,
)


@dataclass
class SpanBuilder:
    """
    Builds and manages the span hierarchy for a benchmark run.

    Hierarchy:
        benchmark (root)
        ├── stage_1
        │   ├── rule_1a
        │   └── rule_1b
        ├── stage_2
        │   ├── rule_2a
        │   └── rule_2b
        └── ...
    """

    trace_id: str = field(default_factory=generate_trace_id)
    service_name: str = "omnibenchmark"

    # Root span for the benchmark
    _root_span: Optional[Span] = None

    # Stage spans keyed by stage_id
    _stage_spans: dict[str, Span] = field(default_factory=dict)

    # Rule spans keyed by rule_name (full snakemake rule name)
    _rule_spans: dict[str, Span] = field(default_factory=dict)

    # Mapping from rule_name to stage_id for hierarchy
    _rule_to_stage: dict[str, str] = field(default_factory=dict)

    def create_benchmark_span(
        self,
        benchmark_name: str,
        benchmark_version: str,
        benchmark_author: str,
        total_rules: int,
        software_backend: str,
        cores: int,
    ) -> Span:
        """Create the root span for the benchmark run."""
        self._root_span = Span(
            trace_id=self.trace_id,
            span_id=generate_span_id(),
            parent_span_id=None,
            name=f"benchmark: {benchmark_name}",
            kind=SpanKind.SERVER,
            start_time_ns=now_ns(),
            status=SpanStatus.UNSET,
            attributes=[
                Attribute("benchmark.name", benchmark_name),
                Attribute("benchmark.version", benchmark_version),
                Attribute("benchmark.author", benchmark_author),
                Attribute("benchmark.total_rules", total_rules),
                Attribute("benchmark.software_backend", software_backend),
                Attribute("benchmark.cores", cores),
                Attribute("service.name", self.service_name),
            ],
        )
        return self._root_span

    def create_stage_span(
        self,
        stage_id: str,
        stage_name: Optional[str] = None,
        module_count: int = 0,
        rule_count: int = 0,
    ) -> Span:
        """Create a span for a benchmark stage."""
        if not self._root_span:
            raise RuntimeError("Must create benchmark span before stage spans")

        span = Span(
            trace_id=self.trace_id,
            span_id=generate_span_id(),
            parent_span_id=self._root_span.span_id,
            name=f"stage: {stage_name or stage_id}",
            kind=SpanKind.INTERNAL,
            status=SpanStatus.UNSET,
            attributes=[
                Attribute("stage.id", stage_id),
                Attribute("stage.name", stage_name or stage_id),
                Attribute("stage.module_count", module_count),
                Attribute("stage.rule_count", rule_count),
            ],
        )
        self._stage_spans[stage_id] = span
        return span

    def create_rule_span(
        self,
        rule_name: str,
        stage_id: str,
        module_id: str,
        node_id: str,
        inputs: dict[str, str],
        outputs: list[str],
        parameters: Optional[str] = None,
    ) -> Span:
        """Create a span for an individual rule/job."""
        parent_span = self._stage_spans.get(stage_id)
        if not parent_span:
            # Fallback to root if stage not found
            parent_span = self._root_span

        if not parent_span:
            raise RuntimeError("Must create benchmark span before rule spans")

        attributes = [
            Attribute("rule.name", rule_name),
            Attribute("rule.stage_id", stage_id),
            Attribute("rule.module_id", module_id),
            Attribute("rule.node_id", node_id),
            Attribute("rule.outputs", outputs),
        ]

        if inputs:
            attributes.append(Attribute("rule.inputs", list(inputs.values())))

        if parameters:
            attributes.append(Attribute("rule.parameters", parameters))

        span = Span(
            trace_id=self.trace_id,
            span_id=generate_span_id(),
            parent_span_id=parent_span.span_id,
            name=f"rule: {rule_name}",
            kind=SpanKind.INTERNAL,
            status=SpanStatus.UNSET,
            attributes=attributes,
        )

        self._rule_spans[rule_name] = span
        self._rule_to_stage[rule_name] = stage_id
        return span

    def get_all_spans(self) -> list[Span]:
        """Get all spans in hierarchy order (root, stages, rules)."""
        spans = []
        if self._root_span:
            spans.append(self._root_span)
        spans.extend(self._stage_spans.values())
        spans.extend(self._rule_spans.values())
        return spans

    def mark_rule_started(self, rule_name: str) -> Optional[Span]:
        """Mark a rule as started (set start time)."""
        span = self._rule_spans.get(rule_name)
        if span:
            span.start_time_ns = now_ns()
            # Also mark parent stage as started if not already
            stage_id = self._rule_to_stage.get(rule_name)
            if stage_id:
                stage_span = self._stage_spans.get(stage_id)
                if stage_span and not stage_span.start_time_ns:
                    stage_span.start_time_ns = now_ns()
        return span

    def mark_rule_completed(
        self,
        rule_name: str,
        stdout: Optional[str] = None,
        stderr: Optional[str] = None,
    ) -> Optional[Span]:
        """Mark a rule as successfully completed."""
        span = self._rule_spans.get(rule_name)
        if span:
            span.end_time_ns = now_ns()
            span.status = SpanStatus.OK

            if stdout:
                span.events.append(
                    SpanEvent(
                        name="stdout",
                        timestamp_ns=span.end_time_ns,
                        attributes=[Attribute("output", stdout)],
                    )
                )

            if stderr:
                span.events.append(
                    SpanEvent(
                        name="stderr",
                        timestamp_ns=span.end_time_ns,
                        attributes=[Attribute("output", stderr)],
                    )
                )
        return span

    def mark_rule_failed(
        self,
        rule_name: str,
        error: str,
        stdout: Optional[str] = None,
        stderr: Optional[str] = None,
    ) -> Optional[Span]:
        """Mark a rule as failed."""
        span = self._rule_spans.get(rule_name)
        if span:
            span.end_time_ns = now_ns()
            span.status = SpanStatus.ERROR
            span.status_message = error

            span.events.append(
                SpanEvent(
                    name="exception",
                    timestamp_ns=span.end_time_ns,
                    attributes=[Attribute("exception.message", error)],
                )
            )

            if stdout:
                span.events.append(
                    SpanEvent(
                        name="stdout",
                        timestamp_ns=span.end_time_ns,
                        attributes=[Attribute("output", stdout)],
                    )
                )

            if stderr:
                span.events.append(
                    SpanEvent(
                        name="stderr",
                        timestamp_ns=span.end_time_ns,
                        attributes=[Attribute("output", stderr)],
                    )
                )
        return span

    def mark_stage_completed(self, stage_id: str) -> Optional[Span]:
        """Mark a stage as completed (all its rules finished)."""
        span = self._stage_spans.get(stage_id)
        if span:
            span.end_time_ns = now_ns()
            # Check if any rules in this stage failed
            has_failures = any(
                self._rule_spans[r].status == SpanStatus.ERROR
                for r, s in self._rule_to_stage.items()
                if s == stage_id and r in self._rule_spans
            )
            span.status = SpanStatus.ERROR if has_failures else SpanStatus.OK
        return span

    def mark_benchmark_completed(
        self, success: bool, message: str = ""
    ) -> Optional[Span]:
        """Mark the benchmark run as completed."""
        if self._root_span:
            self._root_span.end_time_ns = now_ns()
            self._root_span.status = SpanStatus.OK if success else SpanStatus.ERROR
            if message:
                self._root_span.status_message = message
        return self._root_span

    def get_rule_span(self, rule_name: str) -> Optional[Span]:
        """Get a rule span by name."""
        return self._rule_spans.get(rule_name)

    def get_stage_span(self, stage_id: str) -> Optional[Span]:
        """Get a stage span by ID."""
        return self._stage_spans.get(stage_id)
