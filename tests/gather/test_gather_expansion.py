"""Tests for gather stage node expansion (Phase 3).

Tests _build_gather_inputs which enumerates all outputs from provider
stages into a flat input dict for the gather node.
"""

import pytest
from dataclasses import dataclass, field
from typing import List

from omnibenchmark.model.benchmark import (
    GatherInput,
    Stage,
    Module,
    IOFile,
    Repository,
    Benchmark,
    SoftwareEnvironment,
    SoftwareBackendEnum,
)
from omnibenchmark.cli.run import _build_gather_inputs


# ============================================================================
# Minimal stub for ResolvedNode (avoid importing real one which needs
# ResolvedModule with commit validation etc.)
# ============================================================================


@dataclass
class FakeResolvedNode:
    """Minimal stand-in for ResolvedNode with only the fields
    _build_gather_inputs accesses: stage_id and outputs."""

    stage_id: str
    outputs: List[str] = field(default_factory=list)

    # Fields that exist on real ResolvedNode but aren't used by _build_gather_inputs
    id: str = ""
    module_id: str = ""


# ============================================================================
# Factories
# ============================================================================


def _repo(**kwargs):
    defaults = {"url": "https://example.com/repo.git", "commit": "abc123"}
    return Repository(**{**defaults, **kwargs})


def _module(**kwargs):
    defaults = {
        "id": "mod1",
        "software_environment": "host",
        "repository": _repo(),
    }
    return Module(**{**defaults, **kwargs})


def _iofile(**kwargs):
    defaults = {"id": "out1", "path": "output.txt"}
    return IOFile(**{**defaults, **kwargs})


def _stage(**kwargs):
    defaults = {
        "id": "stage1",
        "modules": [_module()],
        "outputs": [_iofile()],
    }
    return Stage(**{**defaults, **kwargs})


def _env(**kwargs):
    defaults = {"id": "host", "description": "Host execution"}
    return SoftwareEnvironment(**{**defaults, **kwargs})


def _benchmark(**kwargs):
    defaults = {
        "id": "test_benchmark",
        "benchmarker": "tester",
        "version": "1.0",
        "software_backend": SoftwareBackendEnum.host,
        "software_environments": [_env()],
        "stages": [],
    }
    return Benchmark(**{**defaults, **kwargs})


# ============================================================================
# Tests for _build_gather_inputs
# ============================================================================


@pytest.mark.short
class TestBuildGatherInputs:
    def test_single_provider_single_output(self):
        """One provider stage, one node, one output."""
        data_stage = _stage(
            id="data",
            provides={"dataset": "data.raw"},
            outputs=[_iofile(id="data.raw", path="data.json")],
        )
        gather_stage = _stage(
            id="summary",
            inputs=[GatherInput(gather="dataset")],
            outputs=[_iofile(id="summary.report", path="report.json")],
        )
        bm = _benchmark(stages=[data_stage, gather_stage])

        output_to_nodes = {
            "data.raw": [("data-D1-.abc123", "data/D1/.abc123/D1_data.json")],
        }

        inputs, name_mapping = _build_gather_inputs(
            ["dataset"], [], bm, output_to_nodes
        )

        assert inputs == {"input_0": "data/D1/.abc123/D1_data.json"}
        assert name_mapping == {"input_0": "dataset"}

    def test_single_provider_multiple_nodes(self):
        """One provider stage with multiple expanded nodes."""
        data_stage = _stage(
            id="data",
            provides={"dataset": "data.raw"},
            outputs=[_iofile(id="data.raw", path="data.json")],
        )
        gather_stage = _stage(
            id="summary",
            inputs=[GatherInput(gather="dataset")],
            outputs=[_iofile(id="summary.report", path="report.json")],
        )
        bm = _benchmark(stages=[data_stage, gather_stage])

        output_to_nodes = {
            "data.raw": [
                ("data-D1-.abc123", "data/D1/.abc123/D1_data.json"),
                ("data-D2-.def456", "data/D2/.def456/D2_data.json"),
            ],
        }

        inputs, name_mapping = _build_gather_inputs(
            ["dataset"], [], bm, output_to_nodes
        )

        assert len(inputs) == 2
        assert inputs["input_0"] == "data/D1/.abc123/D1_data.json"
        assert inputs["input_1"] == "data/D2/.def456/D2_data.json"
        assert all(v == "dataset" for v in name_mapping.values())

    def test_multiple_providers_same_label(self):
        """Two stages both provide "method", gather collects from both."""
        methods_fast = _stage(
            id="methods_fast",
            modules=[_module(id="M1")],
            provides={"method": "methods_fast.result"},
            outputs=[_iofile(id="methods_fast.result", path="fast.json")],
        )
        methods_slow = _stage(
            id="methods_slow",
            modules=[_module(id="M2")],
            provides={"method": "methods_slow.result"},
            outputs=[_iofile(id="methods_slow.result", path="slow.json")],
        )
        gather_stage = _stage(
            id="summary",
            modules=[_module(id="S1")],
            inputs=[GatherInput(gather="method")],
            outputs=[_iofile(id="summary.report", path="report.json")],
        )
        bm = _benchmark(stages=[methods_fast, methods_slow, gather_stage])

        output_to_nodes = {
            "methods_fast.result": [
                ("fast-M1-.abc-D1", "methods_fast/M1/.abc/D1_M1_fast.json"),
                ("fast-M1-.abc-D2", "methods_fast/M1/.abc/D2_M1_fast.json"),
            ],
            "methods_slow.result": [
                ("slow-M2-.def-D1", "methods_slow/M2/.def/D1_M2_slow.json"),
            ],
        }

        inputs, name_mapping = _build_gather_inputs(["method"], [], bm, output_to_nodes)

        assert len(inputs) == 3
        assert inputs["input_0"] == "methods_fast/M1/.abc/D1_M1_fast.json"
        assert inputs["input_1"] == "methods_fast/M1/.abc/D2_M1_fast.json"
        assert inputs["input_2"] == "methods_slow/M2/.def/D1_M2_slow.json"
        assert all(v == "method" for v in name_mapping.values())

    def test_non_provider_output_ids_excluded(self):
        """output_to_nodes entries for non-provided output IDs are not included."""
        data_stage = _stage(
            id="data",
            outputs=[_iofile(id="data.raw", path="data.json")],
        )
        methods_stage = _stage(
            id="methods",
            provides={"method": "methods.result"},
            outputs=[_iofile(id="methods.result", path="result.json")],
        )
        gather_stage = _stage(
            id="summary",
            inputs=[GatherInput(gather="method")],
            outputs=[_iofile(id="summary.report", path="report.json")],
        )
        bm = _benchmark(stages=[data_stage, methods_stage, gather_stage])

        output_to_nodes = {
            "data.raw": [("data-D1-.abc", "data/D1/.abc/D1_data.json")],
            "methods.result": [
                ("methods-M1-.def", "methods/M1/.def/D1_M1_result.json")
            ],
        }

        inputs, name_mapping = _build_gather_inputs(["method"], [], bm, output_to_nodes)

        assert len(inputs) == 1
        assert inputs["input_0"] == "methods/M1/.def/D1_M1_result.json"
        assert "data/D1/.abc/D1_data.json" not in inputs.values()

    def test_stage_with_multiple_outputs_only_declared_one_gathered(self):
        """A provider stage with multiple output files only exposes the declared one."""
        methods_stage = _stage(
            id="methods",
            provides={"method": "methods.result"},
            outputs=[
                _iofile(id="methods.result", path="result.json"),
                _iofile(id="methods.log", path="log.txt"),
            ],
        )
        gather_stage = _stage(
            id="summary",
            inputs=[GatherInput(gather="method")],
            outputs=[_iofile(id="summary.report", path="report.json")],
        )
        bm = _benchmark(stages=[methods_stage, gather_stage])

        output_to_nodes = {
            "methods.result": [("methods-M1-.abc", "methods/M1/.abc/result.json")],
            "methods.log": [("methods-M1-.abc", "methods/M1/.abc/log.txt")],
        }

        inputs, name_mapping = _build_gather_inputs(["method"], [], bm, output_to_nodes)

        # Only methods.result is gathered (declared in provides), not methods.log
        assert len(inputs) == 1
        assert inputs["input_0"] == "methods/M1/.abc/result.json"
        assert "methods/M1/.abc/log.txt" not in inputs.values()

    def test_empty_provider_nodes(self):
        """No provider nodes resolved yet — gather inputs are empty."""
        data_stage = _stage(
            id="data",
            provides={"dataset": "data.raw"},
            outputs=[_iofile(id="data.raw", path="data.json")],
        )
        gather_stage = _stage(
            id="summary",
            inputs=[GatherInput(gather="dataset")],
            outputs=[_iofile(id="summary.report", path="report.json")],
        )
        bm = _benchmark(stages=[data_stage, gather_stage])

        inputs, name_mapping = _build_gather_inputs(["dataset"], [], bm, {})

        assert inputs == {}
        assert name_mapping == {}

    def test_gather_multiple_labels(self):
        """Gathering two different labels collects from both."""
        data_stage = _stage(
            id="data",
            provides={"dataset": "data.raw"},
            outputs=[_iofile(id="data.raw", path="data.json")],
        )
        methods_stage = _stage(
            id="methods",
            provides={"method": "methods.result"},
            outputs=[_iofile(id="methods.result", path="result.json")],
        )
        gather_stage = _stage(
            id="summary",
            inputs=[
                GatherInput(gather="dataset"),
                GatherInput(gather="method"),
            ],
            outputs=[_iofile(id="summary.report", path="report.json")],
        )
        bm = _benchmark(stages=[data_stage, methods_stage, gather_stage])

        output_to_nodes = {
            "data.raw": [("data-D1-.abc", "data/D1/.abc/D1_data.json")],
            "methods.result": [("methods-M1-.def", "methods/M1/.def/M1_result.json")],
        }

        inputs, name_mapping = _build_gather_inputs(
            ["dataset", "method"], [], bm, output_to_nodes
        )

        assert len(inputs) == 2
        assert inputs["input_0"] == "data/D1/.abc/D1_data.json"
        assert inputs["input_1"] == "methods/M1/.def/M1_result.json"
        assert name_mapping["input_0"] == "dataset"
        assert name_mapping["input_1"] == "method"

    def test_keys_are_sequential(self):
        """Input keys are input_0, input_1, ... in order."""
        stage = _stage(
            id="data",
            provides={"dataset": "data.raw"},
            outputs=[_iofile(id="data.raw", path="data.json")],
        )
        gather = _stage(
            id="summary",
            inputs=[GatherInput(gather="dataset")],
            outputs=[_iofile(id="s.out", path="report.json")],
        )
        bm = _benchmark(stages=[stage, gather])

        output_to_nodes = {
            "data.raw": [
                (f"data-D{i}-.default", f"data/D{i}/out.json") for i in range(5)
            ]
        }

        inputs, _ = _build_gather_inputs(["dataset"], [], bm, output_to_nodes)

        assert list(inputs.keys()) == [f"input_{i}" for i in range(5)]
