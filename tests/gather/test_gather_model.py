"""Tests for gather stage model changes (Phase 1).

Tests the GatherInput model, provides field on Stage, input parsing,
and benchmark-level validation of gather semantics.
"""

import pytest
from pydantic import ValidationError

from omnibenchmark.model.benchmark import (
    GatherInput,
    InputCollection,
    Stage,
    Module,
    IOFile,
    Repository,
    Benchmark,
    SoftwareEnvironment,
    SoftwareBackendEnum,
)


# ============================================================================
# Factories (local to gather tests)
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
# GatherInput model
# ============================================================================


@pytest.mark.short
class TestGatherInput:
    def test_basic_creation(self):
        gi = GatherInput(gather="method")
        assert gi.gather == "method"

    def test_empty_label_rejected(self):
        with pytest.raises(ValidationError):
            GatherInput(gather="")

    def test_whitespace_label_rejected(self):
        with pytest.raises(ValidationError):
            GatherInput(gather="   ")

    def test_from_dict(self):
        gi = GatherInput(**{"gather": "dataset"})
        assert gi.gather == "dataset"


# ============================================================================
# Stage.provides field
# ============================================================================


@pytest.mark.short
class TestStageProvides:
    def test_stage_without_provides(self):
        stage = _stage()
        assert stage.provides is None

    def test_stage_with_provides(self):
        stage = _stage(provides=["dataset"])
        assert stage.provides == ["dataset"]

    def test_stage_with_multiple_provides(self):
        stage = _stage(provides=["dataset", "method"])
        assert stage.provides == ["dataset", "method"]

    def test_provides_empty_list(self):
        stage = _stage(provides=[])
        assert stage.provides == []


# ============================================================================
# Stage.is_gather_stage / get_gather_labels
# ============================================================================


@pytest.mark.short
class TestStageGatherHelpers:
    def test_not_gather_when_no_inputs(self):
        stage = _stage()
        assert stage.is_gather_stage() is False
        assert stage.get_gather_labels() == []

    def test_not_gather_with_regular_inputs(self):
        stage = _stage(inputs=[InputCollection(entries=["data.raw"])])
        assert stage.is_gather_stage() is False
        assert stage.get_gather_labels() == []

    def test_is_gather_with_gather_input(self):
        stage = _stage(inputs=[GatherInput(gather="method")])
        assert stage.is_gather_stage() is True
        assert stage.get_gather_labels() == ["method"]

    def test_multiple_gather_labels(self):
        stage = _stage(
            inputs=[GatherInput(gather="method"), GatherInput(gather="dataset")]
        )
        assert stage.is_gather_stage() is True
        assert stage.get_gather_labels() == ["method", "dataset"]


# ============================================================================
# Stage.inputs parsing (YAML-like dicts)
# ============================================================================


@pytest.mark.short
class TestStageInputParsing:
    def test_string_list_becomes_input_collection(self):
        """inputs: [data.raw, data.labels] -> single InputCollection."""
        stage = _stage(inputs=["data.raw", "data.labels"])
        assert len(stage.inputs) == 1
        assert isinstance(stage.inputs[0], InputCollection)
        assert stage.inputs[0].entries == ["data.raw", "data.labels"]

    def test_gather_dict_becomes_gather_input(self):
        """inputs: [{gather: method}] -> GatherInput."""
        stage = _stage(inputs=[{"gather": "method"}])
        assert len(stage.inputs) == 1
        assert isinstance(stage.inputs[0], GatherInput)
        assert stage.inputs[0].gather == "method"

    def test_mixed_inputs_rejected(self):
        """Cannot mix regular and gather inputs in the same stage."""
        with pytest.raises(ValidationError, match="cannot mix"):
            _stage(
                inputs=[
                    {"entries": ["data.raw"]},
                    {"gather": "method"},
                ]
            )

    def test_multiple_gather_inputs_accepted(self):
        """Multiple gather inputs in same stage are fine."""
        stage = _stage(inputs=[{"gather": "method"}, {"gather": "dataset"}])
        assert len(stage.inputs) == 2
        assert all(isinstance(inp, GatherInput) for inp in stage.inputs)

    def test_legacy_entries_format_still_works(self):
        """Legacy {entries: [...]} format still parses (with deprecation warning)."""
        with pytest.warns(FutureWarning, match="entries"):
            stage = _stage(inputs=[{"entries": ["data.raw"]}])
        assert len(stage.inputs) == 1
        assert isinstance(stage.inputs[0], InputCollection)


# ============================================================================
# Benchmark-level validation: provides/gather
# ============================================================================


@pytest.mark.short
class TestBenchmarkGatherValidation:
    def test_gather_with_valid_provider(self):
        """Gather referencing a valid provides label should pass."""
        data_stage = _stage(
            id="data",
            provides=["dataset"],
            outputs=[_iofile(id="data.raw", path="data.json")],
        )
        gather_stage = _stage(
            id="summary",
            inputs=[GatherInput(gather="dataset")],
            outputs=[_iofile(id="summary.report", path="report.json")],
        )
        benchmark = _benchmark(stages=[data_stage, gather_stage])
        # Should not raise
        assert benchmark is not None

    def test_gather_with_no_provider_fails(self):
        """Gather referencing a label that no stage provides should fail."""
        gather_stage = _stage(
            id="summary",
            inputs=[GatherInput(gather="nonexistent")],
            outputs=[_iofile(id="summary.report", path="report.json")],
        )
        from omnibenchmark.model.validation import ValidationError as OBValidationError

        with pytest.raises(OBValidationError, match="no stage provides it"):
            _benchmark(stages=[gather_stage])

    def test_gather_with_template_variable_in_output_fails(self):
        """Gather stage with {something} in output path should fail."""
        data_stage = _stage(
            id="data",
            provides=["dataset"],
            outputs=[_iofile(id="data.raw", path="data.json")],
        )
        gather_stage = _stage(
            id="summary",
            inputs=[GatherInput(gather="dataset")],
            outputs=[_iofile(id="summary.report", path="{dataset}_report.json")],
        )
        from omnibenchmark.model.validation import ValidationError as OBValidationError

        with pytest.raises(OBValidationError, match="template variable"):
            _benchmark(stages=[data_stage, gather_stage])

    def test_gather_with_fixed_output_path_passes(self):
        """Gather stage with no template variables should pass."""
        data_stage = _stage(
            id="data",
            provides=["dataset"],
            outputs=[_iofile(id="data.raw", path="data.json")],
        )
        gather_stage = _stage(
            id="summary",
            inputs=[GatherInput(gather="dataset")],
            outputs=[_iofile(id="summary.report", path="report.json")],
        )
        benchmark = _benchmark(stages=[data_stage, gather_stage])
        assert benchmark is not None

    def test_multiple_providers_for_same_label(self):
        """Multiple stages providing the same label should work."""
        methods1 = _stage(
            id="methods_fast",
            modules=[_module(id="M1")],
            provides=["method"],
            outputs=[_iofile(id="methods_fast.result", path="result.json")],
        )
        methods2 = _stage(
            id="methods_accurate",
            modules=[_module(id="M2")],
            provides=["method"],
            outputs=[_iofile(id="methods_accurate.result", path="result.json")],
        )
        gather_stage = _stage(
            id="summary",
            modules=[_module(id="S1")],
            inputs=[GatherInput(gather="method")],
            outputs=[_iofile(id="summary.report", path="report.json")],
        )
        benchmark = _benchmark(stages=[methods1, methods2, gather_stage])
        assert benchmark is not None
        assert len(benchmark.get_provider_stages("method")) == 2

    def test_provides_without_gather_is_fine(self):
        """A stage can provide a label even if nothing gathers it."""
        stage = _stage(
            id="data",
            provides=["dataset"],
            outputs=[_iofile(id="data.raw", path="data.json")],
        )
        benchmark = _benchmark(stages=[stage])
        assert benchmark is not None
        assert benchmark.get_provider_stages("dataset") == [stage]

    def test_get_provider_stages_returns_empty_for_unknown(self):
        """Querying an unknown label returns empty list."""
        benchmark = _benchmark(stages=[_stage(id="data")])
        assert benchmark.get_provider_stages("nonexistent") == []


# ============================================================================
# Benchmark-level: provides/gather coexists with regular stages
# ============================================================================


@pytest.mark.short
class TestGatherCoexistsWithRegularStages:
    def test_mixed_regular_and_gather_stages(self):
        """Benchmark with both regular map stages and a gather stage."""
        data_stage = _stage(
            id="data",
            provides=["dataset"],
            outputs=[_iofile(id="data.raw", path="data.json")],
        )
        methods_stage = _stage(
            id="methods",
            modules=[_module(id="M1")],
            inputs=[InputCollection(entries=["data.raw"])],
            provides=["method"],
            outputs=[_iofile(id="methods.result", path="result.json")],
        )
        gather_stage = _stage(
            id="summary",
            modules=[_module(id="S1")],
            inputs=[GatherInput(gather="method")],
            outputs=[_iofile(id="summary.report", path="report.json")],
        )
        benchmark = _benchmark(stages=[data_stage, methods_stage, gather_stage])
        assert len(benchmark.stages) == 3
        assert not data_stage.is_gather_stage()
        assert not methods_stage.is_gather_stage()
        assert gather_stage.is_gather_stage()

    def test_regular_inputs_still_validated(self):
        """Regular input references to nonexistent outputs still fail."""
        stage = _stage(
            id="methods",
            inputs=[InputCollection(entries=["nonexistent.output"])],
        )
        from omnibenchmark.model.validation import ValidationError as OBValidationError

        with pytest.raises(OBValidationError, match="not valid"):
            _benchmark(stages=[stage])

    def test_get_stage_implicit_inputs_skips_gather(self):
        """get_stage_implicit_inputs should skip GatherInput items."""
        data_stage = _stage(
            id="data",
            provides=["dataset"],
            outputs=[_iofile(id="data.raw", path="data.json")],
        )
        gather_stage = _stage(
            id="summary",
            inputs=[GatherInput(gather="dataset")],
            outputs=[_iofile(id="summary.report", path="report.json")],
        )
        benchmark = _benchmark(stages=[data_stage, gather_stage])
        # Gather stage has no implicit (regular) inputs
        assert benchmark.get_stage_implicit_inputs("summary") == []
        # Data stage has no inputs at all
        assert benchmark.get_stage_implicit_inputs("data") == []


# ============================================================================
# YAML round-trip: from_yaml with gather syntax
# ============================================================================


@pytest.mark.short
class TestGatherFromYAML:
    def test_parse_gather_from_yaml_string(self):
        yaml_content = """\
id: test_benchmark
benchmarker: tester
version: "1.0"
software_backend: host
software_environments:
  host:
    description: "Host"
stages:
  - id: data
    provides:
      - dataset
    modules:
      - id: D1
        software_environment: host
        repository:
          url: https://example.com/repo.git
          commit: abc123
    outputs:
      - id: data.raw
        path: data.json
  - id: summary
    modules:
      - id: S1
        software_environment: host
        repository:
          url: https://example.com/repo.git
          commit: abc123
    inputs:
      - gather: dataset
    outputs:
      - id: summary.report
        path: report.json
"""
        benchmark = Benchmark.from_yaml(yaml_content)
        assert len(benchmark.stages) == 2

        data_stage = benchmark.stages[0]
        assert data_stage.provides == ["dataset"]
        assert not data_stage.is_gather_stage()

        summary_stage = benchmark.stages[1]
        assert summary_stage.is_gather_stage()
        assert summary_stage.get_gather_labels() == ["dataset"]

    def test_parse_provides_only_from_yaml(self):
        """Stage with provides but nothing gathering it."""
        yaml_content = """\
id: test_benchmark
benchmarker: tester
version: "1.0"
software_backend: host
software_environments:
  host:
    description: "Host"
stages:
  - id: data
    provides:
      - dataset
    modules:
      - id: D1
        software_environment: host
        repository:
          url: https://example.com/repo.git
          commit: abc123
    outputs:
      - id: data.raw
        path: data.json
"""
        benchmark = Benchmark.from_yaml(yaml_content)
        assert benchmark.stages[0].provides == ["dataset"]
        assert benchmark.get_provider_stages("dataset")[0].id == "data"

    def test_parse_multiple_providers_from_yaml(self):
        yaml_content = """\
id: test_benchmark
benchmarker: tester
version: "1.0"
software_backend: host
software_environments:
  host:
    description: "Host"
stages:
  - id: methods_fast
    provides:
      - method
    modules:
      - id: M1
        software_environment: host
        repository:
          url: https://example.com/repo.git
          commit: abc123
    outputs:
      - id: methods_fast.result
        path: fast_result.json
  - id: methods_slow
    provides:
      - method
    modules:
      - id: M2
        software_environment: host
        repository:
          url: https://example.com/repo.git
          commit: abc123
    outputs:
      - id: methods_slow.result
        path: slow_result.json
  - id: summary
    modules:
      - id: S1
        software_environment: host
        repository:
          url: https://example.com/repo.git
          commit: abc123
    inputs:
      - gather: method
    outputs:
      - id: summary.report
        path: report.json
"""
        benchmark = Benchmark.from_yaml(yaml_content)
        providers = benchmark.get_provider_stages("method")
        assert len(providers) == 2
        assert {s.id for s in providers} == {"methods_fast", "methods_slow"}
        assert benchmark.stages[2].is_gather_stage()
