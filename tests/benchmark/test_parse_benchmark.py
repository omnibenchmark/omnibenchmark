from pathlib import Path

import pytest

from omnibenchmark.model import Benchmark


@pytest.mark.short
def test_parse_benchmark():
    benchmark_file = "../data/Benchmark_001.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    try:
        benchmark = Benchmark.from_yaml(benchmark_file_path)
        # Verify basic properties
        assert benchmark is not None
        assert isinstance(benchmark, Benchmark)
        assert hasattr(benchmark, "id")
        assert hasattr(benchmark, "stages")
        assert hasattr(benchmark, "software_environments")
    except Exception as e:
        pytest.fail(f"Parsing benchmark model failed: {e}")


# --- null collections report a YAML hint instead of a bare TypeError (#391) ---


@pytest.mark.short
def test_null_modules_reports_a_yaml_hint(tmp_path):
    """A stage with `modules:` (no value) names the stage and the fix."""
    from omnibenchmark.model.validation import BenchmarkParseError

    yaml_file = tmp_path / "bench.yaml"
    yaml_file.write_text(
        "id: bench\nstages:\n  - id: clustering\n    modules:\n",
        encoding="utf-8",
    )
    with pytest.raises(BenchmarkParseError) as exc_info:
        Benchmark.from_yaml(yaml_file)
    text = str(exc_info.value)
    assert "clustering" in text
    assert "modules" in text
    assert "null" in text
    # The hint shows what to write instead.
    assert "- id: my_module" in text
    # From a file, the error points at the offending line.
    assert exc_info.value.line_number is not None
    assert str(yaml_file) in text


@pytest.mark.short
def test_null_stages_reports_a_yaml_hint():
    """`stages:` with no value names the key and the fix."""
    from omnibenchmark.model.validation import BenchmarkParseError

    with pytest.raises(BenchmarkParseError) as exc_info:
        Benchmark.from_yaml("id: bench\nstages:\n")
    text = str(exc_info.value)
    assert "stages" in text
    assert "null" in text
    assert "- id: my_stage" in text


@pytest.mark.short
def test_null_stage_entry_reports_a_yaml_hint(tmp_path):
    """A `null` entry inside `stages:` names its index."""
    from omnibenchmark.model.validation import BenchmarkParseError

    yaml_file = tmp_path / "bench.yaml"
    yaml_file.write_text("id: bench\nstages:\n  -\n", encoding="utf-8")
    with pytest.raises(BenchmarkParseError) as exc_info:
        Benchmark.from_yaml(yaml_file)
    assert "stages[0]" in str(exc_info.value)


@pytest.mark.short
def test_empty_yaml_reports_a_yaml_hint():
    """An empty YAML document is rejected with a message, not a TypeError."""
    from omnibenchmark.model.validation import BenchmarkParseError

    with pytest.raises(BenchmarkParseError) as exc_info:
        Benchmark.from_yaml("# only a comment, no document\n")
    assert "no document" in str(exc_info.value)
