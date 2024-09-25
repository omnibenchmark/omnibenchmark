import re
from pathlib import Path

import pytest
import yaml

from omni.benchmark import Benchmark


def test_plot_computational_graph():
    benchmark_file = "../data/Benchmark_001.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file
    output_file = "test_plot_computational_graph.png"
    output_file_path = Path(output_file)

    try:
        with open(benchmark_file_path, "r") as file:
            yaml.safe_load(file)
            benchmark = Benchmark(Path(benchmark_file))
            benchmark.plot_computational_graph(output_file=output_file)
            assert (
                output_file_path.exists()
            ), f"Output file {output_file} was not created."
    except Exception as e:
        pytest.fail(f"Plotting benchmark computational graph failed: {e}")
    finally:
        if output_file_path.exists():
            output_file_path.unlink()


def test_export_to_mermaid():
    benchmark_file = "../data/Benchmark_001.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file
    expected_mermaid = """---
    title: Benchmark_001
    ---
    flowchart LR
        \tclassDef param fill:#f96
        \tsubgraph data
            \t\tD1
            \t\tD2
        \tend
        \tsubgraph process
            \t\tP1
            \t\tD1 --> P1
            \t\tD2 --> P1
            \t\tP2
            \t\tD1 --> P2
            \t\tD2 --> P2
        \tend
        \tsubgraph methods
            \t\tM1
            \t\tD1 --> M1
            \t\tD2 --> M1
            \t\tP1 --> M1
            \t\tP2 --> M1
            \t\tM2
            \t\tD1 --> M2
            \t\tD2 --> M2
            \t\tP1 --> M2
            \t\tP2 --> M2
        \tend
        \tsubgraph metrics
            \t\tm1
            \t\tM1 --> m1
            \t\tM2 --> m1
            \t\tm2
            \t\tM1 --> m2
            \t\tM2 --> m2
            \t\tm3
            \t\tM1 --> m3
            \t\tM2 --> m3
        \tend
    """

    try:
        with open(benchmark_file_path, "r") as file:
            yaml.safe_load(file)
            benchmark = Benchmark(Path(benchmark_file))
            mermaid = benchmark.export_to_mermaid()
            assert clean(mermaid).startswith(clean(expected_mermaid))
    except Exception as e:
        pytest.fail(f"Exporting benchmark topology to mermaid failed: {e}")


def clean(output: str) -> str:
    output = output.strip()
    output = output.replace("    ", "")

    # Replace different newline characters with a single '\n'
    normalized_output = re.sub(r"\r\n|\r", "\n", output)

    # Replace multiple spaces and tabs with a single space
    normalized_output = re.sub(r"[ \t]+", " ", normalized_output)

    return normalized_output
