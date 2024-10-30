from pathlib import Path

from omni.benchmark import Benchmark
from tests.cli.cli_setup import OmniCLISetup

benchmark_data = Path("..") / "data"
benchmark_data_path = Path(__file__).parent / benchmark_data


def test_benchmark_computational_plot():
    benchmark = Benchmark(benchmark_data_path / "Benchmark_001.yaml")
    expected_output = benchmark.export_to_dot().to_string()
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "info",
                "computational",
                "-b",
                str(benchmark_data_path / "Benchmark_001.yaml"),
            ],
            input="y",
        )

        assert result.exit_code == 0
        assert result.output.strip("\n") == expected_output.strip("\n")


def test_benchmark_topology_plot():
    benchmark = Benchmark(benchmark_data_path / "Benchmark_001.yaml")
    expected_output = benchmark.export_to_mermaid()
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "info",
                "topology",
                "-b",
                str(benchmark_data_path / "Benchmark_001.yaml"),
            ],
            input="y",
        )

        assert result.exit_code == 0
        assert result.output.strip("\n") == expected_output.strip("\n")
