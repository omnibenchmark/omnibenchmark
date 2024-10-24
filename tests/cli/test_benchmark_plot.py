import os.path
from pathlib import Path

from tests.cli.cli_setup import OmniCLISetup

benchmark_data = Path("..") / "data"
benchmark_data_path = Path(__file__).parent / benchmark_data


def test_benchmark_plot():
    expected_output_plot = Path(os.getcwd()) / "Benchmark_001.png"
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "info",
                "plot",
                "-b",
                str(benchmark_data_path / "Benchmark_001.yaml"),
            ],
            input="y",
        )

        assert result.exit_code == 0
        assert os.path.exists(expected_output_plot)
        os.remove(expected_output_plot)
