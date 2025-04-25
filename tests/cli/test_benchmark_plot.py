from pathlib import Path

from omnibenchmark.benchmark import Benchmark
from tests.cli.cli_setup import OmniCLISetup

benchmark_data = Path("..") / "data"
benchmark_data_path = Path(__file__).parent / benchmark_data


def strip(text):
    return "".join(text.split())


def test_benchmark_computational_plot():
    # XXX what is this testing, exactly? that cli is calling the right method?
    # benchmark = Benchmark(benchmark_data_path / "Benchmark_001.yaml")
    # expected_output = benchmark.export_to_dot().to_string()
    # I get some serious non-determinism, so I'm not sure how this was working before.

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

        assert result.returncode == 0
        assert "process-P1-param_0-after_data" in result.stdout
        assert "methods-M2-param_1-after_data" in result.stdout
        assert "metrics-m3-default-after_methods" in result.stdout


def test_benchmark_topology_plot():
    benchmark = Benchmark(benchmark_data_path / "Benchmark_001.yaml")
    # XXX what is this testing, exactly? that cli is calling the right method?
    # for that we could use a mock object to test the export_to_mermaid method is called
    expected_output = benchmark.export_to_mermaid()

    expected_methods = """
       	subgraph methods
		M1
		D1 --> M1
		D2 --> M1
		P1 --> M1
		P2 --> M1
		M2
		D1 --> M2
		D2 --> M2
		P1 --> M2
		P2 --> M2
	end
	"""

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

        assert result.returncode == 0
        print(result.stdout)
        print(expected_output)
        assert strip(expected_methods) in strip(result.stdout)
        # the output cannot be matched exactly, there's randomness in the layout
        # assert strip(result.stdout) == strip(expected_output)
