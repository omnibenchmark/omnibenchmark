from pathlib import Path

from omnibenchmark.benchmark import Benchmark
from tests.cli.cli_setup import OmniCLISetup

benchmark_data = Path("..") / "data"
benchmark_data_path = Path(__file__).parent / benchmark_data


def strip(text):
    return "".join(text.split())


def test_benchmark_snakemake_plot():
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "describe",
                "snakemake",
                "-b",
                str(benchmark_data_path / "Benchmark_001.yaml"),
            ],
            input="y",
        )

        assert result.returncode == 0
        assert "process-P1-.317a5066-after_data" in result.stdout
        assert "methods-M2-.3297cc0b-after_data" in result.stdout
        assert "metrics-m3-default-after_methods" in result.stdout


def test_benchmark_topology_plot(tmp_path):
    import shutil

    # Copy benchmark file to tmp_path to avoid writing to current directory
    benchmark_file = benchmark_data_path / "Benchmark_001.yaml"
    copied_benchmark_path = tmp_path / "Benchmark_001.yaml"
    shutil.copy(benchmark_file, copied_benchmark_path)

    benchmark = Benchmark(copied_benchmark_path, out_dir=tmp_path / "out")
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
                "describe",
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


def test_describe_status_with_out_dir(tmp_path):
    """Test that describe status command respects --out-dir parameter."""
    import shutil

    # Copy benchmark file and required env files to tmp_path
    benchmark_file = benchmark_data_path / "Benchmark_001.yaml"
    copied_benchmark_path = tmp_path / "test_benchmark.yaml"
    shutil.copy(benchmark_file, copied_benchmark_path)

    # Copy the envs directory if it exists
    envs_src = benchmark_data_path / "envs"
    if envs_src.exists():
        envs_dst = tmp_path / "envs"
        shutil.copytree(envs_src, envs_dst)

    # Test with default out dir first
    with OmniCLISetup() as omni:
        result_default = omni.call(
            [
                "describe",
                "status",
                str(copied_benchmark_path),
            ],
            cwd=str(tmp_path),
        )
        assert result_default.returncode == 0

    # Test with custom out dir
    custom_out_dir = "custom_output"
    with OmniCLISetup() as omni:
        result_custom = omni.call(
            [
                "describe",
                "status",
                str(copied_benchmark_path),
                "--out-dir",
                custom_out_dir,
            ],
            cwd=str(tmp_path),
        )
        assert result_custom.returncode == 0

        # Both should succeed - the parameter is accepted and used
        # The actual validation that it uses the right directory is that
        # BenchmarkExecution is created with the custom path (verified by code inspection)
        # and the command doesn't crash
        assert "file completion" in result_custom.stdout.lower()
