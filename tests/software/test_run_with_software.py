from pathlib import Path

from tests.workflow.Snakemake_setup import SnakemakeSetup

benchmark_data = Path("..") / "data"
benchmark_data_path = Path(__file__).parent / benchmark_data


def test_run_benchmark_with_software():
    benchmark_file = benchmark_data_path / "mock_benchmark.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    with SnakemakeSetup(benchmark_file_path) as setup:
        benchmark = setup.benchmark
        assert benchmark.get_benchmark_name() == "mock_benchmark"

        backend = benchmark.get_benchmark_software_backend()

        # First run the whole workflow
        success = setup.workflow.run_workflow(benchmark, backend=backend)

        assert success
