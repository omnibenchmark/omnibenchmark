from pathlib import Path

import pytest

from tests.workflow.Snakemake_setup import SnakemakeSetup


def test_run_workflow_001():
    benchmark_file = Path("..") / "data" / "Benchmark_001.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    with SnakemakeSetup(benchmark_file_path) as setup:
        benchmark = setup.benchmark
        assert benchmark.get_benchmark_name() == "Benchmark_001"

        # First run the whole workflow
        success = setup.workflow.run_workflow(benchmark)

        assert success


@pytest.mark.skip(reason="Benchmark_002 is not working properly yet")
def test_run_workflow_002():
    benchmark_file = Path("..") / "data" / "Benchmark_002.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    with SnakemakeSetup(benchmark_file_path) as setup:
        benchmark = setup.benchmark
        assert benchmark.get_benchmark_name() == "peiying_with_easyconfigs"

        success = setup.workflow.run_workflow(benchmark)
        assert success
