from pathlib import Path

from omni_schema.datamodel.omni_schema import SoftwareBackendEnum
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


def test_run_workflow_backends_missing():
    benchmark_file = Path("..") / "data" / "benchmark_some_backends_missing.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    with SnakemakeSetup(benchmark_file_path) as setup:
        benchmark = setup.benchmark
        assert benchmark.get_benchmark_name() == "some_backends_missing"

        # First run the whole workflow
        success = setup.workflow.run_workflow(
            benchmark, backend=SoftwareBackendEnum.conda
        )

        assert success
