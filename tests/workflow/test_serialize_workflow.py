import os.path
from pathlib import Path

from tests.workflow.Snakemake_setup import SnakemakeSetup


def test_serialize_workflow_001():
    benchmark_file = Path("..") / "data" / "Benchmark_001.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    with SnakemakeSetup(benchmark_file_path) as setup:
        benchmark = setup.benchmark
        assert benchmark.get_benchmark_name() == "Benchmark_001"

        path = setup.workflow.serialize_workflow(benchmark)
        assert os.path.exists(path)


def test_serialize_workflow_002():
    benchmark_file = Path("..") / "data" / "Benchmark_002.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    with SnakemakeSetup(benchmark_file_path) as setup:
        benchmark = setup.benchmark
        assert benchmark.get_benchmark_name() == "peiying_with_easyconfigs"

        path = setup.workflow.serialize_workflow(benchmark)
        assert os.path.exists(path)
