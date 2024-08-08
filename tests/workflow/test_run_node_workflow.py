import os.path
import shutil
from pathlib import Path

import pytest

from tests.workflow.Snakemake_setup import SnakemakeSetup


def test_run_node_workflow_001():
    benchmark_file = Path("..") / "data" / "Benchmark_001.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    with SnakemakeSetup(benchmark_file_path) as setup:
        benchmark = setup.benchmark
        assert benchmark.get_benchmark_name() == "Benchmark_001"

        benchmark_nodes = benchmark.get_nodes()
        assert len(benchmark_nodes) == 15

        # First run the whole workflow
        assert setup.workflow.run_workflow(benchmark)

        # Delete output of process P1 directory
        process_P1_output = (
            Path(os.getcwd())
            / "out"
            / "data"
            / "D1"
            / "default"
            / "process"
            / "P1"
            / "param_1"
        )
        assert os.path.exists(process_P1_output)
        shutil.rmtree(process_P1_output, onerror=setup._remove_readonly)

        # Then run the workflow for the missing computational node
        benchmark_node_3 = benchmark_nodes[3]
        input_dir = Path(os.getcwd()) / "out" / "data" / "D1" / "default"
        success = setup.workflow.run_node_workflow(
            benchmark_node_3, input_dir=input_dir, dataset="D1"
        )

        assert success == 1


@pytest.mark.skip(reason="Benchmark_002 is not working properly yet")
def test_run_node_workflow_002():
    benchmark_file = Path("..") / "data" / "Benchmark_002.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    with SnakemakeSetup(benchmark_file_path) as setup:
        benchmark = setup.benchmark
        assert benchmark.get_benchmark_name() == "peiying_with_easyconfigs"

        benchmark_nodes = benchmark.get_nodes()
        assert len(benchmark_nodes) == 8

        benchmark_node_3 = benchmark_nodes[3]
        success = setup.workflow.run_node_workflow(benchmark_node_3)
        assert success == 1
