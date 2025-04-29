import os.path
import shutil
from pathlib import Path

from tests.workflow.Snakemake_setup import SnakemakeSetup

from .path import data


def pytest_addoption(parser):
    parser.addoption(
        "--keep-files",
        action="store_true",
        default=False,
        help="Keep temporary files after test execution",
    )
    parser.addoption(
        "--current-dir",
        action="store_true",
        default=False,
        help="Use current directory instead of temporary one",
    )


def test_run_node_workflow_001(snakemake_env, tmp_path):
    # Use current directory if specified, otherwise use tmp_path
    # tmp_path is already provided by pytest fixture if not using current dir
    if snakemake_env["current_dir"]:
        tmp_path = Path(os.getcwd())

    benchmark_file = data / "Benchmark_001.yaml"
    keep_files = snakemake_env["keep_files"]

    with SnakemakeSetup(
        benchmark_file, keep_files=keep_files, cwd=tmp_path.as_posix()
    ) as setup:
        benchmark = setup.benchmark
        assert benchmark.get_benchmark_name() == "Benchmark_001"

        benchmark_nodes = benchmark.get_nodes()
        assert len(benchmark_nodes) == 15

        # First run the whole workflow
        assert setup.workflow.run_workflow(
            benchmark,
            work_dir=tmp_path,
        )

        # Delete output of process P1 directory
        process_P1_output = (
            tmp_path
            / "out"
            / "data"
            / "D1"
            / "default"
            / "process"
            / "P1"
            / ".317a506603d7cb7f079fcc6a38cdf99e3955e1729540d38b9b0f36bd7c16d2a3"
        )

        assert os.path.exists(process_P1_output)

        shutil.rmtree(process_P1_output, onerror=setup._remove_readonly)

        # Then run the workflow for the missing computational node
        benchmark_node_3 = benchmark_nodes[3]

        input_dir = tmp_path / "out" / "data" / "D1" / "default"
        success = setup.workflow.run_node_workflow(
            benchmark_node_3, input_dir=input_dir, dataset="D1"
        )

        assert success is True
