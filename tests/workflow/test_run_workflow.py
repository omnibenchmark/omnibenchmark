from pathlib import Path
import os

from omni_schema.datamodel.omni_schema import SoftwareBackendEnum

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


def test_run_workflow_001(snakemake_env, tmp_path):
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

        # Run the whole workflow
        success = setup.workflow.run_workflow(
            benchmark,
            work_dir=tmp_path,
        )

        assert success


def test_run_workflow_backends_missing(snakemake_env, tmp_path):
    # Use current directory if specified, otherwise use tmp_path
    # tmp_path is already provided by pytest fixture if not using current dir
    if snakemake_env["current_dir"]:
        tmp_path = Path(os.getcwd())

    benchmark_file = data / "benchmark_some_backends_missing.yaml"
    keep_files = snakemake_env["keep_files"]

    with SnakemakeSetup(
        benchmark_file, keep_files=keep_files, cwd=tmp_path.as_posix()
    ) as setup:
        benchmark = setup.benchmark
        assert benchmark.get_benchmark_name() == "some_backends_missing"

        # Run the whole workflow
        success = setup.workflow.run_workflow(
            benchmark,
            backend=SoftwareBackendEnum.conda,
            work_dir=tmp_path,
        )

        assert success is True
