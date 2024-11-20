import os
import shutil
import sys
from pathlib import Path
import os.path as op

import pytest
from omni_schema.datamodel.omni_schema import SoftwareBackendEnum

from omni.benchmark import Benchmark
from omni.workflow.snakemake import SnakemakeEngine

sys.path.insert(0, op.dirname(__file__))

benchmark_data = Path("..") / "data"
benchmark_data_path = Path(__file__).parent / benchmark_data

TO_CLEANUP = [".snakemake", "out", "Snakefile", "snakemake.log"]


def _cleanup_snakemake():
    current_dir = os.getcwd()
    for file in TO_CLEANUP:
        file_path = os.path.join(current_dir, file)
        if os.path.exists(file_path):
            if os.path.isfile(file_path):
                os.remove(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path, onerror=_remove_readonly)


def _remove_readonly(func, path, _):
    """Clear the readonly bit and reattempt the removal"""
    os.chmod(path, os.stat.S_IWRITE)
    func(path)


def test_run_benchmark_with_software_envmodules():
    benchmark_file = benchmark_data_path / "mock_benchmark_with_software.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    if os.path.exists(benchmark_file):
        benchmark = Benchmark(benchmark_file_path)
        workflow = SnakemakeEngine()

        success = workflow.run_workflow(
            benchmark,
            backend=SoftwareBackendEnum.envmodules,
            modulepath=os.path.join(
                os.environ["GITHUB_WORKSPACE"], "tests", "data", "envs"
            ),
        )

        assert success

        _cleanup_snakemake()
    else:
        raise FileNotFoundError(f"Benchmark file {benchmark_file} does not exist.")


def test_run_benchmark_with_software_conda():
    benchmark_file = benchmark_data_path / "mock_benchmark_with_software.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    if os.path.exists(benchmark_file):
        benchmark = Benchmark(benchmark_file_path)
        workflow = SnakemakeEngine()

        success = workflow.run_workflow(
            benchmark,
            backend=SoftwareBackendEnum.conda,
        )

        assert success

        _cleanup_snakemake()
    else:
        raise FileNotFoundError(f"Benchmark file {benchmark_file} does not exist.")


@pytest.mark.skip(reason="Apptainer image is not available yet")
def test_run_benchmark_with_software_apptainer():
    benchmark_file = benchmark_data_path / "mock_benchmark_with_software.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    if os.path.exists(benchmark_file):
        benchmark = Benchmark(benchmark_file_path)
        workflow = SnakemakeEngine()

        success = workflow.run_workflow(
            benchmark,
            backend=SoftwareBackendEnum.apptainer,
        )

        assert success

        _cleanup_snakemake()
    else:
        raise FileNotFoundError(f"Benchmark file {benchmark_file} does not exist.")
