import os
import shutil
import functools

from pathlib import Path

import pytest

from omnibenchmark.model import SoftwareBackendEnum

from omnibenchmark.benchmark import Benchmark
from omnibenchmark.workflow.snakemake import SnakemakeEngine


def skip_if_no_conda(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        if shutil.which("conda") is None:
            pytest.skip("conda not found in system path")
        return func(*args, **kwargs)

    return wrapper


def skip_if_no_envmodules(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        if shutil.which("module") is None:
            pytest.skip("envmodules not found in system path")
        return func(*args, **kwargs)

    return wrapper


benchmark_data = Path("..") / "data"
benchmark_data_path = Path(__file__).parent / benchmark_data


# TODO: use tmp_path, no need to cleanup
TO_CLEANUP = [".snakemake", "out", "Snakefile", "snakemake.log"]

# TODO(ben): reuse flags in tests/io to inject --current-dir and --keep-files


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


@skip_if_no_envmodules
def test_run_benchmark_with_software_envmodules():
    benchmark_file = benchmark_data_path / "Clustering.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    if os.path.exists(benchmark_file):
        benchmark = Benchmark(benchmark_file_path)
        workflow = SnakemakeEngine()

        if os.environ.get("GITHUB_WORKSPACE", None):
            override_module_path = os.path.join(
                os.environ["GITHUB_WORKSPACE"],
                "tests",
                "data",
                "envs",
            )
        else:
            override_module_path = os.environ.get("MODULEPATH", None)

        success = workflow.run_workflow(
            benchmark,
            backend=SoftwareBackendEnum.envmodules,
            module_path=override_module_path,
        )

        assert success

        # TODO: use tmp_path
        _cleanup_snakemake()
    else:
        raise FileNotFoundError(f"Benchmark file {benchmark_file} does not exist.")


@skip_if_no_conda
@pytest.mark.conda
def test_run_benchmark_with_software_conda():
    benchmark_file = benchmark_data_path / "Clustering.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    if os.path.exists(benchmark_file):
        benchmark = Benchmark(benchmark_file_path)
        workflow = SnakemakeEngine()

        success = workflow.run_workflow(
            benchmark,
            backend=SoftwareBackendEnum.conda,
        )

        assert success

        # TODO: use tmp_path
        _cleanup_snakemake()
    else:
        raise FileNotFoundError(f"Benchmark file {benchmark_file} does not exist.")


@pytest.mark.skip(reason="Apptainer image is not available yet")
def test_run_benchmark_with_software_apptainer():
    benchmark_file = benchmark_data_path / "Clustering.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    if os.path.exists(benchmark_file):
        benchmark = Benchmark(benchmark_file_path)
        workflow = SnakemakeEngine()

        success = workflow.run_workflow(
            benchmark,
            backend=SoftwareBackendEnum.apptainer,
        )

        assert success

        # TODO: use tmp_path
        _cleanup_snakemake()
    else:
        raise FileNotFoundError(f"Benchmark file {benchmark_file} does not exist.")
