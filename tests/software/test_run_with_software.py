import os
import shutil
import stat
import sys
from pathlib import Path
import os.path as op

from omni.benchmark import Benchmark
from omni.workflow.snakemake import SnakemakeEngine

sys.path.insert(0, op.dirname(__file__))

benchmark_data = Path("..") / "data"
benchmark_data_path = Path(__file__).parent / benchmark_data


def test_run_benchmark_with_software():
    benchmark_file = benchmark_data_path / "mock_benchmark_with_software.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    if os.path.exists(benchmark_file):
        benchmark = Benchmark(benchmark_file_path)
        workflow = SnakemakeEngine()

        backend = benchmark.get_benchmark_software_backend()

        success = workflow.run_workflow(benchmark, backend=backend)

        assert success

        _cleanup_snakemake()

    else:
        raise FileNotFoundError(f"Benchmark file {benchmark_file} does not exist.")


def _cleanup_snakemake():
    current_dir = os.getcwd()
    for file in [".snakemake", "out", "Snakefile", "snakemake.log"]:
        file_path = os.path.join(current_dir, file)
        if os.path.exists(file_path):
            if os.path.isfile(file_path):
                os.remove(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path, onerror=_remove_readonly)


def _remove_readonly(func, path, _):
    """Clear the readonly bit and reattempt the removal"""
    os.chmod(path, stat.S_IWRITE)
    func(path)
