import os
import sys
from pathlib import Path
import os.path as op

import pytest
from omni_schema.datamodel.omni_schema import SoftwareBackendEnum

from tests.workflow.Snakemake_setup import SnakemakeSetup

sys.path.insert(0, op.dirname(__file__))

benchmark_data = Path("..") / "data"
benchmark_data_path = Path(__file__).parent / benchmark_data


def test_run_benchmark_with_software_envmodules():
    benchmark_file = benchmark_data_path / "mock_benchmark_with_software.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    with SnakemakeSetup(benchmark_file_path) as setup:
        benchmark = setup.benchmark

        success = setup.workflow.run_workflow(
            benchmark,
            backend=SoftwareBackendEnum.envmodules,
            modulepath=os.path.join(
                os.environ["GITHUB_WORKSPACE"], "tests", "data", "envs"
            ),
        )

        assert success


def test_run_benchmark_with_software_conda():
    benchmark_file = benchmark_data_path / "mock_benchmark_with_software.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    with SnakemakeSetup(benchmark_file_path) as setup:
        benchmark = setup.benchmark

        success = setup.workflow.run_workflow(
            benchmark,
            backend=SoftwareBackendEnum.conda,
        )

        assert success


@pytest.mark.skip(reason="Apptainer image is not available yet")
def test_run_benchmark_with_software_apptainer():
    benchmark_file = benchmark_data_path / "mock_benchmark_with_software.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    with SnakemakeSetup(benchmark_file_path) as setup:
        benchmark = setup.benchmark

        success = setup.workflow.run_workflow(
            benchmark,
            backend=SoftwareBackendEnum.apptainer,
        )

        assert success
