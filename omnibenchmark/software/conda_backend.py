from pathlib import Path

import pytest

from snakedeploy import conda as snakeconda

"""
Conda-based software management, per benchmark yaml (and not per omniblock)
"""


@pytest.mark.conda
def pin_conda_envs(benchmark_yaml: Path) -> None:
    """
    Pins a conda env file using snakedeploy
    """
    snakeconda.pin_conda_envs(
        conda_env_paths=[benchmark_yaml], conda_frontend="conda", warn_on_error=True
    )
