#!/usr/bin/env python

from snakedeploy import conda as snakeconda
from pathlib import Path

"""
Conda-based software management, per benchmark yaml (and not per omniblock)
"""


def pin_conda_envs(benchmark_yaml: Path) -> None:
    """
    Pins a conda env file using snakedeploy
    """
    snakeconda.pin_conda_envs(
        conda_env_paths=[benchmark_yaml], conda_frontend="mamba", warn_on_error=False
    )