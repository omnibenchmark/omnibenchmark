#!/usr/bin/env python
##
## Izaskun Mallona
## 22th July 2024

import subprocess
import shutil
from pathlib import Path


def check_cmd_zero_exit(cmd_string):
    completed_process = subprocess.run(cmd_string.split(" "))
    assert completed_process.returncode == 0


def _setup_test_environment(snakefile_path, tmp_path):
    """Copy Snakefile and related files to temporary directory"""
    snakefile_path = Path(snakefile_path)
    test_dir = snakefile_path.parent

    # Copy the entire test directory to tmp_path to preserve structure
    temp_test_dir = tmp_path / test_dir.name
    shutil.copytree(test_dir, temp_test_dir)

    return temp_test_dir / "Snakefile"


def run_snmk_conda(snakefile_path, expected_path, tmp_path):
    temp_snakefile = _setup_test_environment(snakefile_path, tmp_path)
    temp_dir = temp_snakefile.parent

    completed_process = subprocess.run(
        [
            "snakemake",
            "-s",
            temp_snakefile.name,
            "--use-conda",
            "--cores",
            "1",
            "--directory",
            str(temp_dir),
        ]
    )
    assert completed_process.returncode == 0

    # Check output against expected
    produced_file = temp_dir / "test0.out"
    with open(expected_path, "r") as exp_fh, open(produced_file, "r") as prod_fh:
        assert prod_fh.read() == exp_fh.read()


def run_snmk_apptainer(snakefile_path, expected_path, tmp_path):
    temp_snakefile = _setup_test_environment(snakefile_path, tmp_path)
    temp_dir = temp_snakefile.parent

    completed_process = subprocess.run(
        [
            "snakemake",
            "-s",
            temp_snakefile.name,
            "--use-singularity",
            "--cores",
            "1",
            "--directory",
            str(temp_dir),
        ]
    )
    assert completed_process.returncode == 0

    # Check output against expected
    produced_file = temp_dir / "test0.out"
    with open(expected_path, "r") as exp_fh, open(produced_file, "r") as prod_fh:
        assert prod_fh.read() == exp_fh.read()


def run_snmk_envmodules(snakefile_path, expected_path, tmp_path):
    temp_snakefile = _setup_test_environment(snakefile_path, tmp_path)
    temp_dir = temp_snakefile.parent

    completed_process = subprocess.run(
        [
            "snakemake",
            "-s",
            temp_snakefile.name,
            "--use-envmodules",
            "--cores",
            "1",
            "--directory",
            str(temp_dir),
        ]
    )
    assert completed_process.returncode == 0

    # Check output against expected
    produced_file = temp_dir / "test0.out"
    with open(expected_path, "r") as exp_fh, open(produced_file, "r") as prod_fh:
        assert prod_fh.read() == exp_fh.read()


def run(snakefile_path, expected_path, method, tmp_path):
    if method == "conda":
        run_snmk_conda(snakefile_path, expected_path, tmp_path)
    elif method == "apptainer":
        run_snmk_apptainer(snakefile_path, expected_path, tmp_path)
    elif method == "envmodules":
        run_snmk_envmodules(snakefile_path, expected_path, tmp_path)
    else:
        print("not implemented")
