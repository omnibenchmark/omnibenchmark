"""
Archive components handling for different parts of a benchmark.

This module provides functions to prepare different components of a benchmark
for archiving: configuration, code, software environments, and results.
Each component can be included or excluded independently.
"""

import os
from pathlib import Path
from typing import List

from omnibenchmark.benchmark import BenchmarkExecution
from omnibenchmark.model import SoftwareBackendEnum
from omnibenchmark.model.benchmark import _is_environment_url
from omnibenchmark.remote.RemoteStorage import StorageOptions
from omnibenchmark.remote.versioning import get_expected_benchmark_output_files


def prepare_archive_config(benchmark: BenchmarkExecution) -> List[Path]:
    """
    Prepare the configuration files to archive.

    Args:
        benchmark: The benchmark execution object

    Returns:
        List[Path]: The filenames of config files to archive
    """
    config_files = []

    # Always include the main benchmark definition file
    definition_file = benchmark.get_definition_file()
    if definition_file and definition_file.is_file():
        config_files.append(definition_file)

    return config_files


def prepare_archive_code(benchmark: BenchmarkExecution) -> List[Path]:
    """
    Prepare the code files to archive and return list of all filenames.

    Args:
        benchmark: The benchmark execution object

    Returns:
        List[Path]: The filenames of all code to archive
    """
    from omnibenchmark.git.clone import clone_module

    nodes = benchmark.get_nodes()
    repositories = set()
    for node in nodes:
        repo = node.get_repository()
        repositories.add((repo.url, repo.commit))

    files = []

    while repositories:
        repo = repositories.pop()
        try:
            repo_path = clone_module(repo[0], repo[1])
            if repo_path.exists():
                files += list(repo_path.rglob("*"))
        except Exception:
            # Skip repositories that can't be cloned
            continue

    # Filter to only include actual files (not directories)
    return [f for f in files if f.is_file()]


def prepare_archive_software(benchmark: BenchmarkExecution) -> List[Path]:
    """
    Prepare the software files to archive and return list of all filenames.

    Args:
        benchmark: The benchmark execution object

    Returns:
        List[Path]: The filenames of all software to archive
    """
    which_env = benchmark.get_benchmark_software_backend()
    files = []

    if which_env != SoftwareBackendEnum.host:
        if which_env == SoftwareBackendEnum.envmodules:
            files += prepare_archive_software_easyconfig(benchmark)
        elif which_env == SoftwareBackendEnum.apptainer:
            files += prepare_archive_software_apptainer(benchmark)
        elif which_env == SoftwareBackendEnum.conda:
            files += prepare_archive_software_conda(benchmark)
        else:
            # For unknown software backends, return empty list
            pass

    return files


def prepare_archive_software_easyconfig(benchmark: BenchmarkExecution) -> List[Path]:
    """
    Prepare the software easyconfig files to archive.

    Args:
        benchmark: The benchmark execution object

    Returns:
        List[Path]: The filenames of all easyconfig files to archive
    """
    files = []
    softenvs = benchmark.get_benchmark_software_environments()

    for softenv in softenvs.values():
        # Handle envmodule files
        if hasattr(softenv, "envmodule") and softenv.envmodule is not None:
            envmodule_file = (
                benchmark.context.directory / Path(softenv.envmodule)
            ).relative_to(Path(os.getcwd()))
            if envmodule_file.is_file():
                files.append(envmodule_file)

        # Handle easyconfig files
        if hasattr(softenv, "easyconfig") and softenv.easyconfig is not None:
            easyconfig_file = (
                benchmark.context.directory / Path(softenv.easyconfig)
            ).relative_to(Path(os.getcwd()))
            if easyconfig_file.is_file():
                files.append(easyconfig_file)

    return files


def prepare_archive_software_conda(benchmark: BenchmarkExecution) -> List[Path]:
    """
    Prepare the conda environment files to archive.

    Args:
        benchmark: The benchmark execution object

    Returns:
        List[Path]: The filenames of all conda files to archive
    """
    files = []
    softenvs = benchmark.get_benchmark_software_environments()

    for softenv in softenvs.values():
        if hasattr(softenv, "conda") and softenv.conda is not None:
            conda_file = benchmark.context.directory / Path(softenv.conda)
            if conda_file.is_file():
                files.append(conda_file)

    return files


def prepare_archive_software_apptainer(benchmark: BenchmarkExecution) -> List[Path]:
    """
    Prepare the apptainer/singularity container files to archive.

    Args:
        benchmark: The benchmark execution object

    Returns:
        List[Path]: The filenames of all apptainer files to archive
    """
    files = []
    softenvs = benchmark.get_benchmark_software_environments()

    for softenv in softenvs.values():
        if hasattr(softenv, "apptainer") and softenv.apptainer is not None:
            # Skip URL-based images (oras://, http://, https://, docker://)
            if _is_environment_url(softenv.apptainer):
                continue

            apptainer_file = (
                benchmark.context.directory / Path(softenv.apptainer)
            ).relative_to(Path(os.getcwd()))
            if apptainer_file.is_file():
                files.append(apptainer_file)

    return files


def prepare_archive_results(
    benchmark: BenchmarkExecution, results_dir: str, remote_storage: bool = False
) -> List[Path]:
    """
    Prepare the results files to archive.

    This function can work in two modes:
    1. Local mode (remote_storage=False): Only includes existing local files
    2. Remote mode (remote_storage=True): Downloads from remote storage first

    Args:
        benchmark: The benchmark execution object
        results_dir: Directory containing results files
        remote_storage: Whether to use remote storage

    Returns:
        List[Path]: The filenames of all results to archive
    """
    if remote_storage:
        return prepare_archive_results_remote(benchmark, results_dir)
    else:
        return prepare_archive_results_local(benchmark, results_dir)


def prepare_archive_results_local(
    benchmark: BenchmarkExecution, results_dir: str
) -> List[Path]:
    """
    Prepare local results files for archiving (no remote storage required).

    Args:
        benchmark: The benchmark execution object
        results_dir: Directory containing results files

    Returns:
        List[Path]: The filenames of existing local files to archive
    """
    # Get expected files based on benchmark configuration
    storage_options = StorageOptions(out_dir=results_dir)
    expected_files = get_expected_benchmark_output_files(benchmark, storage_options)

    # Filter to only include files that actually exist locally
    existing_files = []
    for file_path in expected_files:
        path_obj = Path(file_path)
        if path_obj.is_file():
            existing_files.append(path_obj)

    # Also include any additional files in the results directory
    results_path = Path(results_dir)
    if results_path.exists() and results_path.is_dir():
        additional_files = [
            f
            for f in results_path.rglob("*")
            if f.is_file() and f not in existing_files
        ]
        existing_files.extend(additional_files)

    return existing_files


def prepare_archive_results_remote(
    benchmark: BenchmarkExecution, results_dir: str
) -> List[Path]:
    """
    Prepare remote results files for archiving (downloads from remote storage).

    Args:
        benchmark: The benchmark execution object
        results_dir: Directory containing results files

    Returns:
        List[Path]: The filenames of results downloaded from remote storage
    """
    # Import here to avoid circular imports
    from omnibenchmark.remote.files import list_files, download_files

    # Get list of files from remote storage
    objectnames, etags = list_files(
        benchmark.get_definition_file().as_posix(),
        type="all",
        stage="",
        module="",
        file_id="",
        remote_storage=True,
        storage_options=StorageOptions(out_dir=results_dir),
    )

    # Download files from remote storage
    download_files(
        benchmark.get_definition_file().as_posix(),
        type="all",
        stage="",
        module="",
        file_id="",
        overwrite=True,
    )

    return [Path(obj) for obj in objectnames]
