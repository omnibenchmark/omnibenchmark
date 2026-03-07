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


def prepare_archive_code(
    benchmark: BenchmarkExecution,
) -> List[tuple[Path, str]]:
    """
    Prepare the code files to archive and return list of (source_path, arcname) tuples.

    Args:
        benchmark: The benchmark execution object

    Returns:
        List of tuples (source_path, arcname) where:
            - source_path: absolute path to the file on disk
            - arcname: relative path to use inside the archive
    """
    from omnibenchmark.git.clone import clone_module
    from omnibenchmark.model.repo import get_repo_hash

    nodes = benchmark.get_nodes()
    repositories = set()
    for node in nodes:
        repo = node.get_repository()
        repositories.add((repo.url, repo.commit))

    files = []

    while repositories:
        repo_url, repo_commit = repositories.pop()
        try:
            repo_path = clone_module(repo_url, repo_commit)
            if repo_path.exists():
                # Create a meaningful name for the repo in the archive
                repo_name = get_repo_hash(repo_url, repo_commit)
                for f in repo_path.rglob("*"):
                    if f.is_file():
                        # Create arcname as modules/{repo_name}/{relative_path}
                        rel_path = f.relative_to(repo_path)
                        arcname = f"modules/{repo_name}/{rel_path}"
                        files.append((f, arcname))
        except Exception:
            # Skip repositories that can't be cloned
            continue

    return files


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

    KNOWN ISSUE: This function currently has a path mismatch problem with Snakemake remote storage.

    The problem:
    - Locally, files are stored in: out/data/D1/.../file.json
    - In S3 (via Snakemake --default-storage-prefix), files are: data/D1/.../file.json
    - The versioning system expects files in tracked_directories=['out', 'versions', 'config', 'software']
    - But S3 has 'data' which isn't in tracked_directories
    - So when 'ob remote version create' runs, it only tags 'config' and 'versions', not 'data'
    - And when archiving, we can't find the data files in the version manifest

    Root cause: Snakemake strips the working directory prefix when uploading to S3.
    If the Snakefile has outputs like "out/data/...", Snakemake uploads them as "data/..."

    Proposed fix:
    1. Update StorageOptions.tracked_directories to include stage names ('data', 'methods', etc.)
       in addition to 'out' when using remote storage
    2. OR: Change --default-storage-prefix to s3://bucket/out to preserve the prefix
    3. OR: Update get_expected_benchmark_output_files to return paths relative to out_dir
       when remote_storage=True

    For now, this function includes a workaround that strips the out/ prefix for matching.

    Args:
        benchmark: The benchmark execution object
        results_dir: Directory containing results files

    Returns:
        List[Path]: The filenames of results downloaded from remote storage
    """
    # Import here to avoid circular imports
    from omnibenchmark.remote.storage import get_storage, remote_storage_args

    # Get expected files based on benchmark configuration
    storage_options = StorageOptions(out_dir=results_dir)
    expected_files = get_expected_benchmark_output_files(benchmark, storage_options)

    # WORKAROUND: Strip out_dir prefix from expected files for S3 comparison
    # S3 keys don't include the out/ prefix when uploaded via Snakemake remote storage
    # because Snakemake uses paths relative to the working directory
    expected_files_normalized = []
    for f in expected_files:
        # Convert to string if it's a Path object
        f_str = str(f) if isinstance(f, Path) else f
        # Strip the results_dir prefix (e.g., "out/") from the path
        if f_str.startswith(results_dir + "/"):
            expected_files_normalized.append(f_str[len(results_dir) + 1 :])
        else:
            expected_files_normalized.append(f_str)

    # Get storage connection with custom storage_options
    auth_options = remote_storage_args(benchmark.model)
    storage_api = benchmark.get_storage_api()
    bucket_name = benchmark.get_storage_bucket_name()

    if storage_api is None:
        raise ValueError("No storage API configured for benchmark")
    if bucket_name is None:
        raise ValueError("No storage bucket configured for benchmark")

    ss = get_storage(
        storage_api,
        auth_options,
        bucket_name,
        storage_options,
    )
    if ss is None:
        raise ValueError("Failed to initialize storage - check credentials")

    # Set version and get objects from the version manifest
    # NOTE: This requires that 'ob remote version create' was run after the benchmark
    # and that it successfully tagged the result files in S3
    ss.set_version(benchmark.get_benchmark_version())
    ss._get_objects()

    # Filter to only expected files that exist in remote storage
    # Use normalized paths (without out/ prefix) for matching against S3 keys
    files_in_storage = {
        k: v for k, v in ss.files.items() if k in expected_files_normalized
    }

    # Download each file to its proper location (with out/ prefix)
    downloaded_files = []
    for objectname in files_in_storage.keys():
        # Download to the local path with out_dir prefix
        filename = f"{results_dir}/{objectname}"
        ss.download_object(objectname, filename)
        downloaded_files.append(Path(filename))

    return downloaded_files
