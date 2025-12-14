"""
Main archive implementation for creating benchmark archives.

This module provides the core archiving functionality that can work with or
without remote storage, separated from the remote storage module to avoid
circular dependencies and architectural confusion.
"""

import zipfile
from pathlib import Path
from typing import List, Optional

from omnibenchmark.benchmark import BenchmarkExecution

from .components import (
    prepare_archive_code,
    prepare_archive_software,
    prepare_archive_results,
    prepare_archive_config,
)


def archive_benchmark(
    benchmark: BenchmarkExecution,
    outdir: Path = Path(),
    config: bool = True,
    code: bool = False,
    software: bool = False,
    results: bool = False,
    results_dir: str = "out",
    compression=zipfile.ZIP_STORED,
    compresslevel: Optional[int] = None,
    dry_run: bool = False,
    remote_storage: bool = False,
    custom_filename: Optional[str] = None,
) -> List[Path]:
    """
    Create an archive of a benchmark and its components.

    Args:
        benchmark: The benchmark execution object
        outdir: Output directory for the archive
        config: Include configuration files (benchmark.yaml)
        code: Include code repositories
        software: Include software environment files
        results: Include result files
        results_dir: Directory containing results (default: "out")
        compression: Compression method for ZIP file
        compresslevel: Compression level
        dry_run: If True, return list of files without creating archive
        remote_storage: Whether to use remote storage for results
        custom_filename: Custom filename for the archive (overrides default naming)

    Returns:
        List of file paths that were/would be archived, or path to created archive
    """
    # Collect all filenames to save
    filenames = []

    # Config (benchmark.yaml) - always include if requested
    if config:
        filenames += prepare_archive_config(benchmark)

    # Code (repository files)
    if code:
        filenames += prepare_archive_code(benchmark)

    # Software (environment files)
    if software:
        filenames += prepare_archive_software(benchmark)

    # Results (output files)
    if results:
        filenames += prepare_archive_results(benchmark, results_dir, remote_storage)

    if dry_run:
        return filenames
    else:
        # Determine file extension based on compression
        match compression:
            case zipfile.ZIP_BZIP2:
                file_extension = ".bz2"
            case zipfile.ZIP_LZMA:
                file_extension = ".xz"
            case _:
                file_extension = ".zip"

        # Create archive filename
        if custom_filename:
            archive_path = outdir / custom_filename
        else:
            benchmark_name = benchmark.get_benchmark_name()
            benchmark_version = benchmark.get_benchmark_version()
            outfile = f"{benchmark_name}_{benchmark_version}{file_extension}"
            archive_path = outdir / outfile

        # Save all files to archive
        with zipfile.ZipFile(
            archive_path, "w", compression=compression, compresslevel=compresslevel
        ) as archive:
            for filename in filenames:
                if Path(filename).is_file():
                    archive.write(filename, filename)

        return [archive_path]


def get_archive_file_extension(compression) -> str:
    """Get appropriate file extension for compression method."""
    match compression:
        case zipfile.ZIP_BZIP2:
            return ".bz2"
        case zipfile.ZIP_LZMA:
            return ".xz"
        case _:
            return ".zip"
