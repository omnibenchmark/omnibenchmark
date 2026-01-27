"""
Main archive implementation for creating benchmark archives.

This module provides the core archiving functionality that can work with or
without remote storage, separated from the remote storage module to avoid
circular dependencies and architectural confusion.
"""

import tarfile
import zipfile
from pathlib import Path
from typing import List, Optional, Union

from omnibenchmark.benchmark import BenchmarkExecution
from omnibenchmark.constants import COMPRESSION_GZIP, DEFAULT_COMPRESSION

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
    compression: Union[int, str] = DEFAULT_COMPRESSION,
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
    # Collect all files to save as (source_path, arcname) tuples
    # Some components return plain paths, others return tuples
    files_to_archive: List[tuple[Path, str]] = []

    # Config (benchmark.yaml) - always include if requested
    if config:
        for f in prepare_archive_config(benchmark):
            files_to_archive.append((f, str(f)))

    # Code (repository files) - returns tuples with proper arcnames
    if code:
        files_to_archive += prepare_archive_code(benchmark)

    # Software (environment files)
    if software:
        for f in prepare_archive_software(benchmark):
            files_to_archive.append((f, str(f)))

    # Results (output files)
    if results:
        for f in prepare_archive_results(benchmark, results_dir, remote_storage):
            files_to_archive.append((f, str(f)))

    if dry_run:
        # Return just the arcnames for display
        return [Path(arcname) for _, arcname in files_to_archive]
    else:
        # Determine file extension based on compression
        file_extension = get_archive_file_extension(compression)

        # Create archive filename
        benchmark_name = benchmark.get_benchmark_name()
        benchmark_version = benchmark.get_benchmark_version()
        # Sanitize: replace whitespace with underscores
        benchmark_name = "_".join(benchmark_name.split())
        benchmark_version = "_".join(benchmark_version.split())
        archive_base_name = f"{benchmark_name}_{benchmark_version}"

        if custom_filename:
            archive_path = outdir / custom_filename
        else:
            outfile = f"{archive_base_name}{file_extension}"
            archive_path = outdir / outfile

        # Save all files to archive
        if compression == COMPRESSION_GZIP:
            # Use tarfile for tar.gz archives
            # Wrap contents in a top-level directory to avoid tar bomb
            with tarfile.open(
                archive_path, "w:gz", compresslevel=compresslevel or 9
            ) as archive:
                for source_path, arcname in files_to_archive:
                    if Path(source_path).is_file():
                        full_arcname = f"{archive_base_name}/{arcname}"
                        archive.add(source_path, arcname=full_arcname)
        else:
            # Use zipfile for ZIP-based archives
            # Wrap contents in a top-level directory to avoid zip bomb
            # At this point compression is guaranteed to be an int (zipfile constant)
            assert isinstance(compression, int)
            with zipfile.ZipFile(
                archive_path, "w", compression=compression, compresslevel=compresslevel
            ) as archive:
                for source_path, arcname in files_to_archive:
                    if Path(source_path).is_file():
                        full_arcname = f"{archive_base_name}/{arcname}"
                        archive.write(source_path, full_arcname)

        return [archive_path]


def get_archive_file_extension(compression: Union[int, str]) -> str:
    """Get appropriate file extension for compression method."""
    if compression == COMPRESSION_GZIP:
        return ".tar.gz"
    match compression:
        case zipfile.ZIP_BZIP2:
            return ".bz2"
        case zipfile.ZIP_LZMA:
            return ".xz"
        case _:
            return ".zip"
