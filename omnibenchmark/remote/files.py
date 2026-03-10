"""Functions to manage files"""

import warnings
from pathlib import Path
from typing import Optional

import tqdm

from omnibenchmark.benchmark import Benchmark
from omnibenchmark.cli.utils.logging import logger
from omnibenchmark.remote.RemoteStorage import StorageOptions
from omnibenchmark.remote.versioning import get_expected_benchmark_output_files


def _setup_remote_storage(
    benchmark_path: str,
    storage_options: Optional[StorageOptions] = None,
):
    """Set up remote storage for a benchmark.

    Returns:
        (storage, benchmark, expected_files) tuple.

    Raises:
        ValueError if storage cannot be configured.
    """
    from .service import StorageService

    service = StorageService.from_path(benchmark_path, storage_options=storage_options)
    service.load_version()
    expected_files = get_expected_benchmark_output_files(
        service._benchmark_model, service.storage.storage_options
    )

    return service.storage, service._benchmark_model, expected_files


def list_files(
    benchmark_path: str,
    type: str,
    stage: str,
    module: str,
    file_id: str,
    remote_storage: bool = True,
    storage_options: Optional[StorageOptions] = None,
):
    """List all available files for a certain benchmark, version and stage"""
    resolved_options = storage_options or StorageOptions(out_dir="out")
    if remote_storage:
        ss, _, expected_files = _setup_remote_storage(benchmark_path, resolved_options)
        files = {k: v for k, v in ss.files.items() if k in expected_files}
        objectnames = list(files.keys())
        etags = [files[objectname]["etag"] for objectname in objectnames]
    else:
        benchmark = Benchmark(Path(benchmark_path))
        expected_files = get_expected_benchmark_output_files(
            benchmark, resolved_options
        )
        file_local_is_local = [Path(f).is_file() for f in expected_files]
        if not all(file_local_is_local):
            logger.warning("Not all expected files are available locally")
            logger.info("Missing files:")
            for i, f in enumerate(expected_files):
                if not file_local_is_local[i]:
                    logger.info(f)
        # make unique
        objectnames = list(set([str(Path(f)) for f in expected_files]))
        etags = []

    return objectnames, etags


def download_files(
    benchmark_path: str,
    type: str,
    stage: str,
    module: str,
    file_id: str,
    verbose: bool = False,
    overwrite: bool = False,
):
    """Download all available files for a certain benchmark, version and stage"""
    from .hash import checksum
    from .sizeof import sizeof_fmt

    ss, _, expected_files = _setup_remote_storage(benchmark_path)

    files = {k: v for k, v in ss.files.items() if k in expected_files}
    objectnames = list(files.keys())
    etags = [files[objectname]["etag"] for objectname in objectnames]
    filenames = objectnames

    logger.debug("Checking if files are already downloaded... ")
    do_download_file = []
    for filename, etag in tqdm.tqdm(
        zip(filenames, etags),
        total=len(filenames),
        desc="Checking local files",
        unit="file",
        delay=1,
        disable=not verbose,
    ):
        if Path(filename).is_file():
            if etag.replace('"', "") == checksum(filename):
                do_download_file.append(False)
            else:
                if overwrite:
                    do_download_file.append(True)
                else:
                    do_download_file.append(False)
        else:
            do_download_file.append(True)

    if verbose:
        size = sum(
            [
                int(files[objectname]["size"])
                for objectname, do_download in zip(objectnames, do_download_file)
                if do_download
            ]
        )
        logger.debug(
            f"Downloading {sum(do_download_file)} files with a total size of {sizeof_fmt(size)} ... ",
        )

    pairs_to_download = [
        (obj, fn) for obj, fn, do in zip(objectnames, filenames, do_download_file) if do
    ]
    for objectname, filename in tqdm.tqdm(
        pairs_to_download,
        total=len(pairs_to_download),
        desc="Downloading",
        unit="file",
        delay=1,
        disable=not verbose,
    ):
        ss.download_object(objectname, filename)
    if verbose:
        logger.debug("Done")

    if verbose:
        logger.debug("Checking MD5 checksums... ")
    for etag, filename in tqdm.tqdm(
        zip(etags, filenames),
        total=len(filenames),
        desc="Verifying checksums",
        unit="file",
        delay=1,
        disable=not verbose,
    ):
        md5sum_val = checksum(filename)
        if not etag.replace('"', "") == md5sum_val:
            warnings.warn(f"MD5 checksum failed for {filename}", Warning)
    if verbose:
        logger.debug("Done")

    return filenames


def checksum_files(
    benchmark: str,
    type: str,
    stage: str,
    module: str,
    file_id: str,
    verbose: bool = False,
):
    """Compare md5 checksums of available files for a certain benchmark, version and stage with local versions"""
    from .hash import checksum

    objectnames, etags = list_files(benchmark, type, stage, module, file_id)

    filenames = objectnames

    failed_checksums = []
    for etag, filename in tqdm.tqdm(
        zip(etags, filenames),
        total=len(filenames),
        desc="Verifying checksums",
        unit="file",
        delay=1,
        disable=not verbose,
    ):
        md5sum_val = checksum(filename)
        if not etag.replace('"', "") == md5sum_val:
            failed_checksums.append(filename)
    return failed_checksums
