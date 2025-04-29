"""Functions to manage files"""

import warnings
from pathlib import Path

import tqdm
import yaml

from omnibenchmark.benchmark import Benchmark
from omnibenchmark.cli.utils.logging import logger
from omnibenchmark.io.versioning import get_expected_benchmark_output_files


def list_files(
    benchmark: str,
    type: str = None,
    stage: str = None,
    module: str = None,
    file_id: str = None,
    local: bool = False,
):
    """List all available files for a certain benchmark, version and stage"""
    from omnibenchmark.io.utils import get_storage, remote_storage_args

    with open(benchmark, "r") as fh:
        yaml.safe_load(fh)
        benchmark = Benchmark(Path(benchmark))

    expected_files = get_expected_benchmark_output_files(benchmark)

    if not local:
        auth_options = remote_storage_args(benchmark)

        ss = get_storage(
            str(benchmark.converter.model.storage_api),
            auth_options,
            str(benchmark.converter.model.storage_bucket_name),
        )
        ss.set_version(benchmark.get_benchmark_version())
        ss._get_objects()
        files = {k: v for k, v in ss.files.items() if k in expected_files}

        # get urls
        objectnames = list(files.keys())
        etags = [files[objectname]["etag"] for objectname in objectnames]

    else:
        file_local_is_local = [Path(f).is_file() for f in expected_files]
        if not all(file_local_is_local):
            logger.warning("Not all expected files are available locally")
            logger.info("Missing files:")
            for i, f in enumerate(expected_files):
                if not file_local_is_local[i]:
                    logger.info(f)
        # make unique
        objectnames = list(set([Path(f) for f in expected_files]))
        etags = []

    return objectnames, etags


def download_files(
    benchmark: str,
    type: str = None,
    stage: str = None,
    module: str = None,
    file_id: str = None,
    verbose: bool = False,
    overwrite: bool = False,
):
    """Download all available files for a certain benchmark, version and stage"""
    from omnibenchmark.io.utils import get_storage, remote_storage_args
    from .utils import md5, sizeof_fmt

    objectnames, etags = list_files(benchmark, type, stage, module, file_id)

    # storage path locally, TODO: maybe add as argument
    filenames = objectnames

    with open(benchmark, "r") as fh:
        yaml.safe_load(fh)
        benchmark = Benchmark(Path(benchmark))

    auth_options = remote_storage_args(benchmark)

    ss = get_storage(
        str(benchmark.converter.model.storage_api),
        auth_options,
        str(benchmark.converter.model.storage_bucket_name),
    )
    ss.set_version(benchmark.get_benchmark_version())
    ss._get_objects()

    logger.debug("Checking if files are already downloaded... ")
    do_download_file = []
    for filename, etag in tqdm.tqdm(
        zip(filenames, etags), delay=5, disable=not verbose
    ):
        if Path(filename).is_file():
            if etag.replace('"', "") == md5(filename):
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
                int(ss.files[objectname]["size"])
                for objectname, do_download in zip(ss.files.keys(), do_download_file)
                if do_download
            ]
        )
        logger.debug(
            f"Downloading {sum(do_download_file)} files with a total size of {sizeof_fmt(size)} ... ",
        )

    for objectname, filename, do_download in tqdm.tqdm(
        zip(objectnames, filenames, do_download_file), delay=5, disable=not verbose
    ):
        if do_download:
            ss.download_object(objectname, filename)
    if verbose:
        logger.debug("Done")

    if verbose:
        logger.debug("Checking MD5 checksums... ")
    for etag, filename in tqdm.tqdm(
        zip(etags, filenames), delay=5, disable=not verbose
    ):
        md5sum_val = md5(filename)
        if not etag.replace('"', "") == md5sum_val:
            warnings.warn(f"MD5 checksum failed for {filename}", Warning)
    if verbose:
        logger.debug("Done")

    return filenames


def checksum_files(
    benchmark: str,
    type: str = None,
    stage: str = None,
    module: str = None,
    file_id: str = None,
    verbose: bool = False,
):
    """Compare md5 checksums of available files for a certain benchmark, version and stage with local versions"""
    from .utils import md5

    objectnames, etags = list_files(benchmark, type, stage, module, file_id)

    filenames = objectnames

    failed_checksums = []
    for etag, filename in tqdm.tqdm(
        zip(etags, filenames), delay=5, disable=not verbose
    ):
        md5sum_val = md5(filename)
        if not etag.replace('"', "") == md5sum_val:
            failed_checksums.append(filename)
    return failed_checksums
