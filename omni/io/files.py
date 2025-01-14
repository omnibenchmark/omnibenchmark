"""Functions to manage files"""

import warnings
from pathlib import Path

import tqdm
import yaml

from omni.benchmark import Benchmark
from omni.cli.utils.logging import logger
from omni.io.utils import get_storage, md5, remote_storage_args, sizeof_fmt


def list_files(
    benchmark: str,
    type: str = None,
    stage: str = None,
    module: str = None,
    file_id: str = None,
):
    """List all available files for a certain benchmark, version and stage"""

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

    all_files = benchmark.get_output_paths()
    expected_files = []
    for file in all_files:
        filter_stage = False
        if stage is not None:
            if file.split("/")[-4] != stage:
                filter_stage = True

        filter_module = False
        if module is not None:
            if file.split("/")[-3] != module:
                filter_module = True

        if not filter_stage and not filter_module:
            expected_files.append(file.replace("{dataset}", file.split("/")[2]))

    files = {k: v for k, v in ss.files.items() if k in expected_files}

    # get urls
    objectnames = list(files.keys())
    etags = [files[objectname]["etag"] for objectname in objectnames]

    return objectnames, etags


def download_files(
    benchmark: str,
    type: str = None,
    stage: str = None,
    module: str = None,
    file_id: str = None,
    verbose: bool = False,
):
    """Download all available files for a certain benchmark, version and stage"""

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

    if verbose:
        size = sum(
            [int(ss.files[objectname]["size"]) for objectname in ss.files.keys()]
        )
        logger.debug(
            f"Downloading {len(ss.files)} files with a total size of {sizeof_fmt(size)} ... ",
        )
    for objectname, filename in tqdm.tqdm(
        zip(objectnames, filenames), delay=5, disable=not verbose
    ):
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
    # benchmark = "tests/data/Benchmark_003.yaml"

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
