"""Functions to manage files"""

import warnings
from typing import Union

import tqdm

from omni.io.utils import get_storage, md5
from omni.sync import get_bench_definition


def list_files(
    benchmark: str,
    type: str = None,
    stage: str = None,
    module: str = None,
    file_id: str = None,
    version: Union[None, str] = None,
    verbose: bool = False,
):
    """List all available files for a certain benchmark, version and stage"""

    # TODO: for testing until get_bench_definition is implemented
    if __name__ == "__main__":
        minio_auth_options_public = {
            "endpoint": "https://omnibenchmark.mls.uzh.ch:9000",
            "secure": False,
        }
        bench_yaml = {
            "auth_options": minio_auth_options_public,
            "storage_type": "minio",
        }
        benchmark = "testversioning"
        version = "0.1"

    ss = get_storage(bench_yaml["storage_type"], bench_yaml["auth_options"], benchmark)
    # set version
    ss.set_version(version)
    # list objects of version
    ss._get_objects()

    # get urls
    objectnames = list(ss.files.keys())
    etags = [ss.files[objectname]["etag"] for objectname in objectnames]

    # TODO: filter by arguments
    # objectnames = ...
    # etags = ...

    return objectnames, etags


def download_files(
    benchmark: str,
    type: str = None,
    stage: str = None,
    module: str = None,
    file_id: str = None,
    version: str = None,
    verbose: bool = False,
):
    """Download all available files for a certain benchmark, version and stage"""
    objectnames, etags = list_files(
        benchmark, type, stage, module, file_id, version, verbose=verbose
    )

    # storage path locally
    filenames = objectnames

    # TODO: for testing until get_bench_definition is implemented
    if __name__ == "__main__":
        minio_auth_options_public = {
            "endpoint": "https://omnibenchmark.mls.uzh.ch:9000",
            "secure": False,
        }
        bench_yaml = {
            "auth_options": minio_auth_options_public,
            "storage_type": "minio",
        }
        benchmark = "testversioning"
        version = "0.1"
    ss = get_storage(bench_yaml["storage_type"], bench_yaml["auth_options"], benchmark)
    ss.set_version(version)

    for objectname, filename in tqdm.tqdm(
        zip(objectnames, filenames), delay=5, disable=not verbose
    ):
        ss.download_object(objectname, filename)

    if verbose:
        print("Checking MD5 checksums...")
    for etag, filename in tqdm.tqdm(
        zip(etags, filenames), delay=5, disable=not verbose
    ):
        md5sum_val = md5(filename)
        if not etag.replace('"', "") == md5sum_val:
            warnings.warn(f"MD5 checksum failed for {filename}", Warning)

    return filenames


def checksum_files(
    benchmark: str,
    type: str,
    stage: str,
    module: str,
    file_id: str,
    version: str,
    verbose: bool = False,
):
    """Compare md5 checksums of available files for a certain benchmark, version and stage with local versions"""

    objectnames, etags = list_files(benchmark, type, stage, module, file_id, version)

    filenames = objectnames

    failed_checksums = []
    for etag, filename in tqdm.tqdm(
        zip(etags, filenames), delay=5, disable=not verbose
    ):
        md5sum_val = md5(filename)
        if not etag.replace('"', "") == md5sum_val:
            failed_checksums.append(filename)
    return failed_checksums
