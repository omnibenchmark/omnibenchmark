"""Functions to manage files"""

import asyncio
import aiohttp
from urllib.parse import urlparse
import re
from pathlib import Path
from omni.io.utils import get_storage, md5
from omni.sync import get_bench_definition
from typing import Union, List
import tqdm
import warnings


def list_files(
    benchmark: str,
    type: str = None,
    stage: str = None,
    module: str = None,
    file_id: str = None,
    version: str = None,
    verbose: bool = False,
):
    """List all available files for a certain benchmark, version and stage"""
    # bench_yaml = get_bench_definition(bench_name, version, stage)

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

    ss = get_storage(bench_yaml["storage_type"], bench_yaml["auth_options"], benchmark)
    if version is None:
        ss.set_current_version()
    else:
        if not version in ss.versions:
            raise ValueError(f"Version {version} not found in {ss.benchmark}")

        # set version
        ss.set_current_version(
            major_version=int(version.split(".")[0]),
            minor_version=int(version.split(".")[1]),
        )
        if not f"{ss.major_version}.{ss.minor_version}" == version:
            raise ValueError("Version specification failed!")

    # list objects of version
    ss._get_objects()

    # get urls
    names = list(ss.files.keys())

    # TODO: filter by arguments
    # names = ...

    # create urls
    urls = {}
    for name in names:
        url = urlparse(
            f"{ss.auth_options['endpoint']}/{ss.benchmark}.{ss.major_version}.{ss.minor_version}/{name}"
        )
        if "secure" in ss.auth_options.keys() and not ss.auth_options["secure"]:
            url = url._replace(scheme="http")
        else:
            url = url._replace(scheme="https")
        urls[name] = {
            "url": url.geturl(),
            "size": int(ss.files[name]["size"]),
            "md5": ss.files[name]["hash"],
        }
    return urls


# adapted from https://realpython.com/python-download-file-from-url/#performing-parallel-file-downloads
async def retrieve_file(url: str):
    async with aiohttp.ClientSession() as session:
        async with session.get(url) as response:
            # to remove schema
            urlp = urlparse(url)
            # to remove benchmark name
            filename = re.sub("^/[a-zA-Z0-9._-]*/", "", urlp.path)
            # create missing directories
            Path(filename).parent.mkdir(parents=True, exist_ok=True)
            # download file
            with open(filename, mode="wb") as file:
                while True:
                    chunk = await response.content.read()
                    if not chunk:
                        break
                    file.write(chunk)
            return response.headers


async def retrieve_files(urls: List[str], verbose: bool = False):
    if verbose:
        print("Downloading files...")
    tasks = [retrieve_file(url) for url in urls]
    from tqdm.asyncio import tqdm_asyncio

    headers = await tqdm_asyncio.gather(*tasks, delay=5, disable=not verbose)
    return headers


def get_file(url: Union[str, list[str]], verbose: bool = False):
    """Download specific file(s) based on its url"""
    if type(url) is str:
        headers = asyncio.run(retrieve_file(url))
        filenames = re.sub("^/[a-zA-Z0-9._-]*/", "", urlparse(url).path)
    elif type(url) is list:
        headers = asyncio.run(retrieve_files(url))
        filenames = [re.sub("^/[a-zA-Z0-9._-]*/", "", urlparse(u).path) for u in url]
    else:
        raise ValueError("url must be a string or a list of strings")

    # md5 checksum
    if verbose:
        print("Checking MD5 checksums...")
    for i in tqdm.tqdm(len(filenames), delay=5, disable=not verbose):
        md5sum_val = md5(filenames[i])
        if not headers[i]["ETag"].replace('"', "") == md5sum_val:
            warnings.warn(f"MD5 checksum failed for {filenames[i]}", Warning)

    return filenames


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
    urls = list_files(benchmark, type, stage, module, file_id, version, verbose=verbose)
    filenames = get_file([urls[i]["url"] for i in urls.keys()], verbose=verbose)
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

    urls = list_files(benchmark, type, stage, module, file_id, version)

    for i in tqdm.tqdm(urls.keys(), delay=5, disable=not verbose):
        # to remove schema
        urlp = urlparse(urls[i]["url"])
        # to remove benchmark name
        filename = re.sub("^/[a-zA-Z0-9._-]*/", "", urlp.path)
        if Path(filename).exists():
            md5_local = md5(filename)
        else:
            md5_local = None
        urls[filename]["md5_local"] = md5_local

    failed_checksums = []
    for i in urls.keys():
        if not urls[i]["md5"] == urls[i]["md5_local"]:
            failed_checksums.append(i)
            if verbose:
                print(f"MD5 checksum failed for {i}")
    return failed_checksums
