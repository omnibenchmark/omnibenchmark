"""Utility functions to manage dataset handling"""

import hashlib

from omni.benchmark import Benchmark
from omni.io.MinIOStorage import MinIOStorage
from omni.io.S3config import S3_access_config_from_env

# from omni_schema.datamodel.omni_schema import StorageAPIEnum


def get_storage(storage_type: str, auth_options: dict, benchmark: str):
    """
    Selects a remote storage type.

    Args:
    - storage_type (str): The type of the remote storage.
    - auth_options (dict): The authentication options.
    - benchmark (str): The benchmark name.

    Returns:
    - RemoteStorage: The remote storage object.
    """
    if storage_type.upper() == "MINIO" or storage_type.upper() == "S3":
        return MinIOStorage(auth_options, benchmark)
    else:
        raise ValueError("Invalid storage type")


# https://stackoverflow.com/a/3431838
def md5(fname: str):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


# from: https://stackoverflow.com/a/1094933
def sizeof_fmt(num: int, suffix: str = "B"):
    if abs(num) < 1024.0:
        return f"{num: 5d}{suffix}"
    for unit in ("", "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi"):
        if abs(num) < 1024.0:
            return f"{num:3.1f}{unit}{suffix}"
        num /= 1024.0
    return f"{num:.1f}Yi{suffix}"


def remote_storage_args(benchmark: Benchmark) -> dict:
    if (
        str(benchmark.converter.model.storage_api).upper() == "MINIO"
        or str(benchmark.converter.model.storage_api).upper() == "S3"
    ):
        auth_options = S3_access_config_from_env()
        return {
            "endpoint": benchmark.converter.model.storage,
            "secure": False,
            "access_key": auth_options["access_key"],
            "secret_key": auth_options["secret_key"],
        }
    else:
        raise ValueError("Invalid storage api")


def remote_storage_snakemake_args(benchmark: Benchmark) -> dict:
    if (
        str(benchmark.converter.model.storage_api).upper() == "MINIO"
        or str(benchmark.converter.model.storage_api).upper() == "S3"
    ):
        auth_options = S3_access_config_from_env()
        return {
            "default-storage-provider": "s3",
            "default-storage-prefix": f"s3://{benchmark.converter.model.id}",
            "storage-s3-endpoint-url": benchmark.converter.model.storage,
            "storage-s3-access-key": auth_options["access_key"],
            "storage-s3-secret-key": auth_options["secret_key"],
        }
    else:
        raise ValueError("Invalid storage api")
