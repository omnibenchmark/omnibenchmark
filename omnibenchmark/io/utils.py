"""Utility functions to manage dataset handling"""

# TODO: rename module to something more descriptive. E.g., storage

import hashlib

from collections import defaultdict

from omnibenchmark.benchmark import Benchmark
from omnibenchmark.io.MinIOStorage import MinIOStorage
from omnibenchmark.io.S3config import S3_access_config_from_env


# XXX revisit this, conceptually. Here we're mixing the storage API with the concrete
# MinIO implementation. We should use a factory pattern to create the appropriate storage object instead,
# assuming we support multiple storage types.
def get_storage(storage_type: str, auth_options: dict, benchmark: str) -> MinIOStorage:
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


def get_storage_from_benchmark(benchmark: Benchmark):
    """
    Selects a remote storage type from a benchmark object.

    Args:
    - benchmark (Benchmark): The benchmark object.

    Returns:
    - RemoteStorage: The remote storage object.
    """
    auth_options = remote_storage_args(benchmark)
    # setup storage
    return get_storage(
        str(benchmark.converter.model.storage_api),
        auth_options,
        str(benchmark.converter.model.storage_bucket_name),
    )


# https://stackoverflow.com/a/3431838
def md5(fname: str):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


# from: https://stackoverflow.com/a/1094933
# TODO(ben): use humanize
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
        base_dic = {
            "endpoint": benchmark.converter.model.storage,
            "secure": (
                True if benchmark.converter.model.storage.startswith("https") else False
            ),
        }
        if "access_key" in auth_options and "secret_key" in auth_options:
            auth_dic = {
                "access_key": auth_options["access_key"],
                "secret_key": auth_options["secret_key"],
            }
        else:
            auth_dic = {}

        return {**base_dic, **auth_dic}
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
            "default-storage-prefix": f"s3://{benchmark.converter.model.storage_bucket_name}",
            "storage-s3-endpoint-url": benchmark.converter.model.storage,
            "storage-s3-access-key": auth_options["access_key"],
            "storage-s3-secret-key": auth_options["secret_key"],
        }
    else:
        raise ValueError("Invalid storage api")


def tree():
    """
    A recursive function that creates a tree structure from a list of file paths.
    """
    return defaultdict(tree)


def tree_string_from_list(file_list):
    file_tree = tree()

    for file_path in file_list:
        parts = file_path.parts
        current_level = file_tree
        for part in parts:
            current_level = current_level[part]

    def print_tree(current_level, indent=""):
        out = ""
        for key, subtree in current_level.items():
            out += f"{indent}{key}/\n"
            out += print_tree(subtree, indent + "    ")
        return out

    return print_tree(file_tree)
