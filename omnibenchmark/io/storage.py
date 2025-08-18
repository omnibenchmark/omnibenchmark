import logging

from typing import Optional

from omnibenchmark.benchmark import Benchmark
from .RemoteStorage import StorageOptions

try:
    from .MinIOStorage import MinIOStorage
    from .S3config import S3_access_config_from_env

    S3_AVAILABLE = True
except ImportError:
    logger = logging.getLogger(__name__)
    logger.warning(
        "S3 storage not available. You might want to install the 'minio' and 'boto3' packages."
    )
    S3_AVAILABLE = False


# XXX ARCHITECTURAL DEBT: Storage Factory Pattern Missing
#
# Current Issue: We're mixing the storage API abstraction with concrete MinIO implementation.
# This creates tight coupling and makes it difficult to add new storage backends.
#
# Recommended Pattern:
# class StorageFactory:
#     @staticmethod
#     def create_storage(storage_config: Storage) -> RemoteStorage:
#         match storage_config.api:
#             case StorageAPIEnum.s3: return MinIOStorage(...)
#             case StorageAPIEnum.gcs: return GCSStorage(...)  # future
#
# This would allow proper dependency injection and easier testing with mock storage.
def get_storage(
    storage_type: str,
    auth_options: dict,
    benchmark: str,
    storage_options: StorageOptions = StorageOptions(out_dir="out"),
) -> Optional["MinIOStorage"]:
    """
    Selects a remote storage type.

    Args:
    - storage_type (str): The type of the remote storage.
    - auth_options (dict): The authentication options.
    - benchmark (str): The benchmark name.

    Returns:
    - Optional[MinIOStorage]: The remote storage object, or None if unavailable.
    """
    if storage_type.upper() == "MINIO" or storage_type.upper() == "S3":
        return MinIOStorage(auth_options, benchmark, storage_options)


def get_storage_from_benchmark(benchmark: Benchmark) -> Optional["MinIOStorage"]:
    """
    Selects a remote storage type from a benchmark object.

    Args:
    - benchmark (Benchmark): The benchmark object.

    Returns:
    - Optional[MinIOStorage]: The remote storage object, or None if unavailable.
    """
    auth_options = remote_storage_args(benchmark)
    # setup storage
    storage_api = benchmark.get_storage_api()
    bucket_name = benchmark.get_storage_bucket_name()

    if storage_api is None or bucket_name is None:
        return None

    return get_storage(
        storage_api, auth_options, bucket_name, StorageOptions(out_dir="out")
    )


def remote_storage_args(benchmark) -> dict:
    storage_api = benchmark.get_storage_api()
    if storage_api and (storage_api.upper() == "MINIO" or storage_api.upper() == "S3"):
        auth_options = S3_access_config_from_env()
        endpoint = benchmark.get_storage_endpoint()
        if endpoint is None:
            return {}
        base_dic = {
            "endpoint": endpoint,
            "secure": (True if endpoint.startswith("https") else False),
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
        return {}


def remote_storage_snakemake_args(benchmark) -> dict:
    storage_api = benchmark.get_storage_api()
    if storage_api and (storage_api.upper() == "MINIO" or storage_api.upper() == "S3"):
        auth_options = S3_access_config_from_env()
        bucket_name = benchmark.get_storage_bucket_name()
        endpoint = benchmark.get_storage_endpoint()
        return {
            "default-storage-provider": "s3",
            "default-storage-prefix": f"s3://{bucket_name}",
            "storage-s3-endpoint-url": endpoint,
            "storage-s3-access-key": auth_options["access_key"],
            "storage-s3-secret-key": auth_options["secret_key"],
        }
    else:
        return {}
