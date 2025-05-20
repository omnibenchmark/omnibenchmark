import logging

from typing import Optional

from omnibenchmark.benchmark import Benchmark

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


# XXX revisit this, conceptually. Here we're mixing the storage API with the concrete
# MinIO implementation. We should use a factory pattern to create the appropriate storage object instead,
# assuming we support multiple storage types.
def get_storage(
    storage_type: str, auth_options: dict, benchmark: str
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
        return MinIOStorage(auth_options, benchmark)


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


def remote_storage_args(benchmark: Benchmark) -> dict:
    if (
        str(benchmark.converter.model.storage_api).upper() == "MINIO"
        or str(benchmark.converter.model.storage_api).upper() == "S3"
    ):
        auth_options = S3_access_config_from_env()
        if benchmark.converter.model.storage is None:
            return {}
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
        return {}


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
        return {}
