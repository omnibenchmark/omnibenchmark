import logging
import warnings

from typing import Optional, TYPE_CHECKING, Any

from omnibenchmark.benchmark import BenchmarkExecution
from .RemoteStorage import StorageOptions

if TYPE_CHECKING:
    from .MinIOStorage import MinIOStorage as MinIOStorageType
else:
    MinIOStorageType = Any

try:
    from .MinIOStorage import MinIOStorage
    from .S3config import S3_access_config_from_env

    S3_AVAILABLE = True
except ImportError:
    S3_AVAILABLE = False
    MinIOStorage = None  # type: ignore
    S3_access_config_from_env = None


def _check_s3_available():
    """Check if S3 storage is available and warn if not."""
    if not S3_AVAILABLE:
        logger = logging.getLogger(__name__)
        logger.warning(
            "S3 storage not available. Install with: pip install omnibenchmark[s3]"
        )
    return S3_AVAILABLE


def _warn_minio_deprecated(storage_type: str) -> None:
    """Emit a deprecation warning if storage_type is MinIO."""
    if storage_type.upper() == "MINIO":
        warnings.warn(
            "Storage type 'MinIO' is deprecated. Use 'S3' with a custom endpoint instead.",
            DeprecationWarning,
            stacklevel=3,
        )


class StorageFactory:
    """Factory for creating RemoteStorage instances.

    Centralises backend selection so that adding a new storage backend
    (e.g. GCS, Azure Blob) only requires a change here, not across all
    call sites.

    Example::

        ss = StorageFactory.create("S3", auth_options, bucket_name)
    """

    @staticmethod
    def create(
        storage_type: str,
        auth_options: dict,
        bucket: str,
        storage_options: StorageOptions = StorageOptions(out_dir="out"),
    ) -> Optional[MinIOStorageType]:
        """Return a RemoteStorage instance for *storage_type*, or None.

        Args:
            storage_type: Storage API identifier ("S3" or the deprecated "MinIO").
            auth_options: Authentication dict (endpoint, access_key, secret_key).
            bucket: S3 bucket name used as the benchmark identifier.
            storage_options: Tracked-directory configuration.

        Returns:
            A RemoteStorage instance, or None when the required extras are not
            installed.
        """
        _warn_minio_deprecated(storage_type)

        if storage_type.upper() in ("MINIO", "S3"):
            if not _check_s3_available() or MinIOStorage is None:
                return None
            return MinIOStorage(auth_options, bucket, storage_options)

        return None


# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------


def get_storage(
    storage_type: str,
    auth_options: dict,
    benchmark: str,
    storage_options: StorageOptions = StorageOptions(out_dir="out"),
) -> Optional[MinIOStorageType]:
    """Return a RemoteStorage instance for *storage_type*.

    Thin wrapper around :class:`StorageFactory` kept for backward
    compatibility. Prefer calling ``StorageFactory.create()`` directly in
    new code.

    Args:
        storage_type: Storage API identifier ("S3" or the deprecated "MinIO").
        auth_options: Authentication dict (endpoint, access_key, secret_key).
        benchmark: S3 bucket / benchmark name.
        storage_options: Tracked-directory configuration.
    """
    return StorageFactory.create(storage_type, auth_options, benchmark, storage_options)


def get_storage_from_benchmark(
    benchmark: BenchmarkExecution,
) -> Optional[MinIOStorageType]:
    """Return a RemoteStorage instance configured from *benchmark*.

    Args:
        benchmark: The benchmark execution object.

    Returns:
        A RemoteStorage instance, or None when storage is not configured.
    """
    auth_options = remote_storage_args(benchmark)
    storage_api = benchmark.get_storage_api()
    bucket_name = benchmark.get_storage_bucket_name()

    if storage_api is None or bucket_name is None:
        return None

    return StorageFactory.create(
        storage_api, auth_options, bucket_name, StorageOptions(out_dir="out")
    )


def remote_storage_args(benchmark, required: bool = False) -> dict:
    """Return authentication args for the benchmark's remote storage.

    Args:
        benchmark: The benchmark object.
        required: If True, require credentials to be present (write operations).

    Returns:
        dict with endpoint, secure, and optionally access_key / secret_key.
    """
    storage_api = benchmark.get_storage_api()
    if storage_api:
        _warn_minio_deprecated(storage_api)

    if storage_api and (storage_api.upper() == "MINIO" or storage_api.upper() == "S3"):
        if not _check_s3_available() or S3_access_config_from_env is None:
            return {}
        auth_options = S3_access_config_from_env(required=required)
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
    """Return Snakemake storage plugin args for *benchmark*'s remote storage."""
    storage_api = benchmark.get_storage_api()
    if storage_api:
        _warn_minio_deprecated(storage_api)

    if storage_api and (storage_api.upper() == "MINIO" or storage_api.upper() == "S3"):
        if not _check_s3_available() or S3_access_config_from_env is None:
            return {}
        auth_options = S3_access_config_from_env()
        bucket_name = benchmark.get_storage_bucket_name()
        endpoint = benchmark.get_storage_endpoint()
        return {
            "default-storage-provider": "s3",
            "default-storage-prefix": f"s3://{bucket_name}",
            "storage-s3-endpoint-url": endpoint,
            "storage-s3-access-key": auth_options["access_key"],
            "storage-s3-secret-key": auth_options["secret_key"],
            "shared-fs-usage": "persistence",  # Keep local copies and use them when available
        }
    else:
        return {}
