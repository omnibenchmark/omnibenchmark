"""High-level service for remote storage operations."""

from typing import Optional

from omnibenchmark.remote.RemoteStorage import StorageOptions


class StorageService:
    """Encapsulates storage setup: configuration validation, authentication, and object loading.

    Accepts any object that implements the benchmark storage interface (i.e. both
    ``BenchmarkExecution`` and the lighter-weight ``Benchmark`` model), so callers
    can avoid building a full DAG when only the storage configuration is needed.

    Usage::

        service = StorageService(benchmark_model)
        service.load_version()          # set current version and fetch object list
        files = service.storage.files   # access the loaded object dict
    """

    def __init__(
        self,
        benchmark_model,
        require_credentials: bool = True,
        storage_options: Optional[StorageOptions] = None,
    ):
        from omnibenchmark.remote.storage import remote_storage_args, StorageFactory

        api = benchmark_model.get_storage_api()
        bucket = benchmark_model.get_storage_bucket_name()

        if api is None:
            raise ValueError(
                "No storage API configured. Set 'storage.api' in your benchmark YAML."
            )
        if bucket is None:
            raise ValueError(
                "No storage bucket configured."
                " Set 'storage.bucket_name' in your benchmark YAML."
            )

        auth_options = remote_storage_args(
            benchmark_model, required=require_credentials
        )
        storage = StorageFactory.create(api, auth_options, bucket, storage_options)

        if storage is None:
            raise ValueError(
                "S3 storage dependencies not installed. Run: pip install omnibenchmark[s3]"
            )

        self.storage = storage
        self._benchmark_model = benchmark_model

    def load_version(self, version: Optional[str] = None) -> None:
        """Set and load a version (defaults to the benchmark's configured version)."""
        resolved = (
            version
            if version is not None
            else self._benchmark_model.get_benchmark_version()
        )
        self.storage.set_version(resolved)
        self.storage.load_objects()
