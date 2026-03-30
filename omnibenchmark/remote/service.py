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
        service.load_version()          # set current version
        files = service.storage.files   # lazily fetches object list on first access
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

        # Derive out_dir from the benchmark's execution context when available,
        # falling back to "out" for lightweight model objects that have no context.
        if storage_options is None:
            from pathlib import Path

            out_dir = getattr(
                getattr(benchmark_model, "context", None), "out_dir", Path("out")
            )
            storage_options = StorageOptions(out_dir=str(out_dir))

        auth_options = remote_storage_args(
            benchmark_model, required=require_credentials
        )
        try:
            storage = StorageFactory.create(api, auth_options, bucket, storage_options)
        except Exception as e:
            # Re-raise S3/storage-level exceptions as ValueError so callers only
            # need to handle one exception type at the service boundary.
            raise ValueError(str(e)) from e

        if storage is None:
            raise ValueError(
                "S3 storage dependencies not installed. Run: pip install omnibenchmark[s3]"
            )

        self.storage = storage
        self._benchmark_model = benchmark_model

    def load_version(self, version: Optional[str] = None) -> None:
        """Set a version (defaults to the benchmark's configured version).

        File metadata is fetched lazily on first access of ``storage.files``.
        """
        resolved = (
            version
            if version is not None
            else self._benchmark_model.get_benchmark_version()
        )
        self.storage.set_version(resolved)

    @classmethod
    def from_path(
        cls,
        benchmark_path: str,
        *,
        full_model: bool = False,
        require_credentials: bool = True,
        storage_options: Optional[StorageOptions] = None,
    ) -> "StorageService":
        """Create a StorageService directly from a benchmark YAML path.

        Args:
            benchmark_path: Path to the benchmark YAML file.
            full_model: When True load a ``BenchmarkExecution`` (full DAG);
                when False load a lightweight ``Benchmark``.
            require_credentials: Whether write credentials are required.
            storage_options: Optional storage configuration override.

        Returns:
            A configured ``StorageService``.  Call ``load_version()`` on the
            result before accessing ``storage.files``.
        """
        from pathlib import Path

        if full_model:
            from omnibenchmark.benchmark import BenchmarkExecution

            benchmark_model = BenchmarkExecution(Path(benchmark_path))
        else:
            from omnibenchmark.benchmark import Benchmark

            benchmark_model = Benchmark(Path(benchmark_path))

        return cls(
            benchmark_model,
            require_credentials=require_credentials,
            storage_options=storage_options,
        )
