"""Base class for remote storage."""

from abc import ABCMeta, abstractmethod
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Optional, Union

if TYPE_CHECKING:
    from omnibenchmark.benchmark import BenchmarkExecution

import packaging.version
from packaging.version import Version

from omnibenchmark.remote.exception import RemoteStorageInvalidInputException


class StorageOptions:
    def __init__(self, out_dir: str):
        self.tracked_directories = [out_dir, "versions", "config", "software"]
        self.results_directories = [out_dir]
        self.extra_files_to_version_not_in_benchmark_yaml = [
            f"{out_dir}/**/parameters.json",
            f"{out_dir}/**/parameters_dict.tsv",
        ]  # glob style

        assert all(
            [
                isinstance(i, str) and i in self.tracked_directories
                for i in self.results_directories
            ]
        )


def is_valid_version(version: str):
    try:
        packaging.version.parse(version)
        return True
    except packaging.version.InvalidVersion:
        return False
    except TypeError:
        return False


class RemoteStorage(metaclass=ABCMeta):
    """
    A class representing a remote storage with version-controlled S3 buckets.

    This class provides an abstract interface for managing benchmark results in S3-compatible
    storage (AWS S3, MinIO) with immutable versioning and retention policies for reproducibility.

    ## S3 Versioning and Object Protection

    When benchmark versions are created, S3 Object Lock (Governance Mode) is automatically
    applied to protect versioned objects from accidental deletion or modification:

    ### Version Creation Process:
    1. **Tag Current Objects**: Latest object versions are tagged with benchmark version
    2. **Apply Retention Policy**: Objects receive Governance Mode retention protection
    3. **Create Version Manifest**: A CSV manifest is stored at `versions/{version}.csv`
    4. **Enable Object Lock**: Bucket-level object locking prevents tampering

    ### WORM Protection (Write Once, Read Many):
    - Tagged objects become **immutable** and cannot be deleted without special permissions
    - Prevents data loss due to accidental deletion or modification
    - Ensures benchmark reproducibility by preserving exact file versions
    - Only users with `s3:BypassGovernanceRetention` permission can override protection

    ### Cleanup Considerations:
    Protected objects require special deletion procedures:
    ```python
    # Standard deletion (will fail for protected objects)
    s3_client.delete_object(Bucket=bucket, Key=key)

    # Required for protected objects
    s3_client.delete_object(
        Bucket=bucket,
        Key=key,
        VersionId=version_id,
        BypassGovernanceRetention=True  # Requires special permission
    )
    ```

    ### Test Environment Cleanup:
    For test environments, cleanup methods automatically handle WORM-protected objects by:
    - Attempting to disable object lock configuration
    - Using `BypassGovernanceRetention=True` for deletions
    - Removing legal holds when present
    - Providing graceful fallbacks for different S3 implementations

    Attributes:
    - version (str): The version of the current benchmark of the format: `major.minor` (for example: `0.1`).
    - versions (list): A list of versions.
    - files (dict): A dictionary of files.
    - benchmark (str): The current benchmark.
    - auth_options (dict): The authentication options.
    - storage_options (dict): The storage options.

    Methods:
    - __init__(auth_options, benchmark): Initializes the RemoteStorage object.
    - _parse_benchmark(benchmark): Parses and validates the benchmark name.
    - _parse_auth_options(auth_options): Parses and validates the authentication options.
    - connect(): Connects to the remote storage.
    - _test_connect(): Tests the connection to the remote storage.
    - _get_containers(): Retrieves the available containers.
    - _get_benchmarks(): Retrieves the available benchmarks.
    - _create_benchmark(benchmark): Creates a new benchmark.
    - _get_versions(): Retrieves the available versions of the current benchmark.
    - _parse_version(version): Parses and validates the current version.
    - set_version(version): Sets the current version of the benchmark.
    - _create_new_version(): Creates a new version of the benchmark.
    - _get_objects(readonly): Retrieves the objects in the storage for the current benchmark.
    - archive_version(version): Archives a specific benchmark version.
    """

    def __init__(
        self,
        auth_options: Dict,
        benchmark: str,
        storage_options: StorageOptions,
    ):
        self.version = None
        self.versions = list()
        self._files: dict = {}
        self._files_loaded: bool = False
        self._parse_benchmark(benchmark)
        self._parse_auth_options(auth_options)
        self._parse_storage_options(storage_options)

    @property
    def files(self) -> dict:
        """File metadata dict; loaded on first access after each set_version() call."""
        if not self._files_loaded:
            self.load_objects()
        return self._files

    @files.setter
    def files(self, value: dict) -> None:
        self._files = value
        self._files_loaded = True

    def _parse_benchmark(self, benchmark: str) -> None:
        if not isinstance(benchmark, str):
            # XXX this is what strict type checking is for :)
            raise RemoteStorageInvalidInputException("benchmark must be a string")
        self.benchmark = benchmark

    def _parse_auth_options(self, auth_options: Dict) -> None:
        if not isinstance(auth_options, dict):
            raise RemoteStorageInvalidInputException(
                "auth_options must be a dictionary"
            )
        self.auth_options = auth_options

    def _parse_storage_options(self, storage_options: StorageOptions) -> None:
        expected_attributes = [t for t in dir(storage_options) if not t.startswith("_")]
        if not all([hasattr(storage_options, t) for t in expected_attributes]):
            raise RemoteStorageInvalidInputException(
                "storage_options must be a StorageOptions object"
            )
        self.storage_options = storage_options

    @abstractmethod
    def connect(self, readonly: bool = False):
        """
        Connects to the storage Service.

        Args:
            readonly: When True, connect without write credentials.

        Returns:
            - Storage Client / Service
        """
        raise NotImplementedError

    @abstractmethod
    def _test_connect(self) -> None:
        """
        This method is used to test the connection to the remote storage.
        It does not perform any actual operations, but checks if the connection is successful.
        """
        raise NotImplementedError

    @abstractmethod
    def _create_benchmark(self, benchmark: str) -> None:
        """
        Creates a new benchmark in the remote storage.

        Args:
            benchmark: The name of the benchmark.
        """
        raise NotImplementedError

    @abstractmethod
    def _get_versions(self) -> None:
        """
        Retrieves the benchmark versions in the remote storage and
        populates ``self.versions``.
        """
        self.versions = []

    def _parse_version(self, version: Union[None, str] = None) -> Version:
        """
        Parses the current version based on the provided major and minor versions.

        Args:
            version (None or str, optional]): The version number. Defaults to None which returns the latest version.

        Returns:
            Version: The parsed version.

        Raises:
            RemoteStorageInvalidInputException: If an invalid version is provided.
        """
        if self.versions == []:
            self._get_versions()
        if version is None:
            if self.versions == []:
                version_out = Version("0.1")
            else:
                max_current_version = max(self.versions)
                version_out = Version(
                    f"{max_current_version.major}.{max_current_version.minor + 1}"
                )
        else:
            version_out = Version(version)
        return version_out

    def set_version(self, version: Union[None, str] = None) -> None:
        """
        Sets the current version of the remote storage.

        Args:
            version (None or str, optional): The version number as a string. Defaults to None which sets the latest version.

        Raises:
            RemoteStorageInvalidInputException: If version does not exist in the available versions.
        """
        self.version = self._parse_version(version)
        self._files_loaded = False

    @abstractmethod
    def create_new_version(
        self, benchmark: "Optional[BenchmarkExecution]" = None
    ) -> None:
        """
        Creates a new version of the benchmark.

        Args:
            benchmark: Optional BenchmarkExecution object. When provided,
                uploads the benchmark YAML and software files and uses
                git-aware version tracking.
        """
        raise NotImplementedError

    @abstractmethod
    def load_objects(self):
        """
        Retrieves the objects in the storage for the current benchmark version.
        """
        raise NotImplementedError

    @abstractmethod
    def download_object(self, object_name: str, local_path: str):
        """
        Downloads an object from the storage.

        Args:
            object_name (str): The name of the object.
            local_path (str): The local path to download the object.
        """
        raise NotImplementedError

    @abstractmethod
    def archive_version(
        self,
        benchmark: "BenchmarkExecution",
        outdir: Path = Path(),
        config: bool = True,
        code: bool = False,
        software: bool = False,
        results: bool = False,
    ):
        """
        Archives/Freezes a specific benchmark version.

        Args:
            benchmark: The benchmark execution context.
            outdir: Local directory for the archive output.
            config: Include config files.
            code: Include code files.
            software: Include software environment files.
            results: Include result files.
        """
        raise NotImplementedError
