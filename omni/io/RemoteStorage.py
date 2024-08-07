"""Base class for remote storage."""

from abc import ABCMeta, abstractmethod
from typing import Dict, Union

import packaging.version
from packaging.version import Version

from omni.io.exception import RemoteStorageInvalidInputException

DEFAULT_STORAGE_OPTIONS = {"tracked_directories": ["out", "versions"]}


def is_valid_version(version: str):
    try:
        packaging.version.parse(version)
        return True
    except packaging.version.InvalidVersion:
        return False


class RemoteStorage(metaclass=ABCMeta):
    """
    A class representing a remote storage.

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
    - delete_version(version): Deletes a specific benchmark version.
    """

    def __init__(self, auth_options: Dict, benchmark: str, storage_options: Dict = {}):
        self.version = None
        self.versions = list()
        self.files = dict()
        self._parse_benchmark(benchmark)
        self._parse_auth_options(auth_options)
        self._parse_storage_options(storage_options)

    def _parse_benchmark(self, benchmark: str) -> None:
        if not type(benchmark) is str:
            raise RemoteStorageInvalidInputException("benchmark must be a string")
        self.benchmark = benchmark

    def _parse_auth_options(self, auth_options: Dict) -> None:
        if not type(auth_options) is dict:
            raise RemoteStorageInvalidInputException(
                "auth_options must be a dictionary"
            )
        self.auth_options = auth_options

    def _parse_storage_options(self, storage_options: Dict) -> None:
        if not type(storage_options) is dict:
            raise RemoteStorageInvalidInputException(
                "storage_options must be a dictionary"
            )
        self.storage_options = {**DEFAULT_STORAGE_OPTIONS, **storage_options}

    @abstractmethod
    def connect(self):
        """
        Connects to the storage Service.

        Returns:
            - Storage Client / Service
        """
        NotImplementedError

    @abstractmethod
    def _test_connect(self) -> None:
        """
        This method is used to test the connection to the remote storage.
        It does not perform any actual operations, but checks if the connection is successful.
        """
        NotImplementedError

    @abstractmethod
    def _create_benchmark(self, benchmark: str, update: bool = True) -> None:
        """
        Creates a new benchmark in the remote storage.

        Args:
            benchmark: The name of the benchmark.
            update: Whether to update the list of benchmarks. Defaults to True.
        """
        NotImplementedError

    @abstractmethod
    def _get_versions(self, update: bool = True, readonly: bool = False) -> None:
        """
        Retrieves the benchmark versions in the remote storage.

        Args:
            update (bool, optional): Whether to update the versions. Defaults to True.
            readonly (bool, optional): Whether to retrieve the versions in read-only mode. Defaults to False.
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

    @abstractmethod
    def create_new_version(self):
        """
        Creates a new version of the benchmark.
        """
        NotImplementedError

    @abstractmethod
    def _get_objects(self):
        """
        Retrieves the objects in the storage for the current benchmark version.

        Args:
            readonly (bool, optional): Whether to retrieve the objects in read-only mode. Defaults to False.
        """
        NotImplementedError

    @abstractmethod
    def download_object(self, object_name: str, local_path: str):
        """
        Downloads an object from the storage.

        Args:
            object_name (str): The name of the object.
            local_path (str): The local path to download the object.
        """
        NotImplementedError

    @abstractmethod
    def archive_version(self, version):
        """
        Archives/Freezes a specific benchmark version.

        Args:
            version (str): The version to archive.
        """
        NotImplementedError

    @abstractmethod
    def delete_version(self, version):
        """
        Deletes a specific benchmark version.

        Args:
            version (str): The version to delete.
        """
        NotImplementedError
