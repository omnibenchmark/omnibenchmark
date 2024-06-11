"""Base class for remote storage."""

from abc import ABCMeta, abstractmethod
from typing import Dict, Union

from packaging.version import Version


class RemoteStorage(metaclass=ABCMeta):
    """
    A class representing a remote storage.

    Attributes:
    - version (str): The version of the current benchmark of the format: `major.minor` (for example: `0.1`).
    - version_new (str): The version of the new benchmark of the format: `major.minor` (for example: `0.1`).
    - benchmarks (list): A list of benchmarks.
    - versions (list): A list of versions.
    - other_versions (list): A list of other (non standard, such as tests) versions.
    - files (dict): A dictionary of files.
    - benchmark (str): The current benchmark.
    - auth_options (dict): The authentication options.

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
    - _parse_current_version(version): Parses and validates the current version.
    - set_current_version(version): Sets the current version of the benchmark.
    - _parse_new_version(version): Parses and validates the new version.
    - set_new_version(version): Sets the new version of the benchmark.
    - _update_overview(): Updates the overview of the benchmark.
    - _create_new_version(): Creates a new version of the benchmark.
    - _get_objects(readonly): Retrieves the objects in the storage for the current benchmark.
    - find_objects_to_copy(reference_time, tagging_type): Finds objects to copy based on a reference time.
    - copy_objects(type): Copies the objects from the current benchmark version to the new benchmark version.
    - create_new_version(version_new, tagging_type, copy_type): Creates a new version of the benchmark and copies the objects.
    - archive_version(version): Archives a specific benchmark version.
    - delete_version(version): Deletes a specific benchmark version.
    """

    def __init__(self, auth_options: Dict, benchmark: str):
        self.version = None
        self.version_new = None
        self.benchmarks = list()
        self.versions = list()
        self.other_versions = list()
        self.files = dict()
        self._parse_benchmark(benchmark)
        self._parse_auth_options(auth_options)

    def _parse_benchmark(self, benchmark: str) -> None:
        if not type(benchmark) is str:
            raise ValueError("benchmark must be a string")
        if not benchmark.isidentifier():
            raise ValueError("benchmark must be a valid identifier")
        self.benchmark = benchmark

    def _parse_auth_options(self, auth_options: Dict) -> None:
        if not type(auth_options) is dict:
            raise ValueError("auth_options must be a dictionary")
        self.auth_options = auth_options

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
    def _get_containers(self) -> None:
        """
        Retrieves the list of containers from the remote storage and saves them in the `containers` attribute.
        """
        NotImplementedError

    @abstractmethod
    def _get_benchmarks(self) -> None:
        """
        Retrieves the benchmarks from the remote storage and saves them in the `benchmarks` attribute.
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

    def _parse_current_version(self, version: Union[None, str] = None) -> Version:
        """
        Parses the current version based on the provided major and minor versions.

        Args:
            version (None or str, optional]): The version number. Defaults to None which returns the latest version.

        Returns:
            Version: The parsed version.

        Raises:
            ValueError: If an invalid version is provided.
        """
        if self.versions == []:
            self._get_versions()
        if version is None:
            version_out = max(self.versions)
        else:
            version_out = Version(version)
            if version_out not in self.versions:
                raise ValueError(f"Version {version_out} not found in {self.benchmark}")
        return version_out

    def set_current_version(self, version: Union[None, str] = None) -> None:
        """
        Sets the current version of the remote storage.

        Args:
            version (None or str, optional): The version number as a string. Defaults to None which sets the latest version.

        Raises:
            ValueError: If version does not exist in the available versions.
        """
        self.version = self._parse_current_version(version)

    def _parse_new_version(self, version: Union[None, str] = None) -> Version:
        """
        Parses the new version based on the provided major and minor versions.

        Args:
            version (None or str, optional]): The version number. Defaults to None which increments the minor version of the current version.

        Returns:
            Version: The parsed version.

        Raises:
            ValueError: If an invalid version is provided.
        """
        if self.versions == []:
            self._get_versions()
        if version is None:
            version_current = max(self.versions)
            version_out = Version(
                f"{version_current.major}.{version_current.minor + 1}"
            )
        else:
            version_out = Version(version)
            if version_out <= max(self.versions):
                raise ValueError(f"Version {version} not newest version")
        return version_out

    def set_new_version(
        self,
        version_new: Union[None, str] = None,
    ) -> None:
        """
        Sets a new version for the remote storage.

        Args:
            version_new (None or str, optional): The version number of the new version. Defaults to None which increments the minor version of the current version.

        Raises:
            ValueError: If no version is provided or if the version already exists.
        """
        self.version_new = self._parse_new_version(version_new)

    @abstractmethod
    def _update_overview(self, cleanup: bool = True):
        """
        Updates the overview of versions of the benchmark.

        Args:
            cleanup (bool, optional): Whether to clean up, i.e. remove not existing, versions of the benchmark. Defaults to True.
        """
        NotImplementedError

    @abstractmethod
    def _create_new_version(self):
        """
        Creates a new version of the benchmark.
        """
        NotImplementedError

    @abstractmethod
    def _get_objects(self, readonly=False):
        """
        Retrieves the objects in the storage for the current benchmark version.

        Args:
            readonly (bool, optional): Whether to retrieve the objects in read-only mode. Defaults to False.
        """
        NotImplementedError

    @abstractmethod
    def find_objects_to_copy(self, reference_time=None, tagging_type="all"):
        """
        Finds objects to copy based on the reference time and tagging type.

        Args:
            reference_time (datetime.datetime, optional): The reference time to compare with the object's time. Defaults to None.
            tagging_type (str, optional): The tagging type to consider. Defaults to "all".

        Raises:
            ValueError: If the tagging type is invalid or the reference time is not a datetime object.
        """
        NotImplementedError

    @abstractmethod
    def copy_objects(self, type="copy"):
        """
        Copy objects from the current version to the new version.

        Args:
            type (str): The type of copying to perform. Valid values are "copy" and "symlink".
        """
        NotImplementedError

    @abstractmethod
    def create_new_version(
        self,
        version_new: Union[None, str] = None,
        tagging_type: str = "all",
        copy_type: str = "copy",
    ):
        """
        Wrapper to create new version and copy objects.

        Args:
            version_new (str or None, optional): The new version to create. Defaults to None which increments the minor version of the current version.
            tagging_type (str, optional): The type of tagging to apply. Defaults to "all".
            copy_type (str, optional): The type of copying to perform. Defaults to "copy".
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
