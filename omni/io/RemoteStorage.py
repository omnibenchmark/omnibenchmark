"""Base class for remote storage."""

from typing import Dict, Union

from packaging.version import Version


class RemoteStorage:
    """
    A class representing a remote storage.

    Attributes:
    - major_version (int): The major version of the storage.
    - major_version_new (int): The new major version of the storage.
    - minor_version (int): The minor version of the storage.
    - minor_version_new (int): The new minor version of the storage.
    - benchmarks (list): A list of benchmarks.
    - versions (list): A list of versions.
    - files (dict): A dictionary of files.
    - benchmark (str): The current benchmark.
    - auth_options (dict): The authentication options.

    Methods:
    - __init__(auth_options, benchmark): Initializes the RemoteStorage object.
    - _parse_benchmark(benchmark): Parses and validates the benchmark name.
    - _parse_auth_options(auth_options): Parses and validates the authentication options.
    - connect(): Connects to the remote storage.
    - _test_connect(): Tests the connection to the remote storage.
    - _get_benchmarks(): Retrieves the available benchmarks.
    - create_benchmark(benchmark): Creates a new benchmark.
    - _get_versions(): Retrieves the available versions of the current benchmark.
    - set_current_version(major_version, minor_version): Sets the current version of the benchmark.
    - _parse_new_version(major_version, minor_version): Parses and validates the new version.
    - set_new_version(major_version_new, minor_version_new): Sets the new version of the benchmark.
    - create_new_version(): Creates a new version of the benchmark.
    - _get_objects(): Retrieves the objects in the storage for the current benchmark.
    - find_objects_to_copy(reference_time): Finds objects to copy based on a reference time.
    - copy_objects(): Copies the objects from the current benchmark version to the new benchmark version.
    - archive_version(major_version, minor_version): Archives a specific benchmark version.
    - delete_version(major_version, minor_version): Deletes a specific benchmark version.
    """

    def __init__(self, auth_options: Dict, benchmark: str):
        self.version = None
        self.version_new = None
        # self.major_version = None
        # self.major_version_new = None
        # self.minor_version = None
        # self.minor_version_new = None
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

    def connect(self):
        """
        Connects to the storage Service.

        Returns:
            - Storage Client / Service
        """
        NotImplementedError

    def _test_connect(self) -> None:
        """
        This method is used to test the connection to the remote storage.
        It does not perform any actual operations, but checks if the connection is successful.
        """
        NotImplementedError

    def _get_container(self) -> None:
        """
        Retrieves the list of containers from the remote storage.

        Returns:
        - A list of container names.
        """
        NotImplementedError

    def _get_benchmarks(self) -> None:
        """
        Retrieves the benchmarks from the remote storage.

        Returns:
            A list of benchmarks.
        """
        NotImplementedError

    def create_benchmark(self, benchmark: str) -> None:
        """
        Creates a new benchmark in the remote storage.

        Args:
            benchmark: The benchmark object to be created.

        Returns:
            None
        """
        NotImplementedError

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
            version (None or str, optional]): The version number. Defaults to None.

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
            version (None or str, optional): The version number as a string. Defaults to None.

        Raises:
            ValueError: If version does not exist in the available versions.
        """
        self.version = self._parse_current_version(version)

    def _parse_new_version(self, version: Union[None, str] = None) -> Version:
        """
        Parses the new version based on the provided major and minor versions.

        Args:
            version (None or str, optional]): The version number. Defaults to None.

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
            version_new (None or str, optional): The version number of the new version. Defaults to None.

        Raises:
            ValueError: If no version is provided or if the version already exists.
        """
        self.version_new = self._parse_new_version(version_new)

    def create_new_version(self):
        NotImplementedError

    def _get_objects(self):
        NotImplementedError

    def find_objects_to_copy(self, reference_time):
        NotImplementedError

    def copy_objects(self):
        NotImplementedError

    def archive_version(self, version):
        NotImplementedError

    def delete_version(self, version):
        NotImplementedError
