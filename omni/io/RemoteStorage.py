"""Base class for remote storage."""
from typing import Union, List, Dict

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
        self.major_version = None
        self.major_version_new = None
        self.minor_version = None
        self.minor_version_new = None
        self.benchmarks = list()
        self.versions = list()
        self.other_versions = list()
        self.files = dict()
        self._parse_benchmark(benchmark)
        self._parse_auth_options(auth_options)

    def _parse_benchmark(self, benchmark: str):
        if not type(benchmark) is str:
            raise ValueError("benchmark must be a string")
        if not benchmark.isidentifier():
            raise ValueError("benchmark must be a valid identifier")
        self.benchmark = benchmark

    def _parse_auth_options(self, auth_options: Dict):
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

    def _test_connect(self):
        """
        This method is used to test the connection to the remote storage.
        It does not perform any actual operations, but checks if the connection is successful.
        """
        NotImplementedError

    def _get_container(self):
        """
        Retrieves the list of containers from the remote storage.

        Returns:
        - A list of container names.
        """
        NotImplementedError

    def _get_benchmarks(self):
            """
            Retrieves the benchmarks from the remote storage.

            Returns:
                A list of benchmarks.
            """
            NotImplementedError

    def create_benchmark(self, benchmark: str):
            """
            Creates a new benchmark in the remote storage.

            Args:
                benchmark: The benchmark object to be created.

            Returns:
                None
            """
            NotImplementedError

    def _get_versions(self, update: bool = True, readonly: bool = False):
            """
            Retrieves the benchmark versions in the remote storage.

            Args:
                update (bool, optional): Whether to update the versions. Defaults to True.
                readonly (bool, optional): Whether to retrieve the versions in read-only mode. Defaults to False.
            """
            self.versions = []

    def set_current_version(self, major_version: Union[None, int] = None, minor_version: Union[None, int] = None):
        """
        Sets the current version of the remote storage.

        Args:
            major_version (None or int, optional): The major version number. Defaults to None.
            minor_version (None or int, optional): The minor version number. Defaults to None.

        Raises:
            ValueError: If major_version is not an integer or None.
            ValueError: If minor_version is not an integer or None.
            ValueError: If major_version does not exist in the available versions.
            ValueError: If minor_version does not exist in the available versions for the specified major version.
        """

        if not type(major_version) is int and not major_version is None:
            raise ValueError("major_version must be an integer")
        if not type(minor_version) is int and not minor_version is None:
            raise ValueError("minor_version must be an integer")

        if major_version is None or minor_version is None:
            if self.versions is None:
                self._get_versions()
        if major_version is None:
            self.major_version = int(sorted(self.versions, key=lambda x: float(x), reverse=True)[0].split(".")[0])
        else:
            if not major_version in [int(v.split(".")[0]) for v in self.versions]:
                raise ValueError("Major version does not exist")
            else:
                self.major_version = major_version

        valid_minor_versions = [int(v.split(".")[1]) for v in self.versions if int(v.split(".")[0]) == int(self.major_version)]
        if minor_version is None:
            self.minor_version = sorted(valid_minor_versions, reverse=True)[0]
        else:
            if not minor_version in valid_minor_versions:
                raise ValueError("Minor version does not exist")
            else:
                self.minor_version = minor_version
 
    def _parse_new_version(self, major_version: Union[None, int] = None, minor_version: Union[None, int] = None):
        """
        Parses the new version based on the provided major and minor versions.

        Args:
            major_version (None or int, optional]): The major version number. Defaults to None.
            minor_version (None or int, optional): The minor version number. Defaults to None.

        Returns:
            Tuple[int, int]: A tuple containing the parsed major and minor versions.
        
        Raises:
            ValueError: If the major_version or minor_version is not of the expected type.
            ValueError: If an invalid version is provided.
        """

        if not (type(major_version) is int or type(major_version) is bool or major_version is None):
            raise ValueError("major_version must be an integer")
        if not (type(minor_version) is int or type(minor_version) is bool or minor_version is None):
            raise ValueError("minor_version must be an integer")

        if major_version is None or minor_version is None or self.major_version is None or self.minor_version is None:
            if self.versions is None:
                self._get_versions()
            if self.major_version is None or self.minor_version is None:
                self.set_current_version()

        if major_version is None and minor_version is None:
            major_version_p = int(self.major_version)
            minor_version_p = int(self.minor_version) + 1
        # increment major version
        elif (type(major_version) is bool and major_version) and (minor_version is None or (type(minor_version) is bool and minor_version)):
            major_version_p = int(self.major_version) + 1
            minor_version_p = 0
        # increment minor version
        elif type(minor_version) is bool and minor_version and major_version is None:
            major_version_p = int(self.major_version)
            minor_version_p = int(self.minor_version) + 1
        # set version
        elif type(major_version) is int and type(minor_version) is int:
            major_version_p = major_version
            minor_version_p = minor_version
        elif type(major_version) is str or type(minor_version) is str or type(major_version) is float or type(minor_version) is float:
            if int(major_version) == major_version:
                major_version_p = int(major_version)
            if int(minor_version) == minor_version:
                minor_version_p = int(minor_version)
        else:
            raise ValueError("Invalid version")
        return major_version_p, minor_version_p
    
    def set_new_version(self, major_version_new: Union[None, int] = None, minor_version_new: Union[None, int] = None):
        """
        Sets a new version for the remote storage.

        Args:
            major_version_new (None or int, optional): The major version number of the new version. Defaults to None.
            minor_version_new (None or int, optional): The minor version number of the new version. Defaults to None.

        Raises:
            ValueError: If no version is provided or if the version already exists.
        """
        major_version_new, minor_version_new = self._parse_new_version(major_version_new, minor_version_new)
        if major_version_new is None or minor_version_new is None:
            raise ValueError("No version provided")
        if f"{major_version_new}.{minor_version_new}" in self.versions:
            raise ValueError("Version already exists")

        self.major_version_new = major_version_new
        self.minor_version_new = minor_version_new


    def create_new_version(self):
        NotImplementedError

    def _get_objects(self):
        NotImplementedError

    def find_objects_to_copy(self, reference_time):
        NotImplementedError

    def copy_objects(self):
        NotImplementedError

    def archive_version(self, major_version, minor_version):
        NotImplementedError

    def delete_version(self, major_version, minor_version):
        NotImplementedError

