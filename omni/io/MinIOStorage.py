"""MinIO class for remote storage."""
import minio
import io
import datetime
import json
import dateutil.parser
import requests
from bs4 import BeautifulSoup
from urllib.parse import urlparse
import re
import logging
from omni.io.RemoteStorage import RemoteStorage
from omni.io.S3config import bucket_readonly_policy

logging.basicConfig(level=logging.ERROR)
logging.getLogger("requests").setLevel(logging.DEBUG)
logging.getLogger("minio").setLevel(logging.DEBUG)
logger = logging.getLogger(__name__)


def get_meta_mtime(preauthurl, containername, objectname):
    """
    Retrieves the metadata modification time and access time of a file from a MinIO storage.

    Args:
        preauthurl (str): The pre-authenticated URL of the MinIO storage.
        containername (str): The name of the container where the file is stored.
        objectname (str): The name of the file.

    Returns:
        tuple: A tuple containing the file's metadata modification time and access time.

    Raises:
        requests.HTTPError: If the HTTP request to retrieve the file fails.
    """
    urlfile = f"{preauthurl}/{containername}/{objectname}"
    response = requests.get(urlfile)
    if response.ok:
        response_headers = response.headers
        if "X-Object-Meta-Mtime" in response_headers.keys():
            file_time = datetime.datetime.fromtimestamp(
                float(response_headers["X-Object-Meta-Mtime"]), datetime.timezone.utc
            )
        elif "X-Object-Meta-Last-Modified" in response_headers.keys():
            file_time = dateutil.parser.parse(
                response_headers["X-Object-Meta-Last-Modified"]
            )
        else:
            file_time = dateutil.parser.parse(response_headers["Last-Modified"])
        access_time = dateutil.parser.parse(response_headers["Date"])
        return file_time, access_time
    else:
        response.raise_for_status()


def set_bucket_public_readonly(client, bucket_name):
    policy = bucket_readonly_policy(bucket_name)
    client.set_bucket_policy(bucket_name, json.dumps(policy))


class MinIOStorage(RemoteStorage):
    def __init__(self, auth_options, benchmark):
        super().__init__(auth_options, benchmark)
        self.containers = list()
        if "access_key" in self.auth_options.keys():
            self.client = self.connect()
            self._test_connect()
            self._get_benchmarks()

            if not benchmark in self.benchmarks:
                logger.warning(
                    f"Benchmark {benchmark} does not exist, creating new benchmark."
                )
                self._create_benchmark(benchmark)
                self._get_containers()
                self._get_benchmarks()

            self._get_versions()
        else:
            self._get_versions(readonly=True)
            # read-only mode
            pass

    def connect(self):
        """
        Connects to the MinIO storage.

        Returns:
        - A MinIO client object.
        """
        if (
            "endpoint" in self.auth_options.keys()
            and "access_key" in self.auth_options.keys()
            and "secret_key" in self.auth_options.keys()
        ):
            try:
                return minio.Minio(**self.auth_options)
            except Exception as e:
                tmp_auth_options = self.auth_options.copy()
                url = urlparse(tmp_auth_options["endpoint"])
                tmp_auth_options["endpoint"] = url.netloc
                return minio.Minio(**tmp_auth_options)
        else:
            raise ValueError("Invalid auth options")

    def _test_connect(self):
        try:
            _ = self.client.list_buckets()
        except Exception as e:
            raise e

    def _get_containers(self):
        """
        Retrieves the list of containers.

        Returns:
        - A list of container names.
        """
        try:
            containers = list()
            for bucket in self.client.list_buckets():
                containers.append(bucket.name)
            self.containers = containers
            return self.containers

        except Exception as e:
            logger.error(e.value)

    def _get_benchmarks(self, update=True):
        if len(self.containers) == 0 or update:
            self._get_containers()
        benchmarks = list()
        for con in self.containers:
            # remove major and minor version from benchmark name
            bm = "".join(con.split(".")[:-2])
            if bm not in benchmarks and bm != "":
                benchmarks.append(bm)
        self.benchmarks = benchmarks
        return self.benchmarks

    def _create_benchmark(self, benchmark, update=True):
        if len(self.benchmarks) == 0 or update:
            self._get_benchmarks()
        if benchmark in self.benchmarks:
            raise ValueError("Benchmark already exists")
        # create new version
        self.client.make_bucket(bucket_name=f"{benchmark}.test.1")
        set_bucket_public_readonly(self.client, f"{benchmark}.test.1")
        if not self.client.bucket_exists(f"{benchmark}.test.1"):
            raise Exception(f"Benchmark creation of {benchmark}.test.1 failed")
        self.client.make_bucket(bucket_name=f"{benchmark}.overview")
        set_bucket_public_readonly(self.client, f"{benchmark}.overview")
        if not self.client.bucket_exists(f"{benchmark}.overview"):
            raise Exception(f"Benchmark creation of {benchmark}.overview failed")

        self.client.make_bucket(bucket_name=f"{benchmark}.0.1")
        set_bucket_public_readonly(self.client, f"{benchmark}.0.1")
        if not self.client.bucket_exists(f"{benchmark}.0.1"):
            raise Exception(f"Benchmark creation of {benchmark}.0.1 failed")
        self._update_overview(cleanup=True)

    def _get_versions(self, update=True, readonly=False):
        if not "secret_key" in self.auth_options.keys() or readonly:
            url = urlparse(f"{self.auth_options['endpoint']}/{self.benchmark}.overview")
            if self.auth_options["secure"]:
                url = url._replace(scheme="https")
            else:
                url = url._replace(scheme="http")
            params = {"format": "xml"}
            response = requests.get(url.geturl(), params=params)
            if response.ok:
                response_text = response.text
            else:
                response.raise_for_status()
            soup = BeautifulSoup(response_text, "xml")
            allversions = [obj.find("Key").text for obj in soup.find_all("Contents")]
            versions = list()
            other_versions = list()
            for version in allversions:
                if re.match(f"(\\d+).(\\d+)", version):
                    versions.append(version)
                elif re.match(f"test.(\\d+)", version):
                    other_versions.append(version)
            self.versions = versions
            self.other_versions = other_versions
            return self.versions, self.other_versions

        else:
            if update:
                self._get_containers()
            versions = list()
            other_versions = list()
            for con in self.containers:
                if re.match(f"{self.benchmark}.(\\d+).(\\d+)", con):
                    versions.append(".".join(con.split(".")[-2:]))
                elif re.match(f"{self.benchmark}.test.(\\d+)", con):
                    other_versions.append(".".join(con.split(".")[-2:]))
            self.versions = versions
            self.other_versions = other_versions
            return self.versions, self.other_versions

    def _update_overview(self, cleanup=True):
        # check versions in overview
        versions_in_overview = []
        for element in self.client.list_objects(f"{self.benchmark}.overview"):
            versions_in_overview.append(element.object_name)
        # available benchmark versions
        self._get_versions()
        # missing versions
        versions_not_in_overview = [
            v
            for v in self.versions + self.other_versions
            if v not in versions_in_overview
        ]
        # create missing versions
        for v in versions_not_in_overview:
            result = self.client.put_object(
                f"{self.benchmark}.overview", v, io.BytesIO(b""), 0
            )

        # remove unavailable versions
        if cleanup:
            versions_in_overview = []
            for element in self.client.list_objects(f"{self.benchmark}.overview"):
                versions_in_overview.append(element.object_name)

            versions_in_overview_but_unavailable = [
                v
                for v in versions_in_overview
                if v not in self.versions + self.other_versions
            ]
            deletes = self.client.remove_objects(
                f"{self.benchmark}.overview", versions_in_overview_but_unavailable
            )
            for delete in deletes:
                raise Exception(f"Deletion failed: {delete}")

    def _create_new_version(self):
        if self.major_version_new is None or self.minor_version_new is None:
            raise ValueError("No version provided")

        # update
        self._get_versions()
        # check if version exists
        if not f"{self.major_version_new}.{self.minor_version_new}" in self.versions:
            # create new version
            self.client.make_bucket(
                bucket_name=f"{self.benchmark}.{self.major_version_new}.{self.minor_version_new}"
            )
            set_bucket_public_readonly(
                self.client,
                f"{self.benchmark}.{self.major_version_new}.{self.minor_version_new}",
            )
            if not self.client.bucket_exists(
                f"{self.benchmark}.{self.major_version_new}.{self.minor_version_new}"
            ):
                raise Exception(
                    f"Benchmark creation of {self.benchmark}.{self.major_version_new}.{self.minor_version_new} failed"
                )
            # update
            self._update_overview()
        else:
            raise ValueError("Version already exists")

        if not f"{self.major_version_new}.{self.minor_version_new}" in self.versions:
            raise ValueError("Version creation failed")

    def _get_objects(self, readonly=False):
        if self.major_version is None or self.minor_version is None:
            raise ValueError("No version provided")
        if not "secret_key" in self.auth_options.keys() or readonly:
            containername = (
                f"{self.benchmark}.{self.major_version}.{self.minor_version}"
            )
            url = urlparse(f"{self.auth_options['endpoint']}/{containername}")
            if self.auth_options["secure"]:
                url = url._replace(scheme="https")
            else:
                url = url._replace(scheme="http")
            params = {"format": "xml"}
            response = requests.get(url.geturl(), params=params)
            if response.ok:
                response_text = response.text
            else:
                response.raise_for_status()
            soup = BeautifulSoup(response_text, "xml")
            # names = [obj.find('Key').text for obj in soup.find_all('Contents')]
            names = soup.find_all("Contents")

            files = dict()
            for obj in names:
                url = urlparse(f"{self.auth_options['endpoint']}")
                if self.auth_options["secure"]:
                    url = url._replace(scheme="https")
                else:
                    url = url._replace(scheme="http")
                mtime, accesstime = get_meta_mtime(
                    url.geturl(), containername, obj.find("Key").text
                )
                files[obj.find("Key").text] = {
                    "hash": obj.find("ETag").text.replace('"', ""),
                    "size": obj.find("Size").text,
                    "last_modified": obj.find("LastModified").text,
                    "symlink_path": obj.find("symlink_path").text
                    if obj.find("symlink_path") is not None
                    else "",
                    "x-object-meta-mtime": mtime,
                    "accesstime": accesstime,
                }
            self.files = files

        else:
            files = dict()
            # get all objects

            for element in self.client.list_objects(
                f"{self.benchmark}.{self.major_version}.{self.minor_version}",
                recursive=True,
            ):
                print(element)
                if type(element.last_modified) is str:
                    files[element.object_name] = {
                        "size": element.size,
                        "last_modified": element.last_modified,
                        "hash": element.etag.replace('"', ""),
                        "symlink_path": "",
                    }
                elif type(element.last_modified) is datetime.datetime:
                    files[element.object_name] = {
                        "size": element.size,
                        "last_modified": element.last_modified.strftime(
                            "%Y-%m-%dT%H:%M:%S.%f"
                        ),
                        "hash": element.etag.replace('"', ""),
                        "symlink_path": "",
                    }
                else:
                    ValueError("Invalid last_modified")

            self.files = files

    def find_objects_to_copy(self, reference_time=None, tagging_type="all"):
        """
        Finds objects to copy based on the reference time and tagging type.

        Args:
            reference_time (datetime.datetime, optional): The reference time to compare with the object's time. Defaults to None.
            tagging_type (str, optional): The tagging type to consider. Defaults to "all".

        Raises:
            ValueError: If the tagging type is invalid or the reference time is not a datetime object.

        Returns:
            None
        """
        if tagging_type not in ["all"]:
            raise ValueError("Invalid tagging type")
        if tagging_type == "all":
            if reference_time is None:
                reference_time = datetime.datetime.now()
            elif not type(reference_time) is datetime.datetime:
                raise ValueError("Invalid reference time, must be datetime object")
            if len(self.files) == 0:
                self._get_objects()

            if len(self.files) == 0:
                return
            else:
                for filename in self.files.keys():
                    print(filename)
                    if "x-object-meta-mtime" in self.files[filename].keys():
                        file_time = datetime.datetime.fromtimestamp(
                            float(self.files[filename]["x-object-meta-mtime"])
                        )
                    elif "x-object-meta-last-modified" in self.files[filename].keys():
                        file_time = dateutil.parser.parse(
                            self.files[filename]["x-object-meta-last-modified"]
                        )
                    else:
                        file_time = dateutil.parser.parse(
                            self.files[filename]["last_modified"]
                        )
                    if file_time < reference_time:
                        self.files[filename]["copy"] = True
                    else:
                        self.files[filename]["copy"] = False

    def copy_objects(self, type="copy"):
        """
        Copy objects from the current version to the new version.

        Args:
            type (str): The type of copying to perform. Valid values are "copy" and "symlink".
        """
        if type not in ["copy", "symlink"]:
            raise ValueError("Invalid type")
        if (
            self.major_version is None
            or self.minor_version is None
            or self.major_version_new is None
            or self.minor_version_new is None
        ):
            raise ValueError("No version provided")

        if type == "copy":
            filenames = [
                filename
                for filename in self.files.keys()
                if self.files[filename]["copy"]
            ]
            for filename in filenames:
                out = self.client.copy_object(
                    f"{self.benchmark}.{self.major_version_new}.{self.minor_version_new}",
                    filename,
                    minio.commonconfig.CopySource(
                        bucket_name=f"{self.benchmark}.{self.major_version}.{self.minor_version}",
                        object_name=filename,
                    ),
                )
                self.files[out.object_name]["copied"] = True
        if type == "symlink":
            raise NotImplementedError("Symlink copying not implemented")

    def create_new_version(
        self, version_type="minor", tagging_type="all", copy_type="copy"
    ):
        """
        Wrapper to create new version and copy objects.

        Args:
            version_type (str, optional): The type of version to create. Can be "minor" or "major". Defaults to "minor".
            tagging_type (str, optional): The type of tagging to apply. Defaults to "all".
            copy_type (str, optional): The type of copying to perform. Defaults to "copy".
        """
        if version_type not in ["minor", "major"]:
            raise ValueError("Invalid version type: 'minor' or 'major'")
        self.set_current_version()
        if version_type == "minor":
            self.set_new_version(None, True)
        elif version_type == "major":
            self.set_new_version(True, None)

        self._create_new_version()
        self._get_objects()
        self.find_objects_to_copy(tagging_type=tagging_type)
        self.copy_objects(copy_type)

    def archive_version(self, major_version, minor_version):
        NotImplementedError
        # self._update_overview(cleanup=True)

    def delete_version(self, major_version, minor_version):
        NotImplementedError
        # self._update_overview(cleanup=True)
