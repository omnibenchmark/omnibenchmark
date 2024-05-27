"""S3 class for remote storage."""
import boto3
import datetime
import json
import dateutil.parser
import requests
from bs4 import BeautifulSoup
from urllib.parse import urlparse
import re
import logging
from omni.io.RemoteStorage import RemoteStorage

logging.basicConfig(level=logging.ERROR)
logging.getLogger("requests").setLevel(logging.DEBUG)
logging.getLogger("minio").setLevel(logging.DEBUG)
logger = logging.getLogger(__name__)


class S3Storage(RemoteStorage):
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
        Connects to the Swift storage.

        Returns:
            SwiftService: The SwiftService object representing the connection to the Swift storage.
        """
        if (
            "endpoint" in self.auth_options.keys()
            and "access_key" in self.auth_options.keys()
            and "secret_key" in self.auth_options.keys()
        ):
            return boto3.resource(
                "s3",
                endpoint_url=self.auth_options["endpoint"],
                aws_access_key_id=self.auth_options["access_key"],
                aws_secret_access_key=self.auth_options["secret_key"],
                aws_session_token=None,
                config=boto3.session.Config(signature_version="s3v4"),
                verify=False,
            )
        else:
            raise ValueError("Invalid auth options")

    def _test_connect(self):
        try:
            _ = list(self.client.buckets.all())
        except Exception as e:
            raise e

    def _get_containers(self):
        """
        Retrieves the list of containers from the Swift storage.

        Returns:
        - A list of container names.
        """
        try:
            containers = list()
            for bucket in self.client.buckets.iterator():
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

    def _set_bucket_public_readonly(self, bucket_name):
        policy = {
            "Version": "2012-10-17",
            "Statement": [
                {
                    "Effect": "Allow",
                    "Principal": {"AWS": "*"},
                    "Action": [
                        "s3:GetBucketLocation",
                        "s3:ListBucket",
                        "s3:ListObjects",
                    ],
                    "Resource": f"arn:aws:s3:::{bucket_name}",
                },
                {
                    "Effect": "Allow",
                    "Principal": {"AWS": "*"},
                    "Action": "s3:GetObject",
                    "Resource": f"arn:aws:s3:::{bucket_name}/*",
                },
            ],
        }
        self.client.Bucket(bucket_name).Policy().put(Policy=json.dumps(policy))

    def _create_benchmark(self, benchmark, update=True):
        if len(self.benchmarks) == 0 or update:
            self._get_benchmarks()
        if benchmark in self.benchmarks:
            raise ValueError("Benchmark already exists")
        # create new version
        self.client.create_bucket(Bucket=f"{benchmark}.test.1")
        self._set_bucket_public_readonly(f"{benchmark}.test.1")
        if self.client.Bucket(f"{benchmark}.test.1").creation_date is None:
            raise Exception(f"Benchmark creation of {benchmark}.test.1 failed")
        self.client.create_bucket(Bucket=f"{benchmark}.overview")
        self._set_bucket_public_readonly(f"{benchmark}.overview")
        if self.client.Bucket(f"{benchmark}.overview").creation_date is None:
            raise Exception(f"Benchmark creation of {benchmark}.overview failed")
        self.client.create_bucket(Bucket=f"{benchmark}.0.1")
        self._set_bucket_public_readonly(f"{benchmark}.0.1")
        if self.client.Bucket(f"{benchmark}.0.1").creation_date is None:
            raise Exception(f"Benchmark creation of {benchmark}.0.1 failed")
        self._update_overview(cleanup=True)
        # except SwiftError as e:
        #     logger.error(e.value)

    def _get_versions(self, update=True, readonly=False):
        if not "secret_key" in self.auth_options.keys() or readonly:
            url = urlparse(f"{self.auth_options['endpoint']}/{self.benchmark}.overview")
            if self.auth_options["secure"]:
                url = url._replace(scheme="https")
            else:
                url = url._replace(scheme="http")
            params = {"format": "xml"}
            response = requests.get(url=url.geturl(), params=params)
            print(response.text)
            print(response.ok)
            if response.ok:
                response_text = response.text
            else:
                response.raise_for_status()
            response_text
            soup = BeautifulSoup(response_text, "xml")
            print(soup)
            allversions = [obj.find("Key").text for obj in soup.find_all("Contents")]
            print(allversions)
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

            # with self.connect() as swift:
            #     versions = list()
            #     for page in swift.list(container=f"{self.benchmark}.overview"):
            #         for element in page['listing']:
            #             versions.append(element['name'])
            # self.versions = versions
            # return self.versions
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
        for element in self.client.Bucket(
            f"{self.benchmark}.overview"
        ).objects.iterator():
            versions_in_overview.append(element.key)
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
            result = self.client.Bucket(f"{self.benchmark}.overview").put_object(
                Key=v, Body=b""
            )

        # remove unavailable versions
        if cleanup:
            versions_in_overview = []
            for element in self.client.Bucket(
                f"{self.benchmark}.overview"
            ).objects.iterator():
                versions_in_overview.append(element.key)

            versions_in_overview_but_unavailable = [
                v
                for v in versions_in_overview
                if v not in self.versions + self.other_versions
            ]
            if len(versions_in_overview_but_unavailable) > 0:
                deletes = self.client.Bucket(
                    f"{self.benchmark}.overview"
                ).delete_objects(
                    Delete={
                        "Objects": [
                            {"Key": v} for v in versions_in_overview_but_unavailable
                        ]
                    }
                )

                if "Deleted" in deletes.keys():
                    deleted = [delete["Key"] for delete in deletes["Deleted"]]
                    for v in versions_in_overview_but_unavailable:
                        if v not in deleted:
                            raise Exception(f"Deletion failed: {v}")

    def _create_new_version(self):
        if self.major_version_new is None or self.minor_version_new is None:
            raise ValueError("No version provided")

        # update
        self._get_versions()
        # check if version exists
        if not f"{self.major_version_new}.{self.minor_version_new}" in self.versions:
            # create new version
            self.client.create_bucket(
                Bucket=f"{self.benchmark}.{self.major_version_new}.{self.minor_version_new}"
            )
            self._set_bucket_public_readonly(
                f"{self.benchmark}.{self.major_version_new}.{self.minor_version_new}"
            )
            if (
                self.client.Bucket(
                    f"{self.benchmark}.{self.major_version_new}.{self.minor_version_new}"
                ).creation_date
                is None
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

    def _get_meta_mtime(self, preauthurl, containername, objectname):
        """
        Retrieves the metadata modification time and access time of a file from a Swift storage.

        Args:
            preauthurl (str): The pre-authenticated URL of the Swift storage.
            containername (str): The name of the container where the file is stored.
            objectname (str): The name of the file.

        Returns:
            tuple: A tuple containing the file's metadata modification time and access time.

        Raises:
            requests.HTTPError: If the HTTP request to retrieve the file fails.
        """
        urlfile = f"{preauthurl}/{containername}/{objectname}"
        urlfile = "http://omnibenchmark.mls.uzh.ch:9000/tb.0.1/test"

        response = requests.get(url=urlfile)
        if response.ok:
            response_headers = response.headers
            if "X-Object-Meta-Mtime" in response_headers.keys():
                file_time = datetime.datetime.fromtimestamp(
                    float(response_headers["X-Object-Meta-Mtime"]),
                    datetime.timezone.utc,
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
            response = requests.get(url=url.geturl(), params=params)
            print(response)
            print(response.ok)
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
                mtime, accesstime = self._get_meta_mtime(
                    url.geturl(), containername, obj.find("Key").text
                )
                files[obj.find("Key").text] = {
                    "hash": obj.find("ETag").text,
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

            for element in self.client.Bucket(
                f"{self.benchmark}.{self.major_version}.{self.minor_version}"
            ).objects.iterator():
                print(element)
                if type(element.last_modified) is str:
                    files[element.key] = {
                        "size": element.size,
                        "last_modified": element.last_modified,
                        "hash": element.e_tag.replace('"', ""),
                        "symlink_path": "",
                    }
                elif type(element.last_modified) is datetime.datetime:
                    files[element.key] = {
                        "size": element.size,
                        "last_modified": element.last_modified.strftime(
                            "%Y-%m-%dT%H:%M:%S.%f"
                        ),
                        "hash": element.e_tag.replace('"', ""),
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
                out = self.client.Object(
                    f"{self.benchmark}.{self.major_version_new}.{self.minor_version_new}",
                    filename,
                ).copy_from(
                    CopySource={
                        "Bucket": f"{self.benchmark}.{self.major_version}.{self.minor_version}",
                        "Key": filename,
                    }
                )
                self.files[filename]["copied"] = True
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
