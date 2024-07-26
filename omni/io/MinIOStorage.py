"""MinIO class for remote storage."""

import datetime
import io
import json
import logging
import os
import re
from typing import Union
from urllib.parse import urlparse

import dateutil.parser
import minio
import minio.commonconfig
import minio.deleteobjects
import minio.retention
import packaging.version
import requests
from bs4 import BeautifulSoup
from packaging.version import Version

from omni.io.RemoteStorage import RemoteStorage
from omni.io.S3config import bucket_readonly_policy

logging.basicConfig(level=logging.ERROR)
logging.getLogger("requests").setLevel(logging.DEBUG)
logging.getLogger("minio").setLevel(logging.DEBUG)
logger = logging.getLogger(__name__)


def valid_version(version: str):
    try:
        packaging.version.parse(version)
        return True
    except packaging.version.InvalidVersion:
        return False


def get_s3_object_versions_and_tags(client, benchmark):
    if not client.bucket_exists(benchmark):
        raise ValueError(f"Benchmark {benchmark} does not exist.")

    object_ls = list(
        client.list_objects(
            benchmark, include_version=True, include_user_meta=True, recursive=True
        )
    )
    object_names = [x.object_name for x in object_ls]
    version_ids = [x.version_id for x in object_ls]
    sizes = [x.size for x in object_ls]
    last_modified = [x.last_modified for x in object_ls]
    uq_names = list(set(object_names))

    di = {}
    for uq in uq_names:
        di[uq] = {}
        # tv = [v for v,n in zip(version_ids,object_names) if n == uq]
        tmpinds = [i for i, n in enumerate(object_names) if n == uq]
        tv = [version_ids[i] for i in tmpinds]
        ts = [sizes[i] for i in tmpinds]
        tlm = [last_modified[i] for i in tmpinds]
        tdl = [object_ls[i].is_delete_marker for i in tmpinds]
        for j, v in enumerate(tv):
            tags_filt = {}
            if not tdl[j]:
                tags = client.get_object_tags(benchmark, uq, version_id=v)
                if not tags is None:
                    for k, w in tags.items():
                        if valid_version(k):
                            tags_filt[k] = w
            di[uq][v] = {}
            di[uq][v]["tags"] = tags_filt
            tmpinds2 = [i for i, tvers in enumerate(tv) if tvers == v]
            di[uq][v]["size"] = ts[tmpinds2[0]]
            di[uq][v]["last_modified"] = tlm[tmpinds2[0]]
            di[uq][v]["is_delete_marker"] = tdl[tmpinds2[0]]
    return di


def get_s3_objects_to_tag(objdic, tracked_dirs=["out", "versions"]):
    object_names = sorted(list(objdic.keys()))
    newest_versions = []
    for uq in object_names:
        subdatetimels = []
        for v in objdic[uq].keys():
            subdatetimels.append(objdic[uq][v]["last_modified"])
        maxind = [k for k, sd in enumerate(subdatetimels) if sd == max(subdatetimels)][
            0
        ]
        newest_versions.append(list(objdic[uq].keys())[maxind])

    newest_version_exists = []
    for uq, v in zip(object_names, newest_versions):
        newest_version_exists.append(not objdic[uq][v]["is_delete_marker"])

    is_tracked_basedirls = []
    for n in object_names:
        while not os.path.split(n)[0] == "":
            n = os.path.split(n)[0]
        if n in tracked_dirs:
            is_tracked_basedirls.append(True)
        else:
            is_tracked_basedirls.append(False)

    do_tag = []
    for nve, itb in zip(newest_version_exists, is_tracked_basedirls):
        if nve and itb:
            do_tag.append(True)
        else:
            do_tag.append(False)
    object_names_to_tag = [n for n, dt in zip(object_names, do_tag) if dt]
    versionid_of_objects_to_tag = [
        newest_versions[i] for i, dt in enumerate(do_tag) if dt
    ]
    return object_names_to_tag, versionid_of_objects_to_tag


def get_single_s3version_from_bmversion(di, object_name, query_version):
    return list(
        filter(
            lambda v: query_version in di[object_name][v]["tags"].keys(),
            di[object_name].keys(),
        )
    )


def get_s3version_from_bmversion(di, query_version):
    vv_ls = []
    for uq in di.keys():
        tmpv = get_single_s3version_from_bmversion(di, uq, query_version)
        if len(tmpv) == 0:
            continue
        elif len(tmpv) > 1:
            raise ValueError("Multiple versions found")
        else:
            vv_ls.append(
                [uq, tmpv[0], di[uq][tmpv[0]]["last_modified"], di[uq][tmpv[0]]["size"]]
            )
    return vv_ls


def prepare_csv_s3version_from_bmversion(vv_ls):
    outstr = "name,version_id,last_modified,size\n"
    for vv in vv_ls:
        outstr += f"{vv[0]},{vv[1]},{vv[2]},{vv[3]}\n"
    return outstr


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


def set_bucket_lifecycle_config(client, bucket_name, noncurrent_days=1):
    lifecycle_config = minio.lifecycleconfig.LifecycleConfig(
        [
            minio.lifecycleconfig.Rule(
                minio.commonconfig.ENABLED,
                rule_filter=minio.lifecycleconfig.Filter(prefix="*"),
                rule_id="rule1",
                noncurrent_version_expiration=minio.lifecycleconfig.NoncurrentVersionExpiration(
                    noncurrent_days=noncurrent_days
                ),
            ),
        ],
    )
    client.set_bucket_lifecycle(bucket_name, lifecycle_config)


class MinIOStorage(RemoteStorage):
    def __init__(self, auth_options, benchmark):
        super().__init__(auth_options, benchmark)
        if "access_key" in self.auth_options.keys():
            self.client = self.connect()
            self._test_connect()

            if not self.client.bucket_exists(benchmark):
                logger.warning(
                    f"Benchmark {benchmark} does not exist, creating new benchmark."
                )
                self._create_benchmark(benchmark)

            self._get_versions()
        else:
            self._get_versions()
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

    def _test_connect(self) -> None:
        try:
            _ = self.client.list_buckets()
        except Exception as e:
            raise e

    def _create_benchmark(self, benchmark):
        if self.client.bucket_exists(benchmark):
            raise ValueError("Benchmark already exists")
        # create new version
        self.client.make_bucket(bucket_name=benchmark, object_lock=True)
        set_bucket_public_readonly(self.client, benchmark)
        set_bucket_lifecycle_config(self.client, benchmark, noncurrent_days=1)
        _ = self.client.put_object(
            self.benchmark, "versions/.ignore", io.BytesIO(b""), 0
        )
        if not self.client.bucket_exists(benchmark):
            raise Exception(f"Bucket creation for benchmark {benchmark} failed")

    def _get_versions(self):
        url = urlparse(f"{self.auth_options['endpoint']}/{self.benchmark}")
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
        allversions = [
            obj.find("Key").text.split("/")[1].replace(".csv", "")
            for obj in soup.find_all("Contents")
            if obj.find("Key").text.split("/")[0] == "versions"
        ]
        versions = list()
        for version in allversions:
            if re.match(f"(\\d+).(\\d+)", version):
                versions.append(Version(version))
        self.versions = versions

    def create_new_version(self):
        if self.version is None:
            raise ValueError(
                "No version provided, set version first with method 'set_version'"
            )

        if not "access_key" in self.auth_options.keys():
            raise ValueError("Read-only mode, cannot create new version")

        # update
        self._get_versions()
        if self.version in self.versions:
            raise ValueError("Version already exists")

        # get all objects
        objdic = get_s3_object_versions_and_tags(self.client, self.benchmark)

        # get objects to tag
        object_names_to_tag, versionid_of_objects_to_tag = get_s3_objects_to_tag(objdic)

        # Tag all objects with current version
        tags = minio.datatypes.Tags.new_object_tags()
        tags[str(self.version)] = "1"
        for n, v in zip(object_names_to_tag, versionid_of_objects_to_tag):
            self.client.set_object_tags(self.benchmark, n, tags, version_id=v)

        # set retention policy of objects
        retention_config = minio.retention.Retention(
            minio.commonconfig.GOVERNANCE,
            datetime.datetime.utcnow() + datetime.timedelta(weeks=1000),
        )
        for n, v in zip(object_names_to_tag, versionid_of_objects_to_tag):
            self.client.set_object_retention(
                self.benchmark, n, config=retention_config, version_id=v
            )

        # result = client.put_object(benchmark, "out/test.txt", io.BytesIO(b"asdf"), 4)
        # client.remove_object(benchmark, "out/test.txt")

        # write versions/{{version}}.csv
        vv_ls = get_s3version_from_bmversion(objdic, str(self.version))
        vv_str = prepare_csv_s3version_from_bmversion(vv_ls)
        version_filename = f"versions/{self.version}.csv"
        _ = self.client.put_object(
            self.benchmark, version_filename, io.BytesIO(vv_str.encode()), len(vv_str)
        )

        # tag and set retention policy of version file
        self.client.set_object_tags(self.benchmark, version_filename, tags)
        self.client.set_object_retention(
            self.benchmark, version_filename, config=retention_config
        )

        # check if version exists
        self._get_versions()
        if not self.version in self.versions:
            raise ValueError("Version creation failed")

    def _get_objects(self, readonly=False):
        if self.version is None:
            raise ValueError("No version provided")
        if not "secret_key" in self.auth_options.keys() or readonly:
            url = urlparse(f"{self.auth_options['endpoint']}/{self.benchmark}")
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
                    url.geturl(), self.benchmark, obj.find("Key").text
                )
                files[obj.find("Key").text] = {
                    "hash": obj.find("ETag").text.replace('"', ""),
                    "size": int(obj.find("Size").text),
                    "last_modified": obj.find("LastModified").text,
                    "x-object-meta-mtime": mtime,
                    "accesstime": accesstime,
                }
            self.files = files

        else:
            files = dict()
            # get all objects

            for element in self.client.list_objects(
                self.benchmark,
                recursive=True,
            ):
                if type(element.last_modified) is str:
                    files[element.object_name] = {
                        "size": element.size,
                        "last_modified": element.last_modified,
                        "hash": element.etag.replace('"', ""),
                    }
                elif type(element.last_modified) is datetime.datetime:
                    files[element.object_name] = {
                        "size": element.size,
                        "last_modified": element.last_modified.strftime(
                            "%Y-%m-%dT%H:%M:%S.%f"
                        ),
                        "hash": element.etag.replace('"', ""),
                    }
                else:
                    ValueError("Invalid last_modified")

            self.files = files

    def archive_version(self, version):
        NotImplementedError
        # self._update_overview(cleanup=True)

    def delete_version(self, version):
        NotImplementedError
        # self._update_overview(cleanup=True)


RemoteStorage.register(MinIOStorage)
