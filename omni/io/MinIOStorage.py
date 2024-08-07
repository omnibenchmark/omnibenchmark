"""MinIO class for remote storage."""

import datetime
import io
import json
import logging
import os
import re
from itertools import groupby
from typing import Dict, List, Union
from urllib.parse import urlparse

import boto3
import minio
import minio.commonconfig
import minio.retention
from packaging.version import Version

from omni.io.exception import (
    MinIOStorageBucketManipulationException,
    MinIOStorageConnectionException,
    MinIOStorageVersioningCorruptionException,
    RemoteStorageInvalidInputException,
)
from omni.io.RemoteStorage import RemoteStorage, is_valid_version
from omni.io.S3config import S3_DEFAULT_STORAGE_OPTIONS, bucket_readonly_policy

logging.basicConfig(level=logging.ERROR)
logging.getLogger("requests").setLevel(logging.DEBUG)
logging.getLogger("minio").setLevel(logging.DEBUG)
logger = logging.getLogger(__name__)


def get_s3_object_versions_and_tags(
    client: minio.Minio, benchmark: str, readonly: bool = False
) -> Dict:
    """
    Retrieve the metadata of all objects in a S3 Bucket.

    Args:
    - client: A MinIO client object.
    - benchmark: The name of the S3 bucket.
    - readonly: A boolean indicating whether the client is read-only. Object tags can only be retrieved if the client is not read-only.

    Returns:
    - A dictionary of objects and associated metadata.
    """
    if not client.bucket_exists(benchmark):
        raise MinIOStorageBucketManipulationException(
            f"Benchmark {benchmark} does not exist."
        )

    object_ls = list(
        client.list_objects(
            benchmark, include_version=True, include_user_meta=True, recursive=True
        )
    )

    # Group listed objects by their name
    grouped_objects = groupby(object_ls, key=lambda o: o.object_name)

    di = {}

    # Iterate over object groups
    for object_name, objects in grouped_objects:
        objects = list(objects)
        di[object_name] = {}

        # Extract attributes for the objects group
        version_ids = [o.version_id for o in objects]
        sizes = [o.size for o in objects]
        last_modified = [o.last_modified for o in objects]
        is_delete_markers = [o.is_delete_marker for o in objects]
        etags = [o.etag for o in objects]

        # Iterate to all versions and create a dictionary entry
        for i, v in enumerate(version_ids):
            tags_filt = {}
            if not readonly and not is_delete_markers[i]:
                tags = client.get_object_tags(benchmark, object_name, version_id=v)
                if tags is not None:
                    tags_filt = {k: w for k, w in tags.items() if is_valid_version(k)}

            di[object_name][v] = {
                "tags": tags_filt,
                "size": sizes[i],
                "last_modified": last_modified[i],
                "is_delete_marker": is_delete_markers[i],
                "etag": etags[i],
            }
    return di


def get_s3_objects_to_tag(
    objdic: Dict, tracked_dirs: List = S3_DEFAULT_STORAGE_OPTIONS["tracked_directories"]
) -> List:
    """
    Get a list of objects that need to be tagged with the current version.

    Args:
    - objdic: A dictionary of objects and associated metadata.
    - tracked_dirs: A list of directories to track.

    Returns:
    - A list of object names and a list of version ids.
    """
    object_names = sorted(list(objdic.keys()))
    object_names_to_tag = []
    versionid_of_objects_to_tag = []
    for object_name in object_names:
        # get newest version
        newest_version = sorted(
            objdic[object_name].items(),
            key=lambda it: it[1]["last_modified"],
            reverse=True,
        )[0][0]
        # get root directory
        object_name_red = object_name
        while not os.path.split(object_name_red)[0] == "":
            object_name_red = os.path.split(object_name_red)[0]
            if object_name_red in tracked_dirs:
                break
        # check if newest version exists and
        # if object is in tracked directories
        if (object_name_red in tracked_dirs) and (
            not objdic[object_name][newest_version]["is_delete_marker"]
        ):
            object_names_to_tag.append(object_name)
            versionid_of_objects_to_tag.append(newest_version)
    return object_names_to_tag, versionid_of_objects_to_tag


def get_single_s3version_from_bmversion(
    di: Dict, object_name: str, query_version: str
) -> Union[str, None]:
    """
    Get the object version where the tag matches the benchmark version.

    Args:
    - di: A dictionary of objects and associated metadata.
    - object_name: The name of the object.
    - query_version: The benchmark version to query.

    Returns:
    - The object version or None where the tag matches the benchmark version.
    """
    version_ls = list(
        filter(
            lambda v: query_version in di[object_name][v]["tags"].keys(),
            di[object_name].keys(),
        )
    )
    if len(version_ls) > 1:
        raise MinIOStorageVersioningCorruptionException(
            f"Multiple versions found for object {object_name}"
        )
    return version_ls[0] if version_ls else None


def get_s3version_from_bmversion(di: Dict, query_version: str) -> List:
    """
    Get the versions of all objects where the tag matches the benchmark version.

    Args:
    - di: A dictionary of objects and associated metadata.
    - query_version: The benchmark version to query.

    Returns:
    - A list of objects and associated versions where the tag matches the benchmark
    """
    summary_ls = []
    for object_name in di.keys():
        version_id = get_single_s3version_from_bmversion(di, object_name, query_version)
        if version_id is not None:
            summary_ls.append(
                [
                    object_name,
                    version_id,
                    di[object_name][version_id]["last_modified"],
                    di[object_name][version_id]["size"],
                    di[object_name][version_id]["etag"],
                ]
            )
    return summary_ls


def prepare_csv_s3version_from_bmversion(summary_ls: List) -> str:
    """
    Prepare a CSV string from a list of objects and associated versions.
    """
    outstr = "name,version_id,last_modified,size,etag\n"
    for element in summary_ls:
        outstr += f"{element[0]},{element[1]},{element[2]},{element[3]},{element[4]}\n"
    return outstr


def set_bucket_public_readonly(client: minio.Minio, bucket_name: str):
    policy = bucket_readonly_policy(bucket_name)
    client.set_bucket_policy(bucket_name, json.dumps(policy))


def set_bucket_lifecycle_config(
    client: minio.Minio, bucket_name: str, noncurrent_days: int = 1
):
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
    def __init__(
        self,
        auth_options: Dict,
        benchmark: str,
        storage_options: Dict = S3_DEFAULT_STORAGE_OPTIONS,
    ):
        super().__init__(auth_options, benchmark, storage_options)
        assert "endpoint" in self.auth_options.keys()
        if "access_key" in self.auth_options.keys():
            self.client = self.connect()
            self.roclient = self.connect(readonly=True)
            self._test_connect()

            if not self.client.bucket_exists(benchmark):
                logger.warning(
                    f"Benchmark {benchmark} does not exist, creating new benchmark."
                )
                self._create_benchmark(benchmark)

            self._get_versions()
        else:
            self.roclient = self.connect(readonly=True)
            self._get_versions()

    def connect(self, readonly=False) -> minio.Minio:
        """
        Connects to the MinIO storage.

        Returns:
        - A MinIO client object.
        """
        if readonly:
            assert "endpoint" in self.auth_options.keys()
            tmp_auth_options = self.auth_options.copy()
            if "access_key" in tmp_auth_options.keys():
                del tmp_auth_options["access_key"]
            if "secret_key" in tmp_auth_options.keys():
                del tmp_auth_options["secret_key"]
        else:
            assert (
                "endpoint" in self.auth_options.keys()
                and "access_key" in self.auth_options.keys()
                and "secret_key" in self.auth_options.keys()
            )
            tmp_auth_options = self.auth_options.copy()
        try:
            return minio.Minio(**tmp_auth_options)
        except Exception as e:
            url = urlparse(tmp_auth_options["endpoint"])
            tmp_auth_options["endpoint"] = url.netloc
            return minio.Minio(**tmp_auth_options)

    def _test_connect(self) -> None:
        try:
            _ = self.client.list_objects(self.benchmark)
        except MinIOStorageConnectionException as e:
            raise e

    def _create_benchmark(self, benchmark) -> None:
        if self.client.bucket_exists(benchmark):
            raise MinIOStorageBucketManipulationException(
                f"Benchmark {benchmark} already exists."
            )
        # create new version
        self.client.make_bucket(bucket_name=benchmark, object_lock=True)
        if self.client._base_url.is_aws_host:
            s3 = boto3.client(
                "s3",
                aws_access_key_id=self.auth_options["access_key"],
                aws_secret_access_key=self.auth_options["secret_key"],
            )
            s3.delete_public_access_block(Bucket=benchmark)
        set_bucket_public_readonly(self.client, benchmark)
        set_bucket_lifecycle_config(self.client, benchmark, noncurrent_days=1)
        # _ = self.client.put_object(
        #     self.benchmark, "versions/.ignore.csv", io.BytesIO(b""), 0
        # )
        if not self.client.bucket_exists(benchmark):
            raise MinIOStorageBucketManipulationException(
                f"Bucket creation for benchmark {benchmark} failed"
            )

    def _get_versions(self) -> None:
        versionobjects = list(
            self.roclient.list_objects(
                self.benchmark, prefix="versions", recursive=True
            )
        )
        allversions = [
            os.path.basename(v.object_name).replace(".csv", "") for v in versionobjects
        ]
        versions = list()
        for version in allversions:
            if re.match(f"(\\d+).(\\d+)", version):
                versions.append(Version(version))
        self.versions = versions

    def create_new_version(self) -> None:
        if self.version is None:
            raise RemoteStorageInvalidInputException(
                "No version provided, set version first with method 'set_version'"
            )

        if not "access_key" in self.auth_options.keys():
            raise RemoteStorageInvalidInputException(
                "Read-only mode, cannot create new version, set access_key and secret_key in auth_options"
            )

        # update
        self._get_versions()
        if self.version in self.versions:
            raise RemoteStorageInvalidInputException(
                "Version already exists, set new version with method 'set_version'"
            )

        # get all objects
        objdic = get_s3_object_versions_and_tags(self.client, self.benchmark)

        # get objects to tag
        object_names_to_tag, versionid_of_objects_to_tag = get_s3_objects_to_tag(
            objdic, tracked_dirs=self.storage_options["tracked_directories"]
        )

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

        # refresh object list
        objdic = get_s3_object_versions_and_tags(self.client, self.benchmark)

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
            raise MinIOStorageBucketManipulationException("Version creation failed")

    def _get_objects(self) -> None:
        if self.version is None:
            raise RemoteStorageInvalidInputException(
                "No version provided, set version first with method 'set_version'"
            )

        if self.version in self.versions:
            # read overview file
            response = self.roclient.get_object(
                self.benchmark, f"versions/{self.version}.csv"
            )
            objls = response.data.decode("utf-8")
            objls = objls.split("\n")
            objls = [obj for obj in objls if obj]
            objdict = {}
            header = objls[0].split(",")
            assert header[0] == "name"
            for obj in objls[1:]:
                tmpsplit = obj.split(",")
                objdict[tmpsplit[0]] = {}
                for i, head in enumerate(header[1:]):
                    objdict[tmpsplit[0]][head] = tmpsplit[i + 1]

            # add overview file to files
            response_headers = response.getheaders()
            objdict[f"versions/{self.version}.csv"] = {
                "version_id": response_headers.get("x-amz-version-id"),
                # some parsing of date to get to consistent format
                "last_modified": datetime.datetime.strptime(
                    response_headers.get("last-modified"), "%a, %d %b %Y %H:%M:%S GMT"
                ).strftime("%Y-%m-%d %H:%M:%S.%f+00:00"),
                "size": response_headers.get("content-length"),
                "etag": response_headers.get("etag").replace('"', ""),
            }
        else:
            # get all objects
            objdic = get_s3_object_versions_and_tags(
                self.roclient, self.benchmark, readonly=True
            )

            # get objects to tag
            object_names_to_tag, versionid_of_objects_to_tag = get_s3_objects_to_tag(
                objdic, tracked_dirs=self.storage_options["tracked_directories"]
            )
            objdict = {}
            for obj, vt in zip(object_names_to_tag, versionid_of_objects_to_tag):
                objdict[obj] = {}
                objdict[obj]["version_id"] = vt
                objdict[obj]["etag"] = objdic[obj][vt]["etag"]
                objdict[obj]["size"] = objdic[obj][vt]["size"]
                objdict[obj]["last_modified"] = objdic[obj][vt]["last_modified"]

        self.files = objdict

    def download_object(self, object_name: str, local_path: str) -> None:
        if self.version is None:
            raise RemoteStorageInvalidInputException(
                "No version provided. Set version first with method 'set_version'"
            )
        if self.files is None:
            self._get_objects()
        if self.files is None:
            raise RemoteStorageInvalidInputException("No objects found")
        if object_name not in self.files.keys():
            raise RemoteStorageInvalidInputException(f"Object {object_name} not found")
        _ = self.roclient.fget_object(
            self.benchmark,
            object_name,
            local_path,
            version_id=self.files[object_name]["version_id"],
        )

    def archive_version(self, version):
        NotImplementedError
        # self._update_overview(cleanup=True)

    def delete_version(self, version):
        NotImplementedError
        # self._update_overview(cleanup=True)


RemoteStorage.register(MinIOStorage)
