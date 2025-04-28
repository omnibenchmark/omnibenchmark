"""MinIO class for remote storage."""

import datetime
import io
import json
import logging
import os
import re

from packaging.version import Version
from pathlib import Path
from typing import Dict, Union
from urllib.parse import urlparse

import boto3
import minio
import minio.commonconfig
import minio.retention

from omnibenchmark.benchmark import Benchmark
from omnibenchmark.io.exception import (
    MinIOStorageBucketManipulationException,
    MinIOStorageConnectionException,
    RemoteStorageInvalidInputException,
)
from omnibenchmark.io.archive import prepare_archive_software
from omnibenchmark.io.RemoteStorage import RemoteStorage
from omnibenchmark.io.S3config import S3_DEFAULT_STORAGE_OPTIONS, bucket_readonly_policy
from omnibenchmark.io.S3versioning import get_s3_object_versions_and_tags
from omnibenchmark.io.versioning import (
    filter_objects_to_tag,
    get_objects_to_tag,
    get_remoteversion_from_bmversion,
    prepare_csv_remoteversion_from_bmversion,
)

# TODO: set if we have the DEBUG flag set
logging.getLogger("requests").setLevel(logging.DEBUG)
logging.getLogger("minio").setLevel(logging.DEBUG)
logger = logging.getLogger(__name__)


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
        # TODO: fix type problems here.
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
        except Exception:
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
            if re.match("(\\d+).(\\d+)", version):
                versions.append(Version(version))
        self.versions = versions

    def create_new_version(self, benchmark: Union[None, Benchmark] = None) -> None:
        if self.version is None:
            raise RemoteStorageInvalidInputException(
                "No version provided, set version first with method 'set_version'"
            )

        if "access_key" not in self.auth_options.keys():
            raise RemoteStorageInvalidInputException(
                "Read-only mode, cannot create new version, set access_key and secret_key in auth_options"
            )

        # update
        self._get_versions()
        if self.version in self.versions:
            raise RemoteStorageInvalidInputException(
                "Version already exists, set new version with method 'set_version'"
            )

        # upload benchmark definition file
        if benchmark is not None:
            # read file
            bmfile = Path(benchmark.get_definition_file())
            with open(bmfile, "r") as fh:
                bmstr = fh.read()

            # upload file
            _ = self.client.put_object(
                self.benchmark,
                "config/benchmark.yaml",
                io.BytesIO(bmstr.encode()),
                len(bmstr),
            )

            # read software files
            software_files = prepare_archive_software(benchmark)

            # upload software files
            for software_file in software_files:
                print(software_file)
                with open(software_file, "r") as fh:
                    softstr = fh.read()

                _ = self.client.put_object(
                    self.benchmark,
                    f"software/{software_file.name}",
                    io.BytesIO(softstr.encode()),
                    len(softstr),
                )

        # get all objects
        objdic = get_s3_object_versions_and_tags(self.client, self.benchmark)

        # get objects to tag
        object_names_to_tag, versionid_of_objects_to_tag = get_objects_to_tag(
            objdic, storage_options=self.storage_options
        )

        # filter objects based on workflow
        object_names_to_tag, versionid_of_objects_to_tag = filter_objects_to_tag(
            object_names_to_tag,
            versionid_of_objects_to_tag,
            benchmark,
            self.storage_options,
        )

        # Tag all objects with current version
        tags = minio.datatypes.Tags.new_object_tags()
        tags[str(self.version)] = "1"
        for n, v in zip(object_names_to_tag, versionid_of_objects_to_tag):
            self.client.set_object_tags(self.benchmark, n, tags, version_id=v)

        # set retention policy of objects
        retention_config = minio.retention.Retention(
            minio.commonconfig.GOVERNANCE,
            datetime.datetime.now(datetime.UTC) + datetime.timedelta(weeks=1000),
        )
        for n, v in zip(object_names_to_tag, versionid_of_objects_to_tag):
            self.client.set_object_retention(
                self.benchmark, n, config=retention_config, version_id=v
            )

        # refresh object list
        objdic = get_s3_object_versions_and_tags(self.client, self.benchmark)

        # write versions/{{version}}.csv
        vv_ls = get_remoteversion_from_bmversion(objdic, str(self.version))
        vv_str = prepare_csv_remoteversion_from_bmversion(vv_ls)
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
        if self.version not in self.versions:
            raise MinIOStorageBucketManipulationException("Version creation failed")

    # TODO(ben): review if we want this public, it's used in tests
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
            response_headers = response.headers
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
            object_names_to_tag, versionid_of_objects_to_tag = get_objects_to_tag(
                objdic, storage_options=self.storage_options
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

    def archive_version(
        self,
        benchmark: Benchmark,
        outdir: Path = Path(),
        config: bool = True,
        code: bool = False,
        software: bool = False,
        results: bool = False,
    ):
        from omnibenchmark.io.archive import archive_version

        # TODO: upload the zip archive
        _ = archive_version(benchmark, outdir, config, code, software, results)

    def delete_version(self, version):
        raise NotImplementedError


RemoteStorage.register(MinIOStorage)
