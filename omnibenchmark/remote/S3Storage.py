"""S3-compatible storage implementation."""

import datetime
import json
import logging
import os
import re

from botocore import UNSIGNED
from botocore.config import Config
from botocore.exceptions import ClientError
from packaging.version import Version
from pathlib import Path
from typing import Dict, Optional

try:
    import boto3
except ImportError:
    boto3 = None  # type: ignore

from omnibenchmark.benchmark import BenchmarkExecution
from omnibenchmark.remote.exception import (
    RemoteStorageInvalidInputException,
    S3StorageConnectionException,
    S3StorageBucketManipulationException,
)
from omnibenchmark.archive.components import prepare_archive_software
from omnibenchmark.remote.RemoteStorage import RemoteStorage, StorageOptions
from omnibenchmark.remote.S3config import bucket_readonly_policy
from omnibenchmark.remote.S3versioning import get_s3_object_versions_and_tags
from omnibenchmark.remote.versioning import (
    filter_objects_to_tag,
    get_objects_to_tag,
    get_remoteversion_from_bmversion,
    prepare_csv_remoteversion_from_bmversion,
)
from omnibenchmark.versioning.git import GitAwareBenchmarkVersionManager

logging.getLogger("requests").setLevel(logging.DEBUG)
logger = logging.getLogger(__name__)


class S3CompatibleStorage(RemoteStorage):
    """
    S3-compatible storage implementation with versioning and object protection.

    Provides a concrete implementation of RemoteStorage for MinIO and AWS S3,
    with automatic versioning, object locking, and retention policies for
    benchmark result reproducibility.

    Uses boto3 exclusively for all S3 operations, which works with both
    AWS S3 and S3-compatible services like MinIO (same wire protocol).

    Args:
        auth_options (Dict): Authentication configuration with 'endpoint',
            and optionally 'access_key', 'secret_key', 'secure'.
        benchmark (str): Benchmark name/identifier used as S3 bucket name.
        storage_options (StorageOptions): Configuration for tracked directories.

    Raises:
        S3StorageConnectionException: When connection to S3 fails.
        S3StorageBucketManipulationException: When bucket operations fail.
        RemoteStorageInvalidInputException: When version operations fail.
    """

    def __init__(
        self,
        auth_options: Dict,
        benchmark: str,
        storage_options: StorageOptions,
    ):
        super().__init__(auth_options, benchmark, storage_options)
        assert "endpoint" in self.auth_options.keys()

        if "access_key" in self.auth_options.keys():
            self.client = self.connect()
            self.roclient = self.connect(readonly=True)

            if not self._bucket_exists(benchmark):
                logger.warning(
                    f"Benchmark {benchmark} does not exist, creating new benchmark."
                )
                self._create_benchmark(benchmark)

            self._test_connect()

            self._get_versions()
        else:
            self.roclient = self.connect(readonly=True)
            self._get_versions()

    # ------------------------------------------------------------------
    # Connection helpers
    # ------------------------------------------------------------------

    def _endpoint_url(self) -> str:
        """Return the endpoint as a full URL (with scheme)."""
        endpoint = self.auth_options.get("endpoint", "")
        if endpoint.startswith(("http://", "https://")):
            return endpoint
        secure = self.auth_options.get("secure", True)
        scheme = "https" if secure else "http"
        return f"{scheme}://{endpoint}"

    @property
    def _is_aws(self) -> bool:
        """True when the endpoint is AWS S3 (no custom endpoint_url needed)."""
        return "amazonaws.com" in self.auth_options.get("endpoint", "")

    def connect(self, readonly: bool = False):
        """
        Connects to the S3-compatible storage using boto3.

        Args:
            readonly: When True, connect anonymously (no credentials).

        Returns:
            A boto3 S3 client.
        """
        assert "endpoint" in self.auth_options

        kwargs: dict = {}
        if not self._is_aws:
            kwargs["endpoint_url"] = self._endpoint_url()
            kwargs["region_name"] = "us-east-1"

        if readonly or "access_key" not in self.auth_options:
            kwargs["config"] = Config(signature_version=UNSIGNED)
        else:
            assert (
                "access_key" in self.auth_options and "secret_key" in self.auth_options
            )
            kwargs["aws_access_key_id"] = self.auth_options["access_key"]
            kwargs["aws_secret_access_key"] = self.auth_options["secret_key"]

        return boto3.client("s3", **kwargs)  # type: ignore[union-attr]

    def _bucket_exists(self, bucket_name: str) -> bool:
        """Return True if the bucket exists and is accessible."""
        try:
            self.client.head_bucket(Bucket=bucket_name)
            return True
        except ClientError:
            return False

    # ------------------------------------------------------------------
    # Abstract method implementations
    # ------------------------------------------------------------------

    def _test_connect(self) -> None:
        try:
            self.client.list_objects_v2(Bucket=self.benchmark, MaxKeys=1)
        except ClientError as e:
            raise S3StorageConnectionException(str(e)) from e

    def _create_benchmark(self, benchmark: str) -> None:
        if self._bucket_exists(benchmark):
            raise S3StorageBucketManipulationException(
                f"Benchmark {benchmark} already exists."
            )

        self.client.create_bucket(
            Bucket=benchmark,
            ObjectLockEnabledForBucket=True,
        )

        if self._is_aws:
            self.client.delete_public_access_block(Bucket=benchmark)

        policy = bucket_readonly_policy(benchmark)
        self.client.put_bucket_policy(Bucket=benchmark, Policy=json.dumps(policy))

        # Lifecycle configuration is best-effort: MinIO requires Content-MD5 for
        # this operation, which recent boto3 no longer sends by default (it uses
        # x-amz-checksum-* instead). The policy is non-essential — it only
        # auto-expires old noncurrent object versions to save storage.
        try:
            self.client.put_bucket_lifecycle_configuration(
                Bucket=benchmark,
                LifecycleConfiguration={
                    "Rules": [
                        {
                            "ID": "rule1",
                            "Filter": {"Prefix": ""},
                            "Status": "Enabled",
                            "NoncurrentVersionExpiration": {"NoncurrentDays": 1},
                        }
                    ]
                },
            )
        except ClientError as e:
            code = e.response["Error"]["Code"]
            if code in ("MissingContentMD5", "NotImplemented"):
                logger.debug("Lifecycle configuration not applied (%s): %s", code, e)
            else:
                raise

        if not self._bucket_exists(benchmark):
            raise S3StorageBucketManipulationException(
                f"Bucket creation for benchmark {benchmark} failed"
            )

    def _get_versions(self) -> None:
        try:
            paginator = self.roclient.get_paginator("list_objects_v2")
            pages = paginator.paginate(Bucket=self.benchmark, Prefix="versions/")
            versionobjects = []
            for page in pages:
                versionobjects.extend(page.get("Contents", []))
        except ClientError as e:
            code = e.response["Error"]["Code"]
            if code in ("AccessDenied", "403"):
                logger.debug("S3 access denied", exc_info=True)
                raise S3StorageConnectionException(
                    f"Access denied to S3 bucket '{self.benchmark}'."
                    " Check your credentials have the correct permissions."
                ) from e
            if code in ("NoSuchBucket", "404"):
                logger.debug("S3 bucket not found", exc_info=True)
                raise S3StorageConnectionException(
                    f"S3 bucket '{self.benchmark}' does not exist."
                ) from e
            raise

        allversions = [
            os.path.basename(v["Key"]).replace(".csv", "") for v in versionobjects
        ]
        versions = []
        for version in allversions:
            if re.match(r"(\d+)\.(\d+)", version):
                versions.append(Version(version))
        self.versions = versions

    def create_new_version(
        self, benchmark: Optional[BenchmarkExecution] = None
    ) -> None:
        self._validate_for_write()
        self._get_versions()

        version_manager = self._build_version_manager(benchmark)
        try:
            self._persist_version(version_manager, benchmark)
            if benchmark is not None:
                self._upload_benchmark_config(benchmark)
            self._tag_and_protect_objects(benchmark)
        except Exception as e:
            logger.error(f"Error creating version {self.version}, rolling back: {e}")
            raise RemoteStorageInvalidInputException(f"Failed to create version: {e}")

        self._write_version_manifest()
        self._get_versions()
        if self.version not in self.versions:
            raise S3StorageBucketManipulationException("Version creation failed")

    def load_objects(self) -> None:
        if self.version is None:
            raise RemoteStorageInvalidInputException(
                "No version provided, set version first with method 'set_version'"
            )

        if self.version in self.versions:
            response = self.roclient.get_object(
                Bucket=self.benchmark, Key=f"versions/{self.version}.csv"
            )
            body = response["Body"].read().decode("utf-8")
            objls = [line for line in body.split("\n") if line]
            objdict: dict = {}
            header = objls[0].split(",")
            assert header[0] == "name"
            for obj in objls[1:]:
                tmpsplit = obj.split(",")
                objdict[tmpsplit[0]] = {}
                for i, head in enumerate(header[1:]):
                    objdict[tmpsplit[0]][head] = tmpsplit[i + 1]

            objdict[f"versions/{self.version}.csv"] = {
                "version_id": response.get("VersionId"),
                "last_modified": response["LastModified"],
                "size": str(response["ContentLength"]),
                "etag": response["ETag"].strip('"'),
            }
        else:
            objdic = get_s3_object_versions_and_tags(
                self.roclient, self.benchmark, readonly=True
            )

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
            self.load_objects()
        if self.files is None:
            raise RemoteStorageInvalidInputException("No objects found")
        if object_name not in self.files.keys():
            raise RemoteStorageInvalidInputException(f"Object {object_name} not found")

        version_id = self.files[object_name].get("version_id")
        extra_args = {"VersionId": version_id} if version_id else {}

        Path(local_path).parent.mkdir(parents=True, exist_ok=True)
        self.roclient.download_file(
            Bucket=self.benchmark,
            Key=object_name,
            Filename=local_path,
            ExtraArgs=extra_args or None,
        )

    def archive_version(
        self,
        benchmark: BenchmarkExecution,
        outdir: Path = Path(),
        config: bool = True,
        code: bool = False,
        software: bool = False,
        results: bool = False,
    ):
        from omnibenchmark.remote.archive import archive_version

        # TODO: outdir is in benchmarkExecution
        # TODO: upload the zip archive
        _ = archive_version(benchmark, outdir, config, code, software, results)

    def upload_version(self, *args, **kwargs):
        raise NotImplementedError("upload_version is not yet implemented")

    def delete_version(self, version: str) -> None:
        """Delete a benchmark version and its WORM-protected S3 object versions.

        Uses BypassGovernanceRetention=True to remove objects protected by
        Governance Mode retention. Requires the s3:BypassGovernanceRetention
        permission on the bucket.

        Note: if the same S3 object version is referenced by multiple benchmark
        versions (i.e. the object was not re-uploaded between versions), deleting
        one benchmark version will also remove that object version from other
        benchmark versions that share it.
        """
        if "access_key" not in self.auth_options:
            raise RemoteStorageInvalidInputException(
                "Read-only mode, cannot delete version,"
                " set access_key and secret_key in auth_options"
            )

        self._get_versions()
        target = Version(version)
        if target not in self.versions:
            raise RemoteStorageInvalidInputException(
                f"Version {version} does not exist"
            )

        self.set_version(version)
        self.load_objects()
        objects_to_delete = dict(self.files)

        failed: list = []
        for obj_name, meta in objects_to_delete.items():
            version_id = meta.get("version_id")
            try:
                delete_kwargs: dict = {
                    "Bucket": self.benchmark,
                    "Key": obj_name,
                    "BypassGovernanceRetention": True,
                }
                if version_id:
                    delete_kwargs["VersionId"] = version_id
                self.client.delete_object(**delete_kwargs)
            except Exception as e:
                logger.warning(
                    f"Could not delete {obj_name} (version {version_id}): {e}"
                )
                failed.append(obj_name)

        if failed:
            raise RemoteStorageInvalidInputException(
                f"delete_version {version} completed with errors;"
                f" {len(failed)} object(s) could not be deleted: {failed}"
            )

        self._get_versions()

    # ------------------------------------------------------------------
    # create_new_version helpers
    # ------------------------------------------------------------------

    def _validate_for_write(self) -> None:
        if self.version is None:
            raise RemoteStorageInvalidInputException(
                "No version provided, set version first with method 'set_version'"
            )
        if "access_key" not in self.auth_options:
            raise RemoteStorageInvalidInputException(
                "Read-only mode, cannot create new version,"
                " set access_key and secret_key in auth_options"
            )

    def _build_version_manager(self, benchmark: Optional[BenchmarkExecution]):
        """Return a version manager scoped to this call (not stored on self)."""
        from omnibenchmark.config import get_temp_dir

        temp_dir = get_temp_dir()
        known = [str(v) for v in self.versions]

        if benchmark is not None:
            try:
                vm = GitAwareBenchmarkVersionManager(
                    benchmark_path=benchmark.get_definition_file(),
                    git_repo_path=benchmark.context.directory,
                    lock_dir=temp_dir / "locks",
                    lock_timeout=30.0,
                )
                vm.set_known_versions(known)
                return vm
            except Exception:
                pass  # fall through to basic version manager

            from omnibenchmark.versioning import BenchmarkVersionManager

            vm = BenchmarkVersionManager(
                benchmark_path=benchmark.get_definition_file(),
                lock_dir=temp_dir / "locks",
                lock_timeout=30.0,
            )
        else:
            from omnibenchmark.versioning import BenchmarkVersionManager

            vm = BenchmarkVersionManager(
                benchmark_path=temp_dir / f"{self.benchmark}.yaml",
                lock_dir=temp_dir / "locks",
                lock_timeout=30.0,
            )

        vm.set_known_versions(known)
        return vm

    def _persist_version(
        self, version_manager, benchmark: Optional[BenchmarkExecution]
    ) -> None:
        """Record the new version via the version manager (git-aware or basic)."""
        if benchmark is not None and hasattr(
            version_manager, "create_version_with_persistence"
        ):
            try:
                version_manager.create_version_with_persistence(  # type: ignore
                    benchmark,
                    str(self.version),
                    f"Create benchmark version {self.version}",
                )
                return
            except Exception:
                # git-based persistence failed; fall back to in-memory tracking
                from omnibenchmark.config import get_temp_dir
                from omnibenchmark.versioning import BenchmarkVersionManager

                temp_dir = get_temp_dir()
                version_manager = BenchmarkVersionManager(
                    benchmark_path=benchmark.get_definition_file(),
                    lock_dir=temp_dir / "locks",
                    lock_timeout=30.0,
                )
                version_manager.set_known_versions([str(v) for v in self.versions])

        version_manager.set_known_versions(
            version_manager.get_versions() + [str(self.version)]
        )

    def _upload_benchmark_config(self, benchmark: BenchmarkExecution) -> None:
        """Upload the benchmark YAML and its software environment files."""
        bmfile = Path(benchmark.get_definition_file())
        with open(bmfile, "r") as fh:
            bmstr = fh.read()
        self.client.put_object(
            Bucket=self.benchmark,
            Key="config/benchmark.yaml",
            Body=bmstr.encode(),
        )
        for software_file in prepare_archive_software(benchmark):
            with open(software_file, "r") as fh:
                softstr = fh.read()
            self.client.put_object(
                Bucket=self.benchmark,
                Key=f"software/{software_file.name}",
                Body=softstr.encode(),
            )

    def _tag_and_protect_objects(self, benchmark: Optional[BenchmarkExecution]) -> None:
        """Tag current objects with the benchmark version and apply WORM retention.

        Rolls back any tags it applied on failure before re-raising.
        """
        objdic = get_s3_object_versions_and_tags(self.client, self.benchmark)
        names, vids = get_objects_to_tag(objdic, storage_options=self.storage_options)
        names, vids = filter_objects_to_tag(
            names, vids, self.storage_options, benchmark
        )

        tag_set = [{"Key": str(self.version), "Value": "1"}]
        tagged: list = []
        try:
            for n, v in zip(names, vids):
                self.client.put_object_tagging(
                    Bucket=self.benchmark,
                    Key=n,
                    VersionId=v,
                    Tagging={"TagSet": tag_set},
                )
                tagged.append((n, v))
            retain_until = datetime.datetime.now(datetime.UTC) + datetime.timedelta(
                weeks=1000
            )
            for n, v in zip(names, vids):
                self.client.put_object_retention(
                    Bucket=self.benchmark,
                    Key=n,
                    VersionId=v,
                    Retention={"Mode": "GOVERNANCE", "RetainUntilDate": retain_until},
                    ChecksumAlgorithm="CRC32",
                )
        except Exception:
            self._rollback_tags(tagged)
            raise

    def _rollback_tags(self, tagged_objects: list) -> None:
        """Best-effort removal of version tags; called on create_new_version failure."""
        for obj_name, version_id in tagged_objects:
            try:
                self.client.delete_object_tagging(
                    Bucket=self.benchmark,
                    Key=obj_name,
                    VersionId=version_id,
                )
            except Exception:
                pass

    def _write_version_manifest(self) -> None:
        """Write versions/{version}.csv and protect it with a tag and WORM retention."""
        objdic = get_s3_object_versions_and_tags(self.client, self.benchmark)
        vv_ls = get_remoteversion_from_bmversion(objdic, str(self.version))
        vv_str = prepare_csv_remoteversion_from_bmversion(vv_ls)
        version_filename = f"versions/{self.version}.csv"
        self.client.put_object(
            Bucket=self.benchmark,
            Key=version_filename,
            Body=vv_str.encode(),
        )
        tag_set = [{"Key": str(self.version), "Value": "1"}]
        retain_until = datetime.datetime.now(datetime.UTC) + datetime.timedelta(
            weeks=1000
        )
        self.client.put_object_tagging(
            Bucket=self.benchmark,
            Key=version_filename,
            Tagging={"TagSet": tag_set},
        )
        self.client.put_object_retention(
            Bucket=self.benchmark,
            Key=version_filename,
            Retention={"Mode": "GOVERNANCE", "RetainUntilDate": retain_until},
            ChecksumAlgorithm="CRC32",
        )


RemoteStorage.register(S3CompatibleStorage)
