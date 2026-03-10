"""S3-compatible storage implementation."""
# pyright: reportCallIssue=false
# TODO: Fix boto3/minio API call type errors
# TODO: Migrate from minio library to boto3-only
#       The minio Python library is not actively maintained as of late 2025.
#       Since MinIO servers are S3-compatible, we should migrate to use boto3
#       exclusively, which is actively maintained by AWS and is the industry
#       standard. All current minio library operations (bucket management,
#       versioning, policies, lifecycle) can be done with boto3.

import datetime
import io
import json
import logging
import os
import re

from packaging.version import Version
from pathlib import Path
from typing import Dict, Optional
from urllib.parse import urlparse

try:
    import boto3
    import minio
    import minio.commonconfig
    import minio.retention
    from minio.lifecycleconfig import LifecycleConfig, Rule, NoncurrentVersionExpiration
    from minio.commonconfig import Filter, Tags
except ImportError:
    # perhaps should be handled at a higher level, and advice the user to install the required packages
    # via pip install minio boto3 or extras['s3'].
    pass

from omnibenchmark.benchmark import BenchmarkExecution
from omnibenchmark.remote.exception import (
    RemoteStorageInvalidInputException,
    MinIOStorageConnectionException,
    MinIOStorageBucketManipulationException,
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

# TODO: set if we have the DEBUG flag set
logging.getLogger("requests").setLevel(logging.DEBUG)
logging.getLogger("minio").setLevel(logging.DEBUG)
logger = logging.getLogger(__name__)


def set_bucket_public_readonly(client: "minio.Minio", bucket_name: str):
    policy = bucket_readonly_policy(bucket_name)
    client.set_bucket_policy(bucket_name, json.dumps(policy))


def set_bucket_lifecycle_config(
    client: "minio.Minio", bucket_name: str, noncurrent_days: int = 1
):
    lifecycle_config = LifecycleConfig(
        [
            Rule(
                minio.commonconfig.ENABLED,
                rule_filter=Filter(prefix="*"),
                rule_id="rule1",
                noncurrent_version_expiration=NoncurrentVersionExpiration(
                    noncurrent_days=noncurrent_days
                ),
            ),
        ],
    )
    client.set_bucket_lifecycle(bucket_name, lifecycle_config)


class S3CompatibleStorage(RemoteStorage):
    """
    MinIO/S3-compatible storage implementation with versioning and object protection.

    This class provides concrete implementation of RemoteStorage for MinIO and AWS S3,
    with automatic versioning, object locking, and retention policies for benchmark
    result reproducibility.

    ## S3 Object Lock and Versioning Implementation

    MinIOStorage automatically enables S3 Object Lock on buckets during creation:
    - **Bucket Creation**: `object_lock=True` enables versioning and object lock
    - **Version Tagging**: Objects are tagged with benchmark versions during `create_new_version()`
    - **Governance Retention**: Tagged objects receive automatic retention policies
    - **WORM Protection**: Objects become immutable once versioned

    ## Version Creation Workflow

    When `create_new_version()` is called:
    1. **Object Discovery**: Scans tracked directories (`out/`, `config/`, `versions/`, `software/`)
    2. **Version Tagging**: Tags current object versions with benchmark version
    3. **Retention Policy**: Applies S3 Object Lock Governance Mode protection
    4. **Manifest Creation**: Stores version manifest at `versions/{version}.csv`
    5. **Git Integration**: Commits version changes to Git repository if available

    ## Object Protection Details

    Protected objects exhibit the following behaviors:
    - **Immutability**: Cannot be overwritten or deleted without bypass permissions
    - **Version Preservation**: All historical versions remain accessible
    - **Error Messages**: Deletion attempts return "Object is WORM protected" errors
    - **Bypass Required**: Only `BypassGovernanceRetention=True` can override protection

    ## Cleanup Considerations for Tests

    Test environments require special cleanup procedures:
    ```python
    # This will fail for versioned objects
    s3_client.delete_object(Bucket=bucket, Key=key)

    # Required approach for cleanup
    versions = s3_client.list_object_versions(Bucket=bucket)
    for version in versions.get('Versions', []):
        s3_client.delete_object(
            Bucket=bucket,
            Key=version['Key'],
            VersionId=version['VersionId'],
            BypassGovernanceRetention=True
        )
    ```

    ## MinIO vs AWS S3 Differences

    - **MinIO**: Object lock cannot be disabled once enabled (`only 'Enabled' value is allowed`)
    - **AWS S3**: Supports more granular object lock configuration changes
    - **Both**: Support governance retention bypass with appropriate permissions
    - **Both**: Maintain version history and support object tagging

    ## Security and Access Control

    Generated IAM policies automatically exclude dangerous permissions:
    - ❌ `s3:BypassGovernanceRetention` - Reserved for administrative operations
    - ❌ `s3:DeleteObjectVersion` - Prevents version history tampering
    - ❌ `s3:DeleteBucket` - Prevents accidental benchmark deletion
    - ✅ Standard read/write operations on current object versions

    Args:
        auth_options (Dict): Authentication configuration with 'endpoint', 'access_key', 'secret_key'
        benchmark (str): Benchmark name/identifier used as S3 bucket name
        storage_options (StorageOptions): Configuration for tracked directories and file patterns

    Raises:
        MinIOStorageConnectionException: When connection to MinIO/S3 fails
        MinIOStorageBucketManipulationException: When bucket operations fail
        RemoteStorageInvalidInputException: When version operations fail due to validation errors

    Example:
        >>> storage = MinIOStorage(
        ...     auth_options={'endpoint': 'https://s3.amazonaws.com', 'access_key': 'key', 'secret_key': 'secret'},
        ...     benchmark='my-benchmark',
        ...     storage_options=StorageOptions(out_dir='out')
        ... )
        >>> storage.set_version('1.0')
        >>> storage.create_new_version()  # Creates WORM-protected version 1.0
        >>> storage.versions  # ['1.0']
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

    def connect(self, readonly=False) -> "minio.Minio":
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
            return minio.Minio(
                endpoint=tmp_auth_options["endpoint"],
                **{k: v for k, v in tmp_auth_options.items() if k != "endpoint"},
            )
        except NameError:
            import click
            import sys

            logger.error(
                click.style("[ERROR]", fg="red", bold=True)
                + " S3/MinIO storage libraries not installed. Install with: pip install omnibenchmark[s3]"
            )
            if logger.isEnabledFor(logging.DEBUG):
                logger.exception("Full traceback:")
            sys.exit(1)
        except Exception:
            url = urlparse(tmp_auth_options["endpoint"])
            tmp_auth_options["endpoint"] = url.netloc
            return minio.Minio(
                endpoint=tmp_auth_options["endpoint"],
                **{k: v for k, v in tmp_auth_options.items() if k != "endpoint"},
            )

    def _test_connect(self) -> None:
        try:
            _ = self.client.list_objects(self.benchmark)
        except MinIOStorageConnectionException as e:
            raise e

    def _create_benchmark(self, benchmark: str) -> None:
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
        try:
            versionobjects = list(
                self.roclient.list_objects(
                    self.benchmark, prefix="versions", recursive=True
                )
            )
        except Exception as e:
            code = getattr(e, "code", None)
            if code == "AccessDenied":
                logger.debug("S3 access denied", exc_info=True)
                raise MinIOStorageConnectionException(
                    f"Access denied to S3 bucket '{self.benchmark}'."
                    " Check your credentials have the correct permissions."
                ) from e
            if code == "NoSuchBucket":
                logger.debug("S3 bucket not found", exc_info=True)
                raise MinIOStorageConnectionException(
                    f"S3 bucket '{self.benchmark}' does not exist."
                ) from e
            raise

        allversions = [
            os.path.basename(v.object_name).replace(".csv", "")
            for v in versionobjects
            if v.object_name is not None
        ]
        versions = list()
        for version in allversions:
            if re.match("(\\d+).(\\d+)", version):
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
            raise MinIOStorageBucketManipulationException("Version creation failed")

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
            self.benchmark,
            "config/benchmark.yaml",
            io.BytesIO(bmstr.encode()),
            len(bmstr),
        )
        for software_file in prepare_archive_software(benchmark):
            with open(software_file, "r") as fh:
                softstr = fh.read()
            self.client.put_object(
                self.benchmark,
                f"software/{software_file.name}",
                io.BytesIO(softstr.encode()),
                len(softstr),
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

        tags = Tags.new_object_tags()
        tags[str(self.version)] = "1"
        tagged: list = []
        try:
            for n, v in zip(names, vids):
                self.client.set_object_tags(self.benchmark, n, tags, version_id=v)
                tagged.append((n, v))
            retention_config = minio.retention.Retention(
                minio.commonconfig.GOVERNANCE,
                datetime.datetime.now(datetime.UTC) + datetime.timedelta(weeks=1000),
            )
            for n, v in zip(names, vids):
                self.client.set_object_retention(
                    self.benchmark, n, config=retention_config, version_id=v
                )
        except Exception:
            self._rollback_tags(tagged)
            raise

    def _rollback_tags(self, tagged_objects: list) -> None:
        """Best-effort removal of version tags; called on create_new_version failure."""
        for obj_name, version_id in tagged_objects:
            try:
                self.client.delete_object_tags(
                    self.benchmark, obj_name, version_id=version_id
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
            self.benchmark, version_filename, io.BytesIO(vv_str.encode()), len(vv_str)
        )
        tags = Tags.new_object_tags()
        tags[str(self.version)] = "1"
        retention_config = minio.retention.Retention(
            minio.commonconfig.GOVERNANCE,
            datetime.datetime.now(datetime.UTC) + datetime.timedelta(weeks=1000),
        )
        self.client.set_object_tags(self.benchmark, version_filename, tags)
        self.client.set_object_retention(
            self.benchmark, version_filename, config=retention_config
        )

    def load_objects(self) -> None:
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
                    response_headers.get("last-modified") or "",
                    "%a, %d %b %Y %H:%M:%S GMT",
                ).strftime("%Y-%m-%d %H:%M:%S.%f+00:00"),
                "size": response_headers.get("content-length"),
                "etag": (response_headers.get("etag") or "").replace('"', ""),
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
            self.load_objects()
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

    def delete_version(self, version):
        raise NotImplementedError


RemoteStorage.register(S3CompatibleStorage)

# Backward-compatibility alias so that existing code importing MinIOStorage
# continues to work without changes during the deprecation period.
MinIOStorage = S3CompatibleStorage
