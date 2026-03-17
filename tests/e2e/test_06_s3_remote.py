"""
S3 remote storage e2e tests.

Tests share a single pipeline run via a module-scoped fixture so the expensive
pipeline only executes once. Each test validates a distinct concern:

1. test_s3_pipeline_execution         - pipeline runs and bucket is created
2. test_s3_bucket_structure           - cartesian product of files is correct
3. test_s3_file_content_and_checksums - files download with correct content
4. test_s3_version_create_and_list    - version 1.0 is created and listed
5. test_s3_benchmark_modification_and_version_2 - v2 workflow end-to-end
"""

import hashlib
import json
import os
import shutil
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import pytest

from tests.cli.cli_setup import OmniCLISetup


# ---------------------------------------------------------------------------
# S3 environment helper (local fixture copy, kept self-contained)
# ---------------------------------------------------------------------------


class S3TestEnvironment:
    """Manages S3 test environment setup for local and remote testing."""

    def __init__(self, use_remote: bool = False):
        self.use_remote = use_remote
        self.bucket_name = self._generate_bucket_name()
        self.endpoint = None
        self.access_key = None
        self.secret_key = None
        self.region = None
        self._s3_client = None

    def _generate_bucket_name(self) -> str:
        timestamp = datetime.now().strftime("%Y%m%d-%H%M%S-%f")[:-3]
        prefix = os.getenv("OB_E2E_BUCKET_PREFIX", "obdata-e2e")
        pr_number = os.getenv("GITHUB_PR_NUMBER") or os.getenv("CI_MERGE_REQUEST_IID")
        if pr_number:
            return f"{prefix}-pr{pr_number}-{timestamp}"
        return f"{prefix}-{timestamp}"

    def setup_local_environment(self) -> bool:
        """Set up local RustFS environment."""
        try:
            import socket

            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            result = sock.connect_ex(("localhost", 9000))
            sock.close()

            if result == 0:
                print("✓ Found existing RustFS on port 9000")
                self.endpoint = "http://localhost:9000"
                self.region = "us-east-1"

                self.access_key = os.getenv("OB_STORAGE_S3_ACCESS_KEY")
                self.secret_key = os.getenv("OB_STORAGE_S3_SECRET_KEY")

                if not self.access_key or not self.secret_key:
                    print(
                        "  Trying to read credentials from /tmp/rustfs-credentials..."
                    )
                    try:
                        with open("/tmp/rustfs-credentials", "r") as f:
                            for line in f:
                                if (
                                    line.startswith("OB_STORAGE_S3_ACCESS_KEY=")
                                    and not self.access_key
                                ):
                                    self.access_key = line.split("=", 1)[1].strip()
                                elif (
                                    line.startswith("OB_STORAGE_S3_SECRET_KEY=")
                                    and not self.secret_key
                                ):
                                    self.secret_key = line.split("=", 1)[1].strip()
                        print("  ✓ Found credentials in file")
                    except FileNotFoundError:
                        print("  No credentials file found, using defaults")

                if not self.access_key:
                    self.access_key = "rustfsadmin"
                if not self.secret_key:
                    self.secret_key = "rustfsadmin"

                print(f"  Using access key: {self.access_key}")
                return True

            print("No RustFS found on port 9000")
            return False

        except Exception as e:
            print(f"Failed to setup local RustFS: {e}")
            return False

    def setup_remote_environment(self) -> bool:
        """Set up remote S3 environment using environment variables."""
        self.access_key = os.getenv("OB_STORAGE_S3_ACCESS_KEY")
        self.secret_key = os.getenv("OB_STORAGE_S3_SECRET_KEY")
        self.endpoint = os.getenv("OB_STORAGE_S3_ENDPOINT_URL")
        self.region = os.getenv("AWS_DEFAULT_REGION", "eu-central-1")

        if not self.access_key or not self.secret_key:
            return False

        if not self.endpoint:
            self.endpoint = "https://s3.eu-central-1.amazonaws.com"

        print(
            f"Remote S3: endpoint={self.endpoint} region={self.region} bucket={self.bucket_name}"
        )
        return True

    def setup(self) -> bool:
        return (
            self.setup_remote_environment()
            if self.use_remote
            else self.setup_local_environment()
        )

    def set_environment_variables(self):
        if self.access_key and self.secret_key:
            os.environ["OB_STORAGE_S3_ACCESS_KEY"] = self.access_key
            os.environ["OB_STORAGE_S3_SECRET_KEY"] = self.secret_key
        if self.endpoint:
            os.environ["OB_STORAGE_S3_ENDPOINT_URL"] = self.endpoint
        if self.region:
            os.environ["AWS_DEFAULT_REGION"] = self.region

    def cleanup_environment_variables(self):
        for var in [
            "OB_STORAGE_S3_ACCESS_KEY",
            "OB_STORAGE_S3_SECRET_KEY",
            "OB_STORAGE_S3_ENDPOINT_URL",
            "AWS_DEFAULT_REGION",
        ]:
            os.environ.pop(var, None)

    def test_bucket_creation_permissions(self) -> bool:
        s3_client = self.get_s3_client()
        if not s3_client:
            return False
        try:
            if self.region == "us-east-1":
                s3_client.create_bucket(Bucket=self.bucket_name)
            else:
                s3_client.create_bucket(
                    Bucket=self.bucket_name,
                    CreateBucketConfiguration={"LocationConstraint": self.region},
                )
            return True
        except Exception as e:
            print(f"❌ Failed to create bucket: {e}")
            return False

    def get_s3_client(self):
        if self._s3_client is None:
            try:
                import boto3

                self._s3_client = boto3.client(
                    "s3",
                    endpoint_url=self.endpoint,
                    aws_access_key_id=self.access_key,
                    aws_secret_access_key=self.secret_key,
                    region_name=self.region,
                    config=boto3.session.Config(s3={"addressing_style": "path"})
                    if "localhost" in self.endpoint
                    else None,
                )
            except ImportError:
                print("boto3 not available")
                return None
        return self._s3_client

    def list_bucket_contents(self) -> List[Dict[str, Any]]:
        s3_client = self.get_s3_client()
        if not s3_client:
            return []
        try:
            response = s3_client.list_objects_v2(Bucket=self.bucket_name)
            return response.get("Contents", [])
        except Exception as e:
            print(f"Failed to list bucket contents: {e}")
            return []

    def download_object(self, key: str) -> Optional[bytes]:
        s3_client = self.get_s3_client()
        if not s3_client:
            return None
        try:
            response = s3_client.get_object(Bucket=self.bucket_name, Key=key)
            return response["Body"].read()
        except Exception as e:
            print(f"Failed to download {key}: {e}")
            return None

    def cleanup_bucket(self) -> bool:
        s3_client = self.get_s3_client()
        if not s3_client:
            return False

        try:
            deleted = 0

            try:
                s3_client.put_object_lock_configuration(
                    Bucket=self.bucket_name, ObjectLockConfiguration={}
                )
            except Exception:
                pass

            try:
                versions_resp = s3_client.list_object_versions(Bucket=self.bucket_name)

                for ver in versions_resp.get("Versions", []):
                    try:
                        try:
                            s3_client.delete_object(
                                Bucket=self.bucket_name,
                                Key=ver["Key"],
                                VersionId=ver["VersionId"],
                                BypassGovernanceRetention=True,
                            )
                        except Exception:
                            s3_client.delete_object(
                                Bucket=self.bucket_name,
                                Key=ver["Key"],
                                VersionId=ver["VersionId"],
                            )
                        deleted += 1
                    except Exception:
                        try:
                            s3_client.put_object_legal_hold(
                                Bucket=self.bucket_name,
                                Key=ver["Key"],
                                VersionId=ver["VersionId"],
                                LegalHold={"Status": "OFF"},
                            )
                            s3_client.delete_object(
                                Bucket=self.bucket_name,
                                Key=ver["Key"],
                                VersionId=ver["VersionId"],
                                BypassGovernanceRetention=True,
                            )
                            deleted += 1
                        except Exception as e2:
                            print(f"  Failed to delete {ver['Key']}: {e2}")

                for marker in versions_resp.get("DeleteMarkers", []):
                    try:
                        s3_client.delete_object(
                            Bucket=self.bucket_name,
                            Key=marker["Key"],
                            VersionId=marker["VersionId"],
                            BypassGovernanceRetention=True,
                        )
                        deleted += 1
                    except Exception as e:
                        print(f"  Failed to delete marker {marker['Key']}: {e}")

            except Exception:
                resp = s3_client.list_objects_v2(Bucket=self.bucket_name)
                for obj in resp.get("Contents", []):
                    try:
                        s3_client.delete_object(
                            Bucket=self.bucket_name,
                            Key=obj["Key"],
                            BypassGovernanceRetention=True,
                        )
                        deleted += 1
                    except Exception as e:
                        print(f"  Failed to delete {obj['Key']}: {e}")

            try:
                mp_resp = s3_client.list_multipart_uploads(Bucket=self.bucket_name)
                for upload in mp_resp.get("Uploads", []):
                    try:
                        s3_client.abort_multipart_upload(
                            Bucket=self.bucket_name,
                            Key=upload["Key"],
                            UploadId=upload["UploadId"],
                        )
                    except Exception:
                        pass
            except Exception:
                pass

            print(f"  ✓ Deleted {deleted} objects/versions")

            for attempt in range(3):
                try:
                    s3_client.delete_bucket(Bucket=self.bucket_name)
                    print(f"✓ Cleaned up bucket: {self.bucket_name}")
                    return True
                except Exception as e:
                    if attempt == 2:
                        print(f"Failed to delete bucket: {e}")
                        return False
                    time.sleep(2)

        except Exception as e:
            print(f"Failed to cleanup bucket {self.bucket_name}: {e}")
            return False


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def _update_config(
    config_path: Path, tmp_path: Path, s3_env: S3TestEnvironment
) -> Path:
    """Write a copy of the config with the test bucket name and endpoint injected."""
    import yaml

    with open(config_path) as f:
        config = yaml.safe_load(f)

    config["storage"]["bucket_name"] = s3_env.bucket_name
    if s3_env.endpoint:
        config["storage"]["endpoint"] = s3_env.endpoint

    out = tmp_path / "06_s3_remote_updated.yaml"
    with open(out, "w") as f:
        yaml.dump(config, f)
    return out


def _checksum(content: bytes) -> str:
    return hashlib.md5(content).hexdigest()


def _wait_for_bucket_objects(
    env: S3TestEnvironment, min_count: int = 1, timeout: int = 30
) -> List[Dict[str, Any]]:
    """Poll until at least *min_count* objects appear in the bucket, then return them.

    Replaces bare time.sleep() calls with a condition-driven wait so tests
    fail fast on quick systems and stay reliable on slow ones.
    """
    deadline = time.time() + timeout
    while time.time() < deadline:
        contents = env.list_bucket_contents()
        if len(contents) >= min_count:
            return contents
        time.sleep(1)
    contents = env.list_bucket_contents()
    raise TimeoutError(
        f"Expected ≥{min_count} objects in bucket after {timeout}s, got {len(contents)}"
    )


# ---------------------------------------------------------------------------
# Module-scoped fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def s3_config_path():
    return Path(__file__).parent / "configs" / "06_s3_remote.yaml"


@pytest.fixture(scope="module")
def s3_environment(s3_config_path):
    """Single bucket that persists across all tests in this module."""
    use_remote = (
        os.getenv("OB_E2E_USE_REMOTE_S3", "false").lower() == "true"
        or os.getenv("CI", "false").lower() == "true"
    )
    env = S3TestEnvironment(use_remote=use_remote)
    if not env.setup():
        pytest.skip(
            "Remote S3 credentials not available"
            if use_remote
            else "Local RustFS setup failed - docker may not be available"
        )
    env.set_environment_variables()

    yield env

    print("\n=== S3 Cleanup ===")
    if env.cleanup_bucket():
        print("✓ S3 bucket cleanup successful")
    else:
        print("⚠ S3 bucket cleanup failed or skipped")
    env.cleanup_environment_variables()


@dataclass
class S3WorkflowState:
    config_path: Path
    tmp_path: Path
    pipeline_returncode: int
    pipeline_stdout: str
    pipeline_stderr: str


@pytest.fixture(scope="module")
def s3_workflow(s3_environment, s3_config_path, tmp_path_factory):
    """Run the S3 pipeline once; all tests in this module share the result."""
    tmp_path = tmp_path_factory.mktemp("s3_workflow")

    # Replicate the bundled_repos fixture at module scope
    bundles_src = Path(__file__).resolve().parents[1] / "data" / "bundles"
    (tmp_path / "bundles").symlink_to(bundles_src, target_is_directory=True)

    if (
        s3_environment.use_remote
        and not s3_environment.test_bucket_creation_permissions()
    ):
        pytest.fail(
            "S3 credentials do not have permission to create buckets. "
            "Ensure the IAM user/role has s3:CreateBucket permission."
        )

    config_path = _update_config(s3_config_path, tmp_path, s3_environment)
    # Copy to the working directory so relative paths resolve correctly
    config_in_cwd = tmp_path / "06_s3_remote.yaml"
    shutil.copy2(config_path, config_in_cwd)

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(config_in_cwd),
                "--use-remote-storage",
                "--continue-on-error",
                "-y",
            ],
            cwd=str(tmp_path),
        )

    yield S3WorkflowState(
        config_path=config_in_cwd,
        tmp_path=tmp_path,
        pipeline_returncode=result.returncode,
        pipeline_stdout=result.stdout,
        pipeline_stderr=result.stderr,
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@pytest.mark.e2e_s3
@pytest.mark.timeout(300)
def test_s3_pipeline_execution(s3_workflow, s3_environment):
    """Pipeline runs successfully and the S3 bucket is created."""
    assert (
        s3_workflow.pipeline_returncode == 0
    ), f"S3 pipeline failed\nSTDOUT: {s3_workflow.pipeline_stdout}\nSTDERR: {s3_workflow.pipeline_stderr}"

    s3_client = s3_environment.get_s3_client()
    if s3_client:
        try:
            s3_client.head_bucket(Bucket=s3_environment.bucket_name)
        except Exception as e:
            raise AssertionError(f"S3 bucket was not created: {e}") from e


@pytest.mark.e2e_s3
@pytest.mark.timeout(60)
def test_s3_bucket_structure(s3_workflow, s3_environment):
    """Bucket contains the expected cartesian product (3 data + 5 method files)."""
    contents = _wait_for_bucket_objects(s3_environment, min_count=8)
    keys = [obj["Key"] for obj in contents]

    # 3 data files
    data_files = [k for k in keys if k.endswith("_data.json")]
    for name in ["D1_data.json", "D2_data.json", "D3_data.json"]:
        assert any(name in k for k in data_files), f"Missing data file: {name}"
    assert len(data_files) == 3, f"Expected 3 data files, got {len(data_files)}"

    # 5 method files — D2 excludes M2
    method_files = [k for k in keys if "_method.json" in k]
    assert (
        len(method_files) == 5
    ), f"Expected 5 method files (D2 excludes M2), got {len(method_files)}: {method_files}"
    d2_methods = [k for k in method_files if "D2" in k]
    assert len(d2_methods) == 1, f"D2 should only have M1, got: {d2_methods}"


@pytest.mark.e2e_s3
@pytest.mark.timeout(120)
def test_s3_file_content_and_checksums(s3_workflow, s3_environment):
    """All files download with the expected content values and valid MD5 checksums."""
    contents = _wait_for_bucket_objects(s3_environment, min_count=8)

    downloaded: Dict[str, bytes] = {}
    for obj in contents:
        key = obj["Key"]
        data = s3_environment.download_object(key)
        assert data is not None, f"Failed to download {key}"
        downloaded[key] = data

    # Expected values from the benchmark config
    expected_data = {"D1": 100, "D2": 200, "D3": 300}
    expected_methods = {
        "D1_M1": 1100,
        "D1_M2": 2100,
        "D2_M1": 1200,
        "D3_M1": 1300,
        "D3_M2": 2300,
    }

    for dataset, value in expected_data.items():
        key = next(k for k in downloaded if f"{dataset}_data.json" in k)
        actual = json.loads(downloaded[key])["result"]
        assert actual == value, f"{dataset}: expected {value}, got {actual}"

    validated = 0
    for combo, value in expected_methods.items():
        ds, method = combo.split("_")
        matches = [
            k for k in downloaded if ds in k and method in k and "_method.json" in k
        ]
        if matches:
            actual = json.loads(downloaded[matches[0]])["result"]
            assert actual == value, f"{combo}: expected {value}, got {actual}"
            validated += 1
    assert validated == 5, f"Expected 5 method validations, got {validated}"

    # Checksums must be valid, non-empty MD5 hashes
    empty_md5 = "d41d8cd98f00b204e9800998ecf8427e"
    for key, content in downloaded.items():
        if key.endswith(".json"):
            cs = _checksum(content)
            assert len(cs) == 32 and all(c in "0123456789abcdef" for c in cs)
            assert cs != empty_md5, f"Empty file detected: {key}"


@pytest.mark.e2e_s3
@pytest.mark.timeout(60)
def test_s3_version_create_and_list(s3_workflow):
    """Version 1.0 is created and appears in the version list."""
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "remote",
                "version",
                "create",
                "--benchmark",
                str(s3_workflow.config_path),
            ],
            cwd=str(s3_workflow.tmp_path),
        )
    assert (
        result.returncode == 0
    ), f"version create failed\nSTDOUT: {result.stdout}\nSTDERR: {result.stderr}"
    assert "Create a new benchmark version" in result.stdout

    with OmniCLISetup() as omni:
        result = omni.call(
            ["remote", "version", "list", "--benchmark", str(s3_workflow.config_path)],
            cwd=str(s3_workflow.tmp_path),
        )
    assert (
        result.returncode == 0
    ), f"version list failed\nSTDOUT: {result.stdout}\nSTDERR: {result.stderr}"
    assert "1.0" in result.stdout, f"Version 1.0 not found: {result.stdout}"


@pytest.mark.e2e_s3
@pytest.mark.timeout(300)
def test_s3_benchmark_modification_and_version_2(s3_workflow):
    """Adding method M3, re-running, and creating version 2.0 all succeed."""
    import yaml

    # Work on a separate copy so other tests are unaffected
    v2_config = s3_workflow.tmp_path / "06_s3_remote_v2.yaml"
    shutil.copy2(s3_workflow.config_path, v2_config)

    with open(v2_config) as f:
        config = yaml.safe_load(f)

    config["version"] = "2.0"
    methods_stage = next(s for s in config["stages"] if s["id"] == "methods")
    methods_stage["modules"].append(
        {
            "id": "M3",
            "name": "Method 3",
            "software_environment": "host",
            "repository": {
                "url": "bundles/dummymodule_4ff8427.bundle",
                "commit": "4ff8427",
            },
            "parameters": [
                {"evaluate": "input+3000", "input": "data.raw", "kind": "method"}
            ],
        }
    )
    with open(v2_config, "w") as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(v2_config),
                "--use-remote-storage",
                "--continue-on-error",
                "-y",
            ],
            cwd=str(s3_workflow.tmp_path),
        )
    assert (
        result.returncode == 0
    ), f"Modified benchmark run failed\nSTDOUT: {result.stdout}\nSTDERR: {result.stderr}"

    with OmniCLISetup() as omni:
        result = omni.call(
            ["remote", "version", "create", "--benchmark", str(v2_config)],
            cwd=str(s3_workflow.tmp_path),
        )
    assert (
        result.returncode == 0
    ), f"Version 2.0 create failed\nSTDOUT: {result.stdout}\nSTDERR: {result.stderr}"

    # List all versions using the original config (same bucket)
    with OmniCLISetup() as omni:
        result = omni.call(
            ["remote", "version", "list", "--benchmark", str(s3_workflow.config_path)],
            cwd=str(s3_workflow.tmp_path),
        )
    assert result.returncode == 0
    assert "1.0" in result.stdout, f"Version 1.0 missing: {result.stdout}"
    assert "2.0" in result.stdout, f"Version 2.0 missing: {result.stdout}"
