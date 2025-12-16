import os
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import pytest

from tests.fixtures import bundled_repos  # noqa: F401 - pytest fixture


def pytest_addoption(parser):
    # Check if --keep-files option already exists to avoid conflicts
    try:
        parser.addoption(
            "--keep-files",
            action="store_true",
            default=False,
            help="Keep temporary files after test execution for inspection",
        )
    except ValueError:
        # Option already exists, skip adding it
        pass


@pytest.fixture
def e2e_env(request):
    """Fixture that provides e2e test environment parameters."""
    # Check if the request has a param attribute (indirect parametrization)
    if hasattr(request, "param"):
        params = request.param.copy()  # Create a copy to avoid modifying the original
    else:
        # Default values if no parametrization
        params = {"keep_files": False, "current_dir": False}

    # Override with command-line options if provided
    if request.config.getoption("--keep-files"):
        params["keep_files"] = True

    return params


@pytest.fixture
def keep_files(request):
    """Simple fixture to check if files should be kept."""
    return request.config.getoption("--keep-files")


# S3 test fixtures and helpers


class S3TestEnvironment:
    """Manages S3 test environment setup for local and remote testing.

    This class is shared across all S3 e2e tests via conftest.py fixtures.
    """

    def __init__(self, use_remote: bool = False):
        self.use_remote = use_remote
        self.bucket_name = self._generate_bucket_name()
        self.endpoint = None
        self.access_key = None
        self.secret_key = None
        self.region = None
        self._s3_client = None

    def _generate_bucket_name(self) -> str:
        """Generate unique bucket name with timestamp."""
        timestamp = datetime.now().strftime("%Y%m%d-%H%M%S-%f")[:-3]
        prefix = os.getenv("OB_E2E_BUCKET_PREFIX", "obdata-e2e")
        pr_number = os.getenv("GITHUB_PR_NUMBER") or os.getenv("CI_MERGE_REQUEST_IID")
        if pr_number:
            return f"{prefix}-pr{pr_number}-{timestamp}"
        else:
            return f"{prefix}-{timestamp}"

    def setup_local_environment(self) -> bool:
        """Set up local MinIO environment."""
        try:
            import socket

            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            result = sock.connect_ex(("localhost", 9000))
            sock.close()

            if result == 0:
                print("✓ Found existing MinIO on port 9000")
                self.endpoint = "http://localhost:9000"
                self.region = "us-east-1"

                self.access_key = os.getenv("OB_STORAGE_S3_ACCESS_KEY")
                self.secret_key = os.getenv("OB_STORAGE_S3_SECRET_KEY")

                if not self.access_key or not self.secret_key:
                    print("  Trying to read credentials from /tmp/minio-credentials...")
                    try:
                        with open("/tmp/minio-credentials", "r") as f:
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
                    self.access_key = "minioadmin"
                if not self.secret_key:
                    self.secret_key = "minioadmin123"

                print(f"  Using access key: {self.access_key}")
                return True
            else:
                print("No MinIO found on port 9000")
                return False
        except Exception as e:
            print(f"Failed to setup local MinIO: {e}")
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

        return True

    def setup(self) -> bool:
        """Setup the appropriate environment based on configuration."""
        if self.use_remote:
            return self.setup_remote_environment()
        else:
            return self.setup_local_environment()

    def set_environment_variables(self):
        """Set environment variables for the test."""
        if self.access_key and self.secret_key:
            os.environ["OB_STORAGE_S3_ACCESS_KEY"] = self.access_key
            os.environ["OB_STORAGE_S3_SECRET_KEY"] = self.secret_key
        if self.endpoint:
            os.environ["OB_STORAGE_S3_ENDPOINT_URL"] = self.endpoint
        if self.region:
            os.environ["AWS_DEFAULT_REGION"] = self.region

    def cleanup_environment_variables(self):
        """Clean up environment variables after test."""
        for var in [
            "OB_STORAGE_S3_ACCESS_KEY",
            "OB_STORAGE_S3_SECRET_KEY",
            "OB_STORAGE_S3_ENDPOINT_URL",
            "AWS_DEFAULT_REGION",
        ]:
            if var in os.environ:
                del os.environ[var]

    def get_s3_client(self):
        """Get boto3 S3 client for bucket operations."""
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
                print("boto3 not available - cannot validate S3 operations")
                return None
        return self._s3_client

    def list_bucket_contents(self) -> List[Dict[str, Any]]:
        """List all objects in the test bucket."""
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
        """Download an object from S3 and return its content."""
        s3_client = self.get_s3_client()
        if not s3_client:
            return None
        try:
            response = s3_client.get_object(Bucket=self.bucket_name, Key=key)
            return response["Body"].read()
        except Exception as e:
            print(f"Failed to download object {key}: {e}")
            return None

    def cleanup_bucket(self) -> bool:
        """Clean up the test bucket and all its contents."""
        s3_client = self.get_s3_client()
        if not s3_client:
            print("No S3 client available for cleanup")
            return False

        try:
            deleted_objects = 0

            # Try to disable object lock
            try:
                s3_client.put_object_lock_configuration(
                    Bucket=self.bucket_name, ObjectLockConfiguration={}
                )
                print("  ✓ Disabled object lock on bucket")
            except Exception as e:
                print(f"  Note: Could not disable object lock: {e}")

            # Delete all versions and delete markers
            try:
                versions_response = s3_client.list_object_versions(
                    Bucket=self.bucket_name
                )

                for version in versions_response.get("Versions", []):
                    try:
                        try:
                            s3_client.delete_object(
                                Bucket=self.bucket_name,
                                Key=version["Key"],
                                VersionId=version["VersionId"],
                                BypassGovernanceRetention=True,
                            )
                        except Exception:
                            s3_client.delete_object(
                                Bucket=self.bucket_name,
                                Key=version["Key"],
                                VersionId=version["VersionId"],
                            )
                        deleted_objects += 1
                    except Exception as e:
                        print(f"  Failed to delete version {version['Key']}: {e}")

                for delete_marker in versions_response.get("DeleteMarkers", []):
                    try:
                        s3_client.delete_object(
                            Bucket=self.bucket_name,
                            Key=delete_marker["Key"],
                            VersionId=delete_marker["VersionId"],
                            BypassGovernanceRetention=True,
                        )
                        deleted_objects += 1
                    except Exception as e:
                        print(f"  Failed to delete marker {delete_marker['Key']}: {e}")

            except Exception:
                response = s3_client.list_objects_v2(Bucket=self.bucket_name)
                if "Contents" in response:
                    for obj in response["Contents"]:
                        try:
                            s3_client.delete_object(
                                Bucket=self.bucket_name, Key=obj["Key"]
                            )
                            deleted_objects += 1
                        except Exception as e:
                            print(f"  Failed to delete {obj['Key']}: {e}")

            print(f"  ✓ Deleted {deleted_objects} objects/versions")

            # Delete bucket
            max_retries = 3
            for attempt in range(max_retries):
                try:
                    s3_client.delete_bucket(Bucket=self.bucket_name)
                    print(f"✓ Successfully cleaned up bucket: {self.bucket_name}")
                    return True
                except Exception as bucket_err:
                    if attempt == max_retries - 1:
                        print(
                            f"Failed to delete bucket {self.bucket_name}: {bucket_err}"
                        )
                        return False
                    else:
                        import time

                        time.sleep(2)

        except Exception as e:
            print(f"Failed to cleanup bucket {self.bucket_name}: {e}")
            return False


@pytest.fixture
def s3_config_path():
    """Get the path to the S3 remote config."""
    return Path(__file__).parent / "configs" / "06_s3_remote.yaml"


@pytest.fixture(scope="session")
def _minio_container_e2e():
    """Fixture to set up and tear down the MinIO test container for E2E tests."""
    import sys

    if sys.platform != "linux":
        pytest.skip(
            "MinIO container tests only work on Linux (GitHub Actions limitation)",
            allow_module_level=True,
        )

    from tests.remote.MinIOStorage_setup import MinIOSetup

    # Initialize a MinIO test container with a lifetime of this test session
    minio = MinIOSetup()

    # Yield the container for use in tests
    yield minio

    # Cleanup is handled by MinIOSetup context manager


@pytest.fixture
def s3_environment(_minio_container_e2e):
    """Set up S3 test environment based on configuration."""
    # Only use remote S3 if explicitly requested via OB_E2E_USE_REMOTE_S3
    # Otherwise use local MinIO via testcontainers
    use_remote = os.getenv("OB_E2E_USE_REMOTE_S3", "false").lower() == "true"

    if use_remote:
        env = S3TestEnvironment(use_remote=True)
        if not env.setup():
            pytest.skip("Remote S3 credentials not available")
    else:
        # Use local MinIO from testcontainers
        env = S3TestEnvironment(use_remote=False)
        # Set up connection to the MinIO container
        env.endpoint = _minio_container_e2e.minio.get_config()["endpoint"].replace(
            "localhost", "http://localhost"
        )
        env.access_key = _minio_container_e2e.minio.access_key
        env.secret_key = _minio_container_e2e.minio.secret_key
        env.region = "us-east-1"
        env.bucket_name = env._generate_bucket_name()

    env.set_environment_variables()

    try:
        yield env
    finally:
        print("\n=== S3 Cleanup ===")
        if env.cleanup_bucket():
            print("✓ S3 bucket cleanup successful")
        else:
            print("⚠ S3 bucket cleanup failed or skipped")
        env.cleanup_environment_variables()
