"""
Comprehensive end-to-end test for S3 remote storage.

This test validates the complete S3 workflow in sequence:
1. Run pipeline with --use-remote-storage
2. List files with boto3
3. Check remote files list
4. Check remote files download
5. Check remote files checksum

Uses the 06_s3_remote.yaml config (cartesian product: 3 datasets × 2 methods with D2 excluding M2)
"""

import os
import hashlib
import json
import time
import shutil
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional
import pytest

from tests.cli.cli_setup import OmniCLISetup


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
        """Generate unique bucket name with timestamp."""
        timestamp = datetime.now().strftime("%Y%m%d-%H%M%S-%f")[:-3]

        # Use custom prefix if provided (useful for CI)
        prefix = os.getenv("OB_E2E_BUCKET_PREFIX", "obdata-e2e")

        # Add PR number if available (GitHub Actions)
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
                self.region = "us-east-1"  # MinIO default

                # Priority order: env vars > credentials file > defaults
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

                # Fallback to default MinIO root credentials
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

        # Use default S3 endpoint if not specified (EU region)
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
                    # Force path style for MinIO compatibility
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

            # Step 1: Delete all object versions and delete markers
            try:
                versions_response = s3_client.list_object_versions(
                    Bucket=self.bucket_name
                )

                # Delete all current versions
                for version in versions_response.get("Versions", []):
                    try:
                        s3_client.delete_object(
                            Bucket=self.bucket_name,
                            Key=version["Key"],
                            VersionId=version["VersionId"],
                        )
                        deleted_objects += 1
                    except Exception as e:
                        print(f"  Failed to delete version {version['Key']}: {e}")

                # Delete all delete markers
                for delete_marker in versions_response.get("DeleteMarkers", []):
                    try:
                        s3_client.delete_object(
                            Bucket=self.bucket_name,
                            Key=delete_marker["Key"],
                            VersionId=delete_marker["VersionId"],
                        )
                        deleted_objects += 1
                    except Exception as e:
                        print(f"  Failed to delete marker {delete_marker['Key']}: {e}")

            except Exception:
                # If versioning API fails, fall back to regular object listing
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

            # Step 2: Clean up any incomplete multipart uploads
            try:
                multipart_response = s3_client.list_multipart_uploads(
                    Bucket=self.bucket_name
                )
                for upload in multipart_response.get("Uploads", []):
                    try:
                        s3_client.abort_multipart_upload(
                            Bucket=self.bucket_name,
                            Key=upload["Key"],
                            UploadId=upload["UploadId"],
                        )
                        print(f"  ✓ Aborted multipart upload for {upload['Key']}")
                    except Exception as e:
                        print(
                            f"  Failed to abort multipart upload {upload['Key']}: {e}"
                        )
            except Exception:
                # Multipart operations may not be supported, continue
                pass

            print(f"  ✓ Deleted {deleted_objects} objects/versions")

            # Step 3: Force delete bucket with retry
            max_retries = 3
            for attempt in range(max_retries):
                try:
                    s3_client.delete_bucket(Bucket=self.bucket_name)
                    print(f"✓ Successfully cleaned up bucket: {self.bucket_name}")
                    return True
                except Exception as bucket_err:
                    if attempt == max_retries - 1:
                        print(
                            f"Failed to delete bucket {self.bucket_name} after {max_retries} attempts: {bucket_err}"
                        )
                        return False
                    else:
                        print(
                            f"  Retry {attempt + 1}/{max_retries} for bucket deletion..."
                        )
                        # Wait a moment for eventual consistency
                        import time

                        time.sleep(2)

        except Exception as e:
            print(f"Failed to cleanup bucket {self.bucket_name}: {e}")
            return False


@pytest.fixture
def s3_config_path():
    """Get the path to the S3 remote config."""
    return Path(__file__).parent / "configs" / "06_s3_remote.yaml"


@pytest.fixture
def s3_environment():
    """Set up S3 test environment based on configuration."""
    # Determine if we should use remote S3 (also check for CI environment)
    use_remote = (
        os.getenv("OB_E2E_USE_REMOTE_S3", "false").lower() == "true"
        or os.getenv("CI", "false").lower() == "true"
    )

    env = S3TestEnvironment(use_remote=use_remote)

    if not env.setup():
        if use_remote:
            pytest.skip("Remote S3 credentials not available")
        else:
            pytest.skip("Local MinIO setup failed - docker may not be available")

    env.set_environment_variables()

    try:
        yield env
    finally:
        # Always attempt bucket cleanup after tests
        print("\n=== S3 Cleanup ===")
        if env.cleanup_bucket():
            print("✓ S3 bucket cleanup successful")
        else:
            print("⚠ S3 bucket cleanup failed or skipped")

        env.cleanup_environment_variables()


def update_config_with_s3_settings(
    config_path: Path, tmp_path: Path, s3_env: S3TestEnvironment
) -> Path:
    """Update the S3 config with dynamic bucket name and endpoint."""
    import yaml

    # Read original config
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    # Update storage settings
    config["storage"]["bucket_name"] = s3_env.bucket_name
    if s3_env.endpoint:
        config["storage"]["endpoint"] = s3_env.endpoint

    # Write updated config to tmp_path
    updated_config_path = tmp_path / "06_s3_remote_updated.yaml"
    with open(updated_config_path, "w") as f:
        yaml.dump(config, f)

    return updated_config_path


def calculate_checksum(content: bytes) -> str:
    """Calculate MD5 checksum of content."""
    return hashlib.md5(content).hexdigest()


def setup_test_environment(config_path: Path, tmp_path: Path) -> Path:
    """Set up test environment with config file and dummy conda environment."""
    # Ensure tmp_path directory exists
    tmp_path.mkdir(parents=True, exist_ok=True)

    # Copy config to tmp_path
    config_file_in_tmp = tmp_path / "06_s3_remote.yaml"
    shutil.copy2(config_path, config_file_in_tmp)

    # Create dummy conda environment file for Snakemake 9.x validation
    # NOTE: This should no longer be needed as the conda issue has been fixed
    # DUMMY_CONDA_ENV_CONTENT = """name: omnibenchmark-dummy
    # channels:
    #   - conda-forge
    # dependencies:
    #   - python=3.12
    # """
    # conda_file_path = tmp_path / "conda_not_provided.yml"
    # conda_file_path.write_text(DUMMY_CONDA_ENV_CONTENT)
    #
    # # Also create in system temp directory for cross-platform compatibility
    # try:
    #     import tempfile
    #     system_tmp_path = Path(tempfile.gettempdir()) / "conda_not_provided.yml"
    #     system_tmp_path.write_text(DUMMY_CONDA_ENV_CONTENT)
    # except (PermissionError, OSError):
    #     pass  # Continue if we can't write to system temp

    return config_file_in_tmp


@pytest.mark.e2e
@pytest.mark.e2e_s3
def test_s3_remote_storage_complete_workflow(
    s3_config_path, s3_environment, tmp_path, bundled_repos, keep_files
):
    """
    Test complete S3 remote storage workflow in sequence:

    1. Run pipeline with --use-remote-storage
    2. List files with boto3
    3. Check remote files list matches expected cartesian product
    4. Download files from S3
    5. Verify file checksums and content

    Expected cartesian product:
    - 3 datasets (D1, D2, D3)
    - 2 methods (M1, M2)
    - D2 excludes M2
    - Results: 3 data files + 5 method files = 8 files total
    """
    print("\n=== S3 Complete Workflow Test ===")
    print(f"Bucket: {s3_environment.bucket_name}")
    print(f"Endpoint: {s3_environment.endpoint}")

    # Update config with dynamic S3 settings
    updated_config = update_config_with_s3_settings(
        s3_config_path, tmp_path, s3_environment
    )

    # Setup test environment
    config_file_in_tmp = setup_test_environment(updated_config, tmp_path)

    # ========================================
    # Step 1: Run pipeline with --use-remote-storage
    # ========================================
    print("\n--- Step 1: Running pipeline with --use-remote-storage ---")

    base_args = [
        "run",
        "benchmark",
        "--benchmark",
        str(config_file_in_tmp),
        "--use-remote-storage",
        "--continue-on-error",
        "-y",
    ]

    with OmniCLISetup() as omni:
        result = omni.call(base_args, cwd=str(tmp_path))

        if keep_files:
            print("\nCLI EXECUTION DEBUG:")
            print(f"Return code: {result.returncode}")
            print(f"STDOUT:\n{result.stdout}")
            print(f"STDERR:\n{result.stderr}")

    # Verify CLI execution succeeded
    assert result.returncode == 0, (
        f"S3 pipeline execution failed\n"
        f"STDOUT: {result.stdout}\n"
        f"STDERR: {result.stderr}"
    )
    print("✓ Pipeline execution completed successfully")

    # ========================================
    # Step 2: List files with boto3
    # ========================================
    print("\n--- Step 2: Listing S3 bucket contents with boto3 ---")

    # Wait a moment for S3 consistency
    time.sleep(2)

    bucket_contents = s3_environment.list_bucket_contents()
    assert len(bucket_contents) > 0, "No objects found in S3 bucket"

    object_keys = [obj["Key"] for obj in bucket_contents]
    print(f"✓ Found {len(bucket_contents)} objects in S3:")
    for obj in sorted(bucket_contents, key=lambda x: x["Key"]):
        print(f"  {obj['Key']} ({obj['Size']} bytes)")

    # ========================================
    # Step 3: Check remote files list matches expected structure
    # ========================================
    print("\n--- Step 3: Validating cartesian product structure ---")

    # Expected data files (3 datasets)
    expected_data_files = ["D1_data.json", "D2_data.json", "D3_data.json"]
    found_data_files = [key for key in object_keys if key.endswith("_data.json")]

    for data_file in expected_data_files:
        matching_files = [key for key in found_data_files if data_file in key]
        assert len(matching_files) >= 1, f"Missing data file pattern: {data_file}"

    print(f"✓ Found {len(found_data_files)} data files (expected 3)")

    # Expected method files (5 combinations: D1-M1, D1-M2, D2-M1, D3-M1, D3-M2)
    # Note: D2-M2 should be excluded
    found_method_files = [key for key in object_keys if "_method.json" in key]

    assert len(found_method_files) == 5, (
        f"Expected 5 method files (D2 excludes M2), found {len(found_method_files)}: "
        f"{found_method_files}"
    )

    # Verify D2-M2 combination is excluded
    d2_method_files = [key for key in found_method_files if "D2" in key]
    assert (
        len(d2_method_files) == 1
    ), f"D2 should only have 1 method file (M1), found: {d2_method_files}"

    print(f"✓ Found {len(found_method_files)} method files (expected 5)")
    print("✓ D2-M2 exclusion rule properly applied")

    # ========================================
    # Step 4: Download files from S3
    # ========================================
    print("\n--- Step 4: Downloading and validating files from S3 ---")

    downloaded_files = {}
    file_checksums = {}

    for obj in bucket_contents:
        key = obj["Key"]
        print(f"  Downloading {key}...")

        content = s3_environment.download_object(key)
        assert content is not None, f"Failed to download {key}"

        # Allow empty files (like performances.tsv) but validate non-empty JSON files
        downloaded_files[key] = content
        file_checksums[key] = calculate_checksum(content)

        # Validate JSON structure for key files (skip empty files)
        if key.endswith(".json") and len(content) > 0:
            try:
                json_content = json.loads(content.decode("utf-8"))
                assert isinstance(
                    json_content, dict
                ), f"File {key} should contain JSON object"

                if "result" in json_content:
                    assert isinstance(
                        json_content["result"], (int, float)
                    ), f"File {key} result should be numeric, got {type(json_content['result'])}"
                    print(f"    {key}: result = {json_content['result']}")

            except json.JSONDecodeError as e:
                pytest.fail(f"File {key} contains invalid JSON: {e}")

    print(f"✓ Successfully downloaded {len(downloaded_files)} files")

    # ========================================
    # Step 5: Verify file checksums and content patterns
    # ========================================
    print("\n--- Step 5: Verifying checksums and content patterns ---")

    # Validate expected content patterns based on cartesian product
    expected_patterns = {
        "data_files": {
            "D1": 100,  # evaluate: "100"
            "D2": 200,  # evaluate: "200"
            "D3": 300,  # evaluate: "300"
        },
        "method_files": {
            "D1_M1": 1100,  # 100 + 1000
            "D1_M2": 2100,  # 100 + 2000
            "D2_M1": 1200,  # 200 + 1000 (M2 excluded)
            "D3_M1": 1300,  # 300 + 1000
            "D3_M2": 2300,  # 300 + 2000
        },
    }

    # Check data file patterns
    for dataset, expected_value in expected_patterns["data_files"].items():
        matching_keys = [
            k for k in downloaded_files.keys() if f"{dataset}_data.json" in k
        ]
        assert len(matching_keys) >= 1, f"Missing data file for {dataset}"

        key = matching_keys[0]
        content = json.loads(downloaded_files[key].decode("utf-8"))
        assert (
            content.get("result") == expected_value
        ), f"Dataset {dataset}: expected result {expected_value}, got {content.get('result')}"

    # Check method file patterns with correct S3 path structure
    validated_methods = 0
    for method_combo, expected_value in expected_patterns["method_files"].items():
        # Parse dataset and method from combo (e.g., "D1_M1" -> "D1", "M1")
        dataset, method = method_combo.split("_")

        # Find matching method files with correct dataset and method in path
        matching_keys = [
            k
            for k in downloaded_files.keys()
            if f"{dataset}" in k and f"{method}" in k and "_method.json" in k
        ]

        if len(matching_keys) >= 1:
            key = matching_keys[0]
            content = json.loads(downloaded_files[key].decode("utf-8"))
            actual_result = content.get("result")
            assert (
                actual_result == expected_value
            ), f"Method {method_combo}: expected result {expected_value}, got {actual_result}"
            print(f"    ✓ {method_combo}: {actual_result} (expected {expected_value})")
            validated_methods += 1
        else:
            print(
                f"    ⚠ No matching files for {method_combo} (pattern: {dataset}.*{method}.*_method.json)"
            )

    assert (
        validated_methods == 5
    ), f"Expected to validate 5 method combinations, got {validated_methods}"

    print("✓ All content patterns validated")
    print(f"✓ All checksums calculated for {len(file_checksums)} files")

    # ========================================
    # Step 6: Validate checksum properties
    # ========================================
    print("\n--- Step 6: Validating checksum properties ---")

    # Validate that all JSON files have valid checksums
    json_files_with_checksums = [
        k for k in file_checksums.keys() if k.endswith(".json")
    ]

    for json_file in json_files_with_checksums:
        checksum = file_checksums[json_file]

        # Validate checksum format (MD5 is 32 hex characters)
        assert (
            len(checksum) == 32
        ), f"Invalid checksum length for {json_file}: {len(checksum)}"
        assert all(
            c in "0123456789abcdef" for c in checksum
        ), f"Invalid checksum format for {json_file}: {checksum}"

        # Validate checksum is not empty or default
        assert (
            checksum != "d41d8cd98f00b204e9800998ecf8427e"
        ), f"Empty file checksum detected for {json_file}"

        print(f"    ✓ {json_file}: {checksum}")

    print(
        f"✓ Validated checksum format for {len(json_files_with_checksums)} JSON files"
    )
    print("✓ All checksums are valid MD5 hashes")

    # ========================================
    # Final Summary
    # ========================================
    print("\n=== S3 Workflow Test Summary ===")
    print("✓ Pipeline executed successfully with --use-remote-storage")
    print(f"✓ boto3 listed {len(bucket_contents)} objects from S3")
    print(
        "✓ Cartesian product structure validated (3 datasets, 2 methods, 1 exclusion)"
    )
    print(f"✓ All {len(downloaded_files)} files downloaded and validated")
    print("✓ Content patterns match expected values")
    print(f"✓ Bucket: {s3_environment.bucket_name}")

    if keep_files:
        print("\nFile checksums:")
        for key, checksum in file_checksums.items():
            print(f"  {key}: {checksum}")

    # The bucket will be automatically cleaned up by the s3_environment fixture
    print("✓ Test completed successfully")
