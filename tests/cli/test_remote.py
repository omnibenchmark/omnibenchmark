import hashlib
import os
import shutil
import json

from pathlib import Path

from typing import Optional

from click.testing import CliRunner

from .cli_setup import OmniCLISetup
from .asserts import assert_startswith, assert_in_output
from .fixtures import minio_storage, _minio_container  # noqa: F401


def do_first_run(clisetup, file: str, cwd: Optional[str] = None):
    run1 = clisetup.call(
        [
            "run",
            file,
            "--use-remote-storage",
        ],
        cwd=cwd,
    )
    assert run1.returncode == 0


def test_create_version(minio_storage):  # noqa: F811
    with OmniCLISetup() as omni:
        # first we run the benchmark to generate data
        do_first_run(omni, str(minio_storage.benchmark_file))

        # then we create a new version
        run2 = omni.call(
            [
                "remote",
                "version",
                "create",
                "--benchmark",
                str(minio_storage.benchmark_file),
            ]
        )

        # Check that the command succeeded
        if run2.returncode != 0:
            print(f"STDOUT: {run2.stdout}")
            print(f"STDERR: {run2.stderr}")
        assert run2.returncode == 0

        store = minio_storage.get_storage_client()
        store.set_version("1.0")
        store._get_objects()

        assert "versions/1.0.csv" in store.files.keys()


def test_version_sorting_no_regression():
    """Test that Version objects can be sorted without TypeError.

    This test prevents regression of the bug where ss.versions.sort(key=Version)
    caused TypeError: expected string or bytes-like object, got 'Version'

    The bug occurred when trying to sort a list of Version objects using
    Version as a key function instead of letting Version objects sort naturally.
    """
    from packaging.version import Version

    # Create a list of Version objects like those in ss.versions
    versions = [
        Version("2.0"),
        Version("1.0"),
        Version("1.5"),
        Version("0.9"),
    ]

    # This should work without throwing TypeError
    # The original bug was: versions.sort(key=Version) which is incorrect
    # Correct approach: versions.sort() since Version objects are comparable
    versions.sort()

    # Verify versions are sorted correctly
    expected_order = [
        Version("0.9"),
        Version("1.0"),
        Version("1.5"),
        Version("2.0"),
    ]

    assert versions == expected_order

    # Also test that the buggy approach would fail
    import pytest

    versions_copy = [Version("2.0"), Version("1.0")]

    # This would be the original buggy code that caused the TypeError
    with pytest.raises(
        TypeError, match="expected string or bytes-like object, got 'Version'"
    ):
        versions_copy.sort(key=Version)


def get_md5_hash(content: str) -> str:
    hash_md5 = hashlib.md5()
    hash_md5.update(content.encode())
    return hash_md5.hexdigest()


def test_list_files(minio_storage):  # noqa: F811
    with OmniCLISetup() as omni:
        # first we run the benchmark to generate data
        do_first_run(omni, str(minio_storage.benchmark_file), cwd=minio_storage.out_dir)

        # file content depends on the full file name, thus create a file with known content
        file_content = f"1. Created dataset file .snakemake/storage/s3/{minio_storage.bucket_name}/out/data/D1/default/D1.meta.json.\n"
        hash = get_md5_hash(file_content)

        # now we create a new version of the benchmark
        run2 = omni.call(
            [
                "remote",
                "version",
                "create",
                "--benchmark",
                str(minio_storage.benchmark_file),
            ],
            cwd=minio_storage.out_dir,
        )

        assert run2.returncode == 0

        # retrieve the stored files
        run3 = omni.call(
            [
                "remote",
                "files",
                "list",
                "--benchmark",
                str(minio_storage.benchmark_file),
            ]
        )
        assert run3.returncode == 0

        # compare the md5 hash
        expected = f"{hash} out/data/D1/default/D1.meta.json"
        assert_in_output(run3.stdout, expected)


def test_download_files(minio_storage):  # noqa: F811
    from pathlib import Path
    import os

    with OmniCLISetup() as omni:
        # First, run the benchmark to generate data
        do_first_run(omni, str(minio_storage.benchmark_file), cwd=minio_storage.out_dir)

        # Now create a version of the benchmark
        run2 = omni.call(
            [
                "remote",
                "version",
                "create",
                "--benchmark",
                str(minio_storage.benchmark_file),
            ],
            cwd=minio_storage.out_dir,
        )
        assert run2.returncode == 0

        # Test the download functionality
        # First, clear the output directory to ensure we're testing actual downloads
        for item in Path(minio_storage.out_dir / "out").glob("*"):
            if item.is_file():
                os.unlink(item)
            elif item.is_dir() and item.name != ".snakemake":
                shutil.rmtree(item)

        # Now download the files
        run3 = omni.call(
            [
                "remote",
                "files",
                "download",
                "--benchmark",
                str(minio_storage.benchmark_file),
            ],
            cwd=minio_storage.out_dir,
        )

        assert run3.returncode == 0

        # Verify downloaded files exist by comparing with remote files
        store = minio_storage.get_storage_client()
        store.set_version("1.0")
        store._get_objects()

        # Get remote files (focusing on out directory files)
        remote_files = [f for f in store.files.keys() if Path(f).parts[0] == "out"]

        # Compare to local files
        for remote_file in remote_files:
            local_path = Path(minio_storage.out_dir) / remote_file
            assert local_path.exists(), f"Downloaded file {local_path} does not exist"


def test_S3_storage_missing_access_key(minio_storage):  # noqa: F811
    with OmniCLISetup() as omni:
        # Save the original value to restore it later
        original_access_key = os.environ.get("OB_STORAGE_S3_ACCESS_KEY")

        try:
            # Remove access_key but keep secret_key
            os.environ.pop("OB_STORAGE_S3_ACCESS_KEY", None)

            # Try to create a version (should fail)
            run = omni.call(
                [
                    "remote",
                    "version",
                    "create",
                    "--benchmark",
                    str(minio_storage.benchmark_file),
                ],
                cwd=minio_storage.out_dir,
            )

            # Verify failure and error message
            assert run.returncode == 1
            expected_output = "[ERROR] Missing S3 credentials. Set OB_STORAGE_S3_ACCESS_KEY and OB_STORAGE_S3_SECRET_KEY, or OB_STORAGE_S3_CONFIG"
            assert_startswith(run.stdout, expected_output)

        finally:
            # Restore the original environment variable
            if original_access_key:
                os.environ["OB_STORAGE_S3_ACCESS_KEY"] = original_access_key
            else:
                os.environ.pop("OB_STORAGE_S3_ACCESS_KEY", None)


def test_S3_storage_credentials_from_file(minio_storage):  # noqa: F811
    with OmniCLISetup() as omni:
        # Save original environment variables to restore later
        original_access_key = os.environ.get("OB_STORAGE_S3_ACCESS_KEY")
        original_secret_key = os.environ.get("OB_STORAGE_S3_SECRET_KEY")
        original_config = os.environ.get("OB_STORAGE_S3_CONFIG")

        try:
            # Create a JSON file with S3 credentials
            auth_options = {
                "access_key": os.environ.get("OB_STORAGE_S3_ACCESS_KEY"),
                "secret_key": os.environ.get("OB_STORAGE_S3_SECRET_KEY"),
                "endpoint_url": os.environ.get(
                    "OB_STORAGE_S3_ENDPOINT_URL", "http://localhost:9000"
                ),
                "region": os.environ.get("OB_STORAGE_S3_REGION", "us-east-1"),
            }

            storage_s3_json = Path(minio_storage.out_dir) / "storage_s3.json"
            with open(storage_s3_json, "w") as f:
                json.dump(auth_options, f)

            # Set up environment to use the config file instead of individual variables
            os.environ["OB_STORAGE_S3_CONFIG"] = str(storage_s3_json)
            os.environ.pop("OB_STORAGE_S3_ACCESS_KEY", None)
            os.environ.pop("OB_STORAGE_S3_SECRET_KEY", None)

            # Run the benchmark
            do_first_run(
                omni, str(minio_storage.benchmark_file), cwd=minio_storage.out_dir
            )

            # Try to create a version (should succeed)
            run = omni.call(
                [
                    "remote",
                    "version",
                    "create",
                    "--benchmark",
                    str(minio_storage.benchmark_file),
                ],
                cwd=minio_storage.out_dir,
            )

            # Verify success and message
            assert run.returncode == 0
            expected_output = "Create a new benchmark version"
            assert_startswith(run.stdout, expected_output)

        finally:
            # Restore original environment variables
            if original_access_key:
                os.environ["OB_STORAGE_S3_ACCESS_KEY"] = original_access_key
            else:
                os.environ.pop("OB_STORAGE_S3_ACCESS_KEY", None)

            if original_secret_key:
                os.environ["OB_STORAGE_S3_SECRET_KEY"] = original_secret_key
            else:
                os.environ.pop("OB_STORAGE_S3_SECRET_KEY", None)

            if original_config:
                os.environ["OB_STORAGE_S3_CONFIG"] = original_config
            else:
                os.environ.pop("OB_STORAGE_S3_CONFIG", None)


def test_create_policy_cli_signature(minio_storage):  # noqa: F811
    """Test that create-policy command correctly receives benchmark parameter."""
    from omnibenchmark.cli.remote import remote

    runner = CliRunner()

    # Test that the command accepts --benchmark option and passes it correctly
    result = runner.invoke(
        remote,
        [
            "policy",
            "create",
            "--benchmark",
            str(minio_storage.benchmark_file),
        ],
    )

    # The command should execute without parameter binding errors
    # It may fail for other reasons (e.g., policy creation), but not due to parameter mismatch
    # If there was a parameter mismatch, we'd see a TypeError about unexpected keyword argument
    assert "unexpected keyword argument" not in result.output.lower()
    assert "takes 0 positional arguments" not in result.output.lower()


def test_missing_S3_storage_credentials_in_config_file(minio_storage):  # noqa: F811
    with OmniCLISetup() as omni:
        # Save original environment variables
        original_access_key = os.environ.get("OB_STORAGE_S3_ACCESS_KEY")
        original_secret_key = os.environ.get("OB_STORAGE_S3_SECRET_KEY")
        original_config = os.environ.get("OB_STORAGE_S3_CONFIG")

        try:
            # First run the benchmark with VALID credentials
            do_first_run(
                omni, str(minio_storage.benchmark_file), cwd=minio_storage.out_dir
            )

            # NOW create a JSON file with incomplete S3 credentials (missing access_key)
            auth_options = {
                # Intentionally omit access_key
                "secret_key": os.environ.get("OB_STORAGE_S3_SECRET_KEY"),
                "endpoint_url": os.environ.get(
                    "OB_STORAGE_S3_ENDPOINT_URL", "http://localhost:9000"
                ),
                "region": os.environ.get("OB_STORAGE_S3_REGION", "us-east-1"),
            }

            storage_s3_json = Path(minio_storage.out_dir) / "storage_s3.json"
            with open(storage_s3_json, "w") as f:
                json.dump(auth_options, f)

            # Set up environment to use the config file instead of individual variables
            os.environ["OB_STORAGE_S3_CONFIG"] = str(storage_s3_json)
            os.environ.pop("OB_STORAGE_S3_ACCESS_KEY", None)
            os.environ.pop("OB_STORAGE_S3_SECRET_KEY", None)

            # Try to create a version (should fail)
            run = omni.call(
                [
                    "remote",
                    "version",
                    "create",
                    "--benchmark",
                    str(minio_storage.benchmark_file),
                ],
                cwd=minio_storage.out_dir,
            )

            # Verify failure and error message
            assert run.returncode == 1
            expected_output = "[ERROR] Missing access_key or secret_key in config file"
            assert_startswith(run.stdout, expected_output)

        finally:
            # Restore original environment variables
            if original_access_key:
                os.environ["OB_STORAGE_S3_ACCESS_KEY"] = original_access_key
            else:
                os.environ.pop("OB_STORAGE_S3_ACCESS_KEY", None)

            if original_secret_key:
                os.environ["OB_STORAGE_S3_SECRET_KEY"] = original_secret_key
            else:
                os.environ.pop("OB_STORAGE_S3_SECRET_KEY", None)

            if original_config:
                os.environ["OB_STORAGE_S3_CONFIG"] = original_config
            else:
                os.environ.pop("OB_STORAGE_S3_CONFIG", None)
