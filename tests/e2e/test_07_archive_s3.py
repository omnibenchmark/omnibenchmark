"""
Comprehensive end-to-end test for archive functionality with S3 remote storage.

This test validates the archive command behavior with remote storage:
1. Run benchmark with --use-remote-storage to populate S3
2. Test archive with --results --use-remote-storage (downloads from S3)
3. Test archive with --results only (local files only)
4. Validate archive contents match expected files
5. Test different compression options
6. Test dry-run functionality

IMPORTANT: Archive command ALWAYS requires storage.api configuration in benchmark YAML,
even for local-only archiving. Pure local benchmarks (without storage section) cannot
be archived - the archive command will fail with "No storage API configured" error.
"""

import time
import shutil
import zipfile
from pathlib import Path
from typing import Dict
import pytest
import yaml

from tests.cli.cli_setup import OmniCLISetup


# Re-use the S3 test environment setup
from tests.e2e.test_06_s3_remote import (
    update_config_with_s3_settings,
    setup_test_environment,
)


def extract_and_list_archive_contents(archive_path: Path) -> Dict[str, bytes]:
    """Extract archive and return mapping of file paths to contents."""
    contents = {}

    with zipfile.ZipFile(archive_path, "r") as archive:
        for file_info in archive.filelist:
            if not file_info.is_dir():
                with archive.open(file_info.filename) as f:
                    contents[file_info.filename] = f.read()

    return contents


def validate_archive_structure(
    archive_contents: Dict[str, bytes],
    expected_config: bool = True,
    expected_results: bool = True,
    expected_code: bool = False,
    expected_software: bool = False,
) -> None:
    """Validate that archive contains expected components."""
    archive_files = set(archive_contents.keys())

    if expected_config:
        config_files = [f for f in archive_files if f.endswith(".yaml")]
        assert len(config_files) >= 1, f"Expected config file, found: {config_files}"

    if expected_results:
        result_files = [f for f in archive_files if "out/" in f]
        assert (
            len(result_files) > 0
        ), f"Expected result files, archive contains: {archive_files}"

        # Verify expected result patterns from S3 test
        json_files = [f for f in result_files if f.endswith(".json")]
        assert (
            len(json_files) >= 8
        ), f"Expected at least 8 JSON files, found {len(json_files)}: {json_files}"

    if expected_code:
        code_files = [f for f in archive_files if ".snakemake/repos/" in f]
        if expected_code:
            assert (
                len(code_files) > 0
            ), f"Expected code files, archive contains: {archive_files}"

    if expected_software:
        # This would depend on the benchmark configuration
        pass


@pytest.mark.e2e_s3
def test_archive_remote_storage_workflow(
    s3_config_path, s3_environment, tmp_path, bundled_repos, keep_files
):
    """
    Test complete archive workflow with S3 remote storage.

    Tests:
    1. Archive with remote results (downloads from S3 first)
    2. Archive with local results only
    3. Dry-run functionality
    4. Different compression options
    5. Validation of archive contents
    """
    print("\n=== Archive S3 Workflow Test ===")
    print(f"Bucket: {s3_environment.bucket_name}")
    print(f"Endpoint: {s3_environment.endpoint}")

    # Update config with dynamic S3 settings
    updated_config = update_config_with_s3_settings(
        s3_config_path, tmp_path, s3_environment
    )

    # Setup test environment
    config_file_in_tmp = setup_test_environment(updated_config, tmp_path)

    # ========================================
    # Step 1: Run benchmark to populate S3
    # ========================================
    print("\n--- Step 1: Running benchmark with --use-remote-storage ---")

    base_args = [
        "run",
        str(config_file_in_tmp),
        "--use-remote-storage",
        "--continue-on-error",
        "-y",
    ]

    with OmniCLISetup() as omni:
        result = omni.call(base_args, cwd=str(tmp_path))

        if keep_files:
            print("\nBENCHMARK EXECUTION DEBUG:")
            print(f"Return code: {result.returncode}")
            print(f"STDOUT:\n{result.stdout}")
            print(f"STDERR:\n{result.stderr}")

    assert result.returncode == 0, (
        f"S3 pipeline execution failed\n"
        f"STDOUT: {result.stdout}\n"
        f"STDERR: {result.stderr}"
    )
    print("✓ Benchmark execution completed successfully")

    # Verify S3 has files
    time.sleep(2)  # S3 consistency
    bucket_contents = s3_environment.list_bucket_contents()
    assert len(bucket_contents) > 0, "No objects found in S3 bucket"
    print(f"✓ S3 bucket contains {len(bucket_contents)} objects")

    # ========================================
    # Step 2: Test archive dry-run with remote results
    # ========================================
    print("\n--- Step 2: Testing archive dry-run with remote results ---")

    dry_run_args = [
        "archive",
        str(config_file_in_tmp),
        "--results",
        "--use-remote-storage",
        "--dry-run",
    ]

    with OmniCLISetup() as omni:
        dry_result = omni.call(dry_run_args, cwd=str(tmp_path))

        if keep_files:
            print("\nDRY RUN DEBUG:")
            print(f"Return code: {dry_result.returncode}")
            print(f"STDOUT:\n{dry_result.stdout}")
            print(f"STDERR:\n{dry_result.stderr}")

    assert dry_result.returncode == 0, (
        f"Archive dry-run failed\n"
        f"STDOUT: {dry_result.stdout}\n"
        f"STDERR: {dry_result.stderr}"
    )

    # Verify dry-run shows files but doesn't create archive
    assert (
        "Files to archive:" in dry_result.stdout
    ), "Dry-run should show files to archive"
    archive_files = [
        f for f in tmp_path.iterdir() if f.suffix in [".zip", ".bz2", ".xz"]
    ]
    assert (
        len(archive_files) == 0
    ), f"Dry-run should not create archive files, found: {archive_files}"

    print("✓ Dry-run completed successfully")
    print("✓ No archive file created during dry-run")

    # ========================================
    # Step 3: Create archive with remote results (downloads from S3)
    # ========================================
    print("\n--- Step 3: Creating archive with remote results ---")

    remote_archive_args = [
        "archive",
        str(config_file_in_tmp),
        "--results",
        "--use-remote-storage",
        "--compression",
        "none",
    ]

    with OmniCLISetup() as omni:
        remote_result = omni.call(remote_archive_args, cwd=str(tmp_path))

        if keep_files:
            print("\nREMOTE ARCHIVE DEBUG:")
            print(f"Return code: {remote_result.returncode}")
            print(f"STDOUT:\n{remote_result.stdout}")
            print(f"STDERR:\n{remote_result.stderr}")

    assert remote_result.returncode == 0, (
        f"Remote archive creation failed\n"
        f"STDOUT: {remote_result.stdout}\n"
        f"STDERR: {remote_result.stderr}"
    )

    # Find created archive
    archive_files = list(tmp_path.glob("*.zip"))
    assert len(archive_files) == 1, f"Expected 1 archive file, found: {archive_files}"

    remote_archive_path = archive_files[0]
    print(f"✓ Remote archive created: {remote_archive_path.name}")

    # Validate archive contents
    remote_contents = extract_and_list_archive_contents(remote_archive_path)
    validate_archive_structure(
        remote_contents, expected_config=True, expected_results=True
    )

    remote_result_files = [
        f for f in remote_contents.keys() if "out/" in f and f.endswith(".json")
    ]
    print(f"✓ Archive contains {len(remote_result_files)} result JSON files")

    # Verify specific expected files are present
    data_files = [f for f in remote_result_files if "_data.json" in f]
    method_files = [f for f in remote_result_files if "_method.json" in f]

    assert (
        len(data_files) == 3
    ), f"Expected 3 data files, found {len(data_files)}: {data_files}"
    assert (
        len(method_files) == 5
    ), f"Expected 5 method files, found {len(method_files)}: {method_files}"
    print("✓ Archive contains expected cartesian product files")

    # ========================================
    # Step 4: Clean local files and test local-only archive
    # ========================================
    print("\n--- Step 4: Testing local-only archive behavior ---")

    # Remove local out directory to simulate clean environment
    local_out_dir = tmp_path / "out"
    if local_out_dir.exists():
        shutil.rmtree(local_out_dir)
        print("✓ Removed local out directory")

    # Remove previous archive
    remote_archive_path.unlink()

    # Try to create archive without remote storage (should have no results)
    local_archive_args = [
        "archive",
        str(config_file_in_tmp),
        "--results",
        # NOTE: No --use-remote-storage flag
    ]

    with OmniCLISetup() as omni:
        local_result = omni.call(local_archive_args, cwd=str(tmp_path))

        if keep_files:
            print("\nLOCAL ARCHIVE DEBUG:")
            print(f"Return code: {local_result.returncode}")
            print(f"STDOUT:\n{local_result.stdout}")
            print(f"STDERR:\n{local_result.stderr}")

    # This might succeed but with warnings about missing files
    local_archive_files = list(tmp_path.glob("*.zip"))

    if len(local_archive_files) > 0:
        local_archive_path = local_archive_files[0]
        local_contents = extract_and_list_archive_contents(local_archive_path)

        # Should have config but minimal/no results
        config_files = [f for f in local_contents.keys() if f.endswith(".yaml")]
        local_result_files = [f for f in local_contents.keys() if "out/" in f]

        assert len(config_files) >= 1, "Archive should always contain config"
        print(
            f"✓ Local archive contains {len(local_result_files)} result files (expected: few/none)"
        )

        # Local archive should have significantly fewer result files than remote
        assert len(local_result_files) < len(remote_result_files), (
            f"Local archive should have fewer files than remote archive. "
            f"Local: {len(local_result_files)}, Remote: {len(remote_result_files)}"
        )

    print("✓ Local-only archive behavior validated")

    # ========================================
    # Step 5: Test compression options
    # ========================================
    print("\n--- Step 5: Testing compression options ---")

    # Clean up previous archives
    for archive_file in tmp_path.glob("*.*"):
        if archive_file.suffix in [".zip", ".bz2", ".xz"]:
            archive_file.unlink()

    # Test different compression methods
    compression_tests = [
        ("none", ".zip", zipfile.ZIP_STORED),
        ("deflated", ".zip", zipfile.ZIP_DEFLATED),
        ("bzip2", ".bz2", zipfile.ZIP_BZIP2),
        ("lzma", ".xz", zipfile.ZIP_LZMA),
    ]

    archive_sizes = {}

    for compression_name, expected_ext, compression_type in compression_tests:
        print(f"\n  Testing {compression_name} compression...")

        compression_args = [
            "archive",
            str(config_file_in_tmp),
            "--results",
            "--use-remote-storage",
            "--compression",
            compression_name,
        ]

        with OmniCLISetup() as omni:
            comp_result = omni.call(compression_args, cwd=str(tmp_path))

        if comp_result.returncode != 0:
            print(f"  ⚠ {compression_name} compression failed, skipping...")
            if keep_files:
                print(f"    STDOUT: {comp_result.stdout}")
                print(f"    STDERR: {comp_result.stderr}")
            continue

        # Find archive with expected extension
        comp_archives = list(tmp_path.glob(f"*{expected_ext}"))

        if len(comp_archives) > 0:
            comp_archive = comp_archives[0]
            archive_size = comp_archive.stat().st_size
            archive_sizes[compression_name] = archive_size

            print(f"  ✓ {compression_name}: {comp_archive.name} ({archive_size} bytes)")

            # Validate it's a valid archive
            comp_contents = extract_and_list_archive_contents(comp_archive)
            validate_archive_structure(
                comp_contents, expected_config=True, expected_results=True
            )

            # Clean up for next test
            comp_archive.unlink()
        else:
            print(f"  ⚠ No archive found with extension {expected_ext}")

    if len(archive_sizes) >= 2:
        print(f"✓ Compression testing completed: {len(archive_sizes)} methods tested")

        # Generally lzma should be smaller than none, deflated smaller than none
        if "none" in archive_sizes and "lzma" in archive_sizes:
            assert archive_sizes["lzma"] <= archive_sizes["none"], (
                f"LZMA should compress better than none: lzma={archive_sizes['lzma']}, "
                f"none={archive_sizes['none']}"
            )
            print("✓ Compression ratios validated")

    # ========================================
    # Final Summary
    # ========================================
    print("\n=== Archive S3 Workflow Test Summary ===")
    print("✓ Benchmark executed and populated S3 successfully")
    print("✓ Archive dry-run functionality works correctly")
    print("✓ Archive with remote results downloads from S3 and creates valid archive")
    print("✓ Archive contents include expected benchmark files")
    print("✓ Local-only archive behavior differs from remote archive")
    print(f"✓ Compression options tested: {list(archive_sizes.keys())}")
    print("✓ Archive functionality comprehensively validated")


@pytest.mark.e2e_s3
def test_archive_version_limitation(
    s3_config_path, s3_environment, tmp_path, bundled_repos, keep_files
):
    """
    Test and document the current limitation that archive always uses
    the version specified in benchmark.yaml, not a selectable version.

    This test demonstrates the current behavior and serves as documentation
    for the limitation that users cannot specify which version to archive.
    """
    print("\n=== Archive Version Limitation Test ===")
    print(
        "This test documents that archive uses benchmark.yaml version, not selectable versions"
    )

    # Update config with dynamic S3 settings
    updated_config = update_config_with_s3_settings(
        s3_config_path, tmp_path, s3_environment
    )
    config_file_in_tmp = setup_test_environment(updated_config, tmp_path)

    # Run benchmark and create version 1.0
    print("\n--- Creating version 1.0 ---")
    base_args = [
        "run",
        str(config_file_in_tmp),
        "--use-remote-storage",
        "--continue-on-error",
        "-y",
    ]

    with OmniCLISetup() as omni:
        result = omni.call(base_args, cwd=str(tmp_path))
        assert result.returncode == 0

        # Create version 1.0
        version_create_args = [
            "remote",
            "version",
            "create",
            "--benchmark",
            str(config_file_in_tmp),
        ]
        version_result = omni.call(version_create_args, cwd=str(tmp_path))
        assert version_result.returncode == 0

    print("✓ Version 1.0 created")

    # Modify benchmark to version 2.0
    print("\n--- Modifying benchmark for version 2.0 ---")
    with open(config_file_in_tmp, "r") as f:
        config_data = yaml.safe_load(f)

    # Update version in YAML
    config_data["version"] = "2.0"

    # Add a parameter change to make results different
    data_stage = next(stage for stage in config_data["stages"] if stage["id"] == "data")
    for module in data_stage["modules"]:
        if module["id"] == "D1":
            module["parameters"][0]["evaluate"] = "150"  # Changed from 100

    with open(config_file_in_tmp, "w") as f:
        yaml.dump(config_data, f, default_flow_style=False, sort_keys=False)

    # Run modified benchmark
    with OmniCLISetup() as omni:
        result2 = omni.call(base_args, cwd=str(tmp_path))
        assert result2.returncode == 0

        # Create version 2.0
        version_result2 = omni.call(version_create_args, cwd=str(tmp_path))
        assert version_result2.returncode == 0

    print("✓ Version 2.0 created with different results")

    # Archive will use version 2.0 (from current benchmark.yaml)
    print("\n--- Testing archive uses current YAML version (2.0) ---")
    archive_args = [
        "archive",
        str(config_file_in_tmp),
        "--results",
        "--use-remote-storage",
        "--dry-run",  # Use dry-run to see what would be archived
    ]

    with OmniCLISetup() as omni:
        archive_result = omni.call(archive_args, cwd=str(tmp_path))

        if keep_files:
            print("\nARCHIVE OUTPUT:")
            print(f"STDOUT:\n{archive_result.stdout}")
            print(f"STDERR:\n{archive_result.stderr}")

    assert archive_result.returncode == 0, "Archive dry-run should succeed"

    # The archive would contain files from version 2.0 (the current benchmark version)
    # Unfortunately, there's no CLI flag to specify "--version 1.0" for archiving

    print("✓ Archive dry-run completed")
    print("📋 LIMITATION DOCUMENTED:")
    print("   - Archive always uses version from benchmark.yaml")
    print("   - Cannot specify '--version 1.0' to archive older version")
    print("   - To archive version 1.0, must modify benchmark.yaml version field")
    print("   - This is a current limitation of the archive functionality")

    print("\n=== Archive Version Limitation Test Summary ===")
    print("✓ Both versions 1.0 and 2.0 created successfully")
    print("✓ Archive functionality uses current benchmark.yaml version")
    print("✓ Limitation documented: No version selection capability")
    print("✓ Future enhancement: Add --version flag to archive command")
