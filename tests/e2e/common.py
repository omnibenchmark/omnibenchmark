"""
Common utilities for omnibenchmark end-to-end tests.

This module provides shared functionality to reduce code duplication
across E2E test files while maintaining testability and flexibility.
"""

import shutil
from pathlib import Path
from typing import List, Optional, Any

from tests.cli.cli_setup import OmniCLISetup
from tests.e2e.validation import validate_pipeline_results, load_expected_results


# Constants
# XXX This is essentially a hack around or own inability to switch off the conda directive even when
# we're not using conda for a benchmark.
DUMMY_CONDA_ENV_CONTENT = """name: omnibenchmark-dummy
channels:
  - conda-forge
dependencies:
  - python=3.12
"""


class E2ETestRunner:
    """
    Encapsulates common E2E test execution patterns.

    This class provides a consistent interface for running omnibenchmark CLI
    commands with proper setup, debugging, and validation.
    """

    def __init__(self, tmp_path: Path, keep_files: bool = False):
        self.tmp_path = tmp_path
        self.keep_files = keep_files
        self.out_dir = tmp_path / "out"

    def setup_test_environment(self, config_path: Path, config_name: str) -> Path:
        """
        Set up the test environment with config file and dummy conda environment.

        Args:
            config_path: Path to the original config file
            config_name: Name for the config file in tmp_path

        Returns:
            Path to the copied config file in tmp_path
        """
        if self.keep_files:
            print("\n=== TEMP PATH FOR INSPECTION ===")
            print(f"Files will be kept at: {self.tmp_path}")
            print("=== END TEMP PATH INFO ===\n")

        # Ensure tmp_path directory exists
        self.tmp_path.mkdir(parents=True, exist_ok=True)

        # Copy config to tmp_path
        config_file_in_tmp = self.tmp_path / config_name
        shutil.copy2(config_path, config_file_in_tmp)

        # Create dummy conda environment file for Snakemake 9.x validation
        # Create in working directory (preferred location)
        conda_file_path = self.tmp_path / "conda_not_provided.yml"
        conda_file_path.write_text(DUMMY_CONDA_ENV_CONTENT)

        # Also create in system temp directory for cross-platform compatibility
        # On Linux CI systems, Snakemake may look in /tmp regardless of working directory
        try:
            import tempfile

            system_tmp_path = Path(tempfile.gettempdir()) / "conda_not_provided.yml"
            system_tmp_path.write_text(DUMMY_CONDA_ENV_CONTENT)
        except (PermissionError, OSError):
            # If we can't write to system temp, continue - working directory version might be sufficient
            pass

        return config_file_in_tmp

    def execute_cli_command(
        self,
        config_file: Path,
        additional_args: Optional[List[str]] = None,
        debug_label: str = "",
    ) -> Any:
        """
        Execute the omnibenchmark CLI command with standard parameters.

        Args:
            config_file: Path to the config file
            additional_args: Additional CLI arguments to include
            debug_label: Label for debug output (e.g., "first run", "second run")

        Returns:
            CLI execution result
        """
        base_args = [
            "run",
            str(config_file),
            "--out-dir",
            str(self.out_dir),
        ]

        if additional_args:
            base_args.extend(additional_args)

        with OmniCLISetup() as omni:
            result = omni.call(base_args, cwd=str(self.tmp_path))

            # Print debug output
            debug_suffix = f" ({debug_label})" if debug_label else ""
            print(f"\n=== CLI EXECUTION DEBUG{debug_suffix} ===")
            print(f"Return code: {result.returncode}")
            print(f"STDOUT:\n{result.stdout}")
            print(f"STDERR:\n{result.stderr}")
            print("=== END DEBUG ===\n")

            # Assert CLI success
            assert result.returncode == 0, (
                f"CLI execution failed with return code {result.returncode}\n"
                f"STDOUT: {result.stdout}\n"
                f"STDERR: {result.stderr}"
            )

            return result

    def debug_output_structure(self) -> None:
        """Print the output directory structure for debugging."""
        if not self.keep_files:
            return

        print("\n=== ACTUAL DIRECTORY STRUCTURE ===")
        print(f"Output directory: {self.out_dir}")
        all_files = list(self.out_dir.rglob("*"))
        for file_path in sorted(all_files):
            if file_path.is_file():
                print(f"FILE: {file_path.relative_to(self.out_dir)}")
            else:
                print(f"DIR:  {file_path.relative_to(self.out_dir)}/")
        print("=== END DIRECTORY STRUCTURE ===\n")

        # Show JSON files specifically
        json_files = list(self.out_dir.rglob("*.json"))
        print(f"Found {len(json_files)} JSON files:")
        for json_file in sorted(json_files):
            print(f"  - {json_file.relative_to(self.out_dir)}")
        print()

    def validate_results(self, test_name: str) -> None:
        """
        Validate pipeline results using the expected results file.

        Args:
            test_name: Name of the test (used to load expected results)
        """
        # Assert output directory exists
        assert self.out_dir.exists(), "Output directory was not created"

        # Debug output structure if requested
        self.debug_output_structure()

        # Load expected results and validate
        expected_results = load_expected_results(test_name)
        validate_pipeline_results(
            self.out_dir, expected_results, verbose=self.keep_files
        )

    def verify_output_file_count(
        self, min_expected: int, file_pattern: str = "*_data.json"
    ) -> None:
        """
        Verify that the expected number of output files were created.

        Args:
            min_expected: Minimum number of files expected
            file_pattern: Glob pattern for files to count (default: "*_data.json")
        """
        output_files = filter_files_excluding_symlinked_dirs(self.out_dir, file_pattern)

        if self.keep_files:
            print(
                f"\nCreated {len(output_files)} total {file_pattern} output files (excluding symlinked dirs):"
            )
            for output_file in output_files:
                print(f"  - {output_file.relative_to(self.out_dir)}")

        assert (
            len(output_files) >= min_expected
        ), f"Expected at least {min_expected} {file_pattern} files, but found {len(output_files)}"


def filter_files_excluding_symlinked_dirs(
    base_dir: Path, file_pattern: str = "*_data.json"
) -> List[Path]:
    """
    Filter files matching a pattern, excluding those in symlinked directories (to avoid double counting)

    Args:
        base_dir: Base directory to search from
        file_pattern: Glob pattern for files to find (default: "*_data.json")

    Returns:
        List of file paths that are not in symlinked directories
    """
    all_files = list(base_dir.rglob(file_pattern))
    actual_files = []
    for file_path in all_files:
        # Check if any parent directory in the path is a symlink
        is_in_symlink = False
        for parent in file_path.parents:
            if parent.is_symlink():
                is_in_symlink = True
                break
            if parent == base_dir:
                break
        if not is_in_symlink:
            actual_files.append(file_path)
    return actual_files


def run_standard_pipeline_test(
    config_path: Path,
    config_filename: str,
    test_name: str,
    tmp_path: Path,
    keep_files: bool,
    min_expected_files: int = 2,
    additional_cli_args: Optional[List[str]] = None,
) -> None:
    """
    Run a standard pipeline test with common setup, execution, and validation.

    This is a high-level helper that encapsulates the most common E2E test pattern.

    Args:
        config_path: Path to the config file
        config_filename: Name for the config file in tmp_path
        test_name: Name of the test (for loading expected results)
        tmp_path: pytest tmp_path fixture
        keep_files: Whether to keep files for debugging
        min_expected_files: Minimum number of JSON files expected
        additional_cli_args: Additional CLI arguments to pass
    """
    runner = E2ETestRunner(tmp_path, keep_files)

    # Setup environment
    config_file_in_tmp = runner.setup_test_environment(config_path, config_filename)

    # Execute CLI with default --continue-on-error flag
    default_args = ["--continue-on-error"]
    if additional_cli_args:
        default_args.extend(additional_cli_args)
    runner.execute_cli_command(config_file_in_tmp, default_args)

    # Validate results
    runner.validate_results(test_name)

    # Verify file count
    runner.verify_output_file_count(min_expected_files)


def compare_pipeline_runs(
    out_dir_first: Path, out_dir_second: Path, file_pattern: str = "*_data.json"
) -> None:
    """
    Compare files from two pipeline runs to ensure idempotency.

    Only compares core data files, excluding performance metrics and other
    non-deterministic files that naturally differ between runs.

    Args:
        out_dir_first: Output directory from first run
        out_dir_second: Output directory from second run
        file_pattern: Glob pattern for files to compare (default: "*_data.json")
    """
    # Get files matching the pattern from both directories
    first_run_files = filter_files_excluding_symlinked_dirs(out_dir_first, file_pattern)
    second_run_files = filter_files_excluding_symlinked_dirs(
        out_dir_second, file_pattern
    )

    # Create dictionaries mapping relative paths to content
    first_files_dict = {}
    for file_path in first_run_files:
        with open(file_path, "r") as f:
            first_files_dict[str(file_path.relative_to(out_dir_first))] = f.read()

    second_files_dict = {}
    for file_path in second_run_files:
        with open(file_path, "r") as f:
            second_files_dict[str(file_path.relative_to(out_dir_second))] = f.read()

    # Compare file sets
    assert (
        set(first_files_dict.keys()) == set(second_files_dict.keys())
    ), f"File sets differ between runs: {set(first_files_dict.keys())} vs {set(second_files_dict.keys())}"

    # Compare file contents
    for file_key in first_files_dict:
        assert (
            first_files_dict[file_key] == second_files_dict[file_key]
        ), f"File {file_key} changed between runs"


def run_idempotent_pipeline_test(
    config_path: Path,
    config_filename: str,
    test_name: str,
    tmp_path: Path,
    keep_files: bool,
    min_expected_files: int = 2,
    additional_cli_args: Optional[List[str]] = None,
) -> None:
    """
    Run an idempotent pipeline test by running the same pipeline twice in the same directory.

    We run the pipeline twice and verify that outputs don't change between runs.

    Args:
        config_path: Path to the config file
        config_filename: Name for the config file in tmp_path
        test_name: Name of the test (for loading expected results)
        tmp_path: pytest tmp_path fixture
        keep_files: Whether to keep files for debugging
        min_expected_files: Minimum number of JSON files expected
        additional_cli_args: Additional CLI arguments to pass
    """
    runner = E2ETestRunner(tmp_path, keep_files)

    # Setup environment
    config_file_in_tmp = runner.setup_test_environment(config_path, config_filename)

    # Prepare CLI args
    default_args = ["--continue-on-error"]
    if additional_cli_args:
        default_args.extend(additional_cli_args)

    # First run
    runner.execute_cli_command(
        config_file_in_tmp, default_args, debug_label="first run"
    )
    runner.validate_results(test_name)
    runner.verify_output_file_count(min_expected_files)

    # Get file checksums after first run
    first_run_checksums = _get_directory_checksums(runner.out_dir)

    # Second run (should be idempotent)
    runner.execute_cli_command(
        config_file_in_tmp, default_args, debug_label="second run"
    )

    # Get file checksums after second run
    second_run_checksums = _get_directory_checksums(runner.out_dir)

    # Compare checksums to ensure idempotency
    assert (
        first_run_checksums == second_run_checksums
    ), "Pipeline is not idempotent - files changed between runs"


def _get_directory_checksums(directory: Path) -> dict:
    """
    Get MD5 checksums for all files in a directory recursively.

    Returns:
        Dictionary mapping relative file paths to their MD5 checksums
    """
    import hashlib

    checksums = {}
    for file_path in directory.rglob("*"):
        if file_path.is_file() and not file_path.is_symlink():
            # Calculate MD5 checksum
            hasher = hashlib.md5()
            with open(file_path, "rb") as f:
                for chunk in iter(lambda: f.read(4096), b""):
                    hasher.update(chunk)

            # Store with relative path as key
            rel_path = str(file_path.relative_to(directory))
            checksums[rel_path] = hasher.hexdigest()

    return checksums


def extract_test_name_from_config(config_filename: str) -> str:
    """
    Extract test name from config filename.

    Args:
        config_filename: Name of config file (e.g., "00_data_modules.yaml")

    Returns:
        Test name for loading expected results (e.g., "00_data_modules")
    """
    return config_filename.replace(".yaml", "").replace(".yml", "")
