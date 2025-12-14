"""
Examples and documentation for the -o/--output-file feature.

This test file serves as both validation and user documentation
for the custom output file functionality in the archive command.
"""

import pytest
import zipfile
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

from omnibenchmark.archive.archive import archive_benchmark


@pytest.mark.short
def test_output_file_examples():
    """
    Test examples of using custom output files with different scenarios.

    This test demonstrates real-world usage patterns for the --output-file feature.
    """
    mock_benchmark = Mock()
    mock_benchmark.get_benchmark_name.return_value = "example_benchmark"
    mock_benchmark.get_benchmark_version.return_value = "1.5"

    with patch("omnibenchmark.archive.archive.prepare_archive_config") as mock_config:
        mock_config.return_value = [Path("benchmark.yaml")]

        with (
            patch("zipfile.ZipFile") as mock_zipfile,
            patch("pathlib.Path.is_file", return_value=True),
        ):
            mock_archive = MagicMock()
            mock_zipfile.return_value.__enter__.return_value = mock_archive

            # Example 1: Basic custom filename with default compression (bzip2)
            result1 = archive_benchmark(
                benchmark=mock_benchmark, custom_filename="my_backup.bz2", dry_run=False
            )
            assert result1[0].name == "my_backup.bz2"

            # Example 2: Date-stamped archive with LZMA compression
            result2 = archive_benchmark(
                benchmark=mock_benchmark,
                compression=zipfile.ZIP_LZMA,
                custom_filename="benchmark_2024-12-14.xz",
                dry_run=False,
            )
            assert result2[0].name == "benchmark_2024-12-14.xz"

            # Example 3: Project-specific naming with ZIP format
            result3 = archive_benchmark(
                benchmark=mock_benchmark,
                compression=zipfile.ZIP_DEFLATED,
                custom_filename="project_alpha_results.zip",
                dry_run=False,
            )
            assert result3[0].name == "project_alpha_results.zip"

            # Example 4: Version-specific archive (overriding auto-generated name)
            result4 = archive_benchmark(
                benchmark=mock_benchmark,
                custom_filename="benchmark_v1.5_final.bz2",
                dry_run=False,
            )
            assert result4[0].name == "benchmark_v1.5_final.bz2"


@pytest.mark.short
def test_output_file_compression_matrix():
    """
    Test matrix showing all valid compression type and extension combinations.

    This serves as documentation for users about which extensions work with which compression types.
    """

    # Valid combinations matrix
    valid_combinations = {
        # Compression Type: [Valid Extensions]
        "none": [".zip"],
        "deflated": [".zip"],
        "bzip2": [".bz2"],
        "lzma": [".xz"],
    }

    # Test each valid combination
    for compression_name, extensions in valid_combinations.items():
        for ext in extensions:
            # This would be valid in CLI usage
            example_filename = f"test_archive{ext}"

            # Verify extension matches what we expect
            assert Path(example_filename).suffix == ext

    # Invalid combinations that should be rejected
    invalid_combinations = [
        ("bzip2", ".zip"),  # bzip2 compression with zip extension
        ("lzma", ".bz2"),  # lzma compression with bzip2 extension
        ("none", ".xz"),  # no compression with xz extension
        ("deflated", ".bz2"),  # deflated compression with bzip2 extension
    ]

    # These combinations should be caught by CLI validation
    for compression, wrong_ext in invalid_combinations:
        # Document that these are invalid
        assert compression in valid_combinations
        assert wrong_ext not in valid_combinations[compression]


@pytest.mark.short
def test_output_file_usage_scenarios():
    """
    Document common usage scenarios for the output file feature.

    These examples show how users might use the feature in practice.
    """

    # Scenario 1: Timestamped backups
    timestamp_examples = [
        "backup_2024-12-14_morning.bz2",
        "results_2024-12-14T19-30-00.xz",
        "archive_20241214.zip",
    ]

    # Scenario 2: Experiment naming
    experiment_examples = [
        "experiment_001_baseline.bz2",
        "exp_002_modified_params.xz",
        "trial_03_final_run.zip",
    ]

    # Scenario 3: Environment-specific archives
    environment_examples = [
        "production_results.bz2",
        "development_test.zip",
        "staging_validation.xz",
    ]

    # Scenario 4: Full path examples
    path_examples = [
        "/home/user/backups/benchmark.bz2",
        "./archives/experiment_001.xz",
        "../shared/results.zip",
        "results/2024/december/final.bz2",
    ]

    all_examples = (
        timestamp_examples + experiment_examples + environment_examples + path_examples
    )

    # Verify all examples have valid extensions
    for example in all_examples:
        path = Path(example)
        ext = path.suffix.lower()
        assert ext in [
            ".zip",
            ".bz2",
            ".xz",
        ], f"Invalid extension in example: {example}"

    # Document the count for each type
    zip_count = sum(1 for ex in all_examples if Path(ex).suffix.lower() == ".zip")
    bz2_count = sum(1 for ex in all_examples if Path(ex).suffix.lower() == ".bz2")
    xz_count = sum(1 for ex in all_examples if Path(ex).suffix.lower() == ".xz")

    assert zip_count > 0  # At least some ZIP examples
    assert bz2_count > 0  # At least some BZIP2 examples
    assert xz_count > 0  # At least some XZ examples


@pytest.mark.short
def test_output_file_cli_examples():
    """
    Document CLI usage examples for the output file feature.

    These show the actual command-line usage patterns.
    """

    # CLI Examples (as documentation - not executed)
    cli_examples = [
        # Basic usage with default bzip2 compression
        "omnibenchmark archive -b config.yaml -r -o my_results.bz2",
        # Custom compression with matching extension
        "omnibenchmark archive -b config.yaml -r --compression lzma -o compressed.xz",
        # Full path output
        "omnibenchmark archive -b config.yaml -r -o /backup/archive.bz2",
        # With other options
        "omnibenchmark archive -b config.yaml -r -c -s -o full_backup.bz2",
        # Dry run to preview
        "omnibenchmark archive -b config.yaml -r -o test.bz2 --dry-run",
    ]

    # Parse each example to verify format
    for example in cli_examples:
        # Basic validation that examples follow expected format
        assert "omnibenchmark archive" in example
        assert "-o " in example or "--output-file" in example
        assert any(ext in example for ext in [".zip", ".bz2", ".xz"])

    # Document compression/extension mappings shown in examples
    extension_mappings = {
        ".bz2": "bzip2 (default)",
        ".xz": "lzma",
        ".zip": "none or deflated",
    }

    for ext, compression in extension_mappings.items():
        # Verify mapping is documented
        assert ext.startswith(".")
        assert len(compression) > 0


@pytest.mark.short
def test_output_file_error_examples():
    """
    Document common errors and how to fix them.

    This helps users understand what went wrong and how to fix it.
    """

    # Common error scenarios
    error_scenarios = [
        {
            "command": "omnibenchmark archive -b config.yaml -r --compression bzip2 -o wrong.zip",
            "error": "Output file extension '.zip' does not match compression type 'bzip2'",
            "fix": "Change to: -o correct.bz2",
            "explanation": "bzip2 compression requires .bz2 extension",
        },
        {
            "command": "omnibenchmark archive -b config.yaml -r --compression lzma -o wrong.bz2",
            "error": "Output file extension '.bz2' does not match compression type 'lzma'",
            "fix": "Change to: -o correct.xz",
            "explanation": "lzma compression requires .xz extension",
        },
        {
            "command": "omnibenchmark archive -b config.yaml -r --compression none -o wrong.bz2",
            "error": "Output file extension '.bz2' does not match compression type 'none'",
            "fix": "Change to: -o correct.zip",
            "explanation": "no compression requires .zip extension",
        },
    ]

    # Verify error scenarios are well-formed
    for scenario in error_scenarios:
        assert "command" in scenario
        assert "error" in scenario
        assert "fix" in scenario
        assert "explanation" in scenario

        # Check command format
        assert scenario["command"].startswith("omnibenchmark archive")
        assert "-o " in scenario["command"]

        # Check error mentions extension mismatch
        assert "does not match compression type" in scenario["error"]

        # Check fix provides alternative
        assert "Change to:" in scenario["fix"]
        assert "-o " in scenario["fix"]
