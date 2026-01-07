"""
Unit tests for archive custom output file functionality.
"""

import pytest
import zipfile
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

from omnibenchmark.archive.archive import archive_benchmark
from omnibenchmark.cli.archive import archive


@pytest.mark.short
def test_custom_output_file_basic():
    """Test that custom output filename is used correctly."""
    mock_benchmark = Mock()
    mock_benchmark.get_benchmark_name.return_value = "test_benchmark"
    mock_benchmark.get_benchmark_version.return_value = "1.0"

    with patch("omnibenchmark.archive.archive.prepare_archive_config") as mock_config:
        mock_config.return_value = [Path("benchmark.yaml")]

        with (
            patch("zipfile.ZipFile") as mock_zipfile,
            patch("pathlib.Path.is_file", return_value=True),
        ):
            mock_archive = MagicMock()
            mock_zipfile.return_value.__enter__.return_value = mock_archive

            # Test with custom filename
            custom_name = "my_custom_archive.bz2"
            result = archive_benchmark(
                benchmark=mock_benchmark,
                outdir=Path("/tmp"),
                compression=zipfile.ZIP_BZIP2,
                custom_filename=custom_name,
                dry_run=False,
            )

            # Should use custom filename instead of default
            assert result[0] == Path("/tmp") / custom_name
            assert str(result[0]).endswith("my_custom_archive.bz2")


@pytest.mark.short
def test_custom_output_file_preserves_extension():
    """Test that custom filename preserves the provided extension."""
    mock_benchmark = Mock()
    mock_benchmark.get_benchmark_name.return_value = "test_benchmark"
    mock_benchmark.get_benchmark_version.return_value = "1.0"

    with patch("omnibenchmark.archive.archive.prepare_archive_config") as mock_config:
        mock_config.return_value = [Path("benchmark.yaml")]

        with (
            patch("zipfile.ZipFile") as mock_zipfile,
            patch("pathlib.Path.is_file", return_value=True),
        ):
            mock_archive = MagicMock()
            mock_zipfile.return_value.__enter__.return_value = mock_archive

            # Test with different extensions
            test_cases = [
                ("archive.zip", zipfile.ZIP_STORED),
                ("backup.bz2", zipfile.ZIP_BZIP2),
                ("compressed.xz", zipfile.ZIP_LZMA),
                ("my_data.zip", zipfile.ZIP_DEFLATED),
            ]

            for custom_name, compression in test_cases:
                result = archive_benchmark(
                    benchmark=mock_benchmark,
                    compression=compression,
                    custom_filename=custom_name,
                    dry_run=False,
                )

                # Should preserve exact filename
                assert result[0].name == custom_name


@pytest.mark.short
def test_cli_output_file_validation_bzip2():
    """Test CLI validation for bzip2 compression and .bz2 extension."""
    from click.testing import CliRunner

    runner = CliRunner()

    # Create a test benchmark file
    with runner.isolated_filesystem():
        with open("test.yaml", "w") as f:
            f.write("""
id: test
description: test
version: "1.0"
benchmarker: "Test"
benchmark_yaml_spec: "0.3"
software_backend: host
software_environments:
  host:
    description: "Host environment"
stages:
  - id: data
    modules:
      - id: D1
        software_environment: "host"
        repository:
          url: test
          commit: abc123
        parameters:
          - evaluate: "100"
    outputs:
      - id: out1
        description: test output
        format: json
        path: "{dataset}_output.json"
""")

        # Test valid bzip2 extension - should not raise error
        result = runner.invoke(
            archive,
            [
                # removed --benchmark option - now positional
                "test.yaml",  # benchmark is now positional
                "--results",
                "--compression",
                "bzip2",
                "--output-file",
                "my_archive.bz2",
                "--dry-run",
            ],
        )

        # Should succeed (exit code 0) for valid extension
        assert result.exit_code == 0


@pytest.mark.short
def test_cli_output_file_validation_mismatch():
    """Test CLI validation catches extension/compression mismatches."""
    from click.testing import CliRunner

    runner = CliRunner()

    # Test mismatched extensions
    invalid_combinations = [
        ("bzip2", "archive.zip"),  # bzip2 compression with .zip extension
        ("lzma", "archive.bz2"),  # lzma compression with .bz2 extension
        ("none", "archive.xz"),  # no compression with .xz extension
        ("deflated", "archive.bz2"),  # deflated compression with .bz2 extension
    ]

    with runner.isolated_filesystem():
        with open("test.yaml", "w") as f:
            f.write("""
id: test
description: test
version: "1.0"
benchmarker: "Test"
benchmark_yaml_spec: "0.3"
software_backend: host
software_environments:
  host:
    description: "Host environment"
stages:
  - id: data
    modules:
      - id: D1
        software_environment: "host"
        repository:
          url: test
          commit: abc123
        parameters:
          - evaluate: "100"
    outputs:
      - id: out1
        description: test output
        format: json
        path: "{dataset}_output.json"
""")

        for compression, output_file in invalid_combinations:
            result = runner.invoke(
                archive,
                [
                    # removed --benchmark option - now positional
                    "test.yaml",  # benchmark is now positional
                    "--results",
                    "--compression",
                    compression,
                    "--output-file",
                    output_file,
                    "--dry-run",
                ],
            )

            # Should fail with exit code 2 (usage error)
            assert result.exit_code == 2
            assert "does not match compression type" in result.output


@pytest.mark.short
def test_cli_output_file_validation_all_valid_combinations():
    """Test all valid compression/extension combinations."""
    from click.testing import CliRunner

    runner = CliRunner()

    valid_combinations = [
        ("none", "archive.zip"),
        ("deflated", "archive.zip"),
        ("bzip2", "archive.bz2"),
        ("lzma", "archive.xz"),
    ]

    with runner.isolated_filesystem():
        with open("test.yaml", "w") as f:
            f.write("""
id: test
description: test
version: "1.0"
benchmarker: "Test"
benchmark_yaml_spec: "0.3"
software_backend: host
software_environments:
  host:
    description: "Host environment"
stages:
  - id: data
    modules:
      - id: D1
        software_environment: "host"
        repository:
          url: test
          commit: abc123
        parameters:
          - evaluate: "100"
    outputs:
      - id: out1
        description: test output
        format: json
        path: "{dataset}_output.json"
""")

        for compression, output_file in valid_combinations:
            result = runner.invoke(
                archive,
                [
                    # removed --benchmark option - now positional
                    "test.yaml",  # benchmark is now positional
                    "--results",
                    "--compression",
                    compression,
                    "--output-file",
                    output_file,
                    "--dry-run",
                ],
            )

            # Should succeed for valid combinations
            assert (
                result.exit_code == 0
            ), f"Valid combination {compression}/{output_file} should succeed"


@pytest.mark.short
def test_output_file_with_directory_path():
    """Test that output file works with full directory paths."""
    mock_benchmark = Mock()
    mock_benchmark.get_benchmark_name.return_value = "test_benchmark"
    mock_benchmark.get_benchmark_version.return_value = "1.0"

    with patch("omnibenchmark.archive.archive.prepare_archive_config") as mock_config:
        mock_config.return_value = [Path("benchmark.yaml")]

        with (
            patch("zipfile.ZipFile") as mock_zipfile,
            patch("pathlib.Path.is_file", return_value=True),
        ):
            mock_archive = MagicMock()
            mock_zipfile.return_value.__enter__.return_value = mock_archive

            # Test with full path
            result = archive_benchmark(
                benchmark=mock_benchmark,
                outdir=Path("/home/user/archives"),
                compression=zipfile.ZIP_BZIP2,
                custom_filename="my_backup.bz2",
                dry_run=False,
            )

            # Should use the full path correctly
            assert result[0] == Path("/home/user/archives/my_backup.bz2")


@pytest.mark.short
def test_default_filename_when_no_custom_file():
    """Test that default filename generation works when no custom file is specified."""
    mock_benchmark = Mock()
    mock_benchmark.get_benchmark_name.return_value = "my_benchmark"
    mock_benchmark.get_benchmark_version.return_value = "2.5"

    with patch("omnibenchmark.archive.archive.prepare_archive_config") as mock_config:
        mock_config.return_value = [Path("benchmark.yaml")]

        with (
            patch("zipfile.ZipFile") as mock_zipfile,
            patch("pathlib.Path.is_file", return_value=True),
        ):
            mock_archive = MagicMock()
            mock_zipfile.return_value.__enter__.return_value = mock_archive

            # Test without custom filename (should use default)
            result = archive_benchmark(
                benchmark=mock_benchmark,
                compression=zipfile.ZIP_BZIP2,
                custom_filename=None,  # No custom filename
                dry_run=False,
            )

            # Should use default naming pattern
            expected_name = "my_benchmark_2.5.bz2"
            assert result[0].name == expected_name


@pytest.mark.short
def test_cli_output_file_case_insensitive_extensions():
    """Test that extension validation is case insensitive."""
    from click.testing import CliRunner

    runner = CliRunner()

    # Test various case combinations
    valid_cases = [
        ("bzip2", "archive.BZ2"),
        ("bzip2", "archive.bz2"),
        ("lzma", "archive.XZ"),
        ("lzma", "archive.xz"),
        ("none", "archive.ZIP"),
        ("deflated", "archive.zip"),
    ]

    with runner.isolated_filesystem():
        with open("test.yaml", "w") as f:
            f.write("""
id: test
description: test
version: "1.0"
benchmarker: "Test"
benchmark_yaml_spec: "0.3"
software_backend: host
software_environments:
  host:
    description: "Host environment"
stages:
  - id: data
    modules:
      - id: D1
        software_environment: "host"
        repository:
          url: test
          commit: abc123
        parameters:
          - evaluate: "100"
    outputs:
      - id: out1
        description: test output
        format: json
        path: "{dataset}_output.json"
""")

        for compression, output_file in valid_cases:
            result = runner.invoke(
                archive,
                [
                    # removed --benchmark option - now positional
                    "test.yaml",  # benchmark is now positional
                    "--results",
                    "--compression",
                    compression,
                    "--output-file",
                    output_file,
                    "--dry-run",
                ],
            )

            # Should succeed regardless of case
            assert (
                result.exit_code == 0
            ), f"Case variation {compression}/{output_file} should be valid"


@pytest.mark.short
def test_cli_output_file_error_message_format():
    """Test that error messages for extension mismatch are helpful."""
    from click.testing import CliRunner

    runner = CliRunner()

    with runner.isolated_filesystem():
        with open("test.yaml", "w") as f:
            f.write("""
id: test
description: test
version: "1.0"
benchmarker: "Test"
benchmark_yaml_spec: "0.3"
software_backend: host
software_environments:
  host:
    description: "Host environment"
stages:
  - id: data
    modules:
      - id: D1
        software_environment: "host"
        repository:
          url: test
          commit: abc123
        parameters:
          - evaluate: "100"
    outputs:
      - id: out1
        description: test output
        format: json
        path: "{dataset}_output.json"
""")

        result = runner.invoke(
            archive,
            [
                # removed --benchmark option - now positional
                "test.yaml",  # benchmark is now positional
                "--results",
                "--compression",
                "bzip2",
                "--output-file",
                "wrong_extension.zip",
                "--dry-run",
            ],
        )

        # Should fail with usage error
        assert result.exit_code == 2

        # Check that error message contains helpful information
        assert "Output file extension" in result.output
        assert ".zip" in result.output
        assert "bzip2" in result.output
        assert "Expected for 'bzip2': .bz2" in result.output
        assert "Valid combinations:" in result.output
        assert "none/deflated → .zip" in result.output
        assert "bzip2 → .bz2" in result.output
        assert "lzma → .xz" in result.output
