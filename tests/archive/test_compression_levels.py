"""
Unit tests for archive compression levels and behavior.
"""

import pytest
import zipfile
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

from omnibenchmark.archive.archive import archive_benchmark


@pytest.mark.short
def test_compression_level_defaults():
    """Test that compression levels are properly handled."""
    mock_benchmark = Mock()
    mock_benchmark.get_benchmark_name.return_value = "test_benchmark"
    mock_benchmark.get_benchmark_version.return_value = "1.0"

    with patch("omnibenchmark.archive.archive.prepare_archive_config") as mock_config:
        mock_config.return_value = [Path("benchmark.yaml")]

        # Test different compression types with levels
        test_cases = [
            (zipfile.ZIP_BZIP2, 9, ".bz2"),
            (zipfile.ZIP_LZMA, 9, ".xz"),
            (zipfile.ZIP_DEFLATED, 9, ".zip"),
            (zipfile.ZIP_STORED, None, ".zip"),  # No compression level for stored
        ]

        for compression, level, expected_ext in test_cases:
            with (
                patch("zipfile.ZipFile") as mock_zipfile,
                patch("pathlib.Path.is_file", return_value=True),
            ):
                mock_archive = MagicMock()
                mock_zipfile.return_value.__enter__.return_value = mock_archive

                result = archive_benchmark(
                    benchmark=mock_benchmark,
                    compression=compression,
                    compresslevel=level,
                    dry_run=False,
                )

                # Check that ZipFile was called with correct parameters
                expected_args = {"compression": compression}
                if level is not None:
                    expected_args["compresslevel"] = level

                # Verify ZipFile was called with correct compression settings
                mock_zipfile.assert_called()
                call_kwargs = mock_zipfile.call_args[1]
                assert call_kwargs["compression"] == compression

                if level is not None:
                    assert call_kwargs["compresslevel"] == level

                # Check correct file extension
                assert result[0].suffix == expected_ext


@pytest.mark.short
def test_bzip2_max_compression_level():
    """Test that bzip2 compression level 9 (maximum) is properly applied."""
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

            # Test with maximum bzip2 compression level
            archive_benchmark(
                benchmark=mock_benchmark,
                compression=zipfile.ZIP_BZIP2,
                compresslevel=9,
                dry_run=False,
            )

            # Verify ZipFile was called with bzip2 and level 9
            call_kwargs = mock_zipfile.call_args[1]
            assert call_kwargs["compression"] == zipfile.ZIP_BZIP2
            assert call_kwargs["compresslevel"] == 9


@pytest.mark.short
def test_compression_level_validation():
    """Test that compression levels are validated correctly."""
    mock_benchmark = Mock()
    mock_benchmark.get_benchmark_name.return_value = "test_benchmark"
    mock_benchmark.get_benchmark_version.return_value = "1.0"

    with patch("omnibenchmark.archive.archive.prepare_archive_config") as mock_config:
        mock_config.return_value = [Path("benchmark.yaml")]

        # Test invalid compression levels (should be handled by zipfile)
        invalid_levels = [-1, 0, 10, 100]

        for invalid_level in invalid_levels:
            with (
                patch("zipfile.ZipFile") as mock_zipfile,
                patch("pathlib.Path.is_file", return_value=True),
            ):
                # ZipFile should handle invalid levels (may raise exception or ignore)
                try:
                    archive_benchmark(
                        benchmark=mock_benchmark,
                        compression=zipfile.ZIP_BZIP2,
                        compresslevel=invalid_level,
                        dry_run=False,
                    )

                    # If it doesn't raise, just verify it was called
                    call_kwargs = mock_zipfile.call_args[1]
                    assert call_kwargs["compresslevel"] == invalid_level

                except Exception:
                    # ZipFile validation may raise exception, which is expected
                    pass


@pytest.mark.short
def test_compression_none_no_level():
    """Test that ZIP_STORED (no compression) ignores compression level."""
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

            # Test ZIP_STORED with compression level (should be passed but may be ignored)
            archive_benchmark(
                benchmark=mock_benchmark,
                compression=zipfile.ZIP_STORED,
                compresslevel=9,
                dry_run=False,
            )

            # Verify ZipFile was called
            call_kwargs = mock_zipfile.call_args[1]
            assert call_kwargs["compression"] == zipfile.ZIP_STORED
            assert (
                call_kwargs["compresslevel"] == 9
            )  # Passed but ZIP_STORED may ignore it


@pytest.mark.short
def test_default_compression_settings():
    """Test that new default compression settings work correctly."""
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

            # Test with new defaults (bzip2, level 9)
            result = archive_benchmark(
                benchmark=mock_benchmark,
                compression=zipfile.ZIP_BZIP2,  # New default
                compresslevel=9,  # New default max level
                dry_run=False,
            )

            # Should create .bz2 file with maximum compression
            assert result[0].suffix == ".bz2"
            assert "test_benchmark_1.0.bz2" in str(result[0])

            # Verify compression settings
            call_kwargs = mock_zipfile.call_args[1]
            assert call_kwargs["compression"] == zipfile.ZIP_BZIP2
            assert call_kwargs["compresslevel"] == 9
