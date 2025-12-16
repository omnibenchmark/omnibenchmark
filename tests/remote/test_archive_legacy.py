"""
Tests for the legacy archive module wrapper functions.

These tests verify that the legacy archive API continues to work
and properly delegates to the new archive module.
All tests are marked as 'short' since they test wrapper functionality.
"""

import pytest
import zipfile
from pathlib import Path
from unittest.mock import Mock, patch

from omnibenchmark.remote.archive import archive_version


@pytest.mark.short
class TestArchiveLegacy:
    """Test the legacy archive wrapper functions."""

    @patch("omnibenchmark.remote.archive.archive_benchmark")
    def test_archive_version_delegates_to_new_function(self, mock_archive):
        """Test that archive_version properly delegates to archive_benchmark."""
        # Setup
        mock_benchmark = Mock()
        mock_outdir = Path("/test/output")
        mock_archive.return_value = [Path("/test/archive.zip")]

        # Call the legacy function
        result = archive_version(
            benchmark=mock_benchmark,
            outdir=mock_outdir,
            config=True,
            code=False,
            software=True,
            results=False,
            results_dir="custom_out",
            compression=zipfile.ZIP_DEFLATED,
            compresslevel=6,
            dry_run=False,
            remote_storage=True,
        )

        # Verify delegation
        mock_archive.assert_called_once_with(
            benchmark=mock_benchmark,
            outdir=mock_outdir,
            config=True,
            code=False,
            software=True,
            results=False,
            results_dir="custom_out",
            compression=zipfile.ZIP_DEFLATED,
            compresslevel=6,
            dry_run=False,
            remote_storage=True,
        )

        # Verify return value handling
        assert result == Path("/test/archive.zip")

    @patch("omnibenchmark.remote.archive.archive_benchmark")
    def test_archive_version_default_parameters(self, mock_archive):
        """Test archive_version with default parameters."""
        mock_benchmark = Mock()
        mock_archive.return_value = [Path("/test/default.zip")]

        result = archive_version(benchmark=mock_benchmark)

        # Verify call with defaults
        mock_archive.assert_called_once_with(
            benchmark=mock_benchmark,
            outdir=Path(),
            config=True,
            code=False,
            software=False,
            results=False,
            results_dir="out",
            compression=zipfile.ZIP_STORED,
            compresslevel=None,
            dry_run=False,
            remote_storage=True,
        )

        assert result == Path("/test/default.zip")

    @patch("omnibenchmark.remote.archive.archive_benchmark")
    def test_archive_version_dry_run_returns_list(self, mock_archive):
        """Test that dry_run=True returns the list result directly."""
        mock_benchmark = Mock()
        expected_result = [Path("/test/file1.zip"), Path("/test/file2.zip")]
        mock_archive.return_value = expected_result

        result = archive_version(benchmark=mock_benchmark, dry_run=True)

        # Verify dry_run returns the full list
        assert result == expected_result
        mock_archive.assert_called_once()

    @patch("omnibenchmark.remote.archive.archive_benchmark")
    def test_archive_version_empty_result_returns_none(self, mock_archive):
        """Test that empty result returns None."""
        mock_benchmark = Mock()
        mock_archive.return_value = []

        result = archive_version(benchmark=mock_benchmark, dry_run=False)

        assert result is None

    @patch("omnibenchmark.remote.archive.archive_benchmark")
    def test_archive_version_none_result_returns_none(self, mock_archive):
        """Test that None result returns None."""
        mock_benchmark = Mock()
        mock_archive.return_value = None

        result = archive_version(benchmark=mock_benchmark, dry_run=False)

        assert result is None

    @patch("omnibenchmark.remote.archive.archive_benchmark")
    def test_archive_version_multiple_results_returns_first(self, mock_archive):
        """Test that multiple results returns the first one."""
        mock_benchmark = Mock()
        expected_files = [Path("/test/first.zip"), Path("/test/second.zip")]
        mock_archive.return_value = expected_files

        result = archive_version(benchmark=mock_benchmark, dry_run=False)

        # Should return only the first result
        assert result == Path("/test/first.zip")

    @patch("omnibenchmark.remote.archive.archive_benchmark")
    def test_archive_version_passes_all_compression_options(self, mock_archive):
        """Test that all compression-related options are passed through."""
        mock_benchmark = Mock()
        mock_archive.return_value = [Path("/test/compressed.zip")]

        archive_version(
            benchmark=mock_benchmark, compression=zipfile.ZIP_LZMA, compresslevel=9
        )

        # Verify compression options are passed
        call_args = mock_archive.call_args
        assert call_args[1]["compression"] == zipfile.ZIP_LZMA
        assert call_args[1]["compresslevel"] == 9

    @patch("omnibenchmark.remote.archive.archive_benchmark")
    def test_archive_version_passes_all_boolean_flags(self, mock_archive):
        """Test that all boolean flags are passed through correctly."""
        mock_benchmark = Mock()
        mock_archive.return_value = [Path("/test/flags.zip")]

        archive_version(
            benchmark=mock_benchmark,
            config=False,
            code=True,
            software=True,
            results=True,
            dry_run=True,
            remote_storage=False,
        )

        # Verify all boolean flags
        call_args = mock_archive.call_args
        assert call_args[1]["config"] is False
        assert call_args[1]["code"] is True
        assert call_args[1]["software"] is True
        assert call_args[1]["results"] is True
        assert call_args[1]["dry_run"] is True
        assert call_args[1]["remote_storage"] is False

    @patch("omnibenchmark.remote.archive.archive_benchmark")
    def test_archive_version_preserves_path_types(self, mock_archive):
        """Test that Path objects are preserved as Path objects."""
        mock_benchmark = Mock()
        test_outdir = Path("/custom/output/directory")
        mock_archive.return_value = [Path("/test/path.zip")]

        archive_version(
            benchmark=mock_benchmark, outdir=test_outdir, results_dir="custom_results"
        )

        call_args = mock_archive.call_args
        assert isinstance(call_args[1]["outdir"], Path)
        assert call_args[1]["outdir"] == test_outdir
        assert call_args[1]["results_dir"] == "custom_results"

    def test_archive_version_imports_available(self):
        """Test that all required imports are available."""
        # Test that the function can be imported
        from omnibenchmark.remote.archive import archive_version

        assert callable(archive_version)

        # Test that zipfile constants are available
        assert hasattr(zipfile, "ZIP_STORED")
        assert hasattr(zipfile, "ZIP_DEFLATED")

        # Test that Path is available
        assert Path is not None

    @patch("omnibenchmark.remote.archive.archive_benchmark")
    def test_archive_version_exception_propagation(self, mock_archive):
        """Test that exceptions from archive_benchmark are propagated."""
        mock_benchmark = Mock()
        mock_archive.side_effect = ValueError("Test error")

        with pytest.raises(ValueError, match="Test error"):
            archive_version(benchmark=mock_benchmark)
