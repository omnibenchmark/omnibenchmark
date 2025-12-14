"""
Unit tests for omnibenchmark archive functionality.
"""

import pytest
import zipfile
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

from omnibenchmark.archive.archive import archive_benchmark
from omnibenchmark.archive.components import (
    prepare_archive_config,
    prepare_archive_code,
    prepare_archive_results_local,
)


@pytest.mark.short
def test_archive_benchmark_dry_run():
    """Test archive_benchmark in dry-run mode returns file list."""
    # Mock benchmark object
    mock_benchmark = Mock()
    mock_benchmark.get_benchmark_name.return_value = "test_benchmark"
    mock_benchmark.get_benchmark_version.return_value = "1.0"
    mock_benchmark.get_definition_file.return_value = Path("benchmark.yaml")

    # Mock the component preparation functions
    with (
        patch("omnibenchmark.archive.archive.prepare_archive_config") as mock_config,
        patch("omnibenchmark.archive.archive.prepare_archive_code") as mock_code,
        patch(
            "omnibenchmark.archive.archive.prepare_archive_software"
        ) as mock_software,
        patch("omnibenchmark.archive.archive.prepare_archive_results") as mock_results,
    ):
        mock_config.return_value = [Path("benchmark.yaml")]
        mock_code.return_value = [Path("code.py")]
        mock_software.return_value = [Path("environment.yml")]
        mock_results.return_value = [Path("results.json")]

        result = archive_benchmark(
            benchmark=mock_benchmark,
            config=True,
            code=True,
            software=True,
            results=True,
            dry_run=True,
        )

        # Should return list of files
        assert isinstance(result, list)
        assert len(result) == 4
        assert Path("benchmark.yaml") in result
        assert Path("code.py") in result
        assert Path("environment.yml") in result
        assert Path("results.json") in result


@pytest.mark.short
def test_archive_benchmark_compression_types():
    """Test archive_benchmark handles different compression types correctly."""
    mock_benchmark = Mock()
    mock_benchmark.get_benchmark_name.return_value = "test_benchmark"
    mock_benchmark.get_benchmark_version.return_value = "1.0"

    with patch("omnibenchmark.archive.archive.prepare_archive_config") as mock_config:
        mock_config.return_value = [Path("benchmark.yaml")]

        # Test different compression types
        test_cases = [
            (zipfile.ZIP_STORED, ".zip"),
            (zipfile.ZIP_DEFLATED, ".zip"),
            (zipfile.ZIP_BZIP2, ".bz2"),
            (zipfile.ZIP_LZMA, ".xz"),
        ]

        for compression, expected_ext in test_cases:
            with (
                patch("zipfile.ZipFile") as mock_zipfile,
                patch("pathlib.Path.is_file", return_value=True),
            ):
                mock_archive = MagicMock()
                mock_zipfile.return_value.__enter__.return_value = mock_archive

                result = archive_benchmark(
                    benchmark=mock_benchmark, compression=compression, dry_run=False
                )

                # Check that correct file extension is used
                assert result[0].suffix == expected_ext
                assert f"test_benchmark_1.0{expected_ext}" in str(result[0])


@pytest.mark.short
def test_prepare_archive_config():
    """Test prepare_archive_config returns definition file."""
    mock_benchmark = Mock()
    mock_definition_file = Mock(spec=Path)
    mock_definition_file.is_file.return_value = True
    mock_benchmark.get_definition_file.return_value = mock_definition_file

    result = prepare_archive_config(mock_benchmark)

    assert len(result) == 1
    assert result[0] == mock_definition_file


@pytest.mark.short
def test_prepare_archive_config_no_file():
    """Test prepare_archive_config handles missing definition file."""
    mock_benchmark = Mock()
    mock_definition_file = Mock(spec=Path)
    mock_definition_file.is_file.return_value = False
    mock_benchmark.get_definition_file.return_value = mock_definition_file

    result = prepare_archive_config(mock_benchmark)

    assert len(result) == 0


@pytest.mark.short
def test_prepare_archive_results_local_existing_files(tmp_path):
    """Test prepare_archive_results_local includes existing files."""
    # Create test files
    results_dir = tmp_path / "out"
    results_dir.mkdir()
    test_file = results_dir / "test.json"
    test_file.write_text("{}")

    mock_benchmark = Mock()

    with patch(
        "omnibenchmark.archive.components.get_expected_benchmark_output_files"
    ) as mock_expected:
        mock_expected.return_value = [str(test_file)]

        result = prepare_archive_results_local(mock_benchmark, str(results_dir))

        assert len(result) == 1
        assert result[0] == test_file


@pytest.mark.short
def test_prepare_archive_results_local_no_files():
    """Test prepare_archive_results_local handles missing results directory."""
    mock_benchmark = Mock()

    with patch(
        "omnibenchmark.archive.components.get_expected_benchmark_output_files"
    ) as mock_expected:
        mock_expected.return_value = []

        result = prepare_archive_results_local(mock_benchmark, "nonexistent")

        assert len(result) == 0


@pytest.mark.short
def test_prepare_archive_code_handles_clone_failures():
    """Test prepare_archive_code handles repository clone failures gracefully."""
    mock_benchmark = Mock()
    mock_node = Mock()
    mock_repo = Mock()
    mock_repo.url = "https://github.com/test/repo.git"
    mock_repo.commit = "abc123"
    mock_node.get_repository.return_value = mock_repo
    mock_benchmark.get_nodes.return_value = [mock_node]

    with patch("omnibenchmark.git.clone.clone_module") as mock_clone:
        # Simulate clone failure
        mock_clone.side_effect = Exception("Clone failed")

        result = prepare_archive_code(mock_benchmark)

        # Should handle exception gracefully and return empty list
        assert result == []


@pytest.mark.short
def test_archive_benchmark_selective_components():
    """Test archive_benchmark respects component selection flags."""
    mock_benchmark = Mock()
    mock_benchmark.get_benchmark_name.return_value = "test_benchmark"
    mock_benchmark.get_benchmark_version.return_value = "1.0"

    with (
        patch("omnibenchmark.archive.archive.prepare_archive_config") as mock_config,
        patch("omnibenchmark.archive.archive.prepare_archive_code") as mock_code,
        patch(
            "omnibenchmark.archive.archive.prepare_archive_software"
        ) as mock_software,
        patch("omnibenchmark.archive.archive.prepare_archive_results") as mock_results,
    ):
        mock_config.return_value = [Path("benchmark.yaml")]
        mock_code.return_value = [Path("code.py")]
        mock_software.return_value = [Path("environment.yml")]
        mock_results.return_value = [Path("results.json")]

        # Test with only config enabled
        result = archive_benchmark(
            benchmark=mock_benchmark,
            config=True,
            code=False,
            software=False,
            results=False,
            dry_run=True,
        )

        assert len(result) == 1
        assert Path("benchmark.yaml") in result

        # Verify only config function was called
        mock_config.assert_called_once()
        mock_code.assert_not_called()
        mock_software.assert_not_called()
        mock_results.assert_not_called()
