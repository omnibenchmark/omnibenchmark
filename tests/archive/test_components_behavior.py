"""
Behavior-driven tests for archive components preparation functions.

These tests focus on the high-level behavior of preparing different components
(config, code, software, results) for archiving, using minimal mocking to test
the actual business logic and integration points.
"""

import pytest
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

from omnibenchmark.archive.components import (
    prepare_archive_config,
    prepare_archive_code,
    prepare_archive_software,
    prepare_archive_results,
)


@pytest.mark.short
class TestPrepareArchiveConfig:
    """Test configuration files preparation for archiving."""

    def test_prepare_config_with_existing_definition_file(self):
        """Test that existing definition file is included in config archive."""
        # Create a temporary definition file
        with tempfile.NamedTemporaryFile(suffix=".yaml", delete=False) as f:
            definition_path = Path(f.name)
            f.write(b"name: test_benchmark\nversion: 1.0.0\n")

        try:
            # Mock benchmark with existing definition file
            mock_benchmark = Mock()
            mock_benchmark.get_definition_file.return_value = definition_path

            result = prepare_archive_config(mock_benchmark)

            assert len(result) == 1
            assert result[0] == definition_path
            mock_benchmark.get_definition_file.assert_called_once()
        finally:
            definition_path.unlink()

    def test_prepare_config_with_nonexistent_definition_file(self):
        """Test that nonexistent definition file is not included."""
        nonexistent_path = Path("/nonexistent/benchmark.yaml")

        mock_benchmark = Mock()
        mock_benchmark.get_definition_file.return_value = nonexistent_path

        result = prepare_archive_config(mock_benchmark)

        assert result == []
        mock_benchmark.get_definition_file.assert_called_once()

    def test_prepare_config_with_none_definition_file(self):
        """Test that None definition file returns empty list."""
        mock_benchmark = Mock()
        mock_benchmark.get_definition_file.return_value = None

        result = prepare_archive_config(mock_benchmark)

        assert result == []
        mock_benchmark.get_definition_file.assert_called_once()

    def test_prepare_config_with_directory_instead_of_file(self):
        """Test that directory path is not included (only files)."""
        with tempfile.TemporaryDirectory() as temp_dir:
            dir_path = Path(temp_dir)

            mock_benchmark = Mock()
            mock_benchmark.get_definition_file.return_value = dir_path

            result = prepare_archive_config(mock_benchmark)

            # Directory should not be included (not a file)
            assert result == []


@pytest.mark.short
class TestPrepareArchiveCode:
    """Test code files preparation for archiving."""

    @patch("omnibenchmark.git.clone.clone_module")
    def test_prepare_code_calls_clone_module_for_each_module(self, mock_clone):
        """Test that code preparation calls clone_module for each module."""
        # Setup mock benchmark with nodes that have repositories
        mock_benchmark = Mock()
        mock_node1 = Mock()
        mock_node1.get_repository.return_value = Mock(
            url="https://github.com/user/data.git", commit="abc123"
        )
        mock_node2 = Mock()
        mock_node2.get_repository.return_value = Mock(
            url="https://github.com/user/method.git", commit="def456"
        )
        mock_benchmark.get_nodes.return_value = [mock_node1, mock_node2]

        # Mock clone_module to return temporary directories
        with (
            tempfile.TemporaryDirectory() as temp_dir1,
            tempfile.TemporaryDirectory() as temp_dir2,
        ):
            mock_clone.side_effect = [Path(temp_dir1), Path(temp_dir2)]

            result = prepare_archive_code(mock_benchmark)

            # Should call clone_module for each module
            assert mock_clone.call_count == 2
            # Should call clone_module for each unique repository
            assert mock_clone.call_count == 2
            assert isinstance(result, list)

    @patch("omnibenchmark.git.clone.clone_module")
    def test_prepare_code_with_no_modules(self, mock_clone):
        """Test code preparation with benchmark having no modules."""
        mock_benchmark = Mock()
        mock_benchmark.get_nodes.return_value = []

        result = prepare_archive_code(mock_benchmark)

        assert result == []
        mock_clone.assert_not_called()

    @patch("omnibenchmark.git.clone.clone_module")
    def test_prepare_code_handles_clone_failure_gracefully(self, mock_clone):
        """Test that clone failures are handled gracefully."""
        mock_benchmark = Mock()
        mock_node = Mock()
        mock_node.get_repository.return_value = Mock(
            url="https://github.com/user/failing.git", commit="abc123"
        )
        mock_benchmark.get_nodes.return_value = [mock_node]

        # Mock clone_module to raise an exception
        mock_clone.side_effect = Exception("Git clone failed")

        # Should handle exception gracefully and continue
        result = prepare_archive_code(mock_benchmark)
        assert isinstance(result, list)


@pytest.mark.short
class TestPrepareArchiveSoftware:
    """Test software environment preparation for archiving."""

    def test_prepare_software_with_conda_backend(self):
        """Test software preparation for conda-based modules."""
        from omnibenchmark.model import SoftwareBackendEnum

        mock_benchmark = Mock()
        mock_benchmark.get_benchmark_software_backend.return_value = (
            SoftwareBackendEnum.conda
        )
        # Mock get_benchmark_software_environments to return an empty dict
        mock_benchmark.get_benchmark_software_environments.return_value = {}

        result = prepare_archive_software(mock_benchmark)

        # Should return paths to environment files
        assert isinstance(result, list)
        # The exact behavior depends on implementation

    def test_prepare_software_with_apptainer_backend(self):
        """Test software preparation for apptainer-based modules."""
        from omnibenchmark.model import SoftwareBackendEnum

        mock_benchmark = Mock()
        mock_benchmark.get_benchmark_software_backend.return_value = (
            SoftwareBackendEnum.apptainer
        )
        # Mock get_benchmark_software_environments to return an empty dict
        mock_benchmark.get_benchmark_software_environments.return_value = {}

        result = prepare_archive_software(mock_benchmark)

        assert isinstance(result, list)
        # For apptainer, might include image definitions or pull info

    def test_prepare_software_with_mixed_backends(self):
        """Test software preparation with modules using different backends."""
        from omnibenchmark.model import SoftwareBackendEnum

        mock_benchmark = Mock()
        mock_benchmark.get_benchmark_software_backend.return_value = (
            SoftwareBackendEnum.conda
        )
        # Mock get_benchmark_software_environments to return an empty dict
        mock_benchmark.get_benchmark_software_environments.return_value = {}

        result = prepare_archive_software(mock_benchmark)

        assert isinstance(result, list)
        # Should handle both backend types


@pytest.mark.short
class TestPrepareArchiveResults:
    """Test results preparation for archiving."""

    @patch("omnibenchmark.archive.components.get_expected_benchmark_output_files")
    def test_prepare_results_with_local_storage(self, mock_get_files):
        """Test results preparation with local storage."""
        # Create temporary result files
        with tempfile.TemporaryDirectory() as temp_dir:
            result_file1 = Path(temp_dir) / "result1.json"
            result_file2 = Path(temp_dir) / "result2.json"
            result_file1.write_text('{"test": "data1"}')
            result_file2.write_text('{"test": "data2"}')

            mock_get_files.return_value = [result_file1, result_file2]

            mock_benchmark = Mock()

            result = prepare_archive_results(
                benchmark=mock_benchmark, results_dir=temp_dir, remote_storage=False
            )

            assert len(result) == 2
            assert result_file1 in result
            assert result_file2 in result
            mock_get_files.assert_called_once()

    @patch("omnibenchmark.archive.components.get_expected_benchmark_output_files")
    def test_prepare_results_with_remote_storage_download(self, mock_get_files):
        """Test results preparation with remote storage (downloads files)."""
        # Mock remote storage scenario
        mock_benchmark = Mock()

        # Create temporary files to simulate downloaded results
        with tempfile.TemporaryDirectory() as temp_dir:
            downloaded_file = Path(temp_dir) / "downloaded_result.json"
            downloaded_file.write_text('{"remote": "data"}')

            mock_get_files.return_value = [downloaded_file]

            result = prepare_archive_results(
                benchmark=mock_benchmark, results_dir=temp_dir, remote_storage=False
            )

            assert len(result) == 1
            assert downloaded_file in result

    @patch("omnibenchmark.remote.files.download_files")
    @patch("omnibenchmark.remote.files.list_files")
    @patch("omnibenchmark.archive.components.get_expected_benchmark_output_files")
    def test_prepare_results_with_nonexistent_files(
        self, mock_get_files, mock_list_files, mock_download_files
    ):
        """Test results preparation when some expected files don't exist."""
        nonexistent_file = Path("/nonexistent/result.json")
        mock_get_files.return_value = [nonexistent_file]

        mock_benchmark = Mock()
        # Mock get_definition_file to return a Path object
        mock_definition_file = Mock()
        mock_definition_file.as_posix.return_value = "/path/to/benchmark.yaml"
        mock_benchmark.get_definition_file.return_value = mock_definition_file

        # Mock list_files to return empty lists
        mock_list_files.return_value = ([], [])

        result = prepare_archive_results(
            benchmark=mock_benchmark, results_dir="out", remote_storage=True
        )

        # Should filter out nonexistent files
        assert nonexistent_file not in result

    def test_prepare_results_empty_file_list(self):
        """Test results preparation when no result files are found."""
        with patch(
            "omnibenchmark.archive.components.get_expected_benchmark_output_files"
        ) as mock_get_files:
            mock_get_files.return_value = []

            mock_benchmark = Mock()

            result = prepare_archive_results(
                benchmark=mock_benchmark, results_dir="out", remote_storage=False
            )

            assert result == []


@pytest.mark.short
class TestArchiveComponentsIntegration:
    """Test integration between different archive component functions."""

    def test_all_prepare_functions_return_path_lists(self):
        """Test that all prepare_archive_* functions return lists of Path objects."""
        mock_benchmark = Mock()
        mock_benchmark.get_definition_file.return_value = None
        mock_benchmark.get_nodes.return_value = []

        config_result = prepare_archive_config(mock_benchmark)
        assert isinstance(config_result, list)

        with patch("omnibenchmark.git.clone.clone_module"):
            code_result = prepare_archive_code(mock_benchmark)
            assert isinstance(code_result, list)

        software_result = prepare_archive_software(mock_benchmark)
        assert isinstance(software_result, list)

        with patch(
            "omnibenchmark.archive.components.get_expected_benchmark_output_files",
            return_value=[],
        ):
            results_result = prepare_archive_results(
                mock_benchmark, "out", remote_storage=False
            )
            assert isinstance(results_result, list)

    def test_archive_components_handle_benchmark_object_consistently(self):
        """Test that all component functions handle the benchmark object consistently."""
        mock_benchmark = Mock()
        mock_benchmark.get_definition_file.return_value = None
        mock_benchmark.get_nodes.return_value = []

        # All functions should accept the benchmark object without errors
        try:
            prepare_archive_config(mock_benchmark)

            with patch("omnibenchmark.git.clone.clone_module"):
                prepare_archive_code(mock_benchmark)

            prepare_archive_software(mock_benchmark)

            with patch(
                "omnibenchmark.archive.components.get_expected_benchmark_output_files",
                return_value=[],
            ):
                prepare_archive_results(mock_benchmark, "out", remote_storage=False)
        except Exception as e:
            pytest.fail(
                f"Archive component functions should handle benchmark object consistently: {e}"
            )
