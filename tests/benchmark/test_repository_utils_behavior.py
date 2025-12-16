"""Behavior-driven tests for repository utilities."""

from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import tempfile

from omnibenchmark.benchmark.repository_utils import (
    RepositoryManager,
    cleanup_temp_repositories,
    get_module_repository_info,
    resolve_module_repository,
)


class TestRepositoryManager:
    """Behavior tests for RepositoryManager."""

    def test_when_created_then_initializes_with_prefix(self):
        """When RepositoryManager is created, it should use the given prefix."""
        manager = RepositoryManager(prefix="test_prefix")

        assert manager.prefix == "test_prefix"
        assert "omnibenchmark" in str(manager.temp_base_dir)
        assert "tmp_repos" in str(manager.temp_base_dir)

    def test_when_used_as_context_manager_then_cleans_up(self):
        """When used as context manager, should clean up on exit."""
        with patch.object(
            RepositoryManager, "cleanup_temp_directories"
        ) as mock_cleanup:
            with RepositoryManager():
                pass  # Do nothing

            mock_cleanup.assert_called_once()

    def test_when_get_local_repo_path_called_then_returns_expected_path(self):
        """When getting local repo path, should return path under .snakemake/repos."""
        manager = RepositoryManager()

        with patch(
            "omnibenchmark.benchmark.repository_utils.get_repo_hash"
        ) as mock_hash:
            mock_hash.return_value = "abc123"

            path = manager.get_local_repo_path(
                "https://github.com/test/repo", "commit123"
            )

            assert ".snakemake/repos" in str(path)
            assert "abc123" in str(path)
            mock_hash.assert_called_once_with(
                "https://github.com/test/repo", "commit123"
            )

    def test_when_clone_to_temp_succeeds_then_returns_path(self):
        """When clone succeeds, should return the cloned path."""
        manager = RepositoryManager(prefix="test")

        with patch(
            "omnibenchmark.benchmark.repository_utils.clone_module"
        ) as mock_clone:
            mock_path = Path("/tmp/test_clone")
            mock_clone.return_value = mock_path

            result = manager.clone_to_temp(
                "https://github.com/test/repo", "commit123", "test_module"
            )

            assert result == mock_path
            mock_clone.assert_called_once()

    def test_when_clone_fails_then_returns_none_and_cleans_up(self):
        """When clone fails, should return None and attempt cleanup."""
        manager = RepositoryManager(prefix="test")

        with patch(
            "omnibenchmark.benchmark.repository_utils.clone_module"
        ) as mock_clone:
            mock_clone.side_effect = Exception("Clone failed")

            result = manager.clone_to_temp(
                "https://github.com/test/repo", "commit123", "test_module"
            )

            assert result is None

    def test_when_get_repository_files_called_then_reads_standard_files(self):
        """When getting repository files, should read citation, license, and config."""
        # Create temporary directory with test files
        with tempfile.TemporaryDirectory() as temp_dir:
            repo_path = Path(temp_dir)

            # Create test files
            (repo_path / "CITATION.cff").write_text("citation content")
            (repo_path / "LICENSE").write_text("license content")
            (repo_path / "omnibenchmark.yaml").write_text("config content")

            manager = RepositoryManager()
            files = manager.get_repository_files(repo_path)

            assert files["citation"] == "citation content"
            assert files["license"] == "license content"
            assert files["omnibenchmark"] == "config content"

    def test_when_files_missing_then_returns_none_for_missing(self):
        """When files are missing, should return None for those files."""
        with tempfile.TemporaryDirectory() as temp_dir:
            repo_path = Path(temp_dir)
            # Don't create any files

            manager = RepositoryManager()
            files = manager.get_repository_files(repo_path)

            assert files["citation"] is None
            assert files["license"] is None
            assert files["omnibenchmark"] is None

    def test_when_get_files_present_called_then_checks_existence(self):
        """When checking files present, should return dict of booleans."""
        with tempfile.TemporaryDirectory() as temp_dir:
            repo_path = Path(temp_dir)

            # Create some files but not others
            (repo_path / "CITATION.cff").write_text("citation")
            (repo_path / "LICENSE").write_text("license")
            # Don't create omnibenchmark.yaml

            manager = RepositoryManager()
            files_present = manager.get_files_present(repo_path)

            assert files_present["CITATION.cff"] is True
            assert files_present["LICENSE"] is True
            assert files_present["omnibenchmark.yaml"] is False
            assert files_present["LICENSE.txt"] is False


class TestCleanupTempRepositories:
    """Behavior tests for cleanup function."""

    def test_when_temp_dir_exists_then_removes_it(self):
        """When temp directory exists, cleanup should remove it."""
        with patch("omnibenchmark.benchmark.repository_utils.Path") as mock_path_class:
            mock_path = MagicMock()
            mock_path.exists.return_value = True
            mock_path_class.return_value = mock_path

            with patch(
                "omnibenchmark.benchmark.repository_utils.shutil.rmtree"
            ) as mock_rmtree:
                cleanup_temp_repositories()

                mock_rmtree.assert_called_once()

    def test_when_temp_dir_missing_then_does_nothing(self):
        """When temp directory doesn't exist, should do nothing."""
        with patch("tempfile.gettempdir", return_value="/tmp"):
            with patch(
                "omnibenchmark.benchmark.repository_utils.Path"
            ) as mock_path_class:
                # Mock the path construction chain: Path(gettempdir()) / "omnibenchmark" / "tmp_repos"
                mock_final_path = MagicMock()
                mock_final_path.exists.return_value = False

                mock_intermediate_path = MagicMock()
                mock_intermediate_path.__truediv__.return_value = mock_final_path

                mock_base_path = MagicMock()
                mock_base_path.__truediv__.return_value = mock_intermediate_path

                mock_path_class.return_value = mock_base_path

                with patch(
                    "omnibenchmark.benchmark.repository_utils.shutil.rmtree"
                ) as mock_rmtree:
                    cleanup_temp_repositories()

                    mock_rmtree.assert_not_called()


class TestGetModuleRepositoryInfo:
    """Behavior tests for extracting module repository information."""

    def test_when_module_has_url_and_commit_then_returns_both(self):
        """When module has URL and commit, should return both."""
        mock_benchmark = Mock()
        mock_module = Mock()
        mock_repo_info = Mock()
        mock_repo_info.url = "https://github.com/test/repo"
        mock_repo_info.commit_hash = "abc123"

        mock_benchmark.model.get_module_repository.return_value = mock_repo_info

        url, commit = get_module_repository_info(mock_benchmark, mock_module)

        assert url == "https://github.com/test/repo"
        assert commit == "abc123"

    def test_when_using_repository_field_then_returns_it(self):
        """When repo info uses 'repository' instead of 'url', should still work."""
        mock_benchmark = Mock()
        mock_module = Mock()
        mock_repo_info = Mock()
        mock_repo_info.url = None
        mock_repo_info.repository = "https://github.com/test/repo"
        mock_repo_info.commit_hash = "abc123"

        mock_benchmark.model.get_module_repository.return_value = mock_repo_info

        url, commit = get_module_repository_info(mock_benchmark, mock_module)

        assert url == "https://github.com/test/repo"

    def test_when_using_version_field_then_returns_it(self):
        """When repo info uses 'version' instead of 'commit_hash', should still work."""
        mock_benchmark = Mock()
        mock_module = Mock()
        mock_repo_info = Mock()
        mock_repo_info.url = "https://github.com/test/repo"
        mock_repo_info.commit_hash = None
        mock_repo_info.commit = None
        mock_repo_info.version = "v1.0.0"

        mock_benchmark.model.get_module_repository.return_value = mock_repo_info

        url, commit = get_module_repository_info(mock_benchmark, mock_module)

        assert commit == "v1.0.0"

    def test_when_info_missing_then_returns_none_tuple(self):
        """When repository info is missing, should return (None, None)."""
        mock_benchmark = Mock()
        mock_module = Mock()
        mock_repo_info = Mock()
        mock_repo_info.url = None
        mock_repo_info.repository = None
        mock_repo_info.commit_hash = None
        mock_repo_info.commit = None
        mock_repo_info.version = None

        mock_benchmark.model.get_module_repository.return_value = mock_repo_info

        url, commit = get_module_repository_info(mock_benchmark, mock_module)

        assert url is None
        assert commit is None


class TestResolveModuleRepository:
    """Behavior tests for resolving module repository paths."""

    def test_when_local_repo_exists_then_returns_local_path(self):
        """When local repository exists, should return it without cloning."""
        mock_benchmark = Mock()
        mock_module = Mock()
        mock_repo_manager = Mock()

        # Mock repository info
        with patch(
            "omnibenchmark.benchmark.repository_utils.get_module_repository_info"
        ) as mock_get_info:
            mock_get_info.return_value = ("https://github.com/test/repo", "abc123")

            # Mock local path exists
            mock_local_path = Mock()
            mock_local_path.exists.return_value = True
            mock_repo_manager.get_local_repo_path.return_value = mock_local_path

            result = resolve_module_repository(
                mock_benchmark, mock_module, "test_module", mock_repo_manager
            )

            assert result == mock_local_path
            mock_repo_manager.clone_to_temp.assert_not_called()

    def test_when_local_repo_missing_then_clones_to_temp(self):
        """When local repository doesn't exist, should clone to temp."""
        mock_benchmark = Mock()
        mock_module = Mock()
        mock_repo_manager = Mock()

        with patch(
            "omnibenchmark.benchmark.repository_utils.get_module_repository_info"
        ) as mock_get_info:
            mock_get_info.return_value = ("https://github.com/test/repo", "abc123")

            # Mock local path doesn't exist
            mock_local_path = Mock()
            mock_local_path.exists.return_value = False
            mock_repo_manager.get_local_repo_path.return_value = mock_local_path

            # Mock temp clone succeeds
            mock_temp_path = Mock()
            mock_repo_manager.clone_to_temp.return_value = mock_temp_path

            result = resolve_module_repository(
                mock_benchmark, mock_module, "test_module", mock_repo_manager
            )

            assert result == mock_temp_path
            mock_repo_manager.clone_to_temp.assert_called_once()

    def test_when_repo_info_missing_then_returns_none(self):
        """When repository info is missing, should return None."""
        mock_benchmark = Mock()
        mock_module = Mock()
        mock_repo_manager = Mock()

        with patch(
            "omnibenchmark.benchmark.repository_utils.get_module_repository_info"
        ) as mock_get_info:
            mock_get_info.return_value = (None, None)

            result = resolve_module_repository(
                mock_benchmark, mock_module, "test_module", mock_repo_manager
            )

            assert result is None

    def test_when_clone_fails_then_returns_none(self):
        """When cloning fails, should return None."""
        mock_benchmark = Mock()
        mock_module = Mock()
        mock_repo_manager = Mock()

        with patch(
            "omnibenchmark.benchmark.repository_utils.get_module_repository_info"
        ) as mock_get_info:
            mock_get_info.return_value = ("https://github.com/test/repo", "abc123")

            # Mock local path doesn't exist
            mock_local_path = Mock()
            mock_local_path.exists.return_value = False
            mock_repo_manager.get_local_repo_path.return_value = mock_local_path

            # Mock clone fails
            mock_repo_manager.clone_to_temp.return_value = None

            result = resolve_module_repository(
                mock_benchmark, mock_module, "test_module", mock_repo_manager
            )

            assert result is None
