"""
Tests for GitAwareBenchmarkVersionManager.

All tests in this file are marked as 'short' since they mock git operations
and don't require actual git repositories or external dependencies.
"""

from unittest.mock import patch, MagicMock

import pytest

from omnibenchmark.versioning.git import GitAwareBenchmarkVersionManager


@pytest.fixture
def mock_git_repo(tmp_path):
    """Create a mock git repository setup."""
    # Create a fake .git directory
    git_dir = tmp_path / ".git"
    git_dir.mkdir()
    return tmp_path


@pytest.mark.short
class TestGitAwareBenchmarkVersionManager:
    """Test the GitAwareBenchmarkVersionManager class."""

    def test_initialization_with_git_repo(self, mock_git_repo):
        """Test initialization with a valid git repository."""
        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            # Mock successful git repo
            mock_repo_instance = MagicMock()
            MockRepo.return_value = mock_repo_instance

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file,
                git_repo_path=mock_git_repo,
                lock_dir=mock_git_repo / "locks",
            )

            assert manager.benchmark_name == "test"
            assert manager.git_repo_path == mock_git_repo
            assert manager.git_available is True

    def test_initialization_without_git(self, tmp_path):
        """Test initialization when git is not available."""
        benchmark_file = tmp_path / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            # Mock git not available
            from git.exc import InvalidGitRepositoryError

            MockRepo.side_effect = InvalidGitRepositoryError()

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file,
                git_repo_path=tmp_path,
                lock_dir=tmp_path / "locks",
            )

            # Should not fail, just mark git as unavailable
            assert manager.git_available is False

    def test_reconstruct_version_history(self, mock_git_repo):
        """Test reconstructing version history from git."""
        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            mock_repo = MagicMock()
            MockRepo.return_value = mock_repo

            # Mock commits
            mock_commit1 = MagicMock()
            mock_commit1.hexsha = "abc123"
            mock_commit2 = MagicMock()
            mock_commit2.hexsha = "def456"

            mock_repo.iter_commits.return_value = [
                mock_commit2,
                mock_commit1,
            ]  # Newest first

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file, git_repo_path=mock_git_repo
            )

            # Mock the file content extraction
            with patch.object(manager, "_get_file_content_at_commit") as mock_get:
                with patch.object(
                    manager, "_extract_version_from_yaml"
                ) as mock_extract:
                    mock_get.side_effect = ["version: 1.0.0\n", "version: 1.1.0\n"]
                    mock_extract.side_effect = ["1.0.0", "1.1.0"]

                    versions = manager.reconstruct_version_history()

                    # Should be in chronological order (oldest first)
                    assert versions == ["1.0.0", "1.1.0"]

    def test_reconstruct_version_history_no_git(self, tmp_path):
        """Test version history when git is not available."""
        benchmark_file = tmp_path / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            from git.exc import InvalidGitRepositoryError

            MockRepo.side_effect = InvalidGitRepositoryError()

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file, git_repo_path=tmp_path
            )

            versions = manager.reconstruct_version_history()
            assert versions == []

    def test_get_current_git_info_complete(self, mock_git_repo):
        """Test getting complete git information."""
        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            mock_repo = MagicMock()
            MockRepo.return_value = mock_repo

            # Mock commit info
            mock_commit = MagicMock()
            mock_commit.hexsha = "abc123def456"
            mock_commit.author.name = "John Doe"
            mock_commit.author.email = "john@example.com"

            import datetime

            mock_commit.committed_datetime = datetime.datetime(2024, 1, 1, 12, 0, 0)

            mock_repo.head.commit = mock_commit
            mock_repo.head.is_detached = False
            mock_repo.active_branch.name = "main"
            mock_repo.is_dirty.return_value = False
            mock_repo.remotes.origin.url = "https://github.com/example/repo.git"

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file, git_repo_path=mock_git_repo
            )

            info = manager.get_current_git_info()

            assert info["commit"] == "abc123de"  # Short hash
            assert info["commit_full"] == "abc123def456"
            assert info["branch"] == "main"
            assert "John Doe" in info["author"]
            assert info["clean"] is True

    def test_get_current_git_info_dirty_repo(self, mock_git_repo):
        """Test getting git info with uncommitted changes."""
        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            mock_repo = MagicMock()
            MockRepo.return_value = mock_repo

            mock_commit = MagicMock()
            mock_commit.hexsha = "abc123"
            mock_repo.head.commit = mock_commit
            mock_repo.is_dirty.return_value = True  # Dirty repo
            mock_repo.remotes = []  # No remotes

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file, git_repo_path=mock_git_repo
            )

            info = manager.get_current_git_info()

            assert info["clean"] is False

    def test_get_current_git_info_no_git(self, tmp_path):
        """Test getting git info when git is not available."""
        benchmark_file = tmp_path / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            from git.exc import InvalidGitRepositoryError

            MockRepo.side_effect = InvalidGitRepositoryError()

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file, git_repo_path=tmp_path
            )

            info = manager.get_current_git_info()

            assert info == {}

    def test_create_version_with_git_tracking(self, mock_git_repo):
        """Test creating a version with git tracking."""
        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            mock_repo = MagicMock()
            MockRepo.return_value = mock_repo

            mock_commit = MagicMock()
            mock_commit.hexsha = "abc123"
            mock_repo.head.commit = mock_commit
            mock_repo.head.is_detached = False
            mock_repo.active_branch.name = "main"

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file,
                git_repo_path=mock_git_repo,
                lock_dir=mock_git_repo / "locks",
            )

            version = manager.create_version_with_git_tracking("1.0.0")

            assert version == "1.0.0"
            assert "1.0.0" in manager.get_versions()

    def test_create_version_with_git_tracking_auto_increment(self, mock_git_repo):
        """Test creating auto-incremented version with git tracking."""
        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            mock_repo = MagicMock()
            MockRepo.return_value = mock_repo

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file,
                git_repo_path=mock_git_repo,
                lock_dir=mock_git_repo / "locks",
            )

            # Create without specifying version
            version = manager.create_version_with_git_tracking()

            assert version == "0.1"  # Default first version

    def test_extract_version_from_yaml(self, mock_git_repo):
        """Test extracting version from YAML content."""
        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            mock_repo = MagicMock()
            MockRepo.return_value = mock_repo

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file, git_repo_path=mock_git_repo
            )

            # Test valid YAML
            yaml_content = """
version: 2.1.0
id: test
benchmarker: test
software_backend: SNAKEMAKE
software_environments: []
stages: []
"""

            with patch("omnibenchmark.model.Benchmark") as MockBenchmark:
                mock_benchmark = MagicMock()
                mock_benchmark.version = "2.1.0"
                MockBenchmark.return_value = mock_benchmark

                version = manager._extract_version_from_yaml(yaml_content)
                assert version == "2.1.0"

            # Test invalid YAML
            invalid_yaml = "invalid: yaml: content:"
            version = manager._extract_version_from_yaml(invalid_yaml)
            assert version is None

    def test_get_version_at_commit(self, mock_git_repo):
        """Test getting version from a specific commit."""
        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            mock_repo = MagicMock()
            MockRepo.return_value = mock_repo

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file, git_repo_path=mock_git_repo
            )

            with patch.object(manager, "_get_file_content_at_commit") as mock_get:
                with patch.object(
                    manager, "_extract_version_from_yaml"
                ) as mock_extract:
                    mock_get.return_value = "version: 1.5.0\n"
                    mock_extract.return_value = "1.5.0"

                    version = manager.get_version_at_commit("abc123")
                    assert version == "1.5.0"

                    # Test when content not found
                    mock_get.return_value = None
                    version = manager.get_version_at_commit("def456")
                    assert version is None

    def test_get_commits_for_version(self, mock_git_repo):
        """Test getting commits for a specific version."""
        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            mock_repo = MagicMock()
            MockRepo.return_value = mock_repo

            # Mock commits
            mock_commit1 = MagicMock()
            mock_commit1.hexsha = "abc123def"
            mock_commit2 = MagicMock()
            mock_commit2.hexsha = "def456ghi"

            mock_repo.iter_commits.return_value = [mock_commit1, mock_commit2]

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file, git_repo_path=mock_git_repo
            )

            with patch.object(manager, "_get_file_content_at_commit") as mock_get:
                with patch.object(
                    manager, "_extract_version_from_yaml"
                ) as mock_extract:
                    mock_get.side_effect = ["version: 1.0.0\n", "version: 1.0.0\n"]
                    mock_extract.side_effect = ["1.0.0", "1.0.0"]

                    commits = manager.get_commits_for_version("1.0.0")
                    assert len(commits) == 2
                    assert "abc123de" in commits
                    assert "def456gh" in commits

    def test_inherited_functionality(self, mock_git_repo):
        """Test that inherited BenchmarkVersionManager functionality works."""
        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            mock_repo = MagicMock()
            MockRepo.return_value = mock_repo

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file,
                git_repo_path=mock_git_repo,
                lock_dir=mock_git_repo / "locks",
            )

            # Test inherited methods
            assert manager.get_versions() == []
            # Don't call get_current_version as it tries to parse incomplete YAML
            assert manager._current_version is None

            manager.create_version("1.0.0")
            assert manager.version_exists("1.0.0")
            assert manager.compare_versions("1.0.0", "2.0.0") == -1

    def test_detached_head_state(self, mock_git_repo):
        """Test handling of detached HEAD state."""
        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            mock_repo = MagicMock()
            MockRepo.return_value = mock_repo

            mock_commit = MagicMock()
            mock_commit.hexsha = "abc123"
            mock_repo.head.commit = mock_commit
            mock_repo.head.is_detached = True  # Detached HEAD

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file, git_repo_path=mock_git_repo
            )

            info = manager.get_current_git_info()

            # Should not have branch key when in detached HEAD
            assert "branch" not in info
            assert info["commit"] == "abc123"[:8]  # Short hash

    def test_initialize_from_git_history(self, mock_git_repo):
        """Test initializing known versions from git history."""
        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            mock_repo = MagicMock()
            MockRepo.return_value = mock_repo

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file, git_repo_path=mock_git_repo
            )

            with patch.object(
                manager, "reconstruct_version_history"
            ) as mock_reconstruct:
                mock_reconstruct.return_value = ["1.0.0", "1.1.0", "2.0.0"]

                manager.initialize_from_git_history()

                assert manager.get_versions() == ["1.0.0", "1.1.0", "2.0.0"]

    def test_update_benchmark_version_and_commit(self, mock_git_repo):
        """Test updating benchmark version and committing to git."""
        from unittest.mock import Mock, patch, mock_open

        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\nid: test\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            mock_repo = MagicMock()
            mock_index = MagicMock()
            mock_commit = MagicMock()
            mock_commit.hexsha = "abc123def456"

            mock_repo.index = mock_index
            mock_index.commit.return_value = mock_commit
            MockRepo.return_value = mock_repo

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file, git_repo_path=mock_git_repo
            )

            # Create mock benchmark execution
            mock_benchmark_execution = Mock()
            mock_model = Mock()
            mock_model.version = "1.0.0"
            mock_model.to_yaml.return_value = "version: 2.0.0\nid: test\n"
            mock_benchmark_execution.model = mock_model

            # Mock file writing
            with patch("builtins.open", mock_open()) as mock_file:
                commit_hash = manager.update_benchmark_version_and_commit(
                    mock_benchmark_execution, "2.0.0", "Custom commit message"
                )

                # Verify version was updated in model
                assert mock_model.version == "2.0.0"

                # Verify YAML was written to file
                mock_file.assert_called_once_with(benchmark_file, "w")
                mock_file().write.assert_called_once_with("version: 2.0.0\nid: test\n")

                # Verify git operations
                mock_index.add.assert_called_once()
                mock_index.commit.assert_called_once_with("Custom commit message")

                # Verify return value
                assert commit_hash == "abc123def456"

    def test_update_benchmark_version_and_commit_default_message(self, mock_git_repo):
        """Test updating benchmark version with default commit message."""
        from unittest.mock import Mock, patch, mock_open

        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\nid: test\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            mock_repo = MagicMock()
            mock_index = MagicMock()
            mock_commit = MagicMock()
            mock_commit.hexsha = "def456ghi789"

            mock_repo.index = mock_index
            mock_index.commit.return_value = mock_commit
            MockRepo.return_value = mock_repo

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file, git_repo_path=mock_git_repo
            )

            # Create mock benchmark execution
            mock_benchmark_execution = Mock()
            mock_model = Mock()
            mock_model.version = "1.0.0"
            mock_model.to_yaml.return_value = "version: 1.5.0\nid: test\n"
            mock_benchmark_execution.model = mock_model

            # Mock file writing
            with patch("builtins.open", mock_open()):
                commit_hash = manager.update_benchmark_version_and_commit(
                    mock_benchmark_execution, "1.5.0"
                )

                # Verify default commit message was used
                mock_index.commit.assert_called_once_with(
                    "Update benchmark version to 1.5.0"
                )
                assert commit_hash == "def456ghi789"

    def test_update_benchmark_version_and_commit_no_git(self, mock_git_repo):
        """Test updating benchmark version when git is not available."""
        from unittest.mock import Mock
        from omnibenchmark.versioning.exceptions import VersioningError

        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\nid: test\n")

        manager = GitAwareBenchmarkVersionManager(
            benchmark_path=benchmark_file, git_repo_path=mock_git_repo
        )

        # Simulate git not available
        manager.git_available = False
        manager.repo = None

        mock_benchmark_execution = Mock()
        mock_model = Mock()
        mock_benchmark_execution.model = mock_model

        with pytest.raises(VersioningError, match="Git repository not available"):
            manager.update_benchmark_version_and_commit(
                mock_benchmark_execution, "2.0.0"
            )

    def test_update_benchmark_version_and_commit_file_error(self, mock_git_repo):
        """Test handling file write errors during benchmark update."""
        from unittest.mock import Mock, patch
        from omnibenchmark.versioning.exceptions import VersioningError

        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\nid: test\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            mock_repo = MagicMock()
            MockRepo.return_value = mock_repo

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file, git_repo_path=mock_git_repo
            )

            mock_benchmark_execution = Mock()
            mock_model = Mock()
            mock_model.to_yaml.side_effect = Exception("YAML serialization error")
            mock_benchmark_execution.model = mock_model

            with pytest.raises(
                VersioningError, match="Failed to update benchmark and commit"
            ):
                manager.update_benchmark_version_and_commit(
                    mock_benchmark_execution, "2.0.0"
                )

    def test_create_version_with_persistence(self, mock_git_repo):
        """Test creating version with persistence (version creation + file update + commit)."""
        from unittest.mock import Mock, patch, mock_open

        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\nid: test\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            mock_repo = MagicMock()
            mock_index = MagicMock()
            mock_commit = MagicMock()
            mock_commit.hexsha = "persistent123"

            mock_repo.index = mock_index
            mock_index.commit.return_value = mock_commit
            MockRepo.return_value = mock_repo

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file,
                git_repo_path=mock_git_repo,
                lock_dir=mock_git_repo / "locks",
            )

            # Mock benchmark execution
            mock_benchmark_execution = Mock()
            mock_model = Mock()
            mock_model.version = "1.0.0"
            mock_model.to_yaml.return_value = "version: 2.1.0\nid: test\n"
            mock_benchmark_execution.model = mock_model

            # Mock file operations and version creation
            with patch("builtins.open", mock_open()):
                with patch.object(manager, "create_version") as mock_create:
                    mock_create.return_value = "2.1.0"

                    created_version = manager.create_version_with_persistence(
                        mock_benchmark_execution, "2.1.0", "Persistent version update"
                    )

                    # Verify version was created through base class
                    mock_create.assert_called_once()
                    call_args = mock_create.call_args
                    assert call_args[1]["version"] == "2.1.0"

                    # Verify hooks were called (pre_create_hook, post_create_hook)
                    assert "pre_create_hook" in call_args[1]
                    assert "post_create_hook" in call_args[1]

                    # Test the post_create_hook by calling it
                    post_hook = call_args[1]["post_create_hook"]
                    post_hook("2.1.0")

                    # Verify git operations happened in the hook
                    mock_index.add.assert_called_once()
                    mock_index.commit.assert_called_once_with(
                        "Persistent version update"
                    )

                    assert created_version == "2.1.0"

    def test_create_version_with_persistence_auto_version(self, mock_git_repo):
        """Test creating version with persistence using auto-increment."""
        from unittest.mock import Mock, patch, mock_open

        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\nid: test\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            mock_repo = MagicMock()
            mock_index = MagicMock()
            mock_commit = MagicMock()
            mock_commit.hexsha = "auto123"

            mock_repo.index = mock_index
            mock_index.commit.return_value = mock_commit
            MockRepo.return_value = mock_repo

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file,
                git_repo_path=mock_git_repo,
                lock_dir=mock_git_repo / "locks",
            )

            # Set up known versions for auto-increment
            manager.set_known_versions(["1.0.0", "1.1.0"])

            mock_benchmark_execution = Mock()
            mock_model = Mock()
            mock_model.to_yaml.return_value = "version: 1.2.0\nid: test\n"
            mock_benchmark_execution.model = mock_model

            with patch("builtins.open", mock_open()):
                with patch.object(manager, "create_version") as mock_create:
                    mock_create.return_value = "1.2.0"

                    created_version = manager.create_version_with_persistence(
                        mock_benchmark_execution
                    )

                    # Verify version was auto-created
                    mock_create.assert_called_once()
                    call_args = mock_create.call_args
                    assert call_args[1]["version"] is None  # Auto-increment

                    assert created_version == "1.2.0"
