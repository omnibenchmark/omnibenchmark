"""
Tests for GitAwareBenchmarkVersionManager.

All tests in this file are marked as 'short' since they mock git operations
and don't require actual git repositories or external dependencies.
"""

import subprocess
from unittest.mock import patch, MagicMock, Mock, mock_open

import pytest

from omnibenchmark.versioning.git import GitAwareBenchmarkVersionManager
from dulwich.errors import NotGitRepository


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _init_git_repo(path):
    """Initialize a real minimal git repo via subprocess."""
    subprocess.run(["git", "init", str(path)], check=True, capture_output=True)
    subprocess.run(
        ["git", "-C", str(path), "config", "user.email", "test@test.com"],
        check=True,
        capture_output=True,
    )
    subprocess.run(
        ["git", "-C", str(path), "config", "user.name", "Test"],
        check=True,
        capture_output=True,
    )


@pytest.fixture
def mock_git_repo(tmp_path):
    """Create a fake .git directory (not a real repo — used for patched tests)."""
    git_dir = tmp_path / ".git"
    git_dir.mkdir()
    return tmp_path


@pytest.fixture
def real_git_repo(tmp_path):
    """Create a real git repo."""
    _init_git_repo(tmp_path)
    return tmp_path


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestGitAwareBenchmarkVersionManager:
    """Test the GitAwareBenchmarkVersionManager class."""

    def test_initialization_with_git_repo(self, mock_git_repo):
        """Test initialization with a valid git repository."""
        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
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
            MockRepo.side_effect = NotGitRepository()

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file,
                git_repo_path=tmp_path,
                lock_dir=tmp_path / "locks",
            )

            assert manager.git_available is False

    def test_reconstruct_version_history(self, mock_git_repo):
        """Test reconstructing version history from git."""
        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            mock_repo = MagicMock()
            MockRepo.return_value = mock_repo

            # Mock walker entries (dulwich API: entry.commit.id)
            entry1 = MagicMock()
            entry1.commit.id = b"abc123"
            entry2 = MagicMock()
            entry2.commit.id = b"def456"
            mock_repo.get_walker.return_value = [entry2, entry1]  # newest first

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file, git_repo_path=mock_git_repo
            )

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
            MockRepo.side_effect = NotGitRepository()

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

            # Dulwich API: repo.head() returns bytes sha; repo[sha] returns commit obj
            mock_repo.head.return_value = b"abc123def456abc1"
            mock_commit = MagicMock()
            mock_commit.id = b"abc123def456abc1"
            mock_commit.author = b"John Doe <john@example.com>"
            mock_commit.author_time = 1704067200
            mock_repo.__getitem__ = MagicMock(return_value=mock_commit)

            # Branch via symrefs
            mock_repo.refs.get_symrefs.return_value = {b"HEAD": b"refs/heads/main"}

            # Config for remote url
            mock_config = MagicMock()
            mock_config.get.return_value = b"https://github.com/example/repo.git"
            mock_repo.get_config.return_value = mock_config

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file, git_repo_path=mock_git_repo
            )

            with patch("subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(stdout="", returncode=0)
                info = manager.get_current_git_info()

            assert info["commit"] == "abc123de"
            assert info["commit_full"] == "abc123def456abc1"
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

            mock_repo.head.return_value = b"abc123"
            mock_commit = MagicMock()
            mock_commit.id = b"abc123"
            mock_commit.author = b"Test <test@test.com>"
            mock_commit.author_time = 0
            mock_repo.__getitem__ = MagicMock(return_value=mock_commit)
            mock_repo.refs.get_symrefs.return_value = {}
            mock_repo.get_config.side_effect = Exception("no config")

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file, git_repo_path=mock_git_repo
            )

            with patch("subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(
                    stdout=" M dirty_file.py", returncode=0
                )
                info = manager.get_current_git_info()

            assert info["clean"] is False

    def test_get_current_git_info_no_git(self, tmp_path):
        """Test getting git info when git is not available."""
        benchmark_file = tmp_path / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            MockRepo.side_effect = NotGitRepository()

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
            MockRepo.return_value = MagicMock()

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file,
                git_repo_path=mock_git_repo,
                lock_dir=mock_git_repo / "locks",
            )

            with patch.object(manager, "get_current_git_info", return_value={}):
                version = manager.create_version_with_git_tracking("1.0.0")

            assert version == "1.0.0"
            assert "1.0.0" in manager.get_versions()

    def test_create_version_with_git_tracking_auto_increment(self, mock_git_repo):
        """Test creating auto-incremented version with git tracking."""
        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            MockRepo.return_value = MagicMock()

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file,
                git_repo_path=mock_git_repo,
                lock_dir=mock_git_repo / "locks",
            )

            with patch.object(manager, "get_current_git_info", return_value={}):
                version = manager.create_version_with_git_tracking()

            assert version == "0.1"

    def test_extract_version_from_yaml(self, mock_git_repo):
        """Test extracting version from YAML content."""
        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            MockRepo.return_value = MagicMock()

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file, git_repo_path=mock_git_repo
            )

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

            # Invalid YAML returns None
            version = manager._extract_version_from_yaml("invalid: yaml: content:")
            assert version is None

    def test_get_version_at_commit(self, mock_git_repo):
        """Test getting version from a specific commit."""
        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            MockRepo.return_value = MagicMock()

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

                    # Content not found → None
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

            entry1 = MagicMock()
            entry1.commit.id = b"abc123def"
            entry2 = MagicMock()
            entry2.commit.id = b"def456ghi"
            mock_repo.get_walker.return_value = [entry1, entry2]

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
            MockRepo.return_value = MagicMock()

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file,
                git_repo_path=mock_git_repo,
                lock_dir=mock_git_repo / "locks",
            )

            assert manager.get_versions() == []
            assert manager._current_version is None

            manager.create_version("1.0.0")
            assert manager.version_exists("1.0.0")
            assert manager.compare_versions("1.0.0", "2.0.0") == -1

    def test_detached_head_state(self, mock_git_repo):
        """Test handling of detached HEAD state (no refs/heads/ prefix)."""
        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            mock_repo = MagicMock()
            MockRepo.return_value = mock_repo

            mock_repo.head.return_value = b"abc123"
            mock_commit = MagicMock()
            mock_commit.id = b"abc123"
            mock_commit.author = b"Test <t@t.com>"
            mock_commit.author_time = 0
            mock_repo.__getitem__ = MagicMock(return_value=mock_commit)
            # HEAD points directly to a commit hash (detached)
            mock_repo.refs.get_symrefs.return_value = {b"HEAD": b"abc123"}
            mock_repo.get_config.side_effect = Exception("no config")

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file, git_repo_path=mock_git_repo
            )

            with patch("subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(stdout="", returncode=0)
                info = manager.get_current_git_info()

            # No branch key when ref doesn't start with refs/heads/
            assert "branch" not in info
            assert info["commit"] == "abc123"[:8]

    def test_initialize_from_git_history(self, mock_git_repo):
        """Test initializing known versions from git history."""
        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            MockRepo.return_value = MagicMock()

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
        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\nid: test\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            MockRepo.return_value = MagicMock()

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file, git_repo_path=mock_git_repo
            )

            mock_benchmark_execution = Mock()
            mock_model = Mock()
            mock_model.version = "1.0.0"
            mock_model.to_yaml.return_value = "version: 2.0.0\nid: test\n"
            mock_benchmark_execution.model = mock_model

            with patch("omnibenchmark.versioning.git.porcelain") as mock_porcelain:
                mock_porcelain.add = MagicMock()
                mock_porcelain.commit = MagicMock(return_value=b"abc123def456")

                with patch("builtins.open", mock_open()) as mock_file:
                    commit_hash = manager.update_benchmark_version_and_commit(
                        mock_benchmark_execution, "2.0.0", "Custom commit message"
                    )

                assert mock_model.version == "2.0.0"
                mock_file.assert_called_once_with(benchmark_file, "w")
                mock_porcelain.add.assert_called_once()
                mock_porcelain.commit.assert_called_once_with(
                    manager.repo, message=b"Custom commit message"
                )
                assert commit_hash == "abc123def456"

    def test_update_benchmark_version_and_commit_default_message(self, mock_git_repo):
        """Test updating benchmark version with default commit message."""
        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\nid: test\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            MockRepo.return_value = MagicMock()

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file, git_repo_path=mock_git_repo
            )

            mock_benchmark_execution = Mock()
            mock_model = Mock()
            mock_model.version = "1.0.0"
            mock_model.to_yaml.return_value = "version: 1.5.0\nid: test\n"
            mock_benchmark_execution.model = mock_model

            with patch("omnibenchmark.versioning.git.porcelain") as mock_porcelain:
                mock_porcelain.add = MagicMock()
                mock_porcelain.commit = MagicMock(return_value=b"def456ghi789")

                with patch("builtins.open", mock_open()):
                    commit_hash = manager.update_benchmark_version_and_commit(
                        mock_benchmark_execution, "1.5.0"
                    )

                mock_porcelain.commit.assert_called_once_with(
                    manager.repo, message=b"Update benchmark version to 1.5.0"
                )
                assert commit_hash == "def456ghi789"

    def test_update_benchmark_version_and_commit_no_git(self, mock_git_repo):
        """Test updating benchmark version when git is not available."""
        from omnibenchmark.versioning.exceptions import VersioningError

        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\nid: test\n")

        manager = GitAwareBenchmarkVersionManager(
            benchmark_path=benchmark_file, git_repo_path=mock_git_repo
        )
        manager.git_available = False
        manager.repo = None

        mock_benchmark_execution = Mock()
        mock_benchmark_execution.model = Mock()

        with pytest.raises(VersioningError, match="Git repository not available"):
            manager.update_benchmark_version_and_commit(
                mock_benchmark_execution, "2.0.0"
            )

    def test_update_benchmark_version_and_commit_file_error(self, mock_git_repo):
        """Test handling file write errors during benchmark update."""
        from omnibenchmark.versioning.exceptions import VersioningError

        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\nid: test\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            MockRepo.return_value = MagicMock()

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
        """Test creating version with persistence."""
        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\nid: test\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            MockRepo.return_value = MagicMock()

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file,
                git_repo_path=mock_git_repo,
                lock_dir=mock_git_repo / "locks",
            )

            mock_benchmark_execution = Mock()
            mock_model = Mock()
            mock_model.version = "1.0.0"
            mock_model.to_yaml.return_value = "version: 2.1.0\nid: test\n"
            mock_benchmark_execution.model = mock_model

            with patch("omnibenchmark.versioning.git.porcelain") as mock_porcelain:
                mock_porcelain.add = MagicMock()
                mock_porcelain.commit = MagicMock(return_value=b"persistent123")

                with patch("builtins.open", mock_open()):
                    with patch.object(manager, "create_version") as mock_create:
                        mock_create.return_value = "2.1.0"

                        created_version = manager.create_version_with_persistence(
                            mock_benchmark_execution,
                            "2.1.0",
                            "Persistent version update",
                        )

                        mock_create.assert_called_once()
                        call_args = mock_create.call_args
                        assert call_args[1]["version"] == "2.1.0"
                        assert "pre_create_hook" in call_args[1]
                        assert "post_create_hook" in call_args[1]

                        # Invoke post_create_hook to verify git calls happen inside it
                        post_hook = call_args[1]["post_create_hook"]
                        post_hook("2.1.0")

                        mock_porcelain.add.assert_called_once()
                        mock_porcelain.commit.assert_called_once_with(
                            manager.repo, message=b"Persistent version update"
                        )
                        assert created_version == "2.1.0"

    def test_create_version_with_persistence_auto_version(self, mock_git_repo):
        """Test creating version with persistence using auto-increment."""
        benchmark_file = mock_git_repo / "test.yaml"
        benchmark_file.write_text("version: 1.0.0\nid: test\n")

        with patch("omnibenchmark.versioning.git.Repo") as MockRepo:
            MockRepo.return_value = MagicMock()

            manager = GitAwareBenchmarkVersionManager(
                benchmark_path=benchmark_file,
                git_repo_path=mock_git_repo,
                lock_dir=mock_git_repo / "locks",
            )
            manager.set_known_versions(["1.0.0", "1.1.0"])

            mock_benchmark_execution = Mock()
            mock_model = Mock()
            mock_model.to_yaml.return_value = "version: 1.2.0\nid: test\n"
            mock_benchmark_execution.model = mock_model

            with patch("omnibenchmark.versioning.git.porcelain"):
                with patch("builtins.open", mock_open()):
                    with patch.object(manager, "create_version") as mock_create:
                        mock_create.return_value = "1.2.0"

                        created_version = manager.create_version_with_persistence(
                            mock_benchmark_execution
                        )

                        call_args = mock_create.call_args
                        assert call_args[1]["version"] is None  # Auto-increment
                        assert created_version == "1.2.0"
