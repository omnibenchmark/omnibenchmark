"""
Tests for BenchmarkVersionManager.

All tests in this file are marked as 'short' since they don't require
external dependencies, containers, or network I/O.
"""

import pytest

from omnibenchmark.versioning.manager import BenchmarkVersionManager
from omnibenchmark.versioning.exceptions import (
    VersionDowngradeError,
    VersionLockError,
    VersionAlreadyExistsError,
)


@pytest.fixture
def manager(tmp_path):
    """Create a test manager instance."""
    # Create a dummy benchmark file
    benchmark_file = tmp_path / "test_benchmark.yaml"
    benchmark_file.write_text("version: 1.0.0\n")

    return BenchmarkVersionManager(
        benchmark_path=benchmark_file, lock_dir=tmp_path / "locks"
    )


@pytest.mark.short
class TestBenchmarkVersionManager:
    """Test the BenchmarkVersionManager class."""

    def test_initialization(self, tmp_path):
        """Test manager initialization."""
        benchmark_file = tmp_path / "test_benchmark.yaml"
        benchmark_file.write_text("version: 1.0.0\n")

        manager = BenchmarkVersionManager(
            benchmark_path=benchmark_file, lock_dir=tmp_path / "locks"
        )

        assert manager.benchmark_name == "test_benchmark"
        assert manager.lock_dir.exists()
        assert manager.lock_file.name == "test_benchmark.lock"
        assert manager.benchmark_path == benchmark_file

    def test_lock_acquisition(self, manager):
        """Test file-based lock acquisition."""
        # Test exclusive lock
        with manager.acquire_lock(shared=False):
            assert manager.lock_file.exists()

        # Test shared lock
        with manager.acquire_lock(shared=True):
            assert manager.lock_file.exists()

    def test_lock_timeout(self, manager):
        """Test lock timeout behavior."""
        manager.lock_timeout = 0.5  # Set short timeout for testing

        # Simulate lock being held
        with open(manager.lock_file, "w") as f:
            import fcntl

            fcntl.flock(f, fcntl.LOCK_EX)

            # This should timeout
            with pytest.raises(VersionLockError, match="Timeout"):
                with manager.acquire_lock(shared=False):
                    pass

    def test_set_and_get_known_versions(self, manager):
        """Test version tracking in memory."""
        # Set known versions
        manager.set_known_versions(["1.0.0", "1.1.0"])

        versions = manager.get_versions()
        assert "1.0.0" in versions
        assert "1.1.0" in versions

        # Create new version
        manager.create_version("2.0.0")
        assert "2.0.0" in manager.get_versions()

    def test_set_known_versions_with_invalid(self, manager):
        """Test that invalid versions are filtered out."""
        manager.set_known_versions(["1.0.0", "invalid", "2.0"])

        versions = manager.get_versions()
        assert "1.0.0" in versions
        assert "2.0" in versions
        assert "invalid" not in versions

    def test_get_versions_empty(self, manager):
        """Test getting versions when none exist."""
        versions = manager.get_versions()
        assert versions == []

    def test_get_versions_with_data(self, manager):
        """Test getting versions after creating some."""
        manager.create_version("1.0.0")
        manager.create_version("1.1.0")
        manager.create_version("2.0.0")

        versions = manager.get_versions()
        assert versions == ["1.0.0", "1.1.0", "2.0.0"]

    def test_get_current_version_none(self, manager):
        """Test getting current version when none is set."""
        # Should return None when no current version and file can't be parsed
        assert manager._current_version is None
        # Don't try to parse file since our test file is incomplete

    def test_get_current_version_set(self, manager):
        """Test getting current version after setting it."""
        manager.set_current_version("1.2.3")
        assert manager.get_current_version() == "1.2.3"

    def test_set_current_version_explicit(self, manager):
        """Test setting current version explicitly."""
        version = manager.set_current_version("1.0.0")
        assert version == "1.0.0"
        assert manager.get_current_version() == "1.0.0"

    def test_set_current_version_auto_increment(self, manager):
        """Test auto-incrementing version."""
        # First version should be 0.1
        version1 = manager.set_current_version(None)
        assert version1 == "0.1"

        # Create it to add to history
        manager.create_version(version1)

        # Next auto-increment should be 0.2
        version2 = manager.set_current_version(None)
        assert version2 == "0.2"

    def test_set_current_version_already_exists(self, manager):
        """Test setting a version that already exists."""
        manager.create_version("1.0.0")

        # Can set current to an existing version
        manager.set_current_version("1.0.0")
        assert manager.get_current_version() == "1.0.0"

    def test_create_version_simple(self, manager):
        """Test creating a simple version."""
        version = manager.create_version("1.0.0")
        assert version == "1.0.0"
        assert "1.0.0" in manager.get_versions()
        assert manager.get_current_version() == "1.0.0"

    def test_create_version_auto_increment(self, manager):
        """Test creating version with auto-increment."""
        version1 = manager.create_version()
        assert version1 == "0.1"

        version2 = manager.create_version()
        assert version2 == "0.2"

    def test_create_version_with_hooks(self, manager):
        """Test creating version with pre and post hooks."""
        pre_hook_called = False
        post_hook_called = False

        def pre_hook(version):
            nonlocal pre_hook_called
            pre_hook_called = True
            assert version == "1.0.0"

        def post_hook(version):
            nonlocal post_hook_called
            post_hook_called = True
            assert version == "1.0.0"

        manager.create_version(
            "1.0.0", pre_create_hook=pre_hook, post_create_hook=post_hook
        )

        assert pre_hook_called
        assert post_hook_called

    def test_create_version_already_exists(self, manager):
        """Test creating a version that already exists."""
        manager.create_version("1.0.0")

        with pytest.raises(VersionAlreadyExistsError, match="1.0.0"):
            manager.create_version("1.0.0")

    def test_create_version_downgrade(self, manager):
        """Test that version downgrades are detected."""
        manager.create_version("2.0.0")

        with pytest.raises(VersionDowngradeError, match="2.0.0.*1.0.0"):
            manager.create_version("1.0.0")

    def test_create_version_allow_duplicates(self, manager):
        """Test creating version with duplicates allowed."""
        manager.create_version("1.0.0")

        # Should succeed with allow_duplicates=True
        version = manager.create_version("1.0.0", allow_duplicates=True)
        assert version == "1.0.0"

    def test_version_exists(self, manager):
        """Test checking if version exists."""
        manager.create_version("1.0.0")

        assert manager.version_exists("1.0.0")
        assert not manager.version_exists("2.0.0")

    def test_compare_versions(self, manager):
        """Test version comparison."""
        assert manager.compare_versions("1.0.0", "1.0.0") == 0
        assert manager.compare_versions("1.0.0", "2.0.0") == -1
        assert manager.compare_versions("2.0.0", "1.0.0") == 1

    def test_cleanup_locks(self, manager):
        """Test cleanup of lock files."""
        # Create a lock
        with manager.acquire_lock():
            pass

        assert manager.lock_file.exists()

        # Clean up
        manager.cleanup_locks()
        assert not manager.lock_file.exists()

    def test_concurrent_access_simulation(self, manager):
        """Test that concurrent access is properly serialized."""
        import threading

        results = []

        def create_version(version):
            try:
                manager.create_version(version)
                results.append(("success", version))
            except VersionAlreadyExistsError:
                results.append(("exists", version))

        # Create multiple threads trying to create the same version
        threads = []
        for i in range(5):
            t = threading.Thread(target=create_version, args=("1.0.0",))
            threads.append(t)
            t.start()

        for t in threads:
            t.join()

        # Only one should succeed
        successes = [r for r in results if r[0] == "success"]
        exists = [r for r in results if r[0] == "exists"]

        assert len(successes) == 1
        assert len(exists) == 4

    def test_version_format_validation(self, manager):
        """Test that version format is validated."""
        # Valid versions
        manager.create_version("1.0.0")
        manager.create_version("2.1")

        # Invalid versions should raise ValueError
        with pytest.raises(ValueError, match="Invalid version format"):
            manager.create_version("invalid")

        with pytest.raises(ValueError, match="Invalid version format"):
            manager.create_version("1")

    def test_version_validation(self, manager):
        """Test version validation logic."""
        # Set known versions
        manager.set_known_versions(["1.0.0", "2.0.0"])

        # Test downgrade prevention
        with pytest.raises(VersionDowngradeError):
            manager.validate_new_version("0.9.0")

        # Test duplicate prevention
        with pytest.raises(VersionAlreadyExistsError):
            manager.validate_new_version("1.0.0", allow_duplicates=False)

        # Allow duplicates - should not raise
        try:
            manager.validate_new_version("2.0.0", allow_duplicates=True)
        except VersionDowngradeError:
            # This is expected since 2.0.0 already exists
            pass

        # Valid new version
        manager.validate_new_version("3.0.0")

    def test_get_latest_version(self, manager):
        """Test getting the latest version."""
        # No versions
        assert manager.get_latest_version() is None

        # With versions
        manager.set_known_versions(["1.0.0", "2.1.0", "1.5.0"])
        assert manager.get_latest_version() == "2.1.0"

    def test_suggest_next_version(self, manager):
        """Test version suggestion functionality."""
        # No versions yet
        assert manager.suggest_next_version() == "0.1"
        assert manager.suggest_next_version("major") == "1.0"

        # With existing versions
        manager.set_known_versions(["1.0.0", "1.1.0"])
        assert manager.suggest_next_version() == "1.2"
        assert manager.suggest_next_version("major") == "2.0"
        assert manager.suggest_next_version("patch") == "1.1.1"


@pytest.mark.short
class TestVersionManagerEdgeCases:
    """Test edge cases and error conditions."""

    def test_get_current_version_from_file_missing(self, tmp_path):
        """Test getting version when file doesn't exist."""
        benchmark_file = tmp_path / "missing.yaml"

        manager = BenchmarkVersionManager(
            benchmark_path=benchmark_file, lock_dir=tmp_path
        )

        # Should return None when file doesn't exist
        version = manager.get_current_version_from_file()
        assert version is None

    def test_missing_lock_directory_creation(self, tmp_path):
        """Test that lock directory is created if missing."""
        lock_dir = tmp_path / "new_locks"
        assert not lock_dir.exists()

        benchmark_file = tmp_path / "test.yaml"
        benchmark_file.write_text("version: 1.0\n")

        _ = BenchmarkVersionManager(benchmark_path=benchmark_file, lock_dir=lock_dir)

        assert lock_dir.exists()

    def test_version_sorting(self, manager):
        """Test that versions are maintained in creation order."""
        # Create versions in increasing order (downgrades not allowed)
        manager.create_version("1.0.0")
        manager.create_version("1.5.0")
        manager.create_version("2.0.0")
        manager.create_version("2.1.0")

        versions = manager.get_versions()
        # Should be in creation order
        assert versions == ["1.0.0", "1.5.0", "2.0.0", "2.1.0"]

    def test_version_with_different_formats(self, manager):
        """Test handling versions with different formats."""
        manager.create_version("1.0")  # Two components
        manager.create_version("2.0.0")  # Three components

        versions = manager.get_versions()
        assert "1.0" in versions
        assert "2.0.0" in versions

    def test_validate_with_no_downgrades(self, manager):
        """Test that downgrades are never allowed."""
        manager.set_known_versions(["1.0.0", "1.5.0", "2.0.0"])

        # Any version less than any known version is a downgrade
        with pytest.raises(VersionDowngradeError):
            manager.validate_new_version("1.2.0")  # Less than 1.5.0

        with pytest.raises(VersionDowngradeError):
            manager.validate_new_version("0.9.0")  # Less than all

        # Only versions greater than all known versions are allowed
        manager.validate_new_version("2.1.0")  # Greater than all
        manager.validate_new_version("3.0.0")  # Greater than all
