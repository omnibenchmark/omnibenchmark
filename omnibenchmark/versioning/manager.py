"""
Benchmark version manager with file-based locking for concurrency control.

This module provides version management without persisting version history.
Version history should be reconstructed from git history of YAML files.
"""

import fcntl
import time
from contextlib import contextmanager
from pathlib import Path
from typing import List, Optional, Callable
import tempfile

from .exceptions import (
    VersioningError,
    VersionDowngradeError,
    VersionLockError,
    VersionAlreadyExistsError,
)
from .version import Version, increment_version


class BenchmarkVersionManager:
    """
    Manages benchmark versions with file-based locking for concurrency control.

    This class provides:
    - Version tracking and validation
    - File-based locking for concurrent access
    - Version comparison and validation
    - Hooks for storage integration

    Note: Version history is NOT persisted by this class. It should be
    reconstructed from git history or other external sources.
    """

    def __init__(
        self,
        benchmark_path: Path,
        lock_dir: Optional[Path] = None,
        lock_timeout: float = 30.0,
    ):
        """
        Initialize the version manager.

        Args:
            benchmark_path: Path to the benchmark YAML file
            lock_dir: Directory for lock files (defaults to system temp)
            lock_timeout: Maximum time to wait for lock acquisition (seconds)
        """
        self.benchmark_path = Path(benchmark_path)
        self.benchmark_name = self.benchmark_path.stem
        self.lock_timeout = lock_timeout

        # Set up lock directory (transient, not persisted)
        if lock_dir is None:
            # Use system temp directory for transient locks
            lock_dir = Path(tempfile.gettempdir()) / "omnibenchmark_locks"
        self.lock_dir = Path(lock_dir)
        self.lock_dir.mkdir(parents=True, exist_ok=True)

        # Lock file path (transient)
        self.lock_file = self.lock_dir / f"{self.benchmark_name}.lock"

        # In-memory state (not persisted)
        self._current_version: Optional[Version] = None
        self._known_versions: List[str] = []

    @contextmanager
    def acquire_lock(self, shared: bool = False):
        """
        Acquire a file-based lock for concurrent access control.

        Args:
            shared: If True, acquire a shared (read) lock. If False, exclusive (write) lock.

        Yields:
            The lock file handle

        Raises:
            VersionLockError: If lock cannot be acquired within timeout
        """
        lock_handle = None
        start_time = time.time()

        try:
            # Open or create lock file
            lock_handle = open(self.lock_file, "a+")

            # Determine lock type
            lock_type = fcntl.LOCK_SH if shared else fcntl.LOCK_EX

            # Try to acquire lock with timeout
            while True:
                try:
                    fcntl.flock(lock_handle, lock_type | fcntl.LOCK_NB)
                    break
                except IOError:
                    if time.time() - start_time > self.lock_timeout:
                        raise VersionLockError(
                            str(self.lock_file),
                            f"Timeout after {self.lock_timeout} seconds",
                        )
                    time.sleep(0.1)

            yield lock_handle

        finally:
            if lock_handle:
                try:
                    fcntl.flock(lock_handle, fcntl.LOCK_UN)
                except Exception:
                    pass
                lock_handle.close()

    def get_current_version_from_file(self) -> Optional[str]:
        """
        Get the current version from the benchmark YAML file.

        Returns:
            Version string from the file or None if not found
        """
        if not self.benchmark_path.exists():
            return None

        try:
            # Import here to avoid circular dependency
            from omnibenchmark.model import Benchmark

            benchmark = Benchmark.from_yaml(self.benchmark_path)
            return benchmark.version
        except Exception as e:
            raise VersioningError(f"Failed to parse benchmark file: {e}")

    def set_known_versions(self, versions: List[str]) -> None:
        """
        Set the list of known versions.

        This should be called with versions reconstructed from external sources
        like git history or storage systems.

        Args:
            versions: List of version strings
        """
        # Validate all versions
        validated = []
        for v in versions:
            try:
                Version(v)  # Validate format
                validated.append(v)
            except ValueError:
                # Skip invalid versions
                continue

        self._known_versions = validated

    def get_versions(self) -> List[str]:
        """
        Get list of known versions.

        Note: This returns the in-memory list. For persistent version history,
        reconstruct from git or other external sources.

        Returns:
            List of version strings
        """
        return self._known_versions.copy()

    def get_current_version(self) -> Optional[str]:
        """
        Get the current version.

        Returns:
            Current version string or None if no version set
        """
        if self._current_version:
            return str(self._current_version)

        # Try to get from file if not set
        return self.get_current_version_from_file()

    def set_current_version(self, version: Optional[str] = None) -> str:
        """
        Set the current version.

        Args:
            version: Version to set. If None, auto-increment from known versions.

        Returns:
            The version that was set
        """
        if version is None:
            # Auto-increment from known versions
            if self._known_versions:
                # Sort versions to get the latest
                sorted_versions = sorted(self._known_versions, key=Version)
                last_version = sorted_versions[-1]
                version = increment_version(last_version, "minor")
            else:
                version = "0.1"

        # Parse to validate format
        v = Version(version)
        version_str = str(v)

        self._current_version = v
        return version_str

    def validate_new_version(
        self, version: str, allow_duplicates: bool = False
    ) -> None:
        """
        Validate that a new version is acceptable.
        No downgrades are allowed.

        Args:
            version: Version to validate
            allow_duplicates: If False, raise error if version already exists

        Raises:
            VersionAlreadyExistsError: If version already exists and duplicates not allowed
            VersionDowngradeError: If version would be a downgrade
            ValueError: If version format is invalid
        """
        # Parse and validate format
        new_version = Version(version)
        version_str = str(new_version)

        # Check if already exists
        if not allow_duplicates and version_str in self._known_versions:
            raise VersionAlreadyExistsError(version_str)

        # Check for downgrade (no downgrades allowed)
        if self._known_versions:
            for known in self._known_versions:
                known_version = Version(known)
                if new_version < known_version:
                    raise VersionDowngradeError(str(known_version), version_str)

    def create_version(
        self,
        version: Optional[str] = None,
        pre_create_hook: Optional[Callable[[str], None]] = None,
        post_create_hook: Optional[Callable[[str], None]] = None,
        allow_duplicates: bool = False,
    ) -> str:
        """
        Create a new version with optional hooks for storage integration.

        Args:
            version: Version to create. If None, auto-increment.
            pre_create_hook: Function to call before creating version
            post_create_hook: Function to call after creating version
            allow_duplicates: If True, allow creating existing versions

        Returns:
            The created version string

        Raises:
            VersionAlreadyExistsError: If version already exists and duplicates not allowed
            VersionDowngradeError: If version would be a downgrade (never allowed)
        """
        with self.acquire_lock(shared=False):
            # Determine version
            if version is None:
                if self._known_versions:
                    sorted_versions = sorted(self._known_versions, key=Version)
                    last_version = sorted_versions[-1]
                    version = increment_version(last_version, "minor")
                else:
                    version = "0.1"

            # Validate the new version
            self.validate_new_version(version, allow_duplicates=allow_duplicates)

            # Parse for consistent formatting
            new_version = Version(version)
            version_str = str(new_version)

            # Call pre-create hook
            if pre_create_hook:
                pre_create_hook(version_str)

            # Add to known versions (in-memory only)
            if version_str not in self._known_versions:
                self._known_versions.append(version_str)

            # Update current version
            self._current_version = new_version

            # Call post-create hook
            if post_create_hook:
                post_create_hook(version_str)

            return version_str

    def version_exists(self, version: str) -> bool:
        """
        Check if a version exists in the known versions.

        Args:
            version: Version to check

        Returns:
            True if version exists, False otherwise
        """
        return version in self._known_versions

    def compare_versions(self, version1: str, version2: str) -> int:
        """
        Compare two versions.

        Args:
            version1: First version
            version2: Second version

        Returns:
            -1 if version1 < version2
             0 if version1 == version2
             1 if version1 > version2
        """
        v1 = Version(version1)
        v2 = Version(version2)

        if v1 < v2:
            return -1
        elif v1 > v2:
            return 1
        else:
            return 0

    def cleanup_locks(self) -> None:
        """
        Clean up lock files for this benchmark.

        This should be called when done with version management
        to clean up transient lock resources.
        """
        try:
            if self.lock_file.exists():
                self.lock_file.unlink()
        except OSError:
            pass  # Lock file may be in use by another process

    def get_latest_version(self) -> Optional[str]:
        """
        Get the latest known version.

        Returns:
            Latest version string or None if no versions known
        """
        if not self._known_versions:
            return None

        sorted_versions = sorted(self._known_versions, key=Version)
        return sorted_versions[-1]

    def suggest_next_version(self, component: str = "minor") -> str:
        """
        Suggest the next version based on known versions.

        Args:
            component: Which component to increment ("major", "minor", or "patch")

        Returns:
            Suggested next version string
        """
        latest = self.get_latest_version()
        if latest:
            return increment_version(latest, component)
        else:
            # Default starting version
            if component == "major":
                return "1.0"
            else:
                return "0.1"
