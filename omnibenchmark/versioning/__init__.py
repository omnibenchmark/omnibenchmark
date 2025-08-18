"""
Versioning module for omnibenchmark.

Provides robust version management for benchmarks with clean separation of concerns:

Classes:
    Version: Semantic version representation with comparison and increment operations
    BenchmarkVersionManager: File-based locking and version validation
    GitAwareBenchmarkVersionManager: Extends with git commit tracking

Key Features:
    - Semantic versioning (x.y.z or x.y format)
    - Concurrency control via file-based locking
    - Optional git integration for version history
    - Hook system for storage integration
    - No model dependencies (architectural separation)

Example:
    Basic version management:
        manager = BenchmarkVersionManager(benchmark_path=Path("bench.yaml"))
        manager.create_version("1.0.0")

    With git tracking:
        git_manager = GitAwareBenchmarkVersionManager(benchmark_path=Path("bench.yaml"))
        git_manager.create_version_with_git_tracking("1.1.0")
"""

from .manager import BenchmarkVersionManager
from .git import GitAwareBenchmarkVersionManager
from .exceptions import (
    VersioningError,
    VersionDowngradeError,
    VersionLockError,
    VersionFormatError,
)
from .version import Version, parse_version, increment_version

__all__ = [
    "BenchmarkVersionManager",
    "GitAwareBenchmarkVersionManager",
    "VersioningError",
    "VersionDowngradeError",
    "VersionLockError",
    "VersionFormatError",
    "Version",
    "parse_version",
    "increment_version",
]
