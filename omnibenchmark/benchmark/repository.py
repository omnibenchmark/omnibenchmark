"""Shared repository utilities for benchmark module validation and citation extraction.

This module provides common functionality for cloning, accessing, and cleaning up
module repositories used by both validation and citation extraction workflows.
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Tuple

from omnibenchmark.git.cache import clone_module_v2, parse_repo_url

logger = logging.getLogger(__name__)


class RepositoryManager:
    """Manages repository access and cleanup for benchmark modules."""

    def __init__(self):
        """Initialize repository manager."""
        pass

    def get_local_repo_path(self, repo_url: str, commit_hash: str) -> Path:
        """Get the expected local repository path.

        Args:
            repo_url: Repository URL
            commit_hash: Git commit hash

        Returns:
            Path to local repository directory
        """
        repo_name = Path(parse_repo_url(repo_url)).name
        commit_prefix = commit_hash[:7] if commit_hash else "unknown"
        return Path(".modules") / repo_name / commit_prefix

    def clone_to_temp(
        self, repo_url: str, commit_hash: str, module_id: str
    ) -> Optional[Path]:
        """Clone repository to work directory using the central cache.

        Args:
            repo_url: Repository URL to clone
            commit_hash: Git commit hash or tag
            module_id: Module identifier for logging

        Returns:
            Path to repository work directory, or None if clone failed
        """
        try:
            repo_name = Path(parse_repo_url(repo_url)).name
            commit_prefix = commit_hash[:7] if commit_hash else "unknown"
            work_dir = Path(".modules") / repo_name / commit_prefix
            cloned_path, _ = clone_module_v2(repo_url, commit_hash, work_dir=work_dir)
            logger.debug(
                f"Successfully cloned {repo_url}@{commit_hash} to {cloned_path}"
            )
            return cloned_path

        except Exception as e:
            logger.warning(f"Failed to clone repository {repo_url}@{commit_hash}: {e}")
            return None

    def get_repository_files(self, repo_path: Path) -> Dict[str, Optional[str]]:
        """Read standard repository files.

        Args:
            repo_path: Path to repository directory

        Returns:
            Dictionary with file contents (None if file doesn't exist)
        """
        files = {}

        # Standard files to check
        file_paths = {
            "citation": repo_path / "CITATION.cff",
            "license": repo_path / "LICENSE",
            "omnibenchmark": repo_path / "omnibenchmark.yaml",
        }

        for file_key, file_path in file_paths.items():
            if file_path.exists():
                try:
                    files[file_key] = file_path.read_text(encoding="utf-8")
                except Exception as e:
                    logger.warning(f"Failed to read {file_path}: {e}")
                    files[file_key] = None
            else:
                files[file_key] = None

        return files

    def get_files_present(self, repo_path: Path) -> Dict[str, bool]:
        """Check which standard files are present in repository.

        Args:
            repo_path: Path to repository directory

        Returns:
            Dictionary indicating file presence
        """
        return {
            "CITATION.cff": (repo_path / "CITATION.cff").exists(),
            "LICENSE": (repo_path / "LICENSE").exists(),
            "LICENSE.txt": (repo_path / "LICENSE.txt").exists(),
            "LICENSE.md": (repo_path / "LICENSE.md").exists(),
            "COPYING": (repo_path / "COPYING").exists(),
            "COPYING.txt": (repo_path / "COPYING.txt").exists(),
            "omnibenchmark.yaml": (repo_path / "omnibenchmark.yaml").exists(),
            "config.cfg": (repo_path / "config.cfg").exists(),
        }

    def cleanup_temp_directories(self):
        """No-op: module work directories under .modules/ are persistent cache."""
        pass

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit with cleanup."""
        self.cleanup_temp_directories()


def cleanup_temp_repositories():
    """No-op: module work directories under .modules/ are persistent cache."""
    pass


def get_module_repository_info(
    benchmark, module
) -> Tuple[Optional[str], Optional[str]]:
    """Extract repository information from a module.

    Args:
        benchmark: Benchmark instance
        module: Module instance

    Returns:
        Tuple of (repository_url, commit_hash) or (None, None) if not found
    """
    repo_info = benchmark.model.get_module_repository(module)

    repo_url = getattr(repo_info, "url", None) or getattr(repo_info, "repository", None)
    commit_hash = (
        getattr(repo_info, "commit_hash", None)
        or getattr(repo_info, "commit", None)
        or getattr(repo_info, "version", None)
    )

    return repo_url, commit_hash


def resolve_module_repository(
    benchmark, module, module_id: str, repo_manager: RepositoryManager
) -> Optional[Path]:
    """Resolve the path to a module's repository, cloning if necessary.

    Args:
        benchmark: Benchmark instance
        module: Module instance
        module_id: Module identifier
        repo_manager: Repository manager instance

    Returns:
        Path to repository directory, or None if not accessible
    """
    # Get repository information
    repo_url, commit_hash = get_module_repository_info(benchmark, module)

    if not repo_url or not commit_hash:
        logger.warning(f"Module {module_id}: Missing repository information")
        return None

    # Check for local repository first
    local_repo_path = repo_manager.get_local_repo_path(repo_url, commit_hash)

    if local_repo_path.exists():
        logger.debug(f"Using local repository: {local_repo_path}")
        return local_repo_path

    # Try to clone to temporary directory
    logger.debug(f"Local repository not found for {module_id}, attempting to clone...")
    temp_repo_path = repo_manager.clone_to_temp(repo_url, commit_hash, module_id)

    if temp_repo_path is None:
        logger.warning(f"Failed to clone repository for {module_id}: {repo_url}")
        return None

    return temp_repo_path
