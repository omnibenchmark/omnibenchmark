"""Shared repository utilities for benchmark module validation and citation extraction.

This module provides common functionality for cloning, accessing, and cleaning up
module repositories used by both validation and citation extraction workflows.
"""

import logging
import shutil
import tempfile
import glob
from pathlib import Path
from typing import Optional, List, Dict, Tuple

from omnibenchmark.io.code import clone_module
from omnibenchmark.model.repo import get_repo_hash

logger = logging.getLogger(__name__)


class RepositoryManager:
    """Manages repository access and cleanup for benchmark modules."""

    def __init__(self, prefix: str = "omnibenchmark"):
        """Initialize repository manager.

        Args:
            prefix: Prefix for temporary directory names
        """
        self.prefix = prefix
        self.temp_directories: List[Path] = []

    def get_local_repo_path(self, repo_url: str, commit_hash: str) -> Path:
        """Get the expected local repository path.

        Args:
            repo_url: Repository URL
            commit_hash: Git commit hash

        Returns:
            Path to local repository directory
        """
        repos_base_dir = Path(".snakemake/repos")
        folder_name = get_repo_hash(repo_url, commit_hash)
        return repos_base_dir / folder_name

    def clone_to_temp(
        self, repo_url: str, commit_hash: str, module_id: str
    ) -> Optional[Path]:
        """Clone repository to temporary directory.

        Args:
            repo_url: Repository URL to clone
            commit_hash: Git commit hash or tag
            module_id: Module identifier for logging

        Returns:
            Path to temporary repository directory, or None if clone failed
        """
        try:
            # Create temporary directory
            temp_dir = Path(tempfile.mkdtemp(prefix=f"{self.prefix}_{module_id}_"))
            logger.debug(f"Created temporary directory: {temp_dir}")

            # Clone the repository
            cloned_path = clone_module(temp_dir, repo_url, commit_hash)
            logger.info(
                f"Successfully cloned {repo_url}@{commit_hash} to temporary location"
            )

            # Track for cleanup
            self.temp_directories.append(temp_dir.parent)
            return cloned_path

        except Exception as e:
            logger.warning(f"Failed to clone repository {repo_url}@{commit_hash}: {e}")
            # Clean up temp directory if it was created
            if "temp_dir" in locals() and temp_dir.exists():
                try:
                    shutil.rmtree(temp_dir)
                except Exception as cleanup_error:
                    logger.debug(
                        f"Failed to cleanup temp directory {temp_dir}: {cleanup_error}"
                    )
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
        }

    def cleanup_temp_directories(self):
        """Clean up tracked temporary directories."""
        for temp_dir in self.temp_directories:
            try:
                if temp_dir.exists():
                    shutil.rmtree(temp_dir)
                    logger.debug(f"Cleaned up temporary directory: {temp_dir}")
            except Exception as e:
                logger.warning(f"Failed to cleanup temporary directory {temp_dir}: {e}")

        self.temp_directories.clear()

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit with cleanup."""
        self.cleanup_temp_directories()


def cleanup_temp_repositories(prefix: str = "omnibenchmark"):
    """Clean up any temporary repositories created during processing.

    This function removes temporary directories that may have been left behind
    due to interrupted execution or exceptions.

    Args:
        prefix: Prefix pattern to match for cleanup
    """
    temp_pattern = f"{tempfile.gettempdir()}/{prefix}_*"

    for temp_path in glob.glob(temp_pattern):
        try:
            temp_path_obj = Path(temp_path)
            if temp_path_obj.is_dir():
                shutil.rmtree(temp_path_obj)
                logger.debug(
                    f"Cleaned up leftover temporary directory: {temp_path_obj}"
                )
        except Exception as e:
            logger.warning(f"Failed to cleanup leftover directory {temp_path}: {e}")


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
    repo_info = benchmark.get_converter().get_module_repository(module)

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
    logger.info(f"Local repository not found for {module_id}, attempting to clone...")
    temp_repo_path = repo_manager.clone_to_temp(repo_url, commit_hash, module_id)

    if temp_repo_path is None:
        logger.warning(f"Failed to clone repository for {module_id}: {repo_url}")
        return None

    return temp_repo_path
