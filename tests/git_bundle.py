import subprocess

from pathlib import Path
from typing import Optional

import shutil

import logging

logger = logging.getLogger(__name__)


class GitBundleManager:
    """Helper class to manage Git bundle operations for testing."""

    def __init__(self, bundles_dir: Path, repos_dir: Path):
        self.bundles_dir = bundles_dir.resolve()
        self.repos_dir = repos_dir.resolve()
        self.repos_dir.mkdir(parents=True, exist_ok=True)

    def extract_bundle(self, bundle_name: str) -> Optional[Path]:
        """Extract a Git bundle to a repository directory.

        Args:
            bundle_name: Name of the bundle file (with or without .bundle extension)

        Returns:
            Path to the extracted repository or None if failed
        """
        # Normalize bundle name
        if not bundle_name.endswith(".bundle"):
            bundle_name += ".bundle"

        bundle_path = self.bundles_dir / bundle_name

        if not bundle_path.exists():
            logger.error(f"Bundle '{bundle_name}' not found")
            return None

        # Extract commit hash from filename (format: name_commit.bundle)
        commit_hash = bundle_path.stem.split("_")[-1]
        repo_name = bundle_path.stem.replace(f"_{commit_hash}", "")
        repo_path = self.repos_dir / f"{repo_name}.git"

        if repo_path.exists():
            shutil.rmtree(repo_path)

        try:
            # Clone from bundle
            subprocess.run(
                ["git", "clone", bundle_path, repo_path],
                check=True,
                capture_output=True,
                text=True,
            )

            # Checkout the specific commit
            subprocess.run(
                ["git", "checkout", commit_hash],
                cwd=repo_path,
                check=True,
                capture_output=True,
                text=True,
            )

            return repo_path
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to extract bundle {bundle_name}: {e.stderr}")
            return None
