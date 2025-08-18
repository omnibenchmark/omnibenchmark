"""
Git-aware benchmark version manager.

This module extends BenchmarkVersionManager with git integration
for reconstructing version history from git commits.

YAML parsing is delegated to the model layer via Benchmark.from_yaml().
"""

from pathlib import Path
from typing import List, Optional, Dict, Any
from git import Repo
from git.exc import InvalidGitRepositoryError

from .manager import BenchmarkVersionManager


class GitAwareBenchmarkVersionManager(BenchmarkVersionManager):
    """
    Extends BenchmarkVersionManager with git integration.

    This class adds:
    - Version history reconstruction from git commits
    - Git commit information tracking
    - Automatic version discovery from git history

    Versions are extracted from the benchmark YAML files in git history
    using the model's parsing capabilities.
    """

    def __init__(
        self, benchmark_path: Path, git_repo_path: Optional[Path] = None, **kwargs
    ):
        """
        Initialize the git-aware version manager.

        Args:
            benchmark_path: Path to the benchmark YAML file
            git_repo_path: Path to git repository (defaults to current directory)
            **kwargs: Additional arguments passed to BenchmarkVersionManager
        """
        super().__init__(benchmark_path, **kwargs)

        # Set git repo path
        if git_repo_path is None:
            git_repo_path = Path.cwd()
        self.git_repo_path = Path(git_repo_path)

        # Check if git is available
        self.git_available = False
        self.repo = None
        try:
            self.repo = Repo(self.git_repo_path)
            self.git_available = True
        except InvalidGitRepositoryError:
            # Not a git repo, that's OK
            pass

        # Calculate relative path for git operations
        try:
            self.relative_yaml_path = self.benchmark_path.relative_to(
                self.git_repo_path
            )
        except ValueError:
            # benchmark_path is not under git_repo_path
            self.relative_yaml_path = self.benchmark_path

    def reconstruct_version_history(self) -> List[str]:
        """
        Reconstruct version history from git commits.

        Returns:
            List of version strings found in git history, in chronological order
        """
        if not self.git_available or not self.repo:
            return []

        versions = []
        seen_versions = set()

        try:
            # Get all commits that touched the benchmark file
            commits = list(self.repo.iter_commits(paths=str(self.relative_yaml_path)))

            # Process commits from oldest to newest
            for commit in reversed(commits):
                try:
                    # Get the file content at this commit
                    yaml_content = self._get_file_content_at_commit(commit.hexsha)
                    if yaml_content:
                        version = self._extract_version_from_yaml(yaml_content)
                        if version and version not in seen_versions:
                            versions.append(version)
                            seen_versions.add(version)
                except Exception:
                    # Skip commits where we can't extract version
                    continue

        except Exception:
            # If we can't get git history, return empty list
            pass

        return versions

    def _get_file_content_at_commit(self, commit_hash: str) -> Optional[str]:
        """
        Get the content of the benchmark file at a specific commit.

        Args:
            commit_hash: Git commit hash

        Returns:
            File content as string or None if not found
        """
        if self.repo is None:
            return None
        try:
            commit = self.repo.commit(commit_hash)
            # Get the file content from the commit tree
            blob = commit.tree / str(self.relative_yaml_path)
            return blob.data_stream.read().decode("utf-8")
        except Exception:
            return None

    def _extract_version_from_yaml(self, yaml_content: str) -> Optional[str]:
        """
        Extract version from YAML content using the model.

        Args:
            yaml_content: YAML file content as string

        Returns:
            Version string or None if not found/invalid
        """
        try:
            # Import here to avoid circular dependency
            from omnibenchmark.model import Benchmark

            # Load from string (not file)
            import yaml

            data = yaml.safe_load(yaml_content)

            # Create Benchmark instance from data
            benchmark = Benchmark(**data)
            return benchmark.version
        except Exception:
            # If parsing fails, return None
            return None

    def initialize_from_git_history(self) -> None:
        """
        Initialize the known versions from git history.

        This method reconstructs the version history and sets it as known versions.
        """
        versions = self.reconstruct_version_history()
        self.set_known_versions(versions)

    def get_current_git_info(self) -> Dict[str, Any]:
        """
        Get current git repository information.

        Returns:
            Dictionary with git info (commit, branch, author, etc.)
        """
        info = {}

        if not self.git_available or not self.repo:
            return info

        try:
            # Get current commit
            head_commit = self.repo.head.commit
            info["commit"] = head_commit.hexsha[:8]
            info["commit_full"] = head_commit.hexsha
            info["author"] = f"{head_commit.author.name} <{head_commit.author.email}>"
            info["timestamp"] = head_commit.committed_datetime.isoformat()

            # Get current branch
            if not self.repo.head.is_detached:
                info["branch"] = self.repo.active_branch.name

            # Check if working directory is clean
            info["clean"] = not self.repo.is_dirty()

            # Get remote URL if any
            if self.repo.remotes:
                info["remote"] = self.repo.remotes.origin.url

        except Exception:
            # Return whatever info we could gather
            pass

        return info

    def create_version_with_git_tracking(
        self,
        version: Optional[str] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> str:
        """
        Create a new version with git commit tracking.

        Args:
            version: Version to create. If None, auto-increment.
            metadata: Optional metadata to store with version

        Returns:
            The created version string
        """
        # Add git info to metadata
        if metadata is None:
            metadata = {}

        # Get and store git information
        git_info = self.get_current_git_info()
        if git_info:
            metadata["git"] = git_info

        # Create version with base class, passing metadata through hooks
        def pre_create_hook(version_str: str):
            # Could validate git state here if needed
            pass

        def post_create_hook(version_str: str):
            # Could create git tag here if desired
            # But based on requirements, versions come from YAML not tags
            pass

        created_version = self.create_version(
            version=version,
            pre_create_hook=pre_create_hook,
            post_create_hook=post_create_hook,
        )

        return created_version

    def update_benchmark_version_and_commit(
        self,
        benchmark_execution,
        new_version: str,
        commit_message: Optional[str] = None,
    ) -> str:
        """
        Update the benchmark YAML file with new version and commit to git.

        Args:
            benchmark_execution: BenchmarkExecution object with model and context
            new_version: Version string to set in the benchmark
            commit_message: Optional custom commit message

        Returns:
            The commit hash of the new commit

        Raises:
            VersioningError: If file update or git commit fails
        """
        from omnibenchmark.versioning.exceptions import VersioningError

        if not self.git_available or not self.repo:
            raise VersioningError("Git repository not available for committing")

        try:
            # Update the benchmark model's version
            benchmark_execution.model.version = new_version

            # Serialize the updated model to YAML
            yaml_content = benchmark_execution.model.to_yaml()

            # Write the updated YAML to the file
            with open(self.benchmark_path, "w") as f:
                f.write(yaml_content)

            # Stage the file for commit
            self.repo.index.add([str(self.relative_yaml_path)])

            # Create commit message
            if commit_message is None:
                commit_message = f"Update benchmark version to {new_version}"

            # Commit the change
            commit = self.repo.index.commit(commit_message)

            return commit.hexsha

        except Exception as e:
            raise VersioningError(f"Failed to update benchmark and commit: {e}")

    def create_version_with_persistence(
        self,
        benchmark_execution,
        version: Optional[str] = None,
        commit_message: Optional[str] = None,
    ) -> str:
        """
        Create a new version, update the benchmark file, and commit to git.

        This method combines version creation with file persistence:
        1. Creates a new version using the base class create_version method
        2. Updates the benchmark YAML file with the new version
        3. Commits the change to git

        Args:
            benchmark_execution: BenchmarkExecution object with model and context
            version: Version to create. If None, auto-increment.
            commit_message: Optional custom commit message

        Returns:
            The created version string

        Raises:
            VersionAlreadyExistsError: If version already exists
            VersionFormatError: If version format is invalid
            VersioningError: If file update or git commit fails
        """
        from omnibenchmark.versioning.exceptions import VersioningError

        if not self.git_available or not self.repo:
            raise VersioningError("Git repository not available for committing")

        def pre_create_hook(version_str: str):
            """Hook called before version is created - could validate git state here."""
            pass

        def post_create_hook(version_str: str):
            """Hook called after version is created - updates file and commits."""
            try:
                # Update the benchmark model's version
                benchmark_execution.model.version = version_str

                # Serialize the updated model to YAML
                yaml_content = benchmark_execution.model.to_yaml()

                # Write the updated YAML to the file
                with open(self.benchmark_path, "w") as f:
                    f.write(yaml_content)

                if self.repo is not None:
                    # Stage the file for commit
                    self.repo.index.add([str(self.relative_yaml_path)])

                # Create commit message
                final_commit_message = (
                    commit_message or f"Update benchmark version to {version_str}"
                )

                # Commit the change
                if self.repo is not None:
                    _commit = self.repo.index.commit(final_commit_message)

            except Exception as e:
                raise VersioningError(
                    f"Failed to update benchmark file and commit: {e}"
                )

        # Create version with hooks for persistence
        created_version = self.create_version(
            version=version,
            pre_create_hook=pre_create_hook,
            post_create_hook=post_create_hook,
        )

        return created_version

    def get_version_at_commit(self, commit_hash: str) -> Optional[str]:
        """
        Get the version from the benchmark file at a specific commit.

        Args:
            commit_hash: Git commit hash

        Returns:
            Version string or None if not found
        """
        if not self.git_available or not self.repo:
            return None

        yaml_content = self._get_file_content_at_commit(commit_hash)
        if yaml_content:
            return self._extract_version_from_yaml(yaml_content)
        return None

    def get_commits_for_version(self, version: str) -> List[str]:
        """
        Get all commits where the benchmark had a specific version.

        Args:
            version: Version to search for

        Returns:
            List of commit hashes (short form)
        """
        if not self.git_available or not self.repo:
            return []

        commits_with_version = []

        try:
            # Get all commits that touched the benchmark file
            commits = list(self.repo.iter_commits(paths=str(self.relative_yaml_path)))

            for commit in commits:
                try:
                    yaml_content = self._get_file_content_at_commit(commit.hexsha)
                    if yaml_content:
                        commit_version = self._extract_version_from_yaml(yaml_content)
                        if commit_version == version:
                            commits_with_version.append(commit.hexsha[:8])
                except Exception:
                    continue

        except Exception:
            pass

        return commits_with_version
