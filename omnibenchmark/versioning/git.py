"""
Git-aware benchmark version manager.

This module extends BenchmarkVersionManager with git integration
for reconstructing version history from git commits.

YAML parsing is delegated to the model layer via Benchmark.from_yaml().
"""

import subprocess
from pathlib import Path
from typing import List, Optional, Dict, Any

from dulwich import porcelain
from dulwich.errors import NotGitRepository
from dulwich.repo import Repo

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
            self.repo = Repo(str(self.git_repo_path))
            self.git_available = True
        except NotGitRepository:
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
            # Walk commits that touched the benchmark file (HEAD..root)
            rel_path_bytes = str(self.relative_yaml_path).encode()
            commits = list(self.repo.get_walker(paths=[rel_path_bytes]))

            # Process commits from oldest to newest
            for entry in reversed(commits):
                try:
                    commit_sha = entry.commit.id.decode()
                    yaml_content = self._get_file_content_at_commit(commit_sha)
                    if yaml_content:
                        version = self._extract_version_from_yaml(yaml_content)
                        if version and version not in seen_versions:
                            versions.append(version)
                            seen_versions.add(version)
                except Exception:
                    continue

        except Exception:
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
            commit_id = (
                commit_hash.encode() if isinstance(commit_hash, str) else commit_hash
            )
            commit = self.repo[commit_id]
            # Walk the tree to find the blob
            tree = self.repo[commit.tree]
            parts = str(self.relative_yaml_path).split("/")
            obj = tree
            for part in parts:
                part_bytes = part.encode()
                # tree items: (mode, name, sha)
                found = False
                for item in obj.items():
                    if item.path == part_bytes:
                        obj = self.repo[item.sha]
                        found = True
                        break
                if not found:
                    return None
            return obj.data.decode("utf-8")
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
            from omnibenchmark.model import Benchmark
            import yaml

            data = yaml.safe_load(yaml_content)
            benchmark = Benchmark(**data)
            return benchmark.version
        except Exception:
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
        info: Dict[str, Any] = {}

        if not self.git_available or not self.repo:
            return info

        try:
            head_commit = self.repo[self.repo.head()]
            sha = head_commit.id.decode()
            info["commit"] = sha[:8]
            info["commit_full"] = sha
            info["author"] = head_commit.author.decode()
            info["timestamp"] = head_commit.author_time

            # Current branch (HEAD may be detached)
            try:
                head_ref = self.repo.refs.get_symrefs().get(b"HEAD")
                if head_ref and head_ref.startswith(b"refs/heads/"):
                    info["branch"] = head_ref[len(b"refs/heads/") :].decode()
            except Exception:
                pass

            # Dirty check via subprocess (dulwich has no simple is_dirty())
            try:
                result = subprocess.run(
                    ["git", "-C", str(self.git_repo_path), "status", "--porcelain"],
                    capture_output=True,
                    text=True,
                )
                info["clean"] = result.stdout.strip() == ""
            except Exception:
                pass

            # Remote URL
            try:
                config = self.repo.get_config()
                url = config.get((b"remote", b"origin"), b"url")
                if url:
                    info["remote"] = url.decode()
            except Exception:
                pass

        except Exception:
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
        if metadata is None:
            metadata = {}

        git_info = self.get_current_git_info()
        if git_info:
            metadata["git"] = git_info

        def pre_create_hook(version_str: str):
            pass

        def post_create_hook(version_str: str):
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
            benchmark_execution.model.version = new_version
            yaml_content = benchmark_execution.model.to_yaml()

            with open(self.benchmark_path, "w") as f:
                f.write(yaml_content)

            if commit_message is None:
                commit_message = f"Update benchmark version to {new_version}"

            porcelain.add(self.repo, paths=[str(self.relative_yaml_path)])
            commit_sha = porcelain.commit(self.repo, message=commit_message.encode())
            return commit_sha.decode() if isinstance(commit_sha, bytes) else commit_sha

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
            pass

        def post_create_hook(version_str: str):
            try:
                benchmark_execution.model.version = version_str
                yaml_content = benchmark_execution.model.to_yaml()

                with open(self.benchmark_path, "w") as f:
                    f.write(yaml_content)

                if self.repo is not None:
                    porcelain.add(self.repo, paths=[str(self.relative_yaml_path)])

                final_commit_message = (
                    commit_message or f"Update benchmark version to {version_str}"
                )

                if self.repo is not None:
                    porcelain.commit(self.repo, message=final_commit_message.encode())

            except Exception as e:
                raise VersioningError(
                    f"Failed to update benchmark file and commit: {e}"
                )

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
            rel_path_bytes = str(self.relative_yaml_path).encode()
            commits = list(self.repo.get_walker(paths=[rel_path_bytes]))

            for entry in commits:
                try:
                    commit_sha = entry.commit.id.decode()
                    yaml_content = self._get_file_content_at_commit(commit_sha)
                    if yaml_content:
                        commit_version = self._extract_version_from_yaml(yaml_content)
                        if commit_version == version:
                            commits_with_version.append(commit_sha[:8])
                except Exception:
                    continue

        except Exception:
            pass

        return commits_with_version
