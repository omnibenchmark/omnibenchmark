"""
Git module caching with Go build cache style structure.

This module implements the NEW two-tier git caching system that replaces the legacy
per-commit cloning approach. It provides efficient repository management for
OmniBenchmark workflow execution.

Architecture Overview:
======================

Two-Tier Caching Strategy:
    1. Cache Layer (Global)
       - Location: ~/.cache/omnibenchmark/git/{host}/{org}/{repo}/
       - One full clone per repository (not per commit)
       - Persistent across benchmark runs
       - Fetched/updated as needed

    2. Work Layer (Per-Execution)
       - Location: .snakemake/repos/{commit_hash}/ or similar work directories
       - Lightweight checkout of specific commit from cache
       - Created on-demand before module execution
       - Can be cleaned up after execution

Key Benefits:
    - Disk efficiency: One clone per repo instead of per commit
    - Network efficiency: Fetch updates instead of re-cloning
    - XDG compliance: Uses standard cache directory (~/.cache)
    - Go-style paths: Human-readable structure (github.com/user/repo)
    - Concurrent-safe: Git handles its own locking via .git/ files

Cache Structure Example:
    ~/.cache/omnibenchmark/git/
    ├── github.com/
    │   ├── omnibenchmark/
    │   │   ├── clustering_example/     # Full clone
    │   │   │   ├── .git/
    │   │   │   └── ...
    │   │   └── method_template/        # Full clone
    │   │       ├── .git/
    │   │       └── ...
    │   └── user/
    │       └── custom_module/          # Full clone
    └── gitlab.com/
        └── org/
            └── project/                # Full clone

Comparison with Legacy System:
    OLD (clone_legacy.py):
        - .omnibenchmark/git/{hash}/     (hash = repo + commit)
        - One clone per commit
        - GitPython (full git implementation)
        - Application-level file locking

    NEW (this module):
        - ~/.cache/omnibenchmark/git/{host}/{org}/{repo}/
        - One clone per repository
        - dulwich (lightweight porcelain)
        - Git's native locking

Usage:
    # Populate cache (typically during dry-run)
    cache_path = get_or_update_cached_repo("https://github.com/user/repo")

    # Checkout to work directory before execution
    work_path, commit_hash = checkout_to_work_dir(
        "https://github.com/user/repo",
        "main",  # or specific commit
        Path(".snakemake/repos/work")
    )

    # Or use the unified interface
    work_path, commit_hash = clone_module_v2(
        "https://github.com/user/repo",
        "v1.0.0"
    )

Configuration:
    Cache directory is configurable via omnibenchmark config:
        config["dirs"]["git_cache"] = "~/.cache/omnibenchmark/git"

    Or via XDG_CACHE_HOME environment variable:
        XDG_CACHE_HOME=/custom/path

Thread Safety:
    This module is safe for concurrent use. Git/dulwich handles locking internally
    via .git/index.lock and similar lock files. No application-level locking needed.

Migration:
    See docs/GIT_CACHE_MIGRATION.md for migration plan from legacy system.
"""

import logging
import re
import shutil
import threading
from pathlib import Path
from typing import Optional, Tuple
from urllib.parse import urlparse
from dulwich import porcelain
from dulwich.client import get_transport_and_path
from dulwich.errors import NotGitRepository

from omnibenchmark.config import get_git_cache_dir

logger = logging.getLogger(__name__)


def is_local_path(url: str) -> bool:
    """
    Check if a repository URL is a local filesystem path rather than a remote URL.

    Local paths include ".", "..", absolute paths, and file:// URLs.
    These should bypass the git cache and be copied directly.

    Args:
        url: Repository URL or path

    Returns:
        True if this is a local filesystem path
    """
    url = url.strip()
    if url in (".", "..") or url.startswith("./") or url.startswith("../"):
        return True
    if url.startswith("/"):
        return True
    if url.startswith("file://"):
        return True
    return False


def resolve_local_path(url: str, benchmark_dir: Optional[Path] = None) -> Path:
    """
    Resolve a local repository URL to an absolute filesystem path.

    Args:
        url: Local path (e.g., ".", "/home/user/repo", "file:///path/to/repo")
        benchmark_dir: Directory containing the benchmark YAML (for relative paths).
                      Defaults to cwd.

    Returns:
        Resolved absolute Path
    """
    if url.startswith("file://"):
        url = url[len("file://") :]
    p = Path(url)
    if not p.is_absolute():
        base = benchmark_dir or Path.cwd()
        p = base / p
    return p.resolve()


def copy_local_to_work_dir(
    local_path: Path,
    work_dir: Path,
) -> Tuple[Path, str]:
    """
    Copy a local repository to a work directory, resolving the current HEAD commit.

    Only git-tracked files are copied. This avoids copying heavy untracked
    directories (like .pixi/, out/, node_modules/) and prevents recursive
    copies when the work_dir is inside the source tree.

    Args:
        local_path: Absolute path to the local repository
        work_dir: Target work directory

    Returns:
        Tuple of (work_dir, resolved_commit_hash)

    Raises:
        RuntimeError: If the path is not a git repository or HEAD can't be resolved
    """
    if not local_path.is_dir():
        raise RuntimeError(f"Local module path does not exist: {local_path}")

    # Resolve HEAD commit from the local repo
    try:
        repo = porcelain.open_repo(str(local_path))
        head_bytes = repo.head()
        if len(head_bytes) == 20:
            commit_hash = head_bytes.hex()
        else:
            commit_hash = head_bytes.decode("ascii")
    except (NotGitRepository, Exception) as e:
        raise RuntimeError(
            f"Local path '{local_path}' is not a valid git repository: {e}"
        )

    # Get list of git-tracked files via the index (includes staged files).
    # We copy their working-tree versions so that uncommitted edits are
    # picked up in --dirty mode, while untracked dirs (out/, .pixi/, …)
    # are excluded — which also prevents recursive copies when the
    # work_dir lives inside the source tree.
    try:
        tracked_files = [path.decode("utf-8") for path in repo.open_index()]
    except Exception as e:
        raise RuntimeError(f"Failed to list tracked files in '{local_path}': {e}")

    # Copy only tracked files to work directory
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    logger.info(
        f"Copying local module {local_path} ({len(tracked_files)} tracked files) to {work_dir}"
    )
    for rel_path in tracked_files:
        src = local_path / rel_path
        dst = work_dir / rel_path
        if src.exists():
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)

    return work_dir, commit_hash


# Global lock manager for repository checkouts
# Each repo URL gets its own lock to prevent concurrent checkout conflicts
_repo_locks = {}
_repo_locks_lock = threading.Lock()


def _get_repo_lock(repo_url: str) -> threading.Lock:
    """
    Get or create a lock for a specific repository URL.

    This ensures that only one thread at a time can checkout from the same
    cached repository, preventing git index.lock conflicts.

    Args:
        repo_url: Git repository URL

    Returns:
        Thread lock for this repository
    """
    with _repo_locks_lock:
        if repo_url not in _repo_locks:
            _repo_locks[repo_url] = threading.Lock()
        return _repo_locks[repo_url]


def parse_repo_url(url: str) -> str:
    """
    Parse a git repository URL into a Go-style cache path.

    Examples:
        https://github.com/user/repo.git -> github.com/user/repo
        git@github.com:user/repo.git -> github.com/user/repo
        https://gitlab.com/group/subgroup/project -> gitlab.com/group/subgroup/project

    Args:
        url: Git repository URL

    Returns:
        Path-like string (e.g., "github.com/user/repo")
    """
    # Remove .git suffix if present
    url = url.rstrip("/")
    if url.endswith(".git"):
        url = url[:-4]

    # Handle SSH URLs (git@host:path)
    ssh_match = re.match(r"^git@([^:]+):(.+)$", url)
    if ssh_match:
        host, path = ssh_match.groups()
        return f"{host}/{path}"

    # Handle HTTPS URLs
    parsed = urlparse(url)
    if parsed.netloc and parsed.path:
        # Remove leading slash from path
        path = parsed.path.lstrip("/")
        return f"{parsed.netloc}/{path}"

    # Fallback: treat as is
    return url.replace(":", "/")


def get_or_update_cached_repo(repo_url: str, cache_dir: Optional[Path] = None) -> Path:
    """
    Get or clone a full repository to the cache.

    Cache structure:
        ~/.cache/omnibenchmark/git/github.com/user/repo/

    This maintains ONE full clone per repository (not per commit).
    Updates are fetched if the repo already exists.

    Args:
        repo_url: Git repository URL
        cache_dir: Base cache directory (defaults to ~/.cache/omnibenchmark/git)

    Returns:
        Path to the cached repository
    """
    if cache_dir is None:
        cache_dir = get_git_cache_dir()

    # Parse repo URL to path
    repo_path = parse_repo_url(repo_url)

    # Full path: cache_dir/github.com/user/repo
    repo_cache_dir = cache_dir / repo_path
    repo_cache_dir.parent.mkdir(parents=True, exist_ok=True)

    # Note: Git/dulwich handles its own locking internally via .git/ lock files
    # No need for application-level file locking

    if repo_cache_dir.exists():
        # Update existing repository
        try:
            repo = porcelain.open_repo(str(repo_cache_dir))
            logging.info(f"Updating cached repository at {repo_cache_dir}")
            # Fetch all updates from remote
            # Note: porcelain.fetch by default fetches all refs
            porcelain.fetch(repo, repo_url)
            return repo_cache_dir
        except (NotGitRepository, Exception) as e:
            logging.warning(f"Cached repo appears corrupt: {e}. Re-cloning.")
            shutil.rmtree(repo_cache_dir, ignore_errors=True)

    # Clone fresh
    logging.info(f"Cloning {repo_url} to cache at {repo_cache_dir}")

    try:
        # Clone as a full repository (not bare, not shallow)
        # Using checkout=True to ensure all refs and objects are fetched
        porcelain.clone(
            source=repo_url,
            target=str(repo_cache_dir),
            checkout=True,  # Checkout default branch to fetch all objects
            bare=False,
        )
        return repo_cache_dir
    except Exception as e:
        logging.error(f"Failed to clone repository: {e}")
        if repo_cache_dir.exists():
            shutil.rmtree(repo_cache_dir, ignore_errors=True)
        raise RuntimeError(f"Failed to clone {repo_url}: {e}")


def checkout_to_work_dir(
    repo_url: str, ref: str, work_dir: Path, cache_dir: Optional[Path] = None
) -> Tuple[Path, str]:
    """
    Checkout a specific commit/branch/tag to a work directory.

    This does a lightweight clone from the cache to the work directory.

    Args:
        repo_url: Git repository URL
        ref: Git reference (branch, tag, or commit hash)
        work_dir: Work directory path (e.g., .omnibenchmark/modules/github.com/user/repo/abc1234)
        cache_dir: Base cache directory (defaults to ~/.cache/omnibenchmark/git)

    Returns:
        Tuple of (work_dir_path, resolved_commit_hash)
    """
    # Ensure the repo is cached
    repo_cache_dir = get_or_update_cached_repo(repo_url, cache_dir)

    # Get a lock for this repository to prevent concurrent checkout conflicts
    # When multiple threads try to checkout from the same cached repo,
    # git operations like reset/checkout will conflict on .git/index.lock
    repo_lock = _get_repo_lock(repo_url)

    with repo_lock:
        # Open cached repo
        _cached_repo = porcelain.open_repo(str(repo_cache_dir))

        # Create work directory
        work_dir.parent.mkdir(parents=True, exist_ok=True)

        # IMPORTANT: We use the cache repo directly instead of cloning to work_dir
        # This avoids issues with dulwich's clone not copying all pack objects
        # We'll checkout the specific ref directly in the cache, then copy files

        # Open the cache repo
        work_repo = porcelain.open_repo(str(repo_cache_dir))

        # Resolve and checkout the ref
        from dulwich.objectspec import parse_commit

        # Try to resolve the ref - handle commit hashes, local branches, tags,
        # and remote tracking branches (refs/remotes/origin/{ref})
        try:
            commit_obj = parse_commit(work_repo, ref.encode("ascii"))
            commit_hash_hex = commit_obj.id.hex()
            logging.debug(f"Resolved {ref} to {commit_hash_hex}")
        except KeyError:
            commit_obj = None

            # Try remote tracking branch: refs/remotes/origin/{ref}
            remote_ref = f"refs/remotes/origin/{ref}".encode("ascii")
            if remote_ref in work_repo.refs:
                sha = work_repo.refs[remote_ref]
                commit_obj = parse_commit(work_repo, sha)
                commit_hash_hex = commit_obj.id.hex()
                logging.debug(
                    f"Resolved {ref} via remote tracking branch to {commit_hash_hex}"
                )

            # If that fails and ref looks like a short hash, try to expand it
            if commit_obj is None and re.match(r"^[0-9a-f]{4,39}$", ref, re.IGNORECASE):
                ref_prefix = ref.lower()
                matching_commits = []

                for obj_id in work_repo.object_store:
                    obj_hex = obj_id.hex()
                    if obj_hex.startswith(ref_prefix):
                        matching_commits.append(obj_hex)

                if len(matching_commits) == 1:
                    full_hash = matching_commits[0]
                    logging.debug(f"Expanded short hash {ref} to {full_hash}")
                    commit_obj = parse_commit(work_repo, full_hash.encode("ascii"))
                    commit_hash_hex = commit_obj.id.hex()
                elif len(matching_commits) > 1:
                    raise RuntimeError(
                        f"Reference '{ref}' is ambiguous (matches {len(matching_commits)} commits) in repository {repo_url}"
                    )

            if commit_obj is None:
                raise RuntimeError(
                    f"Reference '{ref}' not found in repository {repo_url}"
                )

        # Checkout the commit in the cache repo
        logging.info(f"Checking out {repo_url}@{ref} in cache")
        porcelain.reset(work_repo, "hard", commit_obj)

        # Get the actual commit hash
        actual_commit_bytes = work_repo.head()
        if len(actual_commit_bytes) == 20:
            actual_commit = actual_commit_bytes.hex()
        elif len(actual_commit_bytes) == 40:
            actual_commit = actual_commit_bytes.decode("ascii")
        else:
            actual_commit = actual_commit_bytes.decode("ascii")

        # Copy the checked out files to work directory (excluding .git)
        # This gives us a clean working copy without git metadata
        if work_dir.exists():
            shutil.rmtree(work_dir)

        logging.info(f"Copying files to {work_dir}")
        shutil.copytree(
            repo_cache_dir,
            work_dir,
            ignore=shutil.ignore_patterns(".git"),
            symlinks=True,
        )

        logging.info(f"Checked out {repo_url}@{actual_commit[:7]} to {work_dir}")

        return work_dir, actual_commit


def clone_module_v2(
    repository_url: str,
    ref: str,
    work_dir: Optional[Path] = None,
    cache_dir: Optional[Path] = None,
) -> Tuple[Path, str]:
    """
    Clone a git module using the new two-tier caching system.

    This is the main entry point for the new cloning system. It:
    1. Ensures the repo exists in cache (full clone)
    2. Checks out the specific ref to a work directory
    3. Returns both the path and the resolved commit

    Args:
        repository_url: Git repository URL
        ref: Git reference (branch, tag, or commit hash)
        work_dir: Work directory (defaults to temp location based on commit)
        cache_dir: Base cache directory (defaults to ~/.cache/omnibenchmark/git)

    Returns:
        Tuple of (work_dir_path, commit_hash)
    """
    # If no work_dir specified, create a temp one based on repo and ref
    if work_dir is None:
        from tempfile import gettempdir

        repo_path = parse_repo_url(repository_url)
        ref_short = ref[:7] if len(ref) >= 7 else ref
        work_dir = Path(gettempdir()) / "omnibenchmark" / repo_path / ref_short

    return checkout_to_work_dir(repository_url, ref, work_dir, cache_dir)


def describe_cache(
    cache_dir: Optional[Path] = None,
    fetch: bool = False,
    repos_filter: Optional[set] = None,
) -> list:
    """
    Describe the status of cached git repositories.

    Args:
        cache_dir: Base cache directory (defaults to ~/.cache/omnibenchmark/git)
        fetch: If True, fetch from remote to check for updates
        repos_filter: Optional set of repository URLs to filter by

    Returns:
        List of dictionaries with repo information:
        - repo_path: Relative path in cache (e.g., "github.com/user/repo")
        - url: Git repository URL
        - head: Current HEAD commit hash
        - branch: Current branch name (or "detached")
        - remote_head: Remote HEAD commit (if fetch=True)
        - refs_using_branches: List of branch names used as refs (not commit hashes)
    """
    if cache_dir is None:
        cache_dir = get_git_cache_dir()

    if not cache_dir.exists():
        return []

    results = []

    # Walk through cache directory to find git repos
    for repo_dir in cache_dir.rglob(".git"):
        if not repo_dir.is_dir():
            continue

        repo_path = repo_dir.parent
        repo_rel_path = repo_path.relative_to(cache_dir)

        try:
            repo = porcelain.open_repo(str(repo_path))

            # Get remote URL
            config = repo.get_config()
            remote_url = config.get((b"remote", b"origin"), b"url")
            if remote_url:
                remote_url = remote_url.decode("utf-8")
            else:
                remote_url = "unknown"

            # Filter by repos if specified
            if repos_filter and remote_url not in repos_filter:
                continue

            # Get current HEAD
            try:
                head_bytes = repo.head()
                if len(head_bytes) == 20:
                    head_commit = head_bytes.hex()
                elif len(head_bytes) == 40:
                    head_commit = head_bytes.decode("ascii")
                else:
                    head_commit = head_bytes.decode("ascii")
            except Exception:
                head_commit = "unknown"

            # Get current branch
            try:
                refs = repo.get_refs()
                head_ref = refs.get(b"HEAD")

                # Check if HEAD points to a branch
                branch_name = "detached"
                if head_ref:
                    for ref_name, ref_sha in refs.items():
                        if ref_name.startswith(b"refs/heads/") and ref_sha == head_ref:
                            branch_name = ref_name.decode("ascii").replace(
                                "refs/heads/", ""
                            )
                            break
            except Exception:
                branch_name = "unknown"

            repo_info = {
                "repo_path": str(repo_rel_path),
                "url": remote_url,
                "head": head_commit,
                "branch": branch_name,
            }

            # If fetch is requested, get remote HEAD
            if fetch and remote_url != "unknown":
                try:
                    logging.info(f"Fetching updates for {repo_rel_path}...")
                    porcelain.fetch(repo, remote_url)

                    # Get remote HEAD
                    client, path = get_transport_and_path(remote_url)
                    remote_refs = client.get_refs(path)

                    # Get remote HEAD for the current branch
                    if branch_name != "detached" and branch_name != "unknown":
                        remote_ref_name = f"refs/heads/{branch_name}".encode()
                        if remote_ref_name in remote_refs:
                            remote_head = remote_refs[remote_ref_name].hex()
                            repo_info["remote_head"] = remote_head
                except Exception as e:
                    logging.debug(
                        f"Failed to fetch remote info for {repo_rel_path}: {e}"
                    )

            results.append(repo_info)

        except Exception as e:
            logging.debug(f"Failed to read repo at {repo_path}: {e}")
            continue

    return results
