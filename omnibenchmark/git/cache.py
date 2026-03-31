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

import io
import logging
import re
import shutil
import tempfile
from pathlib import Path
from typing import Optional, Tuple, cast
from urllib.parse import urlparse
from dulwich import porcelain
from dulwich.client import get_transport_and_path
from dulwich.errors import NotGitRepository
from dulwich.refs import Ref
from dulwich.repo import Repo
from filelock import FileLock

from omnibenchmark.config import get_git_cache_dir

logger = logging.getLogger(__name__)


# Dulwich writes clone/fetch progress directly to streams, bypassing Python
# logging.  Redirect both stdout and stderr to /dev/null so progress lines
# ("Enumerating objects", "copied N pack entries", etc.) don't leak to the
# user's terminal.
# errstream expects BinaryIO; outstream expects TextIO.
# Use in-memory null streams to avoid leaking open file handles.
class _NullBinaryStream(io.BytesIO):
    def write(self, b) -> int:
        return len(b) if b else 0


class _NullTextStream(io.StringIO):
    def write(self, s: str) -> int:
        return len(s)


_DEVNULL = _NullBinaryStream()
_DEVNULL_TEXT = _NullTextStream()


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


def _build_gitignore_spec(root: Path):
    """
    Build a pathspec matcher from all .gitignore files in the tree.

    Walks the directory tree pruning already-ignored directories so we never
    descend into large dirs like out/ or .pixi/.  Always ignores .git/.
    """
    import pathspec

    patterns = [".git/", ".git"]

    def _read_gitignore(path: Path, rel_dir: Path):
        try:
            prefix = "" if str(rel_dir) == "." else f"{rel_dir}/"
            for line in path.read_text(errors="replace").splitlines():
                line = line.strip()
                if line and not line.startswith("#"):
                    patterns.append(f"{prefix}{line}")
        except Exception:
            pass

    # Seed from root .gitignore first so we can prune immediately
    root_gi = root / ".gitignore"
    if root_gi.is_file():
        _read_gitignore(root_gi, Path("."))

    # Walk manually, pruning dirs that already match
    dirs_to_visit = [root]
    while dirs_to_visit:
        current = dirs_to_visit.pop()
        spec_so_far = pathspec.PathSpec.from_lines("gitignore", patterns)
        try:
            entries = list(current.iterdir())
        except PermissionError:
            continue
        for entry in entries:
            if not entry.is_dir():
                continue
            if entry.name == ".git":
                continue
            try:
                rel = str(entry.relative_to(root)) + "/"
            except ValueError:
                continue
            if spec_so_far.match_file(rel):
                continue  # pruned — don't descend
            gi = entry / ".gitignore"
            if gi.is_file():
                _read_gitignore(gi, entry.relative_to(root))
            dirs_to_visit.append(entry)

    return pathspec.PathSpec.from_lines("gitignore", patterns)


def copy_local_to_work_dir(
    local_path: Path,
    work_dir: Path,
) -> Tuple[Path, str]:
    """
    Copy a local repository to a work directory, resolving the current HEAD commit.

    Copies all files that are NOT excluded by .gitignore (at any level of the
    tree), including untracked files.  This means --dirty mode picks up new
    scripts and edits without requiring a commit, while heavy directories like
    out/, .pixi/, and node_modules/ (which are typically gitignored) are
    skipped automatically.

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
        repo = cast(Repo, porcelain.open_repo(str(local_path)))
        head_bytes = repo.head()
        if len(head_bytes) == 20:
            commit_hash = head_bytes.hex()
        else:
            commit_hash = head_bytes.decode("ascii")
    except (NotGitRepository, Exception) as e:
        raise RuntimeError(
            f"Local path '{local_path}' is not a valid git repository: {e}"
        )

    spec = _build_gitignore_spec(local_path)

    def _ignore(src_dir: str, names: list) -> set:
        """shutil.copytree ignore callback: return names to skip."""
        src_dir_path = Path(src_dir)
        ignored = set()
        for name in names:
            abs_path = src_dir_path / name
            try:
                rel = str(abs_path.relative_to(local_path))
            except ValueError:
                rel = name
            # Append "/" for directories so directory patterns match
            if abs_path.is_dir():
                rel = rel + "/"
            if spec.match_file(rel):
                ignored.add(name)
        return ignored

    if work_dir.exists():
        shutil.rmtree(work_dir)

    shutil.copytree(str(local_path), str(work_dir), ignore=_ignore, symlinks=True)

    logger.info(f"Copied local module {local_path} to {work_dir} (gitignore-filtered)")

    return work_dir, commit_hash


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

    # Use a cross-process file lock to prevent concurrent clones of the same repo.
    # Two pytest-xdist workers (separate processes) would otherwise both see
    # repo_cache_dir missing and race to create it, causing [Errno 17] File exists.
    lock_path = repo_cache_dir.parent / (repo_cache_dir.name + ".lock")
    with FileLock(str(lock_path)):
        if repo_cache_dir.exists():
            # Update existing repository
            try:
                logging.info(f"Updating cached repository at {repo_cache_dir}")
                # Use the repo as a context manager so it is closed (and any
                # in-memory ref caches are released) before callers open their
                # own handle.  Fetch by remote name ("origin") so dulwich
                # applies the configured refspec and advances
                # refs/remotes/origin/*; fetching by URL alone retrieves
                # objects but does not update tracking refs.
                with cast(Repo, porcelain.open_repo(str(repo_cache_dir))) as repo:
                    porcelain.fetch(
                        repo, "origin", errstream=_DEVNULL, outstream=_DEVNULL_TEXT
                    )

                    # Update the working tree to match the remote tracking branch
                    # so the cache directory reflects the latest fetched state
                    try:
                        from dulwich.index import build_index_from_tree
                        from dulwich.objectspec import parse_commit

                        # Get the current branch from HEAD
                        head_ref = repo.refs.read_ref(Ref(b"HEAD"))
                        if head_ref and head_ref.startswith(b"ref: refs/heads/"):
                            local_branch = head_ref[len(b"ref: refs/heads/") :]
                            remote_ref = Ref(b"refs/remotes/origin/" + local_branch)

                            if remote_ref in repo.refs:
                                # Update local branch to match remote tracking branch
                                remote_sha = repo.refs[remote_ref]
                                repo.refs[Ref(b"refs/heads/" + local_branch)] = (
                                    remote_sha
                                )

                                # Update working tree to the new HEAD
                                commit_obj = parse_commit(repo, remote_sha)
                                # Use temp index file since we don't want to keep it
                                with tempfile.TemporaryDirectory() as _tmp:
                                    build_index_from_tree(
                                        str(repo_cache_dir),
                                        str(Path(_tmp) / "index"),
                                        repo.object_store,
                                        commit_obj.tree,
                                    )
                                logging.debug(
                                    f"Updated cache working tree to {remote_sha.hex()[:7]}"
                                )
                    except Exception as e:
                        logging.debug(f"Could not update working tree after fetch: {e}")
                        # Non-fatal: remote refs are updated, which is sufficient

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
                errstream=_DEVNULL,
                outstream=_DEVNULL,
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

    Writes the git tree directly from the cache object store to work_dir
    without mutating the cache repository's working tree.  This makes
    concurrent checkouts of different refs from the same cached repo safe
    across both threads and processes.

    Args:
        repo_url: Git repository URL
        ref: Git reference (branch, tag, or commit hash)
        work_dir: Work directory path
        cache_dir: Base cache directory (defaults to ~/.cache/omnibenchmark/git)

    Returns:
        Tuple of (work_dir_path, resolved_commit_hash)
    """
    from dulwich.index import build_index_from_tree
    from dulwich.objectspec import parse_commit

    # Ensure the repo is cached (fetch updates or fresh clone)
    repo_cache_dir = get_or_update_cached_repo(repo_url, cache_dir)

    # Open cached repo read-only — we never mutate its working tree
    cached_repo = cast(Repo, porcelain.open_repo(str(repo_cache_dir)))

    # Resolve the ref to a commit object.
    # Order: direct parse (full hash / tag) → remote tracking branch → short hash.
    def _id_to_str(sha: bytes) -> str:
        """dulwich ShaFile.id is 40-byte ASCII hex; repo.head() may be 20-byte binary."""
        return sha.hex() if len(sha) == 20 else sha.decode("ascii")

    commit_obj = None

    # Check remote tracking branch first.
    # dulwich's parse_commit(repo, b"main") resolves via refs/heads/main (the
    # local branch), which is NOT updated by fetch — only refs/remotes/origin/*
    # is advanced.  Checking the remote tracking ref first ensures we always
    # return the commit that the remote HEAD points to after a fresh fetch.
    remote_ref = Ref(f"refs/remotes/origin/{ref}".encode("ascii"))
    if remote_ref in cached_repo.refs:
        sha = cached_repo.refs[remote_ref]
        commit_obj = parse_commit(cached_repo, sha)
        logging.debug(
            f"Resolved {ref} via refs/remotes/origin/{ref} to {_id_to_str(commit_obj.id)}"
        )

    # Direct parse: full commit hash, tag, or any other ref dulwich understands.
    if commit_obj is None:
        try:
            commit_obj = parse_commit(cached_repo, ref.encode("ascii"))
            logging.debug(f"Resolved {ref} to {_id_to_str(commit_obj.id)}")
        except KeyError:
            pass

    # Short hash expansion — object_store yields 20-byte binary SHAs
    if commit_obj is None and re.match(r"^[0-9a-f]{4,39}$", ref, re.IGNORECASE):
        ref_prefix = ref.lower()
        matching = [
            obj_id.hex()
            for obj_id in cached_repo.object_store
            if obj_id.hex().startswith(ref_prefix)
        ]
        if len(matching) == 1:
            commit_obj = parse_commit(cached_repo, matching[0].encode("ascii"))
            logging.debug(f"Expanded short hash {ref} to {_id_to_str(commit_obj.id)}")
        elif len(matching) > 1:
            raise RuntimeError(
                f"Reference '{ref}' is ambiguous (matches {len(matching)} commits)"
                f" in repository {repo_url}"
            )

    if commit_obj is None:
        raise RuntimeError(f"Reference '{ref}' not found in repository {repo_url}")

    # Write the commit's tree directly from the object store to work_dir.
    # build_index_from_tree needs an index path — we use a temp file and
    # discard it so no git metadata leaks into the work directory.
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True)

    with tempfile.TemporaryDirectory() as _tmp:
        build_index_from_tree(
            str(work_dir),
            str(Path(_tmp) / "index"),
            cached_repo.object_store,
            commit_obj.tree,
        )

    actual_commit = _id_to_str(commit_obj.id)
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
            repo = cast(Repo, porcelain.open_repo(str(repo_path)))

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
                head_ref = refs.get(Ref(b"HEAD"))

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
                    porcelain.fetch(
                        repo, remote_url, errstream=_DEVNULL, outstream=_DEVNULL_TEXT
                    )

                    # Get remote HEAD
                    client, path = get_transport_and_path(remote_url)
                    path_bytes = path.encode() if isinstance(path, str) else path
                    ls_result = client.get_refs(path_bytes)
                    remote_refs = ls_result.refs

                    # Get remote HEAD for the current branch
                    if branch_name != "detached" and branch_name != "unknown":
                        remote_ref_name = Ref(f"refs/heads/{branch_name}".encode())
                        if remote_ref_name in remote_refs:
                            ref_sha = remote_refs[remote_ref_name]
                            if ref_sha is not None:
                                remote_head = ref_sha.hex()
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
