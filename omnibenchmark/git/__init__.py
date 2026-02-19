"""
Git operations module for omnibenchmark.

This module provides functionality for cloning and managing git repositories
used in benchmark workflows, separated from remote storage concerns.

Architecture:
    Two-tier caching system (NEW, recommended):
    - Cache Layer: Full clones in ~/.cache/omnibenchmark/git/{host}/{org}/{repo}/
    - Work Layer: Lightweight checkouts to work directories before execution

    Legacy system (OLD, being phased out):
    - Direct clones to .omnibenchmark/git/{hash}/ (one per commit)

See docs/GIT_CACHE_MIGRATION.md for migration details.
"""

# Legacy imports (GitPython-based, scheduled for removal)
from .clone_legacy import clone_git_repo, clone_module

# New cache-based system (dulwich-based with Go-style cache)
from .cache import (
    clone_module_v2,
    copy_local_to_work_dir,
    is_local_path,
    parse_repo_url,
    get_or_update_cached_repo,
    checkout_to_work_dir,
    describe_cache,
    resolve_local_path,
)

__all__ = [
    # Legacy (GitPython-based) - DEPRECATED
    "clone_git_repo",
    "clone_module",
    # New (dulwich-based with Go-style cache) - RECOMMENDED
    "clone_module_v2",
    "copy_local_to_work_dir",
    "is_local_path",
    "parse_repo_url",
    "get_or_update_cached_repo",
    "checkout_to_work_dir",
    "describe_cache",
    "resolve_local_path",
]
