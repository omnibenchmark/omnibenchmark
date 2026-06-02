"""
Git operations module for omnibenchmark.

Two-tier caching system:
- Cache Layer: Full clones in ~/.cache/omnibenchmark/git/{host}/{org}/{repo}/
- Work Layer: Lightweight checkouts to .modules/{repo_name}/{commit}/ before execution
"""

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
    "clone_module_v2",
    "copy_local_to_work_dir",
    "is_local_path",
    "parse_repo_url",
    "get_or_update_cached_repo",
    "checkout_to_work_dir",
    "describe_cache",
    "resolve_local_path",
]
