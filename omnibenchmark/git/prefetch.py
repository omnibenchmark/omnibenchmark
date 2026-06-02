"""Pre-fetch git repositories into the two-tier cache before execution."""

import logging
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import TYPE_CHECKING

from omnibenchmark.config import get_git_cache_dir
from omnibenchmark.git.cache import (
    get_or_update_cached_repo,
    is_local_path,
    parse_repo_url,
)

if TYPE_CHECKING:
    from omnibenchmark.benchmark import BenchmarkExecution

logger = logging.getLogger(__name__)


def populate_git_cache(
    benchmark: "BenchmarkExecution", quiet: bool = False, cores: int = 4
) -> None:
    """Populate the git cache with all repositories referenced by a benchmark.

    For each remote repository in the benchmark definition (modules and metric
    collectors), either clones it fresh or fetches updates into the global cache
    at ~/.cache/omnibenchmark/git/.  Already-cached pinned commits are skipped
    to avoid unnecessary network traffic.

    Args:
        benchmark: The benchmark execution object whose modules will be fetched.
        quiet:     When True, show a Rich progress bar instead of log lines.
        cores:     Maximum number of parallel fetch workers.
    """
    from omnibenchmark.progress import ProgressDisplay

    cache_dir = get_git_cache_dir()

    # Collect unique (url, commit) pairs from modules and metric collectors
    repos: dict[str, str | None] = {}
    for stage in benchmark.model.stages:
        for module in stage.modules:
            if hasattr(module, "repository") and module.repository:
                url = module.repository.url
                commit = module.repository.commit
                if url and not is_local_path(url):
                    repos[url] = commit

    if (
        hasattr(benchmark.model, "metric_collectors")
        and benchmark.model.metric_collectors
    ):
        for collector in benchmark.model.metric_collectors:
            if hasattr(collector, "repository") and collector.repository:
                url = collector.repository.url
                commit = collector.repository.commit
                if url and not is_local_path(url):
                    repos[url] = commit

    if not repos:
        if not quiet:
            logger.info("No repositories found in benchmark definition")
        return

    if not quiet:
        logger.info(f"Populating cache with {len(repos)} repositories...")

    progress = ProgressDisplay() if quiet else None
    if progress is not None:
        progress.start_task("Fetching repositories to cache", total=len(repos))

    success_count = 0
    failed: list[tuple[str, str]] = []
    lock = threading.Lock()

    def _fetch_one(repo_url: str, commit: str | None) -> tuple[str, bool, str | None]:
        """Fetch a single repo; skip if it's a pinned commit already in cache."""
        repo_cache_dir = cache_dir / parse_repo_url(repo_url)

        # Skip the network round-trip when the exact pinned commit is already local
        if (
            repo_cache_dir.exists()
            and commit
            and len(commit) == 40
            and all(c in "0123456789abcdef" for c in commit.lower())
        ):
            try:
                from dulwich import porcelain
                from dulwich.repo import Repo as _Repo
                from typing import cast as _cast

                _repo = _cast(_Repo, porcelain.open_repo(str(repo_cache_dir)))
                _repo[commit.encode("ascii")]
                return (repo_url, True, None)
            except Exception:
                pass

        if not quiet:
            logger.info(f"Caching {repo_url}")

        try:
            get_or_update_cached_repo(repo_url, cache_dir)
            return (repo_url, True, None)
        except Exception as e:
            return (repo_url, False, str(e))

    with ThreadPoolExecutor(max_workers=cores) as executor:
        futures = {
            executor.submit(_fetch_one, url, commit): url
            for url, commit in repos.items()
        }

        for future in as_completed(futures):
            repo_url = futures[future]
            try:
                _, success, error = future.result()
                with lock:
                    if success:
                        success_count += 1
                    else:
                        logger.error(f"Failed to cache {repo_url}: {error}")
                        failed.append((repo_url, error or "unknown error"))
            except Exception as e:
                logger.error(f"Exception while caching {repo_url}: {e}")
                with lock:
                    failed.append((repo_url, str(e)))

            if progress is not None:
                progress.update(advance=1)

    if progress is not None:
        progress.finish()
        progress.success(f"Cached {success_count}/{len(repos)} repositories")
    else:
        logger.info(f"Successfully cached {success_count}/{len(repos)} repositories")

    if failed:
        logger.warning(f"Failed to cache {len(failed)} repositories:")
        for url, error in failed:
            logger.warning(f"  {url}: {error}")
