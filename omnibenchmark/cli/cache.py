"""CLI commands for git cache management"""

import sys
from pathlib import Path

import click

from omnibenchmark.benchmark import BenchmarkExecution
from omnibenchmark.cli.utils.logging import logger
from omnibenchmark.cli.error_formatting import pretty_print_parse_error
from omnibenchmark.config import get_git_cache_dir
from omnibenchmark.git import get_or_update_cached_repo
from omnibenchmark.model.validation import BenchmarkParseError


@click.group(name="cache")
def cache():
    """Manage git repository cache."""
    pass


@cache.command("populate")
@click.argument("benchmark_path", type=click.Path(exists=True))
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    help="Show detailed progress information",
)
def populate(benchmark_path: str, verbose: bool):
    """Populate the cache with repositories from a benchmark definition.

    This command reads the benchmark YAML file, extracts all module repositories,
    and ensures they are cached locally without running any workflows.

    Example:

      ob cache populate benchmark.yaml
    """
    cache_dir = get_git_cache_dir()

    # Parse benchmark
    try:
        b = BenchmarkExecution(Path(benchmark_path))
        if verbose:
            logger.info(f"Loaded benchmark from {benchmark_path}")
    except BenchmarkParseError as e:
        formatted_error = pretty_print_parse_error(e)
        logger.error(f"Failed to load benchmark: {formatted_error}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Failed to load benchmark: {e}")
        sys.exit(1)

    # Extract unique repositories from all stages and modules
    repos = set()
    for stage in b.model.stages:
        for module in stage.modules:
            if hasattr(module, "repository") and module.repository:
                repo_url = module.repository.url
                if repo_url:
                    repos.add(repo_url)

    if not repos:
        logger.info("No repositories found in benchmark definition")
        return

    logger.info(f"Found {len(repos)} unique repositories to cache")
    if verbose:
        logger.info(f"Cache directory: {cache_dir}")

    # Populate cache
    success_count = 0
    failed = []

    for repo_url in sorted(repos):
        try:
            if verbose:
                logger.info(f"Caching {repo_url}...")
            else:
                logger.info(f"Caching {repo_url}")

            cached_path = get_or_update_cached_repo(repo_url, cache_dir)

            if verbose:
                logger.info(f"  Cached at {cached_path}")

            success_count += 1

        except Exception as e:
            logger.error(f"Failed to cache {repo_url}: {e}")
            failed.append((repo_url, str(e)))

    # Summary
    logger.info(f"Successfully cached {success_count}/{len(repos)} repositories")

    if failed:
        logger.warning(f"Failed to cache {len(failed)} repositories:")
        for repo_url, error in failed:
            logger.warning(f"  {repo_url}: {error}")
        sys.exit(1)
    else:
        logger.info("Cache population complete")
