"""cli commands related to benchmark infos and stats"""

from pathlib import Path

import click
import yaml
from packaging.version import Version
import sys

from omni.cli.utils.logging import debug_option, logger
from omni.cli.utils.validation import validate_benchmark
from omni.template.repo import create_repo_files, get_missing_repos
from omni.io.code import check_remote_repo_existance


@click.group(name="template")
@click.pass_context
@debug_option
def template(ctx):
    """List benchmarks and/or information about them."""
    ctx.ensure_object(dict)


@template.command("list-missing")
@click.option(
    "--benchmark",
    "-b",
    required=True,
    type=click.Path(exists=True),
    help="Path to benchmark yaml file or benchmark id.",
    envvar="OB_BENCHMARK",
)
@click.option(
    "--reason",
    "-r",
    is_flag=True,
    help="Include reason for missing repository.",
)
@click.pass_context
def list_missing_repos(ctx, benchmark: str, reason: bool = False):
    """List missing repositories in the benchmark."""

    benchmark = validate_benchmark(benchmark, echo=False)
    urls, nodes = get_missing_repos(benchmark)
    for url in urls:
        if reason:
            if check_remote_repo_existance(url):
                reason_str = "missing_commit"
            else:
                reason_str = "missing_repo"
            click.echo(f"{url}\t{reason_str}")
        else:
            click.echo(url)


@template.command("create-missing")
@click.option(
    "--benchmark",
    "-b",
    required=True,
    type=click.Path(exists=True),
    help="Path to benchmark yaml file or benchmark id.",
    envvar="OB_BENCHMARK",
)
@click.option(
    "--url",
    "-u",
    required=True,
    type=str,
    help="URL of git repository, output from 'list-missing'",
)
@click.option(
    "--output_dir",
    "-o",
    required=True,
    type=click.Path(exists=True),
    help="Output directory to store data files.",
)
@click.option(
    "--language",
    "-l",
    required=False,
    type=str,
    help="Programming language of the module.",
)
@click.pass_context
def create_missing_repos(
    ctx, benchmark: str, url: str, output_dir: Path, language: str = None
):
    """List missing repositories in the benchmark."""

    benchmark = validate_benchmark(benchmark, echo=False)
    urls, nodes = get_missing_repos(benchmark)
    selected_node = [n for u, n in zip(urls, nodes) if url in u]
    if len(selected_node) != 1:
        if len(selected_node) > 1:
            logger.error(f"URL {url} found in multiple missing repositories.")
        else:
            logger.error(f"URL {url} not found in missing repositories.")
        sys.exit(1)
    else:
        selected_node = selected_node[0]
        create_repo_files(selected_node, language=language, output_dir=output_dir)
