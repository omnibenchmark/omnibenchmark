"""cli commands related to benchmark infos and stats"""

from packaging.version import Version
from pathlib import Path

import click
import yaml

from datetime import datetime
from difflib import unified_diff

from omnibenchmark.benchmark import Benchmark
from omnibenchmark.cli.utils.logging import logger
from omnibenchmark.cli.utils.validation import validate_benchmark
from omnibenchmark.io.utils import get_storage, remote_storage_args


@click.group(name="info")
@click.pass_context
# @debug_option
def info(ctx):
    """List benchmarks and/or information about them."""
    ctx.ensure_object(dict)


@info.command("diff")
@click.option(
    "--benchmark",
    "-b",
    required=True,
    type=click.Path(exists=True),
    help="Path to benchmark yaml file or benchmark id.",
    envvar="OB_BENCHMARK",
)
@click.option(
    "--version1",
    "-v1",
    required=True,
    type=str,
    help="Reference version.",
)
@click.option(
    "--version2",
    "-v2",
    required=True,
    type=str,
    help="Version to compare with.",
)
@click.pass_context
def diff_benchmark(ctx, benchmark, version1, version2):
    """Show differences between 2 benchmark versions."""
    logger.info(
        f"Found the following differences in {benchmark} for {version1} and {version2}."
    )
    with open(benchmark, "r") as fh:
        yaml.safe_load(fh)
        benchmark = Benchmark(Path(benchmark))

    auth_options = remote_storage_args(benchmark)

    # setup storage
    ss = get_storage(
        str(benchmark.converter.model.storage_api),
        auth_options,
        str(benchmark.converter.model.storage_bucket_name),
    )

    # get objects for first version
    ss.set_version(version1)
    ss._get_objects()
    files_v1 = [
        f"{f[0]}   {f[1]['size']}   {datetime.fromisoformat(f[1]['last_modified']).strftime('%Y-%m-%d %H:%M:%S')}\n"
        for f in ss.files.items()
    ]
    if f"versions/{version1}.csv" in ss.files.keys():
        creation_time_v1 = datetime.fromisoformat(
            ss.files[f"versions/{version1}.csv"]["last_modified"]
        ).strftime("%Y-%m-%d %H:%M:%S")
    else:
        creation_time_v1 = ""

    # get objects for second version
    ss.set_version(version2)
    ss._get_objects()
    files_v2 = [
        f"{f[0]}   {f[1]['size']}   {datetime.fromisoformat(f[1]['last_modified']).strftime('%Y-%m-%d %H:%M:%S')}\n"
        for f in ss.files.items()
    ]
    if f"versions/{version2}.csv" in ss.files.keys():
        creation_time_v2 = datetime.fromisoformat(
            ss.files[f"versions/{version2}.csv"]["last_modified"]
        ).strftime("%Y-%m-%d %H:%M:%S")
    else:
        creation_time_v2 = ""

    # diff the two versions
    click.echo(
        "".join(
            list(
                unified_diff(
                    files_v1,
                    files_v2,
                    fromfile=f"version {version1}",
                    tofile=f"version {version2}",
                    fromfiledate=creation_time_v1,
                    tofiledate=creation_time_v2,
                )
            )
        )
    )


@info.command("list-versions")
@click.option(
    "--benchmark",
    "-b",
    required=True,
    type=click.Path(exists=True),
    help="Path to benchmark yaml file or benchmark id.",
    envvar="OB_BENCHMARK",
)
@click.pass_context
def list_versions(ctx, benchmark):
    """List all available benchmarks versions at a specific endpoint."""
    logger.info(f"Available versions of {benchmark}:")

    with open(benchmark, "r") as fh:
        yaml.safe_load(fh)
        benchmark = Benchmark(Path(benchmark))

    auth_options = remote_storage_args(benchmark)

    # setup storage
    ss = get_storage(
        str(benchmark.converter.model.storage_api),
        auth_options,
        str(benchmark.converter.model.storage_bucket_name),
    )

    if len(ss.versions) > 0:
        if len(ss.versions) > 1:
            ss.versions.sort(key=Version)
        for version in ss.versions:
            click.echo(f"{str(version):>8}")


@info.command("computational")
@click.option(
    "--benchmark",
    "-b",
    required=True,
    type=click.Path(exists=True),
    help="Path to benchmark yaml file or benchmark id.",
    envvar="OB_BENCHMARK",
)
@click.pass_context
def computational_graph(ctx, benchmark: str):
    """Export computational graph to dot format."""

    b = validate_benchmark(benchmark, echo=False)
    if b is not None:
        dot = b.export_to_dot()
        click.echo(dot.to_string())


@info.command("topology")
@click.option(
    "--benchmark",
    "-b",
    required=True,
    type=click.Path(exists=True),
    help="Path to benchmark yaml file or benchmark id.",
    envvar="OB_BENCHMARK",
)
@click.pass_context
def plot_topology(ctx, benchmark: str):
    """Export benchmark topology to mermaid diagram format."""

    b = validate_benchmark(benchmark, echo=False)
    if b is not None:
        mermaid = b.export_to_mermaid()
        click.echo(mermaid)
