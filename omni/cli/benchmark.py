"""cli commands related to benchmark infos and stats"""

import typer
from packaging.version import Version
from typing_extensions import Annotated

from omni.io.utils import get_storage

cli = typer.Typer(add_completion=False)


@cli.command("archive")
def archive_benchmark(
    benchmark: Annotated[
        str,
        typer.Option(
            "--benchmark", "-b", help="Path to benchmark yaml file or benchmark id."
        ),
    ],
):
    """Archive a benchmark (changes to read-only permissions)"""
    typer.echo(f"Archiving benchmark with yaml {benchmark}", err=True)


@cli.command("cite")
def cite_benchmark(
    benchmark: Annotated[
        str,
        typer.Option(
            "--benchmark",
            "-b",
            help="Path to benchmark yaml file or benchmark id.",
        ),
    ],
):
    """Get the citation for a specific benchmark."""
    typer.echo(f"Citation for benchmark: {benchmark}")


@cli.command("diff")
def diff_benchmark(
    benchmark: Annotated[
        str,
        typer.Option(
            "--benchmark",
            "-b",
            help="Path to benchmark yaml file or benchmark id.",
        ),
    ],
    version1: Annotated[
        str,
        typer.Option(
            "--version1",
            "-v1",
            help="reference version",
        ),
    ],
    version2: Annotated[
        str,
        typer.Option(
            "--version2",
            "-v2",
            help="version to compare with.",
        ),
    ],
):
    """Show differences between 2 benchmark versions."""
    typer.echo(
        f"Found the following differences in {benchmark} for {version1} and {version2}."
    )


@cli.command("list versions")
def list_versions(
    benchmark: Annotated[
        str,
        typer.Option(
            "--benchmark",
            "-b",
            help="Path to benchmark yaml file or benchmark id.",
        ),
    ],
    endpoint: Annotated[
        str,
        typer.Option(
            "--endpoint",
            "-e",
            help="remote/object storage.",
        ),
    ],
):
    """List all available benchmarks versions at a specific endpoint."""
    typer.echo(f"Available versions of {benchmark} at {endpoint}:")
    # TODO: for testing until get_bench_definition is implemented
    if __name__ == "__main__":
        minio_auth_options_public = {
            "endpoint": "https://omnibenchmark.mls.uzh.ch:9000",
            "secure": False,
        }
        bench_yaml = {
            "auth_options": minio_auth_options_public,
            "storage_type": "minio",
        }
        benchmark = "testversioning"
        version = "0.1"

    ss = get_storage(bench_yaml["storage_type"], bench_yaml["auth_options"], benchmark)
    if len(ss.versions) > 0:
        ss.versions.sort(key=Version)
        for version in ss.versions:
            typer.echo(f"{version:>8}")
