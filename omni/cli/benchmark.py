"""cli commands related to benchmark infos and stats"""

import datetime
import difflib

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
    if __name__ == "__main__":
        # TODO: for testing until get_bench_definition is implemented
        minio_auth_options_public = {
            "endpoint": "https://omnibenchmark.mls.uzh.ch:9000",
            "secure": False,
        }
        bench_yaml = {
            "auth_options": minio_auth_options_public,
            "storage_type": "minio",
        }
        benchmark = "testversioning2"
        version1 = "0.1"
        version2 = "0.2"
    # setup storage
    ss = get_storage(bench_yaml["storage_type"], bench_yaml["auth_options"], benchmark)

    # get objects for first version
    ss.set_version(version1)
    ss._get_objects()
    files_v1 = [
        f"{f[0]}   {f[1]['size']}   {datetime.datetime.fromisoformat(f[1]['last_modified']).strftime('%Y-%m-%d %H:%M:%S')}\n"
        for f in ss.files.items()
    ]
    creation_time_v1 = datetime.datetime.fromisoformat(
        ss.files[f"versions/{version1}.csv"]["last_modified"]
    ).strftime("%Y-%m-%d %H:%M:%S")

    # get objects for second version
    ss.set_version(version2)
    ss._get_objects()
    files_v2 = [
        f"{f[0]}   {f[1]['size']}   {datetime.datetime.fromisoformat(f[1]['last_modified']).strftime('%Y-%m-%d %H:%M:%S')}\n"
        for f in ss.files.items()
    ]
    creation_time_v2 = datetime.datetime.fromisoformat(
        ss.files[f"versions/{version2}.csv"]["last_modified"]
    ).strftime("%Y-%m-%d %H:%M:%S")

    # diff the two versions
    typer.echo(
        "".join(
            list(
                difflib.unified_diff(
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
