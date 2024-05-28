"""cli commands related to benchmark infos and stats"""

from typing_extensions import Annotated

import typer
from bs4 import BeautifulSoup
from urllib.parse import urlparse
import requests
import re

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


@cli.command("list")
def list_benchmarks(
    endpoint: Annotated[
        str,
        typer.Option(
            "--endpoint",
            "-e",
            help="remote/object storage.",
        ),
    ],
):
    """List all available benchmarks and versions at a specific endpoint"""
    typer.echo(f"Available benchmarks at {endpoint}:")
    url = urlparse(f"{endpoint}/benchmarks")
    response = requests.get(url.geturl(), params={"format": "xml"})
    if response.ok:
        response_text = response.text
    else:
        response.raise_for_status()
    soup = BeautifulSoup(response_text, "xml")
    benchmark_names = [obj.find("Key").text for obj in soup.find_all("Contents")]
    benchmarks = {}
    for benchmark in benchmark_names:
        url = urlparse(f"{endpoint}/{benchmark}.overview")
        response = requests.get(url.geturl(), params={"format": "xml"})
        if response.ok:
            soup = BeautifulSoup(response.text, "xml")
            buckets = [obj.find("Key").text for obj in soup.find_all("Contents")]
            versions = []
            for bucket in buckets:
                if re.search(r"(\d+\.\d+)", bucket):
                    versions.append(bucket)
            from packaging.version import Version

            versions.sort(key=Version)
            if isinstance(versions, list) and len(versions) > 0:
                benchmarks[benchmark] = versions
            elif isinstance(versions, str) and len(versions) > 0:
                benchmarks[benchmark] = [versions]
            else:
                benchmarks[benchmark] = [""]
    for key, value in sorted(benchmarks.items()):
        typer.echo(f"{key:>20} {value[-1]:>8}")


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
