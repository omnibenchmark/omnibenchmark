"""cli commands related to benchmark infos and stats"""

from pathlib import Path

import click

from omnibenchmark.benchmark import BenchmarkExecution


@click.group(name="describe")
@click.pass_context
# @debug_option
def describe(ctx):
    """Describe benchmarks and/or information about them."""
    ctx.ensure_object(dict)


@describe.command("computational")
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
    b = BenchmarkExecution(benchmark_yaml=Path(benchmark))
    if b is None:
        return
    dot = b.export_to_dot()
    click.echo(dot.to_string())


@describe.command("topology")
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
    b = BenchmarkExecution(benchmark_yaml=Path(benchmark))
    if b is None:
        return
    mermaid = b.export_to_mermaid()
    click.echo(mermaid)
