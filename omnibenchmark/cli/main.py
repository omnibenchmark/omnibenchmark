"""Cli implementation of omni-py via click"""

import click

from omnibenchmark import __version__
from omnibenchmark.cli.benchmark import info
from omnibenchmark.cli.io import storage
from omnibenchmark.cli.run import run
from omnibenchmark.cli.soft import software
from omnibenchmark.cli.utils.logging import configure_logging


@click.group()
@click.option(
    "--debug/--no-debug", default=False, help="Enable debug mode for detailed logging."
)
@click.version_option(__version__, prog_name="OmniBenchmark CLI")
@click.pass_context
def cli(ctx, debug):
    """
    OmniBenchmark Command Line Interface (CLI).
    """

    ctx.ensure_object(dict)
    ctx.obj["DEBUG"] = debug

    configure_logging(debug)


cli.add_command(storage)
cli.add_command(software)
cli.add_command(run)
cli.add_command(info)

if __name__ == "__main__":
    cli(obj={})
