"""Cli implementation of omni-py via click"""

import click

from omni import __version__
from omni.cli.benchmark import info
from omni.cli.io import storage
from omni.cli.run import run
from omni.cli.soft import software
from omni.cli.utils.logging import logger, configure_logging


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
