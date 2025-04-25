"""omnibenchmark CLI"""

import click

from omnibenchmark import __version__
from omnibenchmark.cli.benchmark import info
from omnibenchmark.cli.io import storage
from omnibenchmark.cli.run import run
from omnibenchmark.cli.soft import software

from .debug import add_debug_option


@click.group()
@click.version_option(__version__, prog_name="OmniBenchmark CLI")
@click.pass_context
def cli(ctx):
    """
    OmniBenchmark Command Line Interface (CLI).
    """
    ctx.ensure_object(dict)


# Add subcommands to the CLI
cli.add_command(add_debug_option(storage))
cli.add_command(add_debug_option(software))
cli.add_command(add_debug_option(run))
cli.add_command(add_debug_option(info))

add_debug_option(cli)

if __name__ == "__main__":
    cli(obj={})
