"""Cli implementation of omni-py via click"""

import click

from omni.cli.benchmark import info
from omni.cli.io import storage
from omni.cli.run import run
from omni.cli.soft import software


@click.group()
@click.option("--debug/--no-debug", default=False)
@click.pass_context
def cli(ctx, debug):
    ctx.ensure_object(dict)

    ctx.obj["DEBUG"] = debug


cli.add_command(storage)
cli.add_command(software)
cli.add_command(run)
cli.add_command(info)
