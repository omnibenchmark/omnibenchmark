"""omnibenchmark CLI"""

import click

from omnibenchmark import __version__
from omnibenchmark.cli.archive import archive
from omnibenchmark.cli.cite import cite
from omnibenchmark.cli.create import create
from omnibenchmark.cli.dashboard import dashboard
from omnibenchmark.cli.describe import describe
from omnibenchmark.cli.remote import remote
from omnibenchmark.cli.run import run
from omnibenchmark.cli.validate import validate

from .debug import add_debug_option


def format_recursive_help(ctx, param, value):
    """Custom help formatter that shows all subcommands and their sub-subcommands"""
    if not value or ctx.resilient_parsing:
        return

    click.echo("Usage: ob [OPTIONS] COMMAND [ARGS]...")
    click.echo("")
    click.echo("  OmniBenchmark Command Line Interface (CLI).")
    click.echo("")
    click.echo("Options:")
    click.echo("  --debug / --no-debug  Enable debug mode")
    click.echo("  --version             Show the version and exit.")
    click.echo("  --help                Show this message and exit.")
    click.echo("")
    click.echo("Commands:")

    # Get the main CLI group
    main_cli = ctx.find_root().command

    # Format each main command and its subcommands
    for name, command in main_cli.commands.items():
        click.echo(f"  {name:<12} {command.get_short_help_str(50)}")

        # If it's a group, show its subcommands too
        if hasattr(command, "commands"):
            for subname, subcommand in command.commands.items():
                click.echo(
                    f"    {name} {subname:<10} {subcommand.get_short_help_str(45)}"
                )

    ctx.exit()


@click.group()
@click.version_option(__version__, prog_name="OmniBenchmark CLI")
@click.option(
    "--help",
    "-h",
    is_flag=True,
    expose_value=False,
    is_eager=True,
    callback=format_recursive_help,
    help="Show this message and exit.",
)
@click.pass_context
def cli(ctx):
    """
    OmniBenchmark Command Line Interface (CLI).
    """
    ctx.ensure_object(dict)


# Add subcommands to the CLI
cli.add_command(add_debug_option(cite))
cli.add_command(add_debug_option(create))
cli.add_command(add_debug_option(describe))
cli.add_command(add_debug_option(remote))
cli.add_command(add_debug_option(run))
cli.add_command(add_debug_option(validate))
cli.add_command(archive)
cli.add_command(dashboard)

add_debug_option(cli)

if __name__ == "__main__":
    cli(obj={})
