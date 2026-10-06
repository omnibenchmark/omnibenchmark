"""`ob config`: read and write omnibenchmark.cfg, git-style."""

import sys

import click

from omnibenchmark.config import config as cfg


@click.command(name="config")
@click.argument("key", required=False)
@click.argument("value", required=False)
@click.option("--list", "-l", "list_", is_flag=True, help="List all settings.")
@click.option("--unset", is_flag=True, help="Remove KEY from the config file.")
def config(key, value, list_, unset):
    """Get or set a setting, e.g. `ob config user.name "Ada Lovelace"`.

    KEY is `section.name`. With VALUE, the setting is written to the config file;
    without, its current value is printed. `--unset` removes it.
    """
    if list_:
        for section in cfg.sections():
            for name in cfg.options(section):
                click.echo(f"{section}.{name}={cfg.get(section, name)}")
        return
    if not key or "." not in key:
        raise click.UsageError("KEY must be of the form section.name")
    section, name = key.split(".", 1)
    if unset:
        if not cfg.unset(section, name):
            sys.exit(5)  # git's exit code for unsetting a missing key
        cfg.save()
        return
    if value is None:
        current = cfg.get(section, name)
        if current is None:
            sys.exit(1)
        click.echo(current)
        return
    cfg.set(section, name, value)
    cfg.save()
