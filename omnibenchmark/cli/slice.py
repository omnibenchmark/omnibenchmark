"""`ob slice`: housekeeping for the snapshots a run was started from.
`ob archive --snapshot` writes them, `ob run --use-snapshots` reads them."""

from pathlib import Path

import click

from omnibenchmark import slices

_out = click.option(
    "--out-dir",
    default="out",
    show_default=True,
    type=click.Path(file_okay=False, path_type=Path),
)


@click.group(name="slice")
def slice_group():
    """Snapshots in use by an out dir (experimental)."""


@slice_group.command()
@_out
def unmount(out_dir):
    """Unmount the squashfs slices of OUT_DIR."""
    for name in slices.unmount(out_dir):
        click.echo(f"unmounted {name}")
