"""`ob slice`: pack part of a run, seed a new run from it. SPIKE, local paths only."""

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
    """Pack and reuse slices of a run (experimental)."""


@slice_group.command()
@_out
@click.option("--to", "dest", required=True, type=click.Path(path_type=Path))
@click.option("--until", help="Cut after this stage.")
@click.option("--dataset", help="Keep this dataset's subtree.")
@click.option("--squashfs", is_flag=True, help="Write a squashfs image, not tar.gz.")
@click.option("--base-until", help="Leave out what the slice up to this stage holds.")
def pack(out_dir, dest, until, dataset, squashfs, base_until):
    """Pack a slice of OUT_DIR into the directory TO."""
    try:
        click.echo(slices.pack(out_dir, dest, until, dataset, squashfs, base_until))
    except (ValueError, KeyError, OSError, RuntimeError) as e:
        raise click.ClickException(str(e))


@slice_group.command()
@_out
@click.option(
    "--from",
    "source",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--squashfs",
    is_flag=True,
    help="Mount a squashfs image read-only when there is one; else use tar.gz.",
)
def seed(out_dir, source, squashfs):
    """Seed OUT_DIR from matching slices in FROM. Run `ob run --dry` first."""
    used = slices.seed(out_dir, source, squashfs)
    click.echo("\n".join(used) if used else "no matching slice")


@slice_group.command()
@_out
def unmount(out_dir):
    """Unmount the squashfs slices of OUT_DIR."""
    for name in slices.unmount(out_dir):
        click.echo(f"unmounted {name}")
