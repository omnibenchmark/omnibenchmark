from functools import wraps

import click

from .utils.logging import configure_logging


def add_debug_option(cmd):
    """Decorator to add debug option to commands and groups"""
    if isinstance(cmd, click.Command):
        # For existing commands/groups
        if not any(param.name == "debug" for param in cmd.params):
            cmd.params.insert(
                0,
                click.Option(
                    ["--debug/--no-debug"],
                    is_eager=True,
                    expose_value=False,
                    callback=lambda ctx, param, value: _set_debug(ctx, value),
                    help="Enable debug mode",
                ),
            )
        return cmd

    # For functions that will become commands/groups
    @click.option(
        "--debug/--no-debug",
        is_eager=True,
        expose_value=False,
        callback=lambda ctx, param, value: _set_debug(ctx, value),
        help="Enable debug mode",
    )
    @wraps(cmd)
    def wrapper(*args, **kwargs):
        return cmd(*args, **kwargs)

    return wrapper


def _set_debug(ctx, value: bool):
    """Callback function for debug flag"""
    # Get the root context
    root_ctx = ctx.find_root()
    root_ctx.ensure_object(dict)

    # Get command path depth
    cmd_depth = len(ctx.command_path.split())

    # Initialize DEBUG if not set
    if "DEBUG" not in root_ctx.obj:
        root_ctx.obj["DEBUG"] = False

    # Allow setting debug to True from any level
    # Allow setting debug to False only from the level that set it to True
    if value is True or cmd_depth == 1:
        root_ctx.obj["DEBUG"] = value

    configure_logging(root_ctx.obj["DEBUG"])
    return root_ctx.obj["DEBUG"]
