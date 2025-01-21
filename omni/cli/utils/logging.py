import logging
import sys
import click

logger = logging.getLogger("omnibenchmark")


# Function to configure logging
def configure_logging(debug: bool):
    """
    Configures the logging system based on the debug flag.
    """
    handler = logging.StreamHandler(sys.stdout)
    # formatter = logging.Formatter(
    #     "%(asctime)s - [%(levelname)s] - %(name)s - %(message)s"
    # )
    formatter = logging.Formatter("%(message)s")
    handler.setFormatter(formatter)

    log_level = logging.DEBUG if debug else logging.INFO
    logger.setLevel(log_level)

    if not logger.hasHandlers():
        logger.addHandler(handler)


def debug_option(f):
    """
    Decorator to pass debug configuration to subcommands.

    This ensures that subcommands can access the debug flag
    and logging configuration from the parent context.
    """

    @click.pass_context
    def wrapper(ctx: click.Context, *args, **kwargs):
        # Get the debug flag from the parent context, if available
        debug = ctx.parent.obj.get("DEBUG", False) if ctx.parent else False
        # Configure logging based on the debug flag
        configure_logging(debug)

        # Invoke the command function with the context
        return ctx.invoke(f, *args, **kwargs)

    return wrapper
