"""Shared logging setup: the ``omnibenchmark`` logger and its configuration.

A leaf module — depends on nothing in the package, imported by any layer that
needs to log. (Formerly ``omnibenchmark.cli.utils.logging``.)
"""

import logging
import sys


logger = logging.getLogger("omnibenchmark")


def configure_logging(debug: bool):
    """
    Configures the logging system based on the debug flag.
    """
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter("%(message)s")
    handler.setFormatter(formatter)

    log_level = logging.DEBUG if debug else logging.INFO
    logger.setLevel(log_level)

    if not logger.hasHandlers():
        logger.addHandler(handler)
