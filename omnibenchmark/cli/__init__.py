"""
Command-line interface for omnibenchmark.

Represents: the user-facing ``ob`` commands; each module maps to a subcommand
and wires user input to the layers below.
Layer: cli (top — nothing imports it)
Depends on: everything below (core, model, remote, archive, ...).

The top, user-facing layer: each module maps to an ``ob`` subcommand and wires
user input to the lower layers. It depends downward on everything it needs
(:mod:`omnibenchmark.core`, :mod:`omnibenchmark.model`,
:mod:`omnibenchmark.remote`, ...) and nothing imports back up into it.

Subcommands:
- ``main.py``     — entry point and the root ``ob`` command group
- ``run.py``      — execute a benchmark
- ``create.py``   — scaffold a new benchmark
- ``validate.py`` — validate a benchmark definition
- ``describe.py`` / ``dashboard.py`` — inspect and report on a benchmark
- ``cite.py``     — emit citation metadata
- ``archive.py`` / ``remote.py`` — archiving and remote-storage operations
- ``debug.py``    — debug helpers and logging configuration

Helpers ``formatting.py`` (CLI parse-error rendering) live alongside.
"""
