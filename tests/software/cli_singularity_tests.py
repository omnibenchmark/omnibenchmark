import shlex

import pytest
from click.testing import CliRunner

from omni.cli.main import cli

runner = CliRunner()

# largely as in https://pytest-with-eric.com/pytest-advanced/pytest-argparse-typer/
test_cases = [
    (
        "software singularity build --easyconfig nonexistent.eb",
        "\nERROR: easyconfig not found.",
    ),
]


# command = "software"
# expected_output = "Usage: cli software singularity build [OPTIONS]\nTry 'cli software singularity build --help' for help.\n\nError: Invalid value for '-e' / '--easyconfig': Path 'nonexistent.eb' does not exist.\n"
@pytest.mark.parametrize("command, expected_output", test_cases)
def test_click(command, expected_output, capture_logs):
    runner.invoke(cli, shlex.split(command))
    log_output = capture_logs.getvalue()

    assert expected_output in log_output
