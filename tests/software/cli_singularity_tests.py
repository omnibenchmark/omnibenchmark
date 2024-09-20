import shlex

import pytest
from click.testing import CliRunner

import omni.cli.main
from omni.cli.main import cli

runner = CliRunner()

# largely as in https://pytest-with-eric.com/pytest-advanced/pytest-argparse-typer/
test_cases = [
    ("software", ("Manage and install benchmark-specific", 0)),
    (
        "software singularity build --easyconfig nonexistent.eb",
        (
            "Usage: cli software singularity build [OPTIONS]\nTry 'cli software singularity build --help' for help.\n\nError: Invalid value for '-e' / '--easyconfig': Path 'nonexistent.eb' does not exist.\n",
            2,
        ),
    ),
]


# command = "software"
# expected_output = "Usage: cli software singularity build [OPTIONS]\nTry 'cli software singularity build --help' for help.\n\nError: Invalid value for '-e' / '--easyconfig': Path 'nonexistent.eb' does not exist.\n"
@pytest.mark.parametrize("command, expected_output", test_cases)
def test_click(command, expected_output):
    result = runner.invoke(cli, shlex.split(command))
    assert expected_output[0] in result.stdout
    assert result.exit_code == expected_output[1]
