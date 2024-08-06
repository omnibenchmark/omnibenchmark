import pytest
from typer.testing import CliRunner
import shlex
import omni.cli.main

from omni.cli.main import cli

runner = CliRunner()

# largely as in https://pytest-with-eric.com/pytest-advanced/pytest-argparse-typer/
test_cases = [
    ("software conda pin --env data/ripgrep.yml", "DONE: Pinned data/ripgrep.yml"),
]


@pytest.mark.parametrize("command, expected_output", test_cases)
def test_typer(command, expected_output):
    result = runner.invoke(cli, shlex.split(command))
    assert expected_output in result.stdout
    assert result.exit_code == 0
