import shlex

import pytest
from click.testing import CliRunner

from omnibenchmark.cli.main import cli

runner = CliRunner()

# largely as in https://pytest-with-eric.com/pytest-advanced/pytest-argparse-typer/
test_cases = [
    (
        "software module build -e zlib-1.2.10.eb",
        "Installing software for zlib-1.2.10.eb using easybuild. It will take some time",
    ),
]


@pytest.mark.parametrize("command, expected_output", test_cases)
def test_click(command, expected_output, capture_logs):
    result = runner.invoke(cli, shlex.split(command))
    log_output = capture_logs.getvalue()

    assert expected_output in log_output
    assert result.exit_code == 0
