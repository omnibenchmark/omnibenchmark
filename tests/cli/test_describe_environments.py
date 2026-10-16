"""`ob describe environments`: ids, backends, usage and metadata."""

import textwrap

import pytest
from click.testing import CliRunner

from omnibenchmark.cli.describe import describe

MODULES = """\
stages:
  - id: data
    modules:
      - id: D1
        software_environment: py
        repository: {url: https://example.org/d.git, commit: abc1234}
      - id: D2
        software_environment: py
        repository: {url: https://example.org/d.git, commit: abc1234}
    outputs:
      - {id: data_counts, path: "{dataset}.txt"}
  - id: methods
    modules:
      - id: M1
        software_environment: r
        repository: {url: https://example.org/m.git, commit: abc1234}
    inputs: [data_counts]
    outputs:
      - {id: methods_out, path: "{dataset}.out"}
"""

HEADER = """\
id: b
benchmarker: me
version: "1.0"
"""


def _run(path):
    return CliRunner().invoke(describe, ["environments", str(path)])


@pytest.mark.short
def test_v09_metadata_backends_and_usage(tmp_path):
    plan = tmp_path / "plan.yaml"
    plan.write_text(
        HEADER
        + textwrap.dedent(
            """\
            api_version: "0.9.0"
            software_backend: conda
            software_environments:
              py:
                - {type: conda, uri: envs/py.yaml}
                - {type: apptainer, uri: oras://registry/py.sif}
              r:
                - {type: apptainer, uri: oras://registry/r.sif}
            software_environments_metadata:
              py: {name: Python, description: "3.12 + numpy"}
            """
        )
        + MODULES
    )
    result = _run(plan)
    assert result.exit_code == 0, result.output
    lines = result.output.splitlines()
    assert lines[0] == "Software backend: conda"
    py = next(line for line in lines if line.startswith("py "))
    r = next(line for line in lines if line.startswith("r "))
    assert "Python" in py and "3.12 + numpy" in py
    assert "apptainer, conda*" in py  # active backend marked
    assert "D1, D2" in py
    assert "apptainer" in r and "*" not in r  # no recipe for the active backend
    assert "M1" in r
    assert "* active backend" in result.output


@pytest.mark.short
def test_legacy_description_shown(tmp_path):
    plan = tmp_path / "plan.yaml"
    plan.write_text(
        HEADER
        + textwrap.dedent(
            """\
            api_version: "0.7.0"
            software_backend: apptainer
            software_environments:
              py: {description: "legacy py", apptainer: oras://registry/py.sif}
              r: {name: R env, apptainer: oras://registry/r.sif}
            """
        )
        + MODULES
    )
    result = _run(plan)
    assert result.exit_code == 0, result.output
    assert "legacy py" in result.output
    assert "R env" in result.output
    assert "apptainer*" in result.output


@pytest.mark.short
def test_host_without_environments(tmp_path):
    plan = tmp_path / "plan.yaml"
    plan.write_text(
        HEADER
        + 'api_version: "0.9.0"\nsoftware_backend: host\n'
        + MODULES.replace("        software_environment: py\n", "").replace(
            "        software_environment: r\n", ""
        )
    )
    result = _run(plan)
    assert result.exit_code == 0, result.output
    assert "Software backend: host" in result.output
    assert "No software environments declared." in result.output


@pytest.mark.short
def test_parse_error_reported(tmp_path):
    plan = tmp_path / "plan.yaml"
    plan.write_text(
        HEADER
        + 'api_version: "0.9.0"\nsoftware_backend: conda\n'
        + "software_environments:\n  py:\n    - {type: pip, uri: x}\n"
        + MODULES
    )
    result = _run(plan)
    assert result.exit_code == 1
    assert "unknown backend type 'pip'" in result.output
