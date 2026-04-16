"""Regression tests for IDs that start with a digit being rejected by Snakemake.

Snakemake rule names are Python identifiers — they cannot begin with a digit.
When an omnibenchmark stage/module ID starts with a digit the generated Snakefile
will contain an invalid rule name and Snakemake will refuse to parse it.

Structure:
- Root-cause tests (currently pass): prove Snakemake rejects digit-prefixed rule names.
- Failing test (currently fails, xfail): proves omnibenchmark does NOT yet reject
  digit-starting IDs at parse time.  Remove xfail once the fix lands.

Issue: omnibenchmark/omnibenchmark#312
"""

import subprocess
import sys
from pathlib import Path

import pytest

from omnibenchmark.model.benchmark import Benchmark
from omnibenchmark.model.validation import BenchmarkParseError


def _snakemake_bin() -> str:
    """Return the snakemake binary co-located with the current Python interpreter."""
    candidate = Path(sys.executable).parent / "snakemake"
    if candidate.exists():
        return str(candidate)
    import shutil

    return shutil.which("snakemake") or "snakemake"


def _run_snakemake_dryrun(snakefile: Path) -> subprocess.CompletedProcess:
    cmd = [
        _snakemake_bin(),
        "--snakefile",
        str(snakefile),
        "--dryrun",
        "--cores",
        "1",
        "--quiet",
    ]
    return subprocess.run(cmd, capture_output=True, text=True, cwd=snakefile.parent)


# ---------------------------------------------------------------------------
# Snakefile fixtures
# ---------------------------------------------------------------------------

SNAKEFILE_DIGIT_START = """\
rule all:
    input: "out.txt"

rule 1invalid_stage__mod__default:
    output: "out.txt"
    shell: "touch {output}"
"""

SNAKEFILE_VALID = """\
rule all:
    input: "out.txt"

rule valid_stage__mod__default:
    output: "out.txt"
    shell: "touch {output}"
"""


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_snakemake_rejects_rule_name_starting_with_digit(tmp_path):
    """Snakemake must reject a Snakefile whose rule name begins with a digit.

    This is the root-cause check: if omnibenchmark allows a stage/module id
    like '1stage', the generated Snakefile will contain 'rule 1stage__...'
    and Snakemake will fail.  The fix is to reject such IDs at parse time so
    the Snakefile is never generated in the first place.
    """
    snakefile = tmp_path / "Snakefile"
    snakefile.write_text(SNAKEFILE_DIGIT_START)

    result = _run_snakemake_dryrun(snakefile)

    assert result.returncode != 0, (
        "Expected Snakemake to fail on a rule name starting with a digit, "
        f"but it exited with 0.\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    )


def test_snakemake_accepts_rule_name_with_valid_identifier(tmp_path):
    """Control: Snakemake must accept a Snakefile with a valid rule name.

    Ensures the test infrastructure itself is working and the rejection
    in the sibling test is due to the digit prefix, not an environment issue.
    """
    snakefile = tmp_path / "Snakefile"
    snakefile.write_text(SNAKEFILE_VALID)

    result = _run_snakemake_dryrun(snakefile)

    assert result.returncode == 0, (
        "Expected Snakemake to accept a rule name with a valid identifier, "
        f"but it exited with {result.returncode}.\n"
        f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    )


# ---------------------------------------------------------------------------
# Failing test — documents the missing fix (xfail until #312 is resolved)
# ---------------------------------------------------------------------------

# Minimal benchmark YAML whose stage id starts with a digit.  Snakemake will
# choke on the generated rule name; omnibenchmark should catch it first.
BENCHMARK_YAML_DIGIT_STAGE_ID = """\
id: test_benchmark
version: "1.0"
benchmark_yaml_spec: "0.3"
benchmarker: "Test User"
software_backend: host

software_environments:
  host:
    description: "host env"

stages:
  - id: 1data
    modules:
      - id: mod1
        software_environment: host
        repository:
          url: https://github.com/example/dummy.git
          commit: abc1234
    outputs:
      - id: data.out
        path: "out.txt"
"""

BENCHMARK_YAML_DIGIT_MODULE_ID = """\
id: test_benchmark
version: "1.0"
benchmark_yaml_spec: "0.3"
benchmarker: "Test User"
software_backend: host

software_environments:
  host:
    description: "host env"

stages:
  - id: data
    modules:
      - id: 1mod
        software_environment: host
        repository:
          url: https://github.com/example/dummy.git
          commit: abc1234
    outputs:
      - id: data.out
        path: "out.txt"
"""


def test_benchmark_rejects_stage_id_starting_with_digit():
    """Parsing a benchmark whose stage id starts with a digit must raise an error.

    Without this check the invalid id propagates into the generated Snakefile as
    'rule 1data__mod1__...' and Snakemake fails with an unhelpful SyntaxError.
    The fix catches this at parse time inside omnibenchmark itself.
    """
    with pytest.raises(BenchmarkParseError, match="must not start with a digit"):
        Benchmark.from_yaml(BENCHMARK_YAML_DIGIT_STAGE_ID)


def test_benchmark_rejects_module_id_starting_with_digit():
    """Parsing a benchmark whose module id starts with a digit must raise an error."""
    with pytest.raises(BenchmarkParseError, match="must not start with a digit"):
        Benchmark.from_yaml(BENCHMARK_YAML_DIGIT_MODULE_ID)
