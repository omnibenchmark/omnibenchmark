import os
import re
import textwrap

import pytest

from tests.cli.cli_setup import OmniCLISetup

from .asserts import assert_in_output
from .path import data

# TODO: deprecate fixtures in this module
from ..fixtures import rustfs_storage, _rustfs_container, bundled_repos  # noqa: F401


# TODO: mark as integration
def test_remote(rustfs_storage):  # noqa: F811
    # TODO(ben): the technique of expecting YAML validation in the output is a bit brittle, we could
    # check e.g. that output has been produced.
    # But we should be changing the testing strategy in a gradual way
    expected1 = "Benchmark YAML file integrity check passed."
    expected2 = "Benchmark run completed successfully."

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(rustfs_storage.benchmark_file),
            ]
        )

        print(result.stdout)
        assert result.returncode == 0
        assert_in_output(result.stdout, expected1)
        assert_in_output(result.stdout, expected2)


def test_benchmark_not_found():
    expected = """Error: Invalid value for 'BENCHMARK'"""
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(data / "does_not_exist.yaml"),
            ]
        )
        assert result.returncode == 2
        assert_in_output(result.stderr, expected)


def test_benchmark_format_incorrect():
    # TODO(ben): this should better be a config parsing test
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(data / "benchmark_format_incorrect.yaml"),
            ]
        )
        assert result.returncode == 1


def test_benchmark_software_does_not_exist():
    expected = """
    Failed to load benchmark: Software environment with id 'python' does not have a valid backend definition for: 'conda'.
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(data / "benchmark_software_does_not_exist.yaml"),
            ]
        )

        assert result.returncode == 1
        assert_in_output(result.stdout, expected)


def test_local(tmp_path):
    # Check that benchmark runs successfully (may have deprecation warnings)
    expected = "Benchmark run completed successfully."

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(data / "mock_benchmark.yaml"),
            ],
            cwd=tmp_path,
        )

        assert result.returncode == 0
        assert_in_output(result.stdout, expected)


def test_custom_out_dir(tmp_path):
    # Check that benchmark runs successfully with custom output directory
    expected = "Benchmark run completed successfully."

    custom_out_dir = "out_2313_custom"

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(data / "mock_benchmark.yaml"),
                "--out-dir",
                custom_out_dir,
            ],
            cwd=tmp_path,
        )

        assert result.returncode == 0
        assert_in_output(result.stdout, expected)

        assert os.path.exists(tmp_path / custom_out_dir)


def test_local_dry():
    # Dry run should complete successfully (may have warnings)
    expected_output = "Snakefile generated."
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                str(data / "mock_benchmark.yaml"),
                "--dry",
            ]
        )

        assert result.returncode == 0
        assert_in_output(result.stdout, expected_output)


@pytest.mark.skip(
    reason="bundle-based repos not supported by the new explicit-snakefile resolver; "
    "requires migration of test data to real git repos"
)
def test_benchmark_does_fail_if_one_module_fails(bundled_repos, tmp_path):  # noqa: F811
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                (data / "benchmark_failing_module.yaml").as_posix(),
            ],
            cwd=tmp_path,
        )

        failed_msg = "Benchmark run failed"

        assert failed_msg in result.stdout or failed_msg in result.stderr
        assert result.returncode != 0


@pytest.mark.skip(
    reason="bundle-based repos not supported by the new explicit-snakefile resolver; "
    "requires migration of test data to real git repos"
)
def test_benchmark_ok_if_one_module_fails_with_continue(tmp_path, bundled_repos):  # noqa: F811
    """
    This test checks that the benchmark does not fail if one module fails with continue-on-error.
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                (data / "benchmark_failing_module.yaml").as_posix(),
                "--continue-on-error",
            ],
            cwd=tmp_path,
        )

        failed_msg = "Benchmark run failed"

        assert failed_msg not in result.stdout
        assert failed_msg not in result.stderr
        assert result.returncode == 0


def test_params_key_in_output_path_resolved_in_snakefile(tmp_path, bundled_repos):  # noqa: F811
    """
    Integration test: {params.key} in output path templates must be fully
    substituted to concrete values before being written into the Snakefile.
    Snakemake must never see a {params.*} wildcard in an output: block.

    This test verifies the full pipeline:
      YAML parsing → param expansion → TemplateContext.substitute() → Snakefile generation
    """
    out_dir = tmp_path / "out"
    benchmark_yaml = tmp_path / "bench_params_output.yaml"
    benchmark_yaml.write_text(
        textwrap.dedent("""
            id: test_params_in_output
            version: "1.0"
            benchmarker: "test"
            api_version: "0.5.0"
            software_backend: host
            software_environments:
              host:
                description: "host env"
            stages:
              - id: data
                modules:
                  - id: D1
                    software_environment: host
                    repository:
                      url: bundles/dummymodule_4ff8427.bundle
                      commit: 4ff8427
                    parameters:
                      - k: "3"
                      - k: "5"
                outputs:
                  - id: data.result
                    path: "{name}_k{params.k}_result.txt"
        """).strip()
    )

    with OmniCLISetup() as omni:
        result = omni.call(
            ["run", str(benchmark_yaml), "--dry", "--out-dir", str(out_dir)],
            cwd=str(tmp_path),
        )

    assert (
        result.returncode == 0
    ), f"CLI failed (rc={result.returncode}):\nSTDOUT: {result.stdout}\nSTDERR: {result.stderr}"

    snakefile = out_dir / "Snakefile"
    assert snakefile.exists(), "Snakefile was not generated"
    content = snakefile.read_text()

    # Extract all output: blocks and verify no {params.*} placeholders remain —
    # Snakemake would misinterpret them as wildcards.
    output_blocks = re.findall(
        r"\boutput:\s*\n(.*?)(?=\n[ \t]+\w|\nrule |\Z)", content, re.DOTALL
    )
    assert output_blocks, "No output: blocks found in Snakefile"
    for block in output_blocks:
        assert (
            "{params." not in block
        ), f"Unresolved {{params.*}} placeholder in Snakefile output: block:\n{block}"

    # Both parameter variants must produce distinct, concrete output filenames.
    assert (
        "_k3_result.txt" in content
    ), "Expected k=3 output path not found in Snakefile"
    assert (
        "_k5_result.txt" in content
    ), "Expected k=5 output path not found in Snakefile"
