import pytest
from pathlib import Path

from tests.e2e.common import E2ETestRunner, run_standard_pipeline_test


CONFIG_FILENAME = "09_named_outputs_v06.yaml"
TEST_NAME = "09_named_outputs_v06"


@pytest.fixture
def named_outputs_config():
    return Path(__file__).parent / "configs" / CONFIG_FILENAME


@pytest.mark.e2e
def test_named_outputs_pipeline(
    named_outputs_config, tmp_path, bundled_repos, keep_files
):
    """v0.6 named-outputs: module receives --output path and file is produced."""
    run_standard_pipeline_test(
        config_path=named_outputs_config,
        config_filename=CONFIG_FILENAME,
        test_name=TEST_NAME,
        tmp_path=tmp_path,
        keep_files=keep_files,
        min_expected_files=1,
    )


@pytest.mark.e2e
def test_named_outputs_snakefile_syntax(
    named_outputs_config, tmp_path, bundled_repos, keep_files
):
    """v0.6 Snakefile uses named-output syntax (data_raw=) not positional."""
    runner = E2ETestRunner(tmp_path, keep_files)
    config_file = runner.setup_test_environment(named_outputs_config, CONFIG_FILENAME)
    runner.execute_cli_command(config_file, ["--continue-on-error"])

    snakefile = runner.out_dir / "Snakefile"
    assert snakefile.exists(), "Snakefile not generated"

    content = snakefile.read_text()
    assert "data_raw=" in content, (
        "Expected named-output syntax 'data_raw=' in Snakefile output block; "
        f"actual output block:\n{content}"
    )
    # v0.6 shell contract: single output → --output path (no id= prefix)
    assert (
        "--output {output.data_raw}" in content
    ), "Expected '--output {output.data_raw}' in Snakefile shell block"
