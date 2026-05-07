"""E2E tests for `--capability` capability-gated module pruning."""

from pathlib import Path

import pytest


@pytest.fixture
def capabilities_config():
    """Path to the capabilities-gated benchmark config."""
    return Path(__file__).parent / "configs" / "09_capabilities.yaml"


@pytest.mark.e2e
def test_capability_gated_module_pruned_by_default(
    capabilities_config, tmp_path, bundled_repos, keep_files
):
    """Without `--capability`, the GPU-gated module is pruned; only D1 runs."""
    from tests.e2e.common import (
        E2ETestRunner,
        filter_files_excluding_symlinked_dirs,
    )

    runner = E2ETestRunner(tmp_path, keep_files)
    config_in_tmp = runner.setup_test_environment(
        capabilities_config, "09_capabilities.yaml"
    )

    runner.execute_cli_command(config_in_tmp, ["--continue-on-error"])

    data_files = filter_files_excluding_symlinked_dirs(runner.out_dir, "*_data.json")
    assert (
        len(data_files) == 1
    ), f"expected 1 data output, got {len(data_files)}: {data_files}"
    assert any("D1" in str(p) for p in data_files)
    assert not any(
        "D2_gpu" in str(p) for p in data_files
    ), f"D2_gpu output should be pruned: {data_files}"


@pytest.mark.e2e
def test_capability_gated_module_runs_when_capability_provided(
    capabilities_config, tmp_path, bundled_repos, keep_files
):
    """With `--capability gpu`, both D1 and D2_gpu run."""
    from tests.e2e.common import (
        E2ETestRunner,
        filter_files_excluding_symlinked_dirs,
    )

    runner = E2ETestRunner(tmp_path, keep_files)
    config_in_tmp = runner.setup_test_environment(
        capabilities_config, "09_capabilities.yaml"
    )

    runner.execute_cli_command(
        config_in_tmp,
        ["--continue-on-error", "--capability", "gpu"],
    )

    data_files = filter_files_excluding_symlinked_dirs(runner.out_dir, "*_data.json")
    assert (
        len(data_files) == 2
    ), f"expected 2 data outputs, got {len(data_files)}: {data_files}"
    assert any("D1" in str(p) for p in data_files)
    assert any("D2_gpu" in str(p) for p in data_files)


@pytest.mark.e2e
def test_unrelated_capability_does_not_unprune(
    capabilities_config, tmp_path, bundled_repos, keep_files
):
    """`--capability foo` does not satisfy the gpu gate; D2_gpu still pruned."""
    from tests.e2e.common import (
        E2ETestRunner,
        filter_files_excluding_symlinked_dirs,
    )

    runner = E2ETestRunner(tmp_path, keep_files)
    config_in_tmp = runner.setup_test_environment(
        capabilities_config, "09_capabilities.yaml"
    )

    runner.execute_cli_command(
        config_in_tmp,
        ["--continue-on-error", "--capability", "foo"],
    )

    data_files = filter_files_excluding_symlinked_dirs(runner.out_dir, "*_data.json")
    assert (
        len(data_files) == 1
    ), f"expected 1 data output, got {len(data_files)}: {data_files}"
    assert not any("D2_gpu" in str(p) for p in data_files)
