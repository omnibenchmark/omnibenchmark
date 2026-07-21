import re
import pytest
from pathlib import Path

from tests.e2e.common import E2ETestRunner


@pytest.fixture
def capability_config():
    """Path to the capability-gating config (M_gpu requires the gpu capability)."""
    return Path(__file__).parent / "configs" / "10_capability_gating.yaml"


@pytest.mark.e2e
def test_gpu_module_pruned_without_capability(
    capability_config, tmp_path, bundled_repos, keep_files
):
    """Without --with-capability gpu, the M_gpu module must not appear in the
    generated Snakefile; the plain M_cpu module still does."""
    runner = E2ETestRunner(tmp_path, keep_files)
    config_file = runner.setup_test_environment(
        capability_config, "10_capability_gating.yaml"
    )

    runner.execute_cli_command(config_file, ["--dry"])
    snakefile = (runner.out_dir / "Snakefile").read_text()

    assert re.search(r"methods/M_cpu", snakefile), "M_cpu should always run"
    assert not re.search(
        r"methods/M_gpu", snakefile
    ), "M_gpu requires gpu and must be pruned when the host does not declare it"


@pytest.mark.e2e
def test_gpu_module_included_with_capability(
    capability_config, tmp_path, bundled_repos, keep_files
):
    """With --with-capability gpu, both modules appear in the Snakefile."""
    runner = E2ETestRunner(tmp_path, keep_files)
    config_file = runner.setup_test_environment(
        capability_config, "10_capability_gating.yaml"
    )

    runner.execute_cli_command(config_file, ["--dry", "--with-capability", "gpu"])
    snakefile = (runner.out_dir / "Snakefile").read_text()

    assert re.search(r"methods/M_cpu", snakefile), "M_cpu should always run"
    assert re.search(
        r"methods/M_gpu", snakefile
    ), "M_gpu should be included once the gpu capability is declared"
