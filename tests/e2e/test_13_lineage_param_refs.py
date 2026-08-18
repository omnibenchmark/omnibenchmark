import re
import pytest
from pathlib import Path

from tests.e2e.common import E2ETestRunner


@pytest.fixture
def lineage_param_refs_config():
    """Path to the lineage-resolved parameter config."""
    return Path(__file__).parent / "configs" / "13_lineage_param_refs.yaml"


@pytest.mark.e2e
def test_param_value_resolves_from_the_dataset(
    lineage_param_refs_config, tmp_path, bundled_repos, keep_files
):
    """One `k:` declaration, two concrete values — one per dataset lineage.

    M1 declares ``k: "{dataset.params.ideal_components}"`` once; D1 declares 10
    and D2 declares 25, so the generated command must differ per branch. The
    resolved value is also what the parameter hash is taken from, so the two
    nodes must not share an output directory.
    """
    runner = E2ETestRunner(tmp_path, keep_files)
    config_file = runner.setup_test_environment(
        lineage_param_refs_config, "13_lineage_param_refs.yaml"
    )

    runner.execute_cli_command(config_file, ["--dry"])

    snakefile = (runner.out_dir / "Snakefile").read_text()

    # The template must never survive into the generated workflow.
    assert "{dataset.params" not in snakefile

    # cli_args is emitted once per rule; pair each with the rule's dataset.
    rules = re.findall(
        r'rule ([^\s:]+):.*?cli_args="([^"]*)"', snakefile, flags=re.DOTALL
    )
    k_by_dataset = {
        "D1" if "D1" in name else "D2": re.search(r"--k (\S+)", args).group(1)
        for name, args in rules
        if "--k " in args
    }
    assert k_by_dataset == {"D1": "10", "D2": "25"}

    # Resolved before the param hash: distinct values, distinct directories.
    methods_dirs = set(re.findall(r"\bdata/D\d\S*/methods/M1/(\.\w+)", snakefile))
    assert len(methods_dirs) == 2, methods_dirs
