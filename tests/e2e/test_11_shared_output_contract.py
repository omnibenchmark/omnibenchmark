import re
import pytest
from pathlib import Path

from tests.e2e.common import E2ETestRunner


@pytest.fixture
def shared_output_contract_config():
    """Path to the shared-output-id (alternative producers) config."""
    return Path(__file__).parent / "configs" / "11_shared_output_contract.yaml"


@pytest.mark.e2e
@pytest.mark.parametrize("api_version", ["0.4.0", "0.5.0"])
def test_consumer_binds_to_every_producer_of_the_contract(
    shared_output_contract_config, api_version, tmp_path, bundled_repos, keep_files
):
    """A stage consumes an output id from *any* stage producing it.

    ``pca`` and ``cntfct`` both declare ``embedding``; before the fix the
    resolver kept only the deepest, leaving the other a dead branch that
    computed an embedding nothing consumed — silently.

    Parametrised over api versions on purpose: this must not require an api
    bump, since benchmarks stay pinned at 0.4.0 by the ``--name``/``{dataset}``
    workaround their modules depend on.
    """
    src = shared_output_contract_config.read_text().replace(
        'api_version: "0.5.0"', f'api_version: "{api_version}"'
    )
    pinned = tmp_path / f"shared_contract_{api_version.replace('.', '_')}.yaml"
    pinned.write_text(src)

    runner = E2ETestRunner(tmp_path, keep_files)
    config_file = runner.setup_test_environment(
        pinned, "11_shared_output_contract.yaml"
    )

    runner.execute_cli_command(config_file, ["--dry"])

    snakefile = (runner.out_dir / "Snakefile").read_text()

    via_pca = re.findall(r"\bdata/D1\S*/pca/P1\S*/metrics/M1", snakefile)
    via_cntfct = re.findall(r"\bdata/D1\S*/cntfct/C1\S*/metrics/M1", snakefile)

    # Alternatives fan out: one metrics node per producer. Had the resolver
    # mistaken them for a join, a single node would sit under one branch and
    # the other regex would find nothing.
    assert via_pca, "metrics should expand on the pca branch"
    assert via_cntfct, (
        "metrics should expand on the cntfct branch too — both stages declare "
        "the `embedding` contract, so neither may be silently dropped"
    )
