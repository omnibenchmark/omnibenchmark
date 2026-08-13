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
    """A stage consumes an output id from *any* stage that produces it.

    Topology: ``data`` feeds both ``pca`` and ``cntfct``, which declare the same
    output id ``embedding``; ``metrics`` declares ``inputs: [embedding]``.

    Before the fix the resolver collapsed the producers onto the single deepest
    stage, so ``metrics`` expanded under one branch only and the other computed
    an embedding nothing consumed — a dead branch, with no warning. Both
    branches must now carry a ``metrics`` node.

    Parametrised over api versions because this must NOT require an api bump:
    benchmarks stay pinned at 0.4.0 for unrelated reasons (the `--name`
    /`{dataset}` path-position workaround their modules depend on), and asking
    them to migrate that in order to get multi-producer inputs would couple two
    unrelated changes.

    Asserted on the generated Snakefile (``--dry``) so the check does not
    depend on module execution.
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

    assert via_pca, "metrics should expand on the pca branch"
    assert via_cntfct, (
        "metrics should expand on the cntfct branch too — both stages declare "
        "the `embedding` contract, so neither may be silently dropped"
    )


@pytest.mark.e2e
def test_producers_are_not_joined_into_one_node(
    shared_output_contract_config, tmp_path, bundled_repos, keep_files
):
    """Alternatives fan *out*; they are not a fan-in.

    Two producers of the same id are interchangeable, so each yields its own
    downstream node. A rule taking `embedding` from both branches at once would
    mean the resolver mistook a coproduct for a join.
    """
    runner = E2ETestRunner(tmp_path, keep_files)
    config_file = runner.setup_test_environment(
        shared_output_contract_config, "11_shared_output_contract.yaml"
    )

    runner.execute_cli_command(config_file, ["--dry"])

    snakefile = (runner.out_dir / "Snakefile").read_text()

    for block in re.split(r"^rule ", snakefile, flags=re.M)[1:]:
        # `rule all` aggregates every target, so it legitimately names both
        # branches; it is not a job rule.
        if block.startswith("all:") or "/metrics/M1" not in block:
            continue
        inputs = block.split("output:")[0]
        assert not (
            "/pca/P1" in inputs and "/cntfct/C1" in inputs
        ), f"a metrics rule pulled inputs from both branches at once:\n{block[:400]}"
