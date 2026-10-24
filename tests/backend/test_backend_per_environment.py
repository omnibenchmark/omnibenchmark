"""An environment runs under the backend it declares when the plan's choice is
missing from it; a module without an environment runs on the host (design 014)."""

from types import SimpleNamespace

import pytest

from omnibenchmark.backend.resolver import _backend_for
from omnibenchmark.cli.run import _backend_flags
from omnibenchmark.model.benchmark import SoftwareBackendEnum as B


@pytest.mark.short
def test_backend_for_prefers_the_plan_then_falls_back():
    env = SimpleNamespace(
        conda="envs/x.yml", apptainer=None, docker=None, envmodule=None
    )
    assert _backend_for(env, B.conda) is B.conda
    assert _backend_for(env, B.apptainer) is B.conda
    none = SimpleNamespace(conda=None, apptainer=None, docker=None, envmodule=None)
    assert _backend_for(none, B.apptainer) is None


@pytest.mark.short
def test_flags_cover_every_backend_in_use():
    conda_node = SimpleNamespace(
        module=SimpleNamespace(
            resolved_environment=SimpleNamespace(backend_type=B.conda)
        )
    )
    host_node = SimpleNamespace(module=SimpleNamespace(resolved_environment=None))
    assert _backend_flags(B.apptainer, [conda_node, host_node]) == [
        "--use-singularity",
        "--use-conda",
    ]
    assert _backend_flags(B.host, [host_node]) == []
