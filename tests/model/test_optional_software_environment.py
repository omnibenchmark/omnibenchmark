"""Module-level `software_environment` is optional when unambiguous (#127).

- exactly one environment declared: an omitted reference resolves to it;
- `software_backend: host`: environments are irrelevant, so an omitted
  reference (or a missing `software_environments` block) is fine;
- otherwise an omitted reference is an error naming the module.

An explicit reference is always checked, so typos still fail.
"""

from typing import Optional

import pytest
import yaml
from pydantic import ValidationError

from omnibenchmark.model import Benchmark
from omnibenchmark.model.validation import (
    ValidationError as OmnibenchmarkValidationError,
)

ERRORS = (ValidationError, OmnibenchmarkValidationError)


def _plan(
    backend: str = "conda",
    envs: Optional[list] = ("py",),
    module_env: Optional[str] = None,
    collector_env: Optional[str] = None,
    with_collector: bool = False,
) -> Benchmark:
    module = {
        "id": "D1",
        "repository": {"url": "https://example.org/d.git", "commit": "abc1234"},
    }
    if module_env:
        module["software_environment"] = module_env
    doc = {
        "id": "b",
        "benchmarker": "me",
        "version": "1.0",
        "api_version": "0.7.0",
        "software_backend": backend,
        "stages": [
            {
                "id": "data",
                "modules": [module],
                "outputs": [{"id": "data_counts", "path": "{dataset}.txt"}],
            }
        ],
    }
    if envs is not None:
        doc["software_environments"] = {
            e: {"description": e, "conda": f"envs/{e}.yaml"} for e in envs
        }
    if with_collector:
        collector = {
            "id": "c",
            "name": "c",
            "repository": {"url": "https://example.org/c.git", "commit": "abc1234"},
            "inputs": ["data_counts"],
            "outputs": [{"id": "plot", "path": "plot.html"}],
        }
        if collector_env:
            collector["software_environment"] = collector_env
        doc["metric_collectors"] = [collector]
    return Benchmark.from_yaml(yaml.safe_dump(doc))


def _module(b: Benchmark):
    return b.stages[0].modules[0]


@pytest.mark.short
class TestSoleEnvironment:
    def test_omitted_module_env_resolves_to_sole_env(self):
        b = _plan(envs=["py"])
        assert _module(b).software_environment == "py"
        assert b.get_module_environment("D1") == "py"

    def test_omitted_collector_env_resolves_to_sole_env(self):
        b = _plan(envs=["py"], with_collector=True)
        assert b.metric_collectors[0].software_environment == "py"

    def test_explicit_ref_kept(self):
        assert _module(_plan(envs=["py"], module_env="py")).software_environment == "py"

    def test_explicit_undeclared_ref_still_fails(self):
        with pytest.raises(ERRORS, match="'pyy' for module 'D1' is not declared"):
            _plan(envs=["py"], module_env="pyy")


@pytest.mark.short
class TestSeveralEnvironments:
    def test_omitted_module_env_means_host(self):
        # design 014: no environment = the module runs on the host, whatever the backend
        assert _module(_plan(envs=["py", "r"])).software_environment is None

    def test_omitted_collector_env_means_host(self):
        b = _plan(envs=["py", "r"], module_env="py", with_collector=True)
        assert b.metric_collectors[0].software_environment is None

    def test_explicit_refs_valid(self):
        b = _plan(
            envs=["py", "r"], module_env="r", with_collector=True, collector_env="py"
        )
        assert _module(b).software_environment == "r"
        assert b.metric_collectors[0].software_environment == "py"

    def test_omitted_with_no_environments_block_means_host(self):
        assert _module(_plan(envs=None)).software_environment is None


@pytest.mark.short
class TestHostBackend:
    def test_no_environments_block_needed(self):
        b = _plan(backend="host", envs=None, with_collector=True)
        assert b.software_environments == []
        assert _module(b).software_environment is None

    def test_omitted_ref_with_several_envs(self):
        b = _plan(backend="host", envs=["py", "r"])
        assert _module(b).software_environment is None

    def test_sole_env_still_auto_fills(self):
        # Keeps the plan portable to a non-host backend.
        assert _module(_plan(backend="host", envs=["py"])).software_environment == "py"

    def test_explicit_undeclared_ref_still_fails(self):
        with pytest.raises(ERRORS, match="'pyy' for module 'D1' is not declared"):
            _plan(backend="host", envs=["py", "r"], module_env="pyy")


@pytest.mark.short
def test_host_resolver_skips_env_lookup(tmp_path, caplog):
    """An omitted reference under host resolves to no environment, silently."""
    from omnibenchmark.backend.resolver import ModuleResolver
    from omnibenchmark.model import SoftwareBackendEnum

    resolver = ModuleResolver(
        output_dir=tmp_path, software_backend=SoftwareBackendEnum.host
    )
    with caplog.at_level("WARNING"):
        assert resolver._resolve_environment("", "D1") is None
    assert caplog.records == []
