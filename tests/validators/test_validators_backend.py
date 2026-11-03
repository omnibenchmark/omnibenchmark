"""Unit tests for omnibenchmark.backend.validators: env directive and cleanup."""

from types import SimpleNamespace

import pytest

from omnibenchmark.backend import validators as bv
from omnibenchmark.model.benchmark import SoftwareBackendEnum
from omnibenchmark.model.resolved import ResolvedEnvironment


def _model(backend):
    return SimpleNamespace(
        get_software_backend=lambda: backend,
        get_software_environments=lambda: {},
    )


def _resolve_to(monkeypatch, resolved):
    monkeypatch.setattr(
        "omnibenchmark.backend.resolver.ModuleResolver.resolve_environment",
        lambda self, env_id, module_id: resolved,
    )


@pytest.mark.short
def test_env_directive_keeps_image_refs(tmp_path, monkeypatch):
    env = ResolvedEnvironment(SoftwareBackendEnum.docker, "ubuntu:22.04")
    _resolve_to(monkeypatch, env)
    assert bv._resolve_env_directive(
        _model(SoftwareBackendEnum.docker), tmp_path, tmp_path, "e"
    ) == ('container: "ubuntu:22.04"', "--use-singularity")


@pytest.mark.short
def test_env_directive_makes_local_env_files_absolute(tmp_path, monkeypatch):
    (tmp_path / ".envs").mkdir()
    (tmp_path / ".envs" / "e.yaml").write_text("")
    env = ResolvedEnvironment(SoftwareBackendEnum.conda, ".envs/e.yaml")
    _resolve_to(monkeypatch, env)
    directive, flag = bv._resolve_env_directive(
        _model(SoftwareBackendEnum.conda), tmp_path, tmp_path, "e"
    )
    assert directive == f'conda: "{tmp_path / ".envs" / "e.yaml"}"'
    assert flag == "--use-conda"


@pytest.mark.short
def test_env_directive_refuses_host_fallback(tmp_path, monkeypatch):
    _resolve_to(monkeypatch, None)
    with pytest.raises(ValueError, match="could not resolve"):
        bv._resolve_env_directive(
            _model(SoftwareBackendEnum.conda), tmp_path, tmp_path, "e"
        )


@pytest.mark.short
def test_drop_orphan_results(tmp_path):
    names = ["a.txt.0.json", "a.txt.0.log", "a.txt.1.json", "a.txt.1.log"]
    for name in names + ["snakemake.log", "Snakefile"]:
        (tmp_path / name).write_text("")
    bv._drop_orphan_results(tmp_path, {"a.txt.0": []})
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        "Snakefile",
        "a.txt.0.json",
        "a.txt.0.log",
        "snakemake.log",
    ]


@pytest.mark.short
def test_env_directive_host_backend_needs_none(tmp_path):
    assert bv._resolve_env_directive(
        _model(SoftwareBackendEnum.host), tmp_path, tmp_path, "e"
    ) == ("", None)
