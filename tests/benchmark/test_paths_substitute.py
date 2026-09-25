"""The legacy path layer fills the same template variables as the run path."""

from pathlib import Path
from types import SimpleNamespace

import pytest

from omnibenchmark.core._paths import normalize_dataset_path, substitute_node_vars
from omnibenchmark.model import IOFile, expand_output_path


def _node(parameters=None, name=None):
    return SimpleNamespace(
        stage_id="methods",
        module_id="M1",
        param_id=".abc123",
        module=SimpleNamespace(name=name),
        parameters=parameters,
    )


@pytest.mark.short
def test_substitutes_run_path_variables():
    tmpl = "{input}/{stage}/{module}/{params}/{name}-{module.id}-{module.stage}.txt"
    assert (
        substitute_node_vars(tmpl, _node(), Path("out"))
        == "out/methods/M1/.abc123/M1-M1-methods.txt"
    )


@pytest.mark.short
def test_module_name_falls_back_to_id():
    assert substitute_node_vars("{module.name}", _node(), Path("out")) == "M1"
    assert (
        substitute_node_vars("{module.name}", _node(name="Fast"), Path("o")) == "Fast"
    )


@pytest.mark.short
def test_params_key_and_unknown_left_alone():
    node = _node(parameters={"k": 10})
    assert (
        substitute_node_vars("{params.k}_{params.x}_{dataset}", node, Path("out"))
        == "10_{params.x}_{dataset}"
    )


@pytest.mark.short
def test_dataset_resolved_without_format():
    # str.format would raise on the leftover braces
    path = "out/data/D1/.abc/{dataset}_{params.x}.txt"
    assert (
        normalize_dataset_path(path, Path("out"))
        == "out/data/D1/.abc/D1_{params.x}.txt"
    )


@pytest.mark.short
def test_expand_output_path_keeps_filename():
    assert expand_output_path(IOFile(id="o", path="{module.id}.txt")) == (
        "{input}/{stage}/{module}/{params}/{module.id}.txt"
    )
