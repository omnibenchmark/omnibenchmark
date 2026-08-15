"""Parameter values resolved from the lineage — `{label.params.key}`.

Drives `_expand_scatter_stage` with lightweight stand-ins (same approach as
`test_gather_expand.py`) so the test needs neither module resolution nor a real
benchmark. Asserts the whole point of the feature: one downstream module
declaration expands to per-lineage parameter *values*, with the param hash
following along.
"""

from types import SimpleNamespace

import pytest

from omnibenchmark.cli.run import _expand_scatter_stage
from omnibenchmark.model.benchmark import Parameter


def _output(output_id, path):
    return SimpleNamespace(id=output_id, path=path)


def _module(module_id, params=None):
    return SimpleNamespace(
        id=module_id,
        name=module_id,
        parameters=[Parameter(id="p", params=params)] if params else None,
        requires=None,
        requires_capabilities=None,
        resources=None,
        exclude=None,
    )


def _stage(stage_id, modules, outputs, inputs=None):
    return SimpleNamespace(
        id=stage_id,
        modules=modules,
        inputs=inputs,
        outputs=outputs,
        resources=None,
        provides=None,
    )


def _expand(stages):
    """Expand `stages` in order, returning nodes grouped by stage id."""
    benchmark = SimpleNamespace(
        model=SimpleNamespace(
            stages=stages,
            get_name=lambda: "bench",
            get_version=lambda: "1.0",
            get_author=lambda: "me",
        )
    )
    resolved_modules_cache = {
        (s.id, m.id): SimpleNamespace(id=m.id) for s in stages for m in s.modules
    }
    resolved_nodes, nodes_by_id, output_to_nodes = [], {}, {}
    dag_errors, previous, by_stage = [], [], {}

    for stage in stages:
        previous = _expand_scatter_stage(
            stage=stage,
            benchmark=benchmark,
            resolved_modules_cache=resolved_modules_cache,
            resolved_nodes=resolved_nodes,
            nodes_by_id=nodes_by_id,
            output_to_nodes=output_to_nodes,
            previous_stage_nodes=previous,
            stages_to_expand=stages,
            path_exclusions={},
            nesting_strategy="nested",
            module_filter=None,
            target_stage=None,
            dag_errors=dag_errors,
            quiet=True,
        )
        by_stage[stage.id] = previous

    return by_stage, dag_errors


def _pipeline(pca_params):
    """Two datasets with different `ideal_components`, one PCA module."""
    return [
        _stage(
            "data",
            [
                _module("D1", {"ideal_components": 10}),
                _module("D2", {"ideal_components": 25}),
            ],
            [_output("data.raw", "{name}_data.json")],
        ),
        _stage(
            "pca",
            [_module("PCA", pca_params)],
            [_output("pca.embedding", "{name}_embedding.json")],
            inputs=[SimpleNamespace(entries=["data.raw"])],
        ),
    ]


@pytest.mark.short
def test_param_resolves_per_lineage():
    by_stage, dag_errors = _expand(
        _pipeline({"k": "{dataset.params.ideal_components}"})
    )
    assert not dag_errors

    k_by_dataset = {
        node.template_context.provides["dataset"]: node.parameters["k"]
        for node in by_stage["pca"]
    }
    assert k_by_dataset == {"D1": 10, "D2": 25}


@pytest.mark.short
def test_resolved_value_keeps_native_type():
    by_stage, _ = _expand(_pipeline({"k": "{dataset.params.ideal_components}"}))
    assert len(by_stage["pca"]) == 2
    assert all(isinstance(n.parameters["k"], int) for n in by_stage["pca"])


@pytest.mark.short
def test_param_hash_follows_the_resolved_value():
    """The hash is taken after resolution, so the two nodes never share a
    param directory — which under the `flat` strategy would be one directory."""
    by_stage, _ = _expand(_pipeline({"k": "{dataset.params.ideal_components}"}))
    assert len({n.param_id for n in by_stage["pca"]}) == 2
    assert all(n.param_id in n.id for n in by_stage["pca"])


@pytest.mark.short
def test_literal_params_still_share_one_hash():
    """Regression guard: without a reference, nothing changes."""
    by_stage, _ = _expand(_pipeline({"k": 3}))
    assert len({n.param_id for n in by_stage["pca"]}) == 1
    assert all(n.parameters["k"] == 3 for n in by_stage["pca"])


@pytest.mark.short
def test_missing_key_is_a_dag_error():
    _, dag_errors = _expand(_pipeline({"k": "{dataset.params.nonexistent}"}))
    assert dag_errors
    assert "declares no parameter 'nonexistent'" in dag_errors[0][2]


@pytest.mark.short
def test_unknown_label_is_a_dag_error():
    _, dag_errors = _expand(_pipeline({"k": "{treatment.params.dose}"}))
    assert dag_errors
    assert "Unknown lineage label 'treatment'" in dag_errors[0][2]


@pytest.mark.short
def test_reference_chains_through_an_intermediate_stage():
    """`dataset` is inherited, so a ref works at any depth, not just depth 1."""
    stages = _pipeline({"k": "{dataset.params.ideal_components}"})
    stages.append(
        _stage(
            "metrics",
            [_module("M1", {"n": "{dataset.params.ideal_components}"})],
            [_output("metrics.score", "{name}_score.json")],
            inputs=[SimpleNamespace(entries=["pca.embedding"])],
        )
    )
    by_stage, dag_errors = _expand(stages)
    assert not dag_errors

    n_by_dataset = {
        node.template_context.provides["dataset"]: node.parameters["n"]
        for node in by_stage["metrics"]
    }
    assert n_by_dataset == {"D1": 10, "D2": 25}
