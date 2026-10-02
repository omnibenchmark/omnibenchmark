"""A stage's inputs resolve from their producers, never from expansion order.

Regression: `expand_scatter_stage` resolved declared inputs only when the stage
expanded just before it had nodes. Topological ties are broken by set (hash)
order, so a stage with no modules (or every module pruned) could land between a
producer and its consumer, and the consumer then expanded as a root: no inputs,
no parent, nothing passed on the command line. Which plans broke depended on
PYTHONHASHSEED.
"""

from types import SimpleNamespace

import pytest

from omnibenchmark.core._expand import expand_scatter_stage


def _module(module_id):
    return SimpleNamespace(
        id=module_id,
        name=module_id,
        parameters=None,
        provides=None,
        requires=None,
        requires_capabilities=None,
        resources=None,
    )


def _stage(stage_id, inputs, output_id):
    return SimpleNamespace(
        id=stage_id,
        inputs=[SimpleNamespace(entries=inputs)] if inputs else None,
        outputs=[SimpleNamespace(id=output_id, path="{name}_" + output_id)],
        modules=[_module(f"{stage_id.lower()}-m")],
        resources=None,
        provides=None,
        gather=None,
    )


def _expand(stage, benchmark, state, previous_stage_nodes):
    return expand_scatter_stage(
        stage=stage,
        benchmark=benchmark,
        resolved_modules_cache={(stage.id, m.id): object() for m in stage.modules},
        previous_stage_nodes=previous_stage_nodes,
        stages_to_expand=benchmark.stages,
        path_exclusions={},
        nesting_strategy="nested",
        module_filter=None,
        target_stage=None,
        dag_errors=state["errors"],
        prune_counts={"requires": 0, "exclude": 0},
        quiet=True,
        resolved_nodes=state["resolved"],
        nodes_by_id=state["by_id"],
        output_to_nodes=state["outputs"],
    )


@pytest.mark.short
def test_consumer_after_empty_stage_still_resolves_its_inputs():
    # data -> pca, with an unrelated module-less stage sorted between them.
    data = _stage("data", None, "counts")
    empty = _stage("empty", ["counts"], "other")
    empty.modules = []
    pca = _stage("pca", ["counts"], "embedding")
    stages = [data, empty, pca]
    model = SimpleNamespace(
        get_name=lambda: "b",
        get_version=lambda: "1.0",
        get_author=lambda: "me",
        get_stage=lambda sid: next((s for s in stages if s.id == sid), None),
    )
    benchmark = SimpleNamespace(model=model, stages=stages)
    state = {"errors": [], "resolved": [], "by_id": {}, "outputs": {}}

    data_nodes = _expand(data, benchmark, state, [])
    empty_nodes = _expand(empty, benchmark, state, data_nodes)
    assert data_nodes and empty_nodes == []

    pca_nodes = _expand(pca, benchmark, state, empty_nodes)

    assert not state["errors"]
    assert len(pca_nodes) == 1
    (node,) = pca_nodes
    assert node.parent_id == data_nodes[0].id
    assert set(node.inputs) == {"counts"}
