"""Lean-MVP check for gather expansion (design 010).

Exercises `_expand_gather_stage` with lightweight stand-ins so the test needs
neither module resolution nor a real benchmark. Asserts the fan-in shape:
members grouped by the ancestor module of a `group_by` stage, chain cut,
outputs registered downstream.
"""

from types import SimpleNamespace

from omnibenchmark.cli.run import _expand_gather_stage
from omnibenchmark.model.benchmark import GatherSpec, Stage


def _fake_benchmark():
    model = SimpleNamespace(
        get_name=lambda: "bench",
        get_version=lambda: "1.0",
        get_author=lambda: "me",
    )
    return SimpleNamespace(model=model)


def _member(node_id, parent_id, stage_id, module_id):
    return SimpleNamespace(
        id=node_id, parent_id=parent_id, stage_id=stage_id, module_id=module_id
    )


def test_group_by_stage_partitions_members():
    # Two datasets (d1, d2), each with two clustering methods → 4 producers.
    # group_by dataset ⇒ 2 gather nodes, one per dataset.
    nodes_by_id = {
        "d1.default": _member("d1.default", None, "data", "d1"),
        "d2.default": _member("d2.default", None, "data", "d2"),
        "d1.default-clu-ma.default": _member(
            "d1.default-clu-ma.default", "d1.default", "clu", "ma"
        ),
        "d1.default-clu-mb.default": _member(
            "d1.default-clu-mb.default", "d1.default", "clu", "mb"
        ),
        "d2.default-clu-ma.default": _member(
            "d2.default-clu-ma.default", "d2.default", "clu", "ma"
        ),
        "d2.default-clu-mb.default": _member(
            "d2.default-clu-mb.default", "d2.default", "clu", "mb"
        ),
    }
    output_to_nodes = {
        "clustering": [
            ("d1.default-clu-ma.default", "d1/clu/ma/a.tsv"),
            ("d1.default-clu-mb.default", "d1/clu/mb/b.tsv"),
            ("d2.default-clu-ma.default", "d2/clu/ma/a.tsv"),
            ("d2.default-clu-mb.default", "d2/clu/mb/b.tsv"),
        ]
    }
    stage = SimpleNamespace(
        id="metrics",
        prefix="agg",
        gather=[SimpleNamespace(from_="clustering", group_by="data")],
        modules=[
            SimpleNamespace(
                id="summ", name="summ", parameters=None, provides=None, resources=None
            )
        ],
        # Output template references the group label {data} (the group_by stage).
        outputs=[SimpleNamespace(id="metrics.summary", path="{data}_summary.tsv")],
        resources=None,
    )
    cache = {("metrics", "summ"): object()}

    nodes = _expand_gather_stage(
        stage=stage,
        benchmark=_fake_benchmark(),
        resolved_modules_cache=cache,
        output_to_nodes=output_to_nodes,
        nodes_by_id=nodes_by_id,
    )

    assert len(nodes) == 2
    by_group = {n.id: n for n in nodes}
    assert set(by_group) == {"metrics-summ-d1.default", "metrics-summ-d2.default"}

    d1 = by_group["metrics-summ-d1.default"]
    assert d1.is_gather is True
    assert d1.parent_id is None
    # Only d1's two methods land in the d1 group.
    assert d1.gathered_from == [
        "d1.default-clu-ma.default",
        "d1.default-clu-mb.default",
    ]
    assert set(d1.inputs.values()) == {"d1/clu/ma/a.tsv", "d1/clu/mb/b.tsv"}
    assert set(d1.input_name_mapping.values()) == {"clustering"}
    # Group value bound to the template + baked into the prefix path.
    assert d1.outputs == ["agg/d1/metrics/summ/.default/d1_summary.tsv"]
    # Registered downstream so a later stage can consume it.
    assert ("metrics-summ-d1.default", d1.outputs[0]) in output_to_nodes[
        "metrics.summary"
    ]


def test_gather_collects_across_multiple_stages():
    # Two DIFFERENT stages (method_a, method_b) both produce output id
    # `clustering`, each descending from the same `data` datasets. A single
    # `from: clustering` must collect producers from BOTH stages, grouped by
    # dataset. This is the shared-output-id contract (design 010 §3.1).
    nodes_by_id = {
        "d1.default": _member("d1.default", None, "data", "d1"),
        "d2.default": _member("d2.default", None, "data", "d2"),
        # stage method_a
        "d1.default-ma.default": _member(
            "d1.default-ma.default", "d1.default", "method_a", "ma"
        ),
        "d2.default-ma.default": _member(
            "d2.default-ma.default", "d2.default", "method_a", "ma"
        ),
        # stage method_b (a different stage, same output id)
        "d1.default-mb.default": _member(
            "d1.default-mb.default", "d1.default", "method_b", "mb"
        ),
        "d2.default-mb.default": _member(
            "d2.default-mb.default", "d2.default", "method_b", "mb"
        ),
    }
    output_to_nodes = {
        "clustering": [
            ("d1.default-ma.default", "d1/ma/a.tsv"),  # from method_a
            ("d2.default-ma.default", "d2/ma/a.tsv"),
            ("d1.default-mb.default", "d1/mb/b.tsv"),  # from method_b
            ("d2.default-mb.default", "d2/mb/b.tsv"),
        ]
    }
    stage = SimpleNamespace(
        id="metrics",
        prefix="agg",
        gather=[SimpleNamespace(from_="clustering", group_by="data")],
        modules=[
            SimpleNamespace(
                id="summ", name="summ", parameters=None, provides=None, resources=None
            )
        ],
        outputs=[SimpleNamespace(id="metrics.summary", path="{data}.tsv")],
        resources=None,
    )

    nodes = _expand_gather_stage(
        stage=stage,
        benchmark=_fake_benchmark(),
        resolved_modules_cache={("metrics", "summ"): object()},
        output_to_nodes=output_to_nodes,
        nodes_by_id=nodes_by_id,
    )

    assert len(nodes) == 2
    d1 = next(n for n in nodes if n.id == "metrics-summ-d1.default")
    # d1's group pulls one member from method_a AND one from method_b.
    assert d1.gathered_from == ["d1.default-ma.default", "d1.default-mb.default"]
    assert set(d1.inputs.values()) == {"d1/ma/a.tsv", "d1/mb/b.tsv"}


def test_zero_producer_from_is_plan_time_error():
    stage = SimpleNamespace(
        id="metrics",
        prefix="agg",
        gather=[SimpleNamespace(from_="nonexistent", group_by="data")],
        modules=[],
        outputs=[],
        resources=None,
    )
    try:
        _expand_gather_stage(
            stage=stage,
            benchmark=_fake_benchmark(),
            resolved_modules_cache={},
            output_to_nodes={},
            nodes_by_id={},
        )
    except ValueError as e:
        assert "nonexistent" in str(e)
    else:
        raise AssertionError("expected ValueError for zero-producer gather.from")


def test_gather_requires_prefix():
    """Model validation: gather without prefix is rejected at parse time."""
    try:
        Stage(
            id="metrics",
            modules=[],
            outputs=[],
            gather=[GatherSpec(from_="clustering", group_by="data")],
        )
    except Exception as e:
        assert "prefix" in str(e)
    else:
        raise AssertionError("expected validation error for gather without prefix")


_GATHER_YAML = """
id: t
description: t
version: '1.0'
benchmarker: me
api_version: {api}
software_backend: host
software_environments:
  env: {{description: e, easyconfig: e.eb}}
stages:
  - id: data
    outputs: [{{id: clustering, path: c.tsv}}]
    modules: [{{id: d1, repository: {{url: 'http://x', commit: abc}}, software_environment: env}}]
  - id: metrics
    prefix: agg
    gather: [{{from: clustering, group_by: {group_by}}}]
    modules: [{{id: s, repository: {{url: 'http://x', commit: abc}}, software_environment: env}}]
    outputs: [{{id: metrics.summary, path: summary.tsv}}]
"""


def test_gather_benchmark_parses_and_from_alias():
    from omnibenchmark.model.benchmark import Benchmark

    bench = Benchmark.from_yaml(_GATHER_YAML.format(api="0.6.0", group_by="data"))
    gather_stage = bench.stages[1]
    assert gather_stage.gather[0].from_ == "clustering"
    assert gather_stage.gather[0].group_by == "data"


def test_gather_gated_on_api_0_6():
    from omnibenchmark.model.benchmark import Benchmark

    try:
        Benchmark.from_yaml(_GATHER_YAML.format(api="0.5.0", group_by="data"))
    except Exception as e:
        assert "gather" in str(e) and "0.6.0" in str(e)
    else:
        raise AssertionError("expected gather to be gated on api >= 0.6.0")


def test_shared_output_id_across_stages_is_legal():
    """Two stages may declare the same output id (the gather contract, §3.1)."""
    from omnibenchmark.model.benchmark import Benchmark

    y = """
id: t
description: t
version: '1.0'
benchmarker: me
api_version: 0.6.0
software_backend: host
software_environments:
  env: {description: e, easyconfig: e.eb}
stages:
  - id: method_a
    outputs: [{id: clustering, path: a.tsv}]
    modules: [{id: ma, repository: {url: 'http://x', commit: abc}, software_environment: env}]
  - id: method_b
    outputs: [{id: clustering, path: b.tsv}]
    modules: [{id: mb, repository: {url: 'http://x', commit: abc}, software_environment: env}]
"""
    bench = Benchmark.from_yaml(y)
    assert [s.id for s in bench.stages] == ["method_a", "method_b"]


def test_gather_entries_must_share_one_group_by_axis():
    """Multiple gather entries with differing group_by are rejected; a shared
    axis is accepted (Stage.validate_gather — makes the [0] axis sound)."""
    # Same axis across entries: OK.
    Stage(
        id="metrics",
        prefix="agg",
        modules=[],
        outputs=[],
        gather=[
            GatherSpec(from_="clustering", group_by="data"),
            GatherSpec(from_="embedding", group_by="data"),
        ],
    )
    # Differing axes: rejected.
    try:
        Stage(
            id="metrics",
            prefix="agg",
            modules=[],
            outputs=[],
            gather=[
                GatherSpec(from_="clustering", group_by="data"),
                GatherSpec(from_="embedding", group_by="method"),
            ],
        )
    except Exception as e:
        assert "group_by" in str(e)
    else:
        raise AssertionError("expected rejection of differing group_by axes")


def test_get_stages_by_output_returns_all_producers_in_order():
    """One-to-many output→stage lookup: every producer, declaration order."""
    from omnibenchmark.model.benchmark import Benchmark

    y = """
id: t
description: t
version: '1.0'
benchmarker: me
api_version: 0.6.0
software_backend: host
software_environments:
  env: {description: e, easyconfig: e.eb}
stages:
  - id: method_a
    outputs: [{id: clustering, path: a.tsv}]
    modules: [{id: ma, repository: {url: 'http://x', commit: abc}, software_environment: env}]
  - id: method_b
    outputs: [{id: clustering, path: b.tsv}]
    modules: [{id: mb, repository: {url: 'http://x', commit: abc}, software_environment: env}]
  - id: solo
    outputs: [{id: embedding, path: e.tsv}]
    modules: [{id: ms, repository: {url: 'http://x', commit: abc}, software_environment: env}]
"""
    bench = Benchmark.from_yaml(y)
    shared = bench.get_stages_by_output("clustering")
    assert [s.id for s in shared] == ["method_a", "method_b"]  # both, in order
    assert [s.id for s in bench.get_stages_by_output("embedding")] == ["solo"]
    assert bench.get_stages_by_output("nope") == []


def test_group_by_must_name_a_stage():
    from omnibenchmark.model.benchmark import Benchmark

    try:
        Benchmark.from_yaml(_GATHER_YAML.format(api="0.6.0", group_by="nosuch"))
    except Exception as e:
        assert "not a known stage" in str(e)
    else:
        raise AssertionError("expected group_by to be validated against stage ids")
