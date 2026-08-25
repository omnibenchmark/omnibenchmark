"""Lean-MVP check for gather expansion (design 010).

Exercises `expand_gather_stage` with lightweight stand-ins so the test needs
neither module resolution nor a real benchmark. Asserts the fan-in shape:
members grouped by the ancestor module of a `group_by` stage, chain cut,
outputs registered downstream.
"""

from types import SimpleNamespace

import pytest

from omnibenchmark.backend.snakemake import _human_link_name
from omnibenchmark.core._expand import expand_gather_stage
from omnibenchmark.core._lineage import expansion_segment, select_input_bundles
from omnibenchmark.core._paths import make_human_name as _make_human_name
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

    nodes = expand_gather_stage(
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

    nodes = expand_gather_stage(
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
        expand_gather_stage(
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

    bench = Benchmark.from_yaml(_GATHER_YAML.format(api="0.7.0", group_by="data"))
    gather_stage = bench.stages[1]
    assert gather_stage.gather[0].from_ == "clustering"
    assert gather_stage.gather[0].group_by == "data"


def test_gather_gated_on_api_0_7():
    from omnibenchmark.model.benchmark import Benchmark

    try:
        Benchmark.from_yaml(_GATHER_YAML.format(api="0.5.0", group_by="data"))
    except Exception as e:
        assert "gather" in str(e) and "0.7.0" in str(e)
    else:
        raise AssertionError("expected gather to be gated on api >= 0.7.0")


def test_shared_output_id_across_stages_is_legal():
    """Two stages may declare the same output id (the gather contract, §3.1)."""
    from omnibenchmark.model.benchmark import Benchmark

    y = """
id: t
description: t
version: '1.0'
benchmarker: me
api_version: 0.7.0
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
api_version: 0.7.0
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
        Benchmark.from_yaml(_GATHER_YAML.format(api="0.7.0", group_by="nosuch"))
    except Exception as e:
        assert "not a known stage" in str(e)
    else:
        raise AssertionError("expected group_by to be validated against stage ids")


def test_select_input_bundles_pairs_diamond_branches_by_root():
    """The fan-in join (design 010 §3.9, #289): a stage declaring inputs from
    two divergent branches gets one bundle per (anchor, partner) pair, and
    partners are only drawn from the SAME lineage root — no cross-dataset
    joins. Linear anchors stay 1-tuples (fast path, no behaviour change)."""
    nodes = [
        _member("a1.default", None, "A", "a1"),
        _member("a2.default", None, "A", "a2"),
        _member("a1.default-C1-m.default", "a1.default", "C1", "m"),
        _member("a2.default-C1-m.default", "a2.default", "C1", "m"),
        _member("a1.default-C2-m.default", "a1.default", "C2", "m"),
        _member("a2.default-C2-m.default", "a2.default", "C2", "m"),
    ]
    nodes_by_id = {n.id: n for n in nodes}
    output_to_nodes = {
        "c1.out": [
            ("a1.default-C1-m.default", "p1"),
            ("a2.default-C1-m.default", "p2"),
        ],
        "c2.out": [
            ("a1.default-C2-m.default", "q1"),
            ("a2.default-C2-m.default", "q2"),
        ],
    }

    bundles = select_input_bundles(
        declared_input_ids=["c1.out", "c2.out"],
        output_to_nodes=output_to_nodes,
        resolved_nodes=nodes,
        stage_ids_in_order=["A", "C1", "C2", "E"],
        previous_stage_nodes=[n for n in nodes if n.stage_id == "C2"],
        nodes_by_id=nodes_by_id,
    )

    assert sorted(tuple(n.id for n in b) for b in bundles) == [
        ("a1.default-C2-m.default", "a1.default-C1-m.default"),
        ("a2.default-C2-m.default", "a2.default-C1-m.default"),
    ]

    # Linear case: input covered by the anchor's own lineage → 1-tuples.
    linear = select_input_bundles(
        declared_input_ids=["c2.out"],
        output_to_nodes=output_to_nodes,
        resolved_nodes=nodes,
        stage_ids_in_order=["A", "C1", "C2", "E"],
        previous_stage_nodes=[n for n in nodes if n.stage_id == "C2"],
        nodes_by_id=nodes_by_id,
    )
    assert all(len(b) == 1 for b in linear)
    assert {b[0].stage_id for b in linear} == {"C2"}


def test_gather_stage_is_ordered_after_its_producers():
    """A gather stage declares no `inputs:`, so build_stage_dag must derive its
    edges from `gather.from` — otherwise topological expansion (#289/#367)
    could expand the gather before its members exist. Declare the gather stage
    FIRST to prove declaration order is irrelevant."""
    from omnibenchmark.core._graph import build_stage_dag, compute_stage_order
    from omnibenchmark.model.benchmark import Benchmark

    y = """
id: t
description: t
version: '1.0'
benchmarker: me
api_version: 0.7.0
software_backend: host
software_environments:
  env: {description: e, easyconfig: e.eb}
stages:
  - id: metrics
    prefix: agg
    gather: [{from: clustering, group_by: data}]
    modules: [{id: s, repository: {url: 'http://x', commit: abc}, software_environment: env}]
    outputs: [{id: metrics.summary, path: summary.tsv}]
  - id: data
    outputs: [{id: clustering, path: c.tsv}]
    modules: [{id: d1, repository: {url: 'http://x', commit: abc}, software_environment: env}]
"""
    bench = Benchmark.from_yaml(y)
    dag = build_stage_dag(bench)
    assert ("data", "metrics") in dag.edges
    order = compute_stage_order(dag)
    assert order.index("data") < order.index("metrics")


def test_gather_context_does_not_bind_builtin_dataset():
    """A gather node binds exactly its group key (+ name) — the builtin
    `dataset` label must NOT leak in as the gather module's own id, or it
    shadows the real dataset for every downstream template/requires match."""
    nodes_by_id = {
        "d1.default": _member("d1.default", None, "data", "d1"),
        "d1.default-clu-ma.default": _member(
            "d1.default-clu-ma.default", "d1.default", "clu", "ma"
        ),
    }
    output_to_nodes = {"clustering": [("d1.default-clu-ma.default", "p1")]}
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
    (node,) = expand_gather_stage(
        stage=stage,
        benchmark=_fake_benchmark(),
        resolved_modules_cache={("metrics", "summ"): object()},
        output_to_nodes=output_to_nodes,
        nodes_by_id=nodes_by_id,
    )
    assert "dataset" not in node.template_context.provides
    assert node.template_context.provides["data"] == "d1"


def test_select_input_bundles_rejects_cross_lineage_partners():
    """Sharing a root is not enough: partners must agree with the anchor at
    EVERY shared stage. data A -> process B (b1, b2) -> divergent C1, C2 ->
    join: the b1-chain anchor must never pair with a b2-chain producer."""
    nodes = [_member("a.default", None, "A", "a")]
    for b in ("b1", "b2"):
        bid = f"a.default-B-{b}.default"
        nodes.append(_member(bid, "a.default", "B", b))
        for c_stage, m in (("C1", "m1"), ("C2", "m2")):
            nodes.append(_member(f"{bid}-{c_stage}-{m}.default", bid, c_stage, m))
    nodes_by_id = {n.id: n for n in nodes}
    output_to_nodes = {
        "c1.out": [(n.id, f"{n.id}/p") for n in nodes if n.stage_id == "C1"],
        "c2.out": [(n.id, f"{n.id}/q") for n in nodes if n.stage_id == "C2"],
    }

    bundles = select_input_bundles(
        declared_input_ids=["c1.out", "c2.out"],
        output_to_nodes=output_to_nodes,
        resolved_nodes=nodes,
        stage_ids_in_order=["A", "B", "C1", "C2", "E"],
        previous_stage_nodes=[n for n in nodes if n.stage_id == "C2"],
        nodes_by_id=nodes_by_id,
    )

    assert sorted(tuple(n.id for n in b) for b in bundles) == [
        (
            "a.default-B-b1.default-C2-m2.default",
            "a.default-B-b1.default-C1-m1.default",
        ),
        (
            "a.default-B-b2.default-C2-m2.default",
            "a.default-B-b2.default-C1-m1.default",
        ),
    ]


def test_select_input_bundles_join_anchor_sees_parents_ancestry():
    """An anchor that is itself a fan-in (hash id, ancestry only in .parents)
    must not classify inputs covered by its true ancestry as missing."""
    a = _member("a.default", None, "A", "a")
    b = _member("a.default-B-mb.default", "a.default", "B", "mb")
    c = _member("a.default-C-mc.default", "a.default", "C", "mc")
    # Join node: hash id (no prefix chain), parent_id = primary, parents = both.
    j = _member("D-md-cafe1234.default", b.id, "D", "md")
    j.parents = [b.id, c.id]
    nodes = [a, b, c, j]
    nodes_by_id = {n.id: n for n in nodes}
    output_to_nodes = {
        "d.out": [(j.id, "d/p")],
        "b.out": [(b.id, "b/p")],
    }

    bundles = select_input_bundles(
        declared_input_ids=["d.out", "b.out"],
        output_to_nodes=output_to_nodes,
        resolved_nodes=nodes,
        stage_ids_in_order=["A", "B", "C", "D", "E"],
        previous_stage_nodes=[j],
        nodes_by_id=nodes_by_id,
    )

    # b.out is covered by the anchor's parents-ancestry: one linear 1-tuple,
    # no spurious re-pairing with B producers.
    assert bundles == [(j,)]


def test_gather_respects_module_filter():
    """`ob run -m X`: gather expands only the target module and one combo."""
    nodes_by_id = {
        "d1.default": _member("d1.default", None, "data", "d1"),
        "d2.default": _member("d2.default", None, "data", "d2"),
        "d1.default-clu-ma.default": _member(
            "d1.default-clu-ma.default", "d1.default", "clu", "ma"
        ),
        "d2.default-clu-ma.default": _member(
            "d2.default-clu-ma.default", "d2.default", "clu", "ma"
        ),
    }
    output_to_nodes = {
        "clustering": [
            ("d1.default-clu-ma.default", "d1/p"),
            ("d2.default-clu-ma.default", "d2/p"),
        ]
    }
    stage = SimpleNamespace(
        id="metrics",
        prefix="agg",
        gather=[SimpleNamespace(from_="clustering", group_by="data")],
        modules=[
            SimpleNamespace(
                id="s1", name="s1", parameters=None, provides=None, resources=None
            ),
            SimpleNamespace(
                id="s2", name="s2", parameters=None, provides=None, resources=None
            ),
        ],
        outputs=[SimpleNamespace(id="metrics.summary", path="{data}.tsv")],
        resources=None,
    )
    cache = {("metrics", "s1"): object(), ("metrics", "s2"): object()}

    nodes = expand_gather_stage(
        stage=stage,
        benchmark=_fake_benchmark(),
        resolved_modules_cache=cache,
        output_to_nodes=dict(output_to_nodes),
        nodes_by_id=nodes_by_id,
        module_filter="s2",
        target_stage=stage,
    )
    # Target stage: only the named module, first combo only.
    assert [n.module_id for n in nodes] == ["s2"]

    nodes = expand_gather_stage(
        stage=stage,
        benchmark=_fake_benchmark(),
        resolved_modules_cache=cache,
        output_to_nodes=dict(output_to_nodes),
        nodes_by_id=nodes_by_id,
        module_filter="other",
        target_stage=None,
    )
    # Non-target stage under -m: first module, first combo only.
    assert [(n.module_id, len(n.gathered_from)) for n in nodes] == [("s1", 1)]


def test_gather_applies_exclusions_at_member_level():
    """An exclude pairing the gather module with a member's lineage drops that
    member from that module's gather — it must not poison the whole group."""
    nodes_by_id = {
        "d1.default": _member("d1.default", None, "data", "d1"),
        "d1.default-clu-ma.default": _member(
            "d1.default-clu-ma.default", "d1.default", "clu", "ma"
        ),
        "d1.default-clu-mb.default": _member(
            "d1.default-clu-mb.default", "d1.default", "clu", "mb"
        ),
    }
    output_to_nodes = {
        "clustering": [
            ("d1.default-clu-ma.default", "d1/ma/p"),
            ("d1.default-clu-mb.default", "d1/mb/p"),
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

    (node,) = expand_gather_stage(
        stage=stage,
        benchmark=_fake_benchmark(),
        resolved_modules_cache={("metrics", "summ"): object()},
        output_to_nodes=output_to_nodes,
        nodes_by_id=nodes_by_id,
        path_exclusions={"summ": ["ma"]},
    )
    assert node.gathered_from == ["d1.default-clu-mb.default"]
    assert set(node.inputs.values()) == {"d1/mb/p"}


def test_gather_groups_through_join_partner_branch():
    """Grouping must see ancestors on a join's PARTNER branch (parents edges),
    not just the parent_id spine."""
    a = _member("d1.default", None, "data", "d1")
    b = _member("d1.default-B-mb.default", "d1.default", "B", "mb")
    c = _member("d1.default-C-mc.default", "d1.default", "C", "mc")
    j = _member("J-mj-beef0001.default", b.id, "J", "mj")
    j.parents = [b.id, c.id]
    nodes_by_id = {n.id: n for n in (a, b, c, j)}
    stage = SimpleNamespace(
        id="metrics",
        prefix="agg",
        # group_by C: reachable only via the join's partner branch.
        gather=[SimpleNamespace(from_="j.out", group_by="C")],
        modules=[
            SimpleNamespace(
                id="summ", name="summ", parameters=None, provides=None, resources=None
            )
        ],
        outputs=[SimpleNamespace(id="metrics.summary", path="{C}.tsv")],
        resources=None,
    )

    (node,) = expand_gather_stage(
        stage=stage,
        benchmark=_fake_benchmark(),
        resolved_modules_cache={("metrics", "summ"): object()},
        output_to_nodes={"j.out": [(j.id, "j/p")]},
        nodes_by_id=nodes_by_id,
    )
    assert node.template_context.provides["C"] == "mc"


def test_build_stage_dag_skips_self_edge_for_shared_output_id():
    """A stage re-declaring an output id it consumes must not get a self-edge
    (it would make the topological sort reject a valid plan)."""
    from omnibenchmark.core._graph import build_stage_dag, compute_stage_order
    from omnibenchmark.model.benchmark import Benchmark

    y = """
id: t
description: t
version: '1.0'
benchmarker: me
api_version: 0.7.0
software_backend: host
software_environments:
  env: {description: e, easyconfig: e.eb}
stages:
  - id: raw
    outputs: [{id: counts, path: c.tsv}]
    modules: [{id: r, repository: {url: 'http://x', commit: abc}, software_environment: env}]
  - id: refine
    inputs: [{entries: [counts]}]
    outputs: [{id: counts, path: c2.tsv}]
    modules: [{id: f, repository: {url: 'http://x', commit: abc}, software_environment: env}]
"""
    bench = Benchmark.from_yaml(y)
    dag = build_stage_dag(bench)
    assert ("refine", "refine") not in dag.edges
    order = compute_stage_order(dag)
    assert order.index("raw") < order.index("refine")


def test_stage_adjacency_includes_gather_edges():
    """Topology viz (mermaid/dot/obeditor via stage_adjacency) must show the
    gather dependency even though gather declares no `inputs:`."""
    from omnibenchmark.core._graph import stage_adjacency
    from omnibenchmark.model.benchmark import Benchmark

    bench = Benchmark.from_yaml(_GATHER_YAML.format(api="0.7.0", group_by="data"))
    assert ("data", "metrics", ["clustering"]) in stage_adjacency(bench)


# ---------------------------------------------------------------------------
# Fan-in output paths must distinguish parent sets.
#
# A join's path prefix (`base_path`) is the parent of its DEEPEST input, i.e.
# one branch only. Two joins sharing that branch but differing in another
# parent therefore landed on the same output path, and Snakemake rejected the
# pair with AmbiguousRuleException. Reproduced with an asymmetric diamond:
# root -> shallow (two modules) and root -> mid -> deep, joined together.
# ---------------------------------------------------------------------------


@pytest.mark.short
def test_join_expansion_segment_distinguishes_parent_sets():
    """Same anchor, different partner => different directory segment."""
    deep = _member("root-R1-mid-M1-deep-D1", "root-R1-mid-M1", "deep", "D1")
    s1 = _member("root-R1-shallow-S1", "root-R1", "shallow", "S1")
    s2 = _member("root-R1-shallow-S2", "root-R1", "shallow", "S2")

    seg1 = expansion_segment(".abc12345", (deep, s1))
    seg2 = expansion_segment(".abc12345", (deep, s2))

    assert seg1 != seg2, (
        "two joins sharing their deepest input must not share an output "
        "directory; base_path is identical for both, so the segment is the "
        "only thing that can separate them"
    )
    assert seg1.startswith(".abc12345-") and seg2.startswith(".abc12345-")
    # keyed on the parent *set*: bundle order must not leak into the path,
    # or the same join would land in two directories across runs.
    assert seg1 == expansion_segment(".abc12345", (s1, deep))


@pytest.mark.short
def test_linear_expansion_segment_is_just_the_parameter_hash():
    """Non-join nodes keep the plain parameter segment — their ancestry is
    already carried by the path prefix, so nothing may change for them."""
    only = _member("root-R1-stage-M1", "root-R1", "stage", "M1")

    assert expansion_segment(".abc12345", (only,)) == ".abc12345"
    assert expansion_segment(".default", (only,)) == ".default"
    assert expansion_segment(".default", ()) == ".default"


@pytest.mark.short
def test_human_link_name_is_unique_per_join():
    """The readable sibling symlink must not collide either.

    `ln -sfn` overwrites, so two joins sharing a module and parameters — they
    differ only in parents — would leave a single link pointing at whichever
    job happened to run last, silently hiding the other branch.
    """
    params = SimpleNamespace(
        items=lambda: {"evaluate": "input+10"}.items(),
        hash_short=lambda: "abc12345",
    )
    left = SimpleNamespace(parameters=params, parents=["deep-D1", "shallow-S1"])
    right = SimpleNamespace(parameters=params, parents=["deep-D1", "shallow-S2"])
    linear = SimpleNamespace(parameters=params, parents=[])

    assert _human_link_name(left) != _human_link_name(right)
    # Non-join nodes keep the plain readable name.
    assert _human_link_name(linear) == _make_human_name(params)
