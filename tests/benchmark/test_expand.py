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


@pytest.mark.short
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
    # Group value bound to the template + baked into the path, which roots at
    # the stage id.
    assert d1.outputs == ["metrics/d1/summ/.default/d1_summary.tsv"]
    # Registered downstream so a later stage can consume it.
    assert ("metrics-summ-d1.default", d1.outputs[0]) in output_to_nodes[
        "metrics.summary"
    ]


@pytest.mark.short
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


@pytest.mark.short
def test_zero_producer_from_is_plan_time_error():
    stage = SimpleNamespace(
        id="metrics",
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
    gather: [{{from: clustering, group_by: {group_by}}}]
    modules: [{{id: s, repository: {{url: 'http://x', commit: abc}}, software_environment: env}}]
    outputs: [{{id: metrics.summary, path: summary.tsv}}]
"""


@pytest.mark.short
def test_gather_benchmark_parses_and_from_alias():
    from omnibenchmark.model.benchmark import Benchmark

    bench = Benchmark.from_yaml(_GATHER_YAML.format(api="0.7.0", group_by="data"))
    gather_stage = bench.stages[1]
    assert gather_stage.gather[0].from_ == "clustering"
    assert gather_stage.gather[0].group_by == "data"


@pytest.mark.short
def test_gather_gated_on_api_0_7():
    from omnibenchmark.model.benchmark import Benchmark

    try:
        Benchmark.from_yaml(_GATHER_YAML.format(api="0.5.0", group_by="data"))
    except Exception as e:
        assert "gather" in str(e) and "0.7.0" in str(e)
    else:
        raise AssertionError("expected gather to be gated on api >= 0.7.0")


@pytest.mark.short
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


@pytest.mark.short
def test_gather_entries_must_share_one_group_by_axis():
    """Multiple gather entries with differing group_by are rejected; a shared
    axis is accepted (Stage.validate_gather — makes the [0] axis sound)."""
    # Same axis across entries: OK.
    Stage(
        id="metrics",
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


@pytest.mark.short
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


@pytest.mark.short
def test_group_by_must_name_a_stage():
    from omnibenchmark.model.benchmark import Benchmark

    try:
        Benchmark.from_yaml(_GATHER_YAML.format(api="0.7.0", group_by="nosuch"))
    except Exception as e:
        assert "not a known stage" in str(e)
    else:
        raise AssertionError("expected group_by to be validated against stage ids")


_HEAD = """
id: t
description: t
version: '1.0'
benchmarker: me
api_version: 0.7.0
software_backend: host
software_environments:
  env: {description: e, easyconfig: e.eb}
stages:
"""

_REPO = "repository: {url: 'http://x', commit: abc}, software_environment: env"


def _parse_error(yaml_text):
    """`from_yaml` wraps validation in BenchmarkParseError, which is not a
    ValueError — catch broadly and assert on the message, as the rest of this
    file does."""
    from omnibenchmark.model.benchmark import Benchmark

    try:
        Benchmark.from_yaml(yaml_text)
    except Exception as e:
        return str(e)
    raise AssertionError("expected the benchmark to be rejected")


@pytest.mark.short
def test_gather_stage_is_not_an_initial_stage():
    """A gather consumes `gather.from` from real producers, so it is not a root.

    `is_initial` looked only at `inputs:`, which a gather replaces, so every
    gather read as a root and the initial-stage `requires` check fired on it.
    """
    from omnibenchmark.model.benchmark import Benchmark

    bench = Benchmark.from_yaml(_GATHER_YAML.format(api="0.7.0", group_by="data"))
    data, metrics = bench.stages
    assert bench.get_stage_implicit_inputs(metrics) == [["clustering"]]
    assert bench.get_stage_implicit_inputs(data) == []


@pytest.mark.short
def test_requires_on_gather_module_rejected():
    """A gather cuts the chain, so a member-lineage gate has nothing to match.

    `expand_gather_stage` never evaluates `requires`, so accepting it would make
    the gate silently inert (008 §2). The message must name the gather rather
    than claim the stage has no upstream lineage — it has plenty, just not
    reachable across the cut.
    """
    message = _parse_error(
        _HEAD
        + f"""  - id: data
    outputs: [{{id: clustering, path: c.tsv}}]
    modules: [{{id: d1, {_REPO}}}]
  - id: metrics
    gather: [{{from: clustering, group_by: data}}]
    modules: [{{id: s, {_REPO}, requires: {{data: d1}}}}]
    outputs: [{{id: metrics.summary, path: summary.tsv}}]
"""
    )
    assert "gather stage" in message
    assert "initial stage" not in message


@pytest.mark.short
def test_requires_on_input_less_stage_still_rejected():
    """The original rule survives: a stage with neither `inputs` nor `gather`."""
    message = _parse_error(
        _HEAD
        + f"""  - id: data
    outputs: [{{id: clustering, path: c.tsv}}]
    modules: [{{id: d1, {_REPO}, requires: {{size: lg}}}}]
"""
    )
    assert "initial stage" in message


@pytest.mark.short
def test_provides_label_on_the_group_by_stage_rejected():
    """The two values the group key could take must never both exist.

    Before the cut a node carries `Module.provides['data']`; after it the gather
    binds `data` to the ancestor module id. A downstream `requires: {data: custom}`
    would match upstream and prune downstream with no diagnostic. Rejecting the
    declaration is what keeps that unreachable.
    """
    message = _parse_error(
        _HEAD
        + f"""  - id: data
    provides: [data]
    outputs: [{{id: clustering, path: c.tsv}}]
    modules: [{{id: d1, {_REPO}, provides: {{data: custom}}}}]
  - id: metrics
    gather: [{{from: clustering, group_by: data}}]
    modules: [{{id: s, {_REPO}}}]
    outputs: [{{id: metrics.summary, path: summary.tsv}}]
"""
    )
    assert "also a stage id" in message


@pytest.mark.short
def test_group_by_reserved_stage_name_rejected():
    """A stage called `dataset` cannot be a `group_by` target: the group key is
    bound under the stage id and would clobber the runtime builtin."""
    message = _parse_error(
        _HEAD
        + f"""  - id: dataset
    outputs: [{{id: clustering, path: c.tsv}}]
    modules: [{{id: d1, {_REPO}}}]
  - id: metrics
    gather: [{{from: clustering, group_by: dataset}}]
    modules: [{{id: s, {_REPO}}}]
    outputs: [{{id: metrics.summary, path: summary.tsv}}]
"""
    )
    assert "reserved builtin label" in message


@pytest.mark.short
def test_stage_named_dataset_without_gather_accepted():
    """Only the `group_by` target is constrained. A stage merely *called*
    `dataset` is legal in released 0.6 specs and stays legal."""
    from omnibenchmark.model.benchmark import Benchmark

    bench = Benchmark.from_yaml(
        _HEAD
        + f"""  - id: dataset
    outputs: [{{id: clustering, path: c.tsv}}]
    modules: [{{id: d1, {_REPO}}}]
  - id: metrics
    inputs: [clustering]
    modules: [{{id: s, {_REPO}}}]
    outputs: [{{id: metrics.summary, path: summary.tsv}}]
"""
    )
    assert [s.id for s in bench.stages] == ["dataset", "metrics"]


@pytest.mark.short
def test_select_input_bundles_pairs_diamond_branches_by_root():
    """The fan-in join (design 010 §5.2, #289): a stage declaring inputs from
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


@pytest.mark.short
def test_join_partners_drop_shadowed_producers():
    """A join partner is chosen with the same shadowing rule as an anchor
    (design 010 §3.1): B declares `x.out` at `b1` and again at its descendant
    `b2`, so `b1` is shadowed and only the `b2` bundle exists. Without the
    filter the cartesian product emits one bundle per producer, and the extra
    node runs the module against the upstream file."""
    nodes = [
        _member("root.default", None, "root", "r"),
        _member("b1.default", "root.default", "b1", "mb1"),
        _member("b1.default-b2.default", "b1.default", "b2", "mb2"),
        _member("a.default", "root.default", "a", "ma"),
    ]
    nodes_by_id = {n.id: n for n in nodes}
    output_to_nodes = {
        "x.out": [("b1.default", "b1/x.tsv"), ("b1.default-b2.default", "b2/x.tsv")],
        "a.out": [("a.default", "a/a.tsv")],
    }

    bundles = select_input_bundles(
        declared_input_ids=["a.out", "x.out"],
        output_to_nodes=output_to_nodes,
        resolved_nodes=nodes,
        stage_ids_in_order=["root", "b1", "b2", "a"],
        previous_stage_nodes=[],
        nodes_by_id=nodes_by_id,
    )

    assert [tuple(n.id for n in b) for b in bundles] == [
        ("a.default", "b1.default-b2.default")
    ]


@pytest.mark.short
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


@pytest.mark.short
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


@pytest.mark.short
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


@pytest.mark.short
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


@pytest.mark.short
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


@pytest.mark.short
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


@pytest.mark.short
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


@pytest.mark.short
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


@pytest.mark.short
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


@pytest.mark.short
def test_global_gather_collects_every_producer_into_one_node():
    """No `group_by` is the `metric_collector` shape (design 010 §3.6): one
    node over every producer, and no group segment in the id or the path."""
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
            ("d1.default-clu-ma.default", "d1/clu/ma/a.tsv"),
            ("d2.default-clu-ma.default", "d2/clu/ma/a.tsv"),
        ]
    }
    stage = SimpleNamespace(
        id="report",
        gather=[SimpleNamespace(from_="clustering", group_by=None)],
        modules=[
            SimpleNamespace(
                id="R", name="R", parameters=None, provides=None, resources=None
            )
        ],
        outputs=[SimpleNamespace(id="report.html", path="report.html")],
        resources=None,
    )

    nodes = expand_gather_stage(
        stage=stage,
        benchmark=_fake_benchmark(),
        resolved_modules_cache={("report", "R"): object()},
        output_to_nodes=output_to_nodes,
        nodes_by_id=nodes_by_id,
    )

    assert len(nodes) == 1
    node = nodes[0]
    assert node.id == "report-R.default"
    assert node.outputs == ["report/R/.default/report.html"]
    # Both datasets contribute — grouping is what a global gather forgoes.
    assert sorted(node.gathered_from) == [
        "d1.default-clu-ma.default",
        "d2.default-clu-ma.default",
    ]
    # No group label, and the `dataset` builtin must not leak in its place.
    assert node.template_context.provides == {"name": "R"}


@pytest.mark.short
def test_global_and_grouped_gather_entries_cannot_mix():
    """One axis per stage, and "no axis" is one of them."""
    with pytest.raises(ValueError, match="differing group_by"):
        Stage(
            id="report",
            modules=[],
            outputs=[],
            gather=[
                GatherSpec(from_="clustering", group_by="data"),
                GatherSpec(from_="metrics"),
            ],
        )


@pytest.mark.short
def test_partially_populated_group_is_an_error():
    """Two gather entries, one dataset producing only one of the two ids.

    The node would carry a single `--methods.result` flag and no
    `--methods.other`, failing inside the module's argparse. Reject at plan
    time (design 010 §3.2).
    """
    nodes_by_id = {
        "d1.default": _member("d1.default", None, "data", "d1"),
        "d2.default": _member("d2.default", None, "data", "d2"),
        "d1.default-other.default": _member(
            "d1.default-other.default", "d1.default", "methods", "other"
        ),
        "d1.default-res.default": _member(
            "d1.default-res.default", "d1.default", "methods", "res"
        ),
        # d2 has no `methods.other` producer — the branch is excluded there.
        "d2.default-res.default": _member(
            "d2.default-res.default", "d2.default", "methods", "res"
        ),
    }
    output_to_nodes = {
        "methods.other": [("d1.default-other.default", "d1/other.tsv")],
        "methods.result": [
            ("d1.default-res.default", "d1/res.tsv"),
            ("d2.default-res.default", "d2/res.tsv"),
        ],
    }
    stage = SimpleNamespace(
        id="metrics",
        gather=[
            SimpleNamespace(from_="methods.other", group_by="data"),
            SimpleNamespace(from_="methods.result", group_by="data"),
        ],
        modules=[
            SimpleNamespace(
                id="MC", name="MC", parameters=None, provides=None, resources=None
            )
        ],
        outputs=[SimpleNamespace(id="metrics.out", path="{data}_out.tsv")],
        resources=None,
    )

    with pytest.raises(ValueError, match="methods.other"):
        expand_gather_stage(
            stage=stage,
            benchmark=_fake_benchmark(),
            resolved_modules_cache={("metrics", "MC"): object()},
            output_to_nodes=output_to_nodes,
            nodes_by_id=nodes_by_id,
        )


# --- the scatter path -------------------------------------------------------
#
# `expand_scatter_stage` is the other half of the planner and the one every
# benchmark goes through, but nothing unit-tested it — the gather tests above
# exercise its sibling, and only e2e reached this. These drive it over a real
# parsed model, the way `cli/run.py` does, with the module cache faked so no
# repository is cloned.


def _plan(yaml_text, nesting_strategy="nested"):
    """Expand every stage of a parsed benchmark, mirroring `cli/run.py`'s loop.

    Returns `(nodes, prune_counts, dag_errors)`. Module resolution is faked:
    the cache maps every (stage, module) to a sentinel, which is all the
    expander does with it (stashes it on the node).
    """
    from omnibenchmark.core._expand import expand_scatter_stage
    from omnibenchmark.core._graph import build_stage_dag, compute_stage_order
    from omnibenchmark.core._paths import collect_path_exclusions
    from omnibenchmark.model.benchmark import Benchmark

    bench = Benchmark.from_yaml(yaml_text)
    benchmark = SimpleNamespace(model=bench)
    by_id = {s.id: s for s in bench.stages}
    stages_to_expand = [
        by_id[sid] for sid in compute_stage_order(build_stage_dag(bench))
    ]
    cache = {(s.id, m.id): object() for s in bench.stages for m in s.modules}

    resolved_nodes, nodes_by_id, output_to_nodes = [], {}, {}
    prune_counts = {"requires": 0, "exclude": 0}
    dag_errors = []
    previous = []
    for stage in stages_to_expand:
        if stage.gather:
            previous = expand_gather_stage(
                stage=stage,
                benchmark=benchmark,
                resolved_modules_cache=cache,
                output_to_nodes=output_to_nodes,
                nodes_by_id=nodes_by_id,
            )
            resolved_nodes.extend(previous)
            nodes_by_id.update({n.id: n for n in previous})
        else:
            previous = expand_scatter_stage(
                stage=stage,
                benchmark=benchmark,
                resolved_modules_cache=cache,
                resolved_nodes=resolved_nodes,
                nodes_by_id=nodes_by_id,
                output_to_nodes=output_to_nodes,
                previous_stage_nodes=previous,
                stages_to_expand=stages_to_expand,
                path_exclusions=collect_path_exclusions(bench),
                nesting_strategy=nesting_strategy,
                module_filter=None,
                target_stage=None,
                dag_errors=dag_errors,
                prune_counts=prune_counts,
                quiet=True,
            )
    return resolved_nodes, prune_counts, dag_errors


_BLOCK_REPO = """        repository: {url: 'http://x', commit: abc}
        software_environment: env"""


def _line(text, indent):
    """One optional extra YAML line at `indent` spaces, or nothing."""
    return f"\n{' ' * indent}{text}" if text else ""


def _chain_yaml(data_extra="", d1_extra="", method_module_extra=""):
    """data (D1, D2) -> method (M1). Every test below varies one knob on it."""
    return f"""
id: t
description: t
version: '1.0'
benchmarker: me
api_version: 0.7.0
software_backend: host
software_environments:
  env: {{description: e, easyconfig: e.eb}}
stages:
  - id: data{_line(data_extra, 4)}
    outputs: [{{id: data.out, path: "{{name}}_d.txt"}}]
    modules:
      - id: D1
{_BLOCK_REPO}{_line(d1_extra, 8)}
      - id: D2
{_BLOCK_REPO}
  - id: method
    inputs: [data.out]
    outputs: [{{id: method.out, path: "{{name}}_m.txt"}}]
    modules:
      - id: M1
{_BLOCK_REPO}{_line(method_module_extra, 8)}
"""


@pytest.mark.short
def test_scatter_chains_each_module_onto_every_upstream_node():
    """The base case: 2 data modules x 1 method = 2 method nodes, each parented
    to its own upstream node and nested under its directory."""
    nodes, prune, errors = _plan(_chain_yaml())
    assert errors == []

    data = {n.id: n for n in nodes if n.stage_id == "data"}
    method = {n.id: n for n in nodes if n.stage_id == "method"}
    assert len(data) == 2 and len(method) == 2

    # A root node has no parent and roots its own tree; that is the "layout"
    # half of what `is_initial` used to answer.
    for node in data.values():
        assert node.parent_id is None and node.parents == []
    assert sorted(n.outputs[0] for n in data.values()) == [
        "data/D1/.default/D1_d.txt",
        "data/D2/.default/D2_d.txt",
    ]

    # Each method node extends exactly one data node's directory.
    for node in method.values():
        assert node.parent_id in data
        assert node.outputs[0].startswith(
            data[node.parent_id].outputs[0].rsplit("/", 1)[0]
        )
        assert node.inputs == {"data_out": data[node.parent_id].outputs[0]}
        assert node.input_name_mapping == {"data_out": "data.out"}


@pytest.mark.short
def test_scatter_expands_one_node_per_parameter_set():
    """Parameters are a fan-out axis, and the param hash separates the nodes."""
    nodes, _, errors = _plan(
        _chain_yaml(
            method_module_extra="parameters: [{values: ['-k', '1']}, {values: ['-k', '2']}]"
        )
    )
    assert errors == []
    method = [n for n in nodes if n.stage_id == "method"]
    # 2 data nodes x 2 parameter sets.
    assert len(method) == 4
    assert all(n.param_id != ".default" for n in method)
    # Distinct paths, so Snakemake sees four rules rather than two collisions.
    assert len({n.outputs[0] for n in method}) == 4


@pytest.mark.short
def test_scatter_prunes_excluded_lineages():
    """`exclude` drops the combination and is counted, not silently skipped."""
    nodes, prune, errors = _plan(_chain_yaml(method_module_extra="exclude: [D2]"))
    assert errors == []
    method = [n for n in nodes if n.stage_id == "method"]
    assert len(method) == 1, [n.id for n in method]
    assert "D1" in method[0].id
    assert prune["exclude"] == 1


@pytest.mark.short
def test_scatter_prunes_by_requires_against_upstream_labels():
    """`requires` matches the resolved label, so `Module.provides` decides.

    D1 binds `size: lg`, D2 falls through to the module-id default, so only
    D1's lineage satisfies the gate.
    """
    nodes, prune, errors = _plan(
        _chain_yaml(
            data_extra="provides: [size]",
            d1_extra="provides: {size: lg}",
            method_module_extra="requires: {size: lg}",
        )
    )
    assert errors == []
    method = [n for n in nodes if n.stage_id == "method"]
    assert len(method) == 1, [n.id for n in method]
    assert "D1" in method[0].id
    assert prune["requires"] == 1


@pytest.mark.short
def test_flat_nesting_drops_the_parent_prefix():
    """`flat` roots every stage at its own id instead of extending the parent."""
    nodes, _, errors = _plan(_chain_yaml(), nesting_strategy="flat")
    assert errors == []
    method = [n for n in nodes if n.stage_id == "method"]
    assert {n.outputs[0] for n in method} == {"method/M1/.default/M1_m.txt"}


@pytest.mark.short
def test_unknown_nesting_strategy_is_reported_as_a_dag_error():
    """A bad strategy raises inside the module loop, which records it against
    the stage rather than crashing the whole plan."""
    nodes, _, errors = _plan(_chain_yaml(), nesting_strategy="sideways")
    assert errors, "expected the ValueError to be captured"
    assert any("sideways" in msg for _, _, msg in errors)


@pytest.mark.short
def test_module_missing_from_the_resolution_cache_is_skipped():
    """A module whose repository failed to resolve produces no nodes, and the
    stage warns rather than cascading an empty set silently."""
    from omnibenchmark.core._expand import expand_scatter_stage
    from omnibenchmark.model.benchmark import Benchmark

    bench = Benchmark.from_yaml(_chain_yaml())
    data_stage = bench.stages[0]
    nodes = expand_scatter_stage(
        stage=data_stage,
        benchmark=SimpleNamespace(model=bench),
        resolved_modules_cache={},  # nothing resolved
        resolved_nodes=[],
        nodes_by_id={},
        output_to_nodes={},
        previous_stage_nodes=[],
        stages_to_expand=[data_stage],
        path_exclusions={},
        nesting_strategy="nested",
        module_filter=None,
        target_stage=None,
        dag_errors=[],
        prune_counts={"requires": 0, "exclude": 0},
        quiet=True,
    )
    assert nodes == []


@pytest.mark.short
def test_scatter_joins_divergent_branches_into_one_node():
    """The fan-in join (#289): a stage consuming two output ids from divergent
    branches gets ONE node holding both parents, not one node per branch.

    The id cannot prefix-compose off a single chain, so it is a readable stem
    plus a hash of the parent set; `parents` carries the real edges.
    """
    nodes, _, errors = _plan(f"""
id: t
description: t
version: '1.0'
benchmarker: me
api_version: 0.7.0
software_backend: host
software_environments:
  env: {{description: e, easyconfig: e.eb}}
stages:
  - id: root
    outputs: [{{id: root.out, path: r.txt}}]
    modules:
      - id: R1
{_BLOCK_REPO}
  - id: left
    inputs: [root.out]
    outputs: [{{id: left.out, path: l.txt}}]
    modules:
      - id: L1
{_BLOCK_REPO}
  - id: right
    inputs: [root.out]
    outputs: [{{id: right.out, path: rt.txt}}]
    modules:
      - id: RT1
{_BLOCK_REPO}
  - id: join
    inputs: [left.out, right.out]
    outputs: [{{id: join.out, path: j.txt}}]
    modules:
      - id: J1
{_BLOCK_REPO}
""")
    assert errors == []
    joins = [n for n in nodes if n.stage_id == "join"]
    assert len(joins) == 1, [n.id for n in joins]
    node = joins[0]

    # Both branches feed it, and both are recorded as parents.
    assert len(node.parents) == 2
    assert set(node.inputs) == {"left_out", "right_out"}
    assert node.input_name_mapping == {"left_out": "left.out", "right_out": "right.out"}

    # `parent_id` is the anchor (one branch), so the id cannot be the chain —
    # it carries a digest of the parent set instead.
    assert node.parent_id in node.parents
    assert node.id.startswith("join-J1-")
    assert node.id != f"{node.parent_id}-join-J1{node.param_id}"


@pytest.mark.short
def test_module_resources_win_over_stage_resources():
    """Resources fall back stage -> module, most specific first."""
    nodes, _, errors = _plan(
        _chain_yaml(method_module_extra="resources: {mem_mb: 512}")
    )
    assert errors == []
    method = [n for n in nodes if n.stage_id == "method"]
    assert all(n.resources is not None for n in method)
    assert all(n.resources.mem_mb == 512 for n in method)
    # The data stage declares none, so its nodes carry none.
    assert all(n.resources is None for n in nodes if n.stage_id == "data")


@pytest.mark.short
def test_gather_drops_a_member_with_no_ancestor_in_the_group_by_stage(caplog):
    """Membership is silent-absence-free (§2): a member that cannot be placed
    in any group is dropped with a warning, not into an arbitrary bucket."""
    import logging

    nodes_by_id = {
        "d1.default": _member("d1.default", None, "data", "d1"),
        "d1.default-clu-ma.default": _member(
            "d1.default-clu-ma.default", "d1.default", "clu", "ma"
        ),
        # No `data` ancestor: this one descends from nothing.
        "orphan.default": _member("orphan.default", None, "clu", "mb"),
    }
    output_to_nodes = {
        "clustering": [
            ("d1.default-clu-ma.default", "d1/clu/ma/a.tsv"),
            ("orphan.default", "orphan/b.tsv"),
        ]
    }
    stage = SimpleNamespace(
        id="metrics",
        gather=[SimpleNamespace(from_="clustering", group_by="data")],
        modules=[
            SimpleNamespace(
                id="summ", name="summ", parameters=None, provides=None, resources=None
            )
        ],
        outputs=[SimpleNamespace(id="metrics.summary", path="{data}.tsv")],
        resources=None,
    )

    with caplog.at_level(logging.WARNING):
        nodes = expand_gather_stage(
            stage=stage,
            benchmark=_fake_benchmark(),
            resolved_modules_cache={("metrics", "summ"): object()},
            output_to_nodes=output_to_nodes,
            nodes_by_id=nodes_by_id,
        )

    assert len(nodes) == 1
    assert nodes[0].gathered_from == ["d1.default-clu-ma.default"]
    assert "orphan.default" in caplog.text and "dropped" in caplog.text


@pytest.mark.short
def test_gather_skips_a_module_missing_from_the_resolution_cache():
    """Same contract as the scatter path: an unresolved module makes no nodes."""
    nodes_by_id = {
        "d1.default": _member("d1.default", None, "data", "d1"),
        "d1.default-clu-ma.default": _member(
            "d1.default-clu-ma.default", "d1.default", "clu", "ma"
        ),
    }
    nodes = expand_gather_stage(
        stage=SimpleNamespace(
            id="metrics",
            gather=[SimpleNamespace(from_="clustering", group_by="data")],
            modules=[
                SimpleNamespace(
                    id="summ",
                    name="summ",
                    parameters=None,
                    provides=None,
                    resources=None,
                )
            ],
            outputs=[SimpleNamespace(id="metrics.summary", path="{data}.tsv")],
            resources=None,
        ),
        benchmark=_fake_benchmark(),
        resolved_modules_cache={},  # nothing resolved
        output_to_nodes={
            "clustering": [("d1.default-clu-ma.default", "d1/clu/ma/a.tsv")]
        },
        nodes_by_id=nodes_by_id,
    )
    assert nodes == []


# ---------------------------------------------------------------------------
# {label.params.key}: parameter values resolved from the lineage
# ---------------------------------------------------------------------------


def _param_ref_yaml(method_params, data_extra=""):
    """data (D1: ideal_components=10, D2: 25) -> method (M1 with `method_params`)."""
    return f"""
id: t
description: t
version: '1.0'
benchmarker: me
api_version: 0.7.0
software_backend: host
software_environments:
  env: {{description: e, easyconfig: e.eb}}
stages:
  - id: data{_line(data_extra, 4)}
    outputs: [{{id: data.out, path: "{{name}}_d.txt"}}]
    modules:
      - id: D1
{_BLOCK_REPO}
        parameters: [{{ideal_components: 10}}]
      - id: D2
{_BLOCK_REPO}
        parameters: [{{ideal_components: 25}}]
  - id: method
    inputs: [data.out]
    outputs: [{{id: method.out, path: "{{name}}_m.txt"}}]
    modules:
      - id: M1
{_BLOCK_REPO}
        parameters: [{method_params}]
"""


def _k_by_root(nodes, key="k"):
    method = [n for n in nodes if n.stage_id == "method"]
    return {n.template_context.provides["dataset"]: n.parameters[key] for n in method}


@pytest.mark.short
def test_param_ref_resolves_per_lineage_keeping_native_type():
    nodes, _, errors = _plan(
        _param_ref_yaml('{k: "{dataset.params.ideal_components}"}')
    )
    assert errors == []
    assert _k_by_root(nodes) == {"D1": 10, "D2": 25}
    method = [n for n in nodes if n.stage_id == "method"]
    assert all(isinstance(n.parameters["k"], int) for n in method)


@pytest.mark.short
def test_param_ref_resolves_a_declared_provides_label():
    nodes, _, errors = _plan(
        _param_ref_yaml(
            '{k: "{treatment.params.ideal_components}"}',
            data_extra="provides: [treatment]",
        )
    )
    assert errors == []
    assert _k_by_root(nodes) == {"D1": 10, "D2": 25}


@pytest.mark.short
def test_param_hash_follows_the_resolved_value():
    """The hash is taken after resolution, so the two nodes never share a
    param directory — which under the `flat` strategy would be one directory."""
    nodes, _, _ = _plan(
        _param_ref_yaml('{k: "{dataset.params.ideal_components}"}'), "flat"
    )
    method = [n for n in nodes if n.stage_id == "method"]
    assert len({n.param_id for n in method}) == 2
    assert len({n.outputs[0] for n in method}) == 2


@pytest.mark.short
def test_literal_params_still_share_one_hash():
    nodes, _, _ = _plan(_param_ref_yaml("{k: 3}"))
    method = [n for n in nodes if n.stage_id == "method"]
    assert len({n.param_id for n in method}) == 1


@pytest.mark.short
@pytest.mark.parametrize(
    "ref, message",
    [
        ("{dataset.params.nonexistent}", "declares no parameter 'nonexistent'"),
        ("{treatment.params.dose}", "Unknown lineage label 'treatment'"),
    ],
)
def test_bad_param_ref_is_a_dag_error(ref, message):
    _, _, errors = _plan(_param_ref_yaml(f'{{k: "{ref}"}}'))
    assert errors and message in errors[0][2]
