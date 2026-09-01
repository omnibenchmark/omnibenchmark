"""Pure lineage helpers over resolved nodes.

Ancestry walks over a node's ``parent_id`` chain and ``parents`` fan-in edges
(see ``model.resolved.ResolvedNode``), the input-selection rules built on them
(design 010 §3.1/§5.2), and the lineage labels a node carries
(``TemplateContext``, design 008). Together they back the stage-expansion loop
in ``cli/run.py``; they live here because they are pure and unit-testable in
isolation.
"""

import hashlib
from collections import deque
from itertools import product
from typing import Optional

from omnibenchmark.logging import logger
from omnibenchmark.model.resolved import TemplateContext


def join_hash(node_ids) -> str:
    """Order-independent digest of a fan-in node's parent *set*, so bundle
    order never leaks into a name or path. Shared by resolver and backend."""
    return hashlib.sha1("|".join(sorted(node_ids)).encode()).hexdigest()[:8]


def iter_ancestors(node, nodes_by_id, through_gather: bool = True):
    """Every distinct ancestor of `node`, nearest-first: the `parent_id` chain
    plus explicit `parents` edges, so a node downstream of a join sees every
    branch (010 §5.2).

    `through_gather=False` is the GATING walk — it stops at a gather, whose
    partition deliberately forgets (010 §3.3), so one excluded member among
    hundreds must not poison everything downstream.
    """
    seen = {node.id}
    queue: deque = deque()

    def _push(n):
        if not through_gather and getattr(n, "is_gather", False):
            return
        parent_id = getattr(n, "parent_id", None)
        if parent_id:
            queue.append(parent_id)
        for pid in getattr(n, "parents", None) or []:
            queue.append(pid)

    _push(node)
    while queue:
        node_id = queue.popleft()
        if node_id in seen:
            continue
        seen.add(node_id)
        ancestor = nodes_by_id.get(node_id)
        if ancestor is None:
            continue
        yield ancestor
        _push(ancestor)


def _maximal_producer_stages(producing_stage_ids: set, resolved_nodes, nodes_by_id):
    """Producing stages that are not ancestors of another producing stage.

    Shadowed producers — those with another producer downstream of them — drop
    out here, which is what reproduces nearest-ancestor-wins (design 010 §3.1,
    "shadowed").
    """
    if len(producing_stage_ids) < 2:
        return set(producing_stage_ids)

    dominated = {
        ancestor.stage_id
        for node in resolved_nodes
        if node.stage_id in producing_stage_ids
        for ancestor in iter_ancestors(node, nodes_by_id)
        if ancestor.stage_id in producing_stage_ids
        and ancestor.stage_id != node.stage_id
    }
    return set(producing_stage_ids) - dominated


def _alternative_producer_stages(
    declared_input_ids: list,
    output_to_nodes: dict,
    nodes_by_id: dict,
    maximal_stages: set,
    deepest_stage_id: str,
):
    """The maximal producers that are *alternatives* to ``deepest_stage_id``.

    Several maximal producers mean one of two things, told apart by whether
    they produce the **same** declared input id:

    * same id — alternatives (010 §3.1): each is its own expansion base, so
      the consumer runs once per producer.
    * different ids — fan-in partners (010 §5.2): the consumer needs both at
      once. Left to ``select_input_bundles``; returning them here would anchor
      the same join once per branch and emit it twice.
    """
    alternatives = {deepest_stage_id}
    for input_id in declared_input_ids:
        stages_for_id = {
            node.stage_id
            for node_id, _path in output_to_nodes.get(input_id, [])
            if (node := nodes_by_id.get(node_id)) is not None
        } & maximal_stages
        if deepest_stage_id in stages_for_id:
            alternatives |= stages_for_id
    return alternatives


def select_input_nodes(
    declared_input_ids: list[str],
    output_to_nodes: dict,
    resolved_nodes: list,
    stage_ids_in_order: list[str],
    previous_stage_nodes: list,
    nodes_by_id: dict | None = None,
) -> list:
    """Return the node list to use as the cartesian expansion base for a stage.

    Per design 010 §3.1: a producer with another producer downstream of it is
    shadowed and drops out, so producers on one chain collapse to the deepest;
    producers on parallel branches declaring the *same* id are alternatives,
    each its own expansion base; parallel producers of *different* ids are a
    fan-in, left to ``select_input_bundles``.

    Ungated on api_version, unlike gather and fan-in: with one producer per id
    the maximal set is a singleton, so existing plans are byte-identical, and a
    gate would force benchmarks pinned at 0.4.0 to migrate modules for a
    resolver fix. To gate later, thread a flag from the caller.
    """
    if not declared_input_ids:
        return previous_stage_nodes

    if nodes_by_id is None:
        nodes_by_id = {n.id: n for n in resolved_nodes}

    providing_stage_id_to_depth: dict[str, int] = {}
    for input_id in declared_input_ids:
        for node_id, _path in output_to_nodes.get(input_id, []):
            node_obj = nodes_by_id.get(
                node_id
            )  # O(1); was an O(N) scan of resolved_nodes
            if node_obj is not None:
                depth = (
                    stage_ids_in_order.index(node_obj.stage_id)
                    if node_obj.stage_id in stage_ids_in_order
                    else -1
                )
                existing = providing_stage_id_to_depth.get(node_obj.stage_id, -1)
                if depth > existing:
                    providing_stage_id_to_depth[node_obj.stage_id] = depth

    if not providing_stage_id_to_depth:
        return previous_stage_nodes

    deepest_stage_id = max(
        providing_stage_id_to_depth,
        key=providing_stage_id_to_depth.__getitem__,
    )

    # The deepest producer is always maximal — a dominator would be one of its
    # descendants, hence a higher topo index — so a singleton maximal set is
    # exactly {deepest} and falls out of the general path below.
    selected = _alternative_producer_stages(
        declared_input_ids,
        output_to_nodes,
        nodes_by_id,
        _maximal_producer_stages(
            set(providing_stage_id_to_depth), resolved_nodes, nodes_by_id
        ),
        deepest_stage_id,
    )
    return [n for n in resolved_nodes if n.stage_id in selected]


def inherited_provides(input_nodes) -> dict:
    """The lineage labels a node inherits from its producer bundle.

    The **union** over every member, not just the anchor: a node downstream of
    a fan-in sees every branch (design 010 §5.2), the same rule `exclude`
    already follows. Branches cannot disagree on a value because a label is
    owned by exactly one stage (008 §3.5, enforced at parse time); the only
    overlap is the builtin `name`, which every node overwrites with its own id.
    """
    provides: dict = {}
    for member in input_nodes:
        if member.template_context is not None:
            provides.update(member.template_context.provides)
    return provides


def satisfies_requires(requires: dict, provides: dict) -> bool:
    """Return True if every `requires` label matches the upstream lineage's.

    `provides` is the upstream label set — for a fan-in, the union over every
    branch (`inherited_provides`), so a gate can name a label carried by any
    of them.
    """
    return all(provides.get(label) == value for label, value in requires.items())


def lineage_module_ids(input_node, nodes_by_id) -> set:
    """Return the module_ids along input_node's GATING lineage, inclusive.

    Walks the explicit ``parent_id`` chain plus join ``parents`` edges (so
    exclusions see every branch of a diamond, #289), but does NOT cross a
    gather cut — the partition forgets (design 010 §3.3), and one excluded
    member among hundreds gathered must not poison every downstream node.
    """
    lineage = {input_node.module_id}
    for ancestor in iter_ancestors(input_node, nodes_by_id, through_gather=False):
        lineage.add(ancestor.module_id)
    return lineage


def ancestor_module_at_stage(member_id, group_stage, nodes_by_id):
    """The module id of the `group_stage` node a member descends from (the
    group value). Parents-aware: sees through joins and prior gathers. Returns
    None when the member has no ancestor in that stage, or when the ancestry
    crosses a fan-in that yields several distinct modules of that stage —
    grouping would be ambiguous."""
    node = nodes_by_id.get(member_id)
    if node is None:
        return None
    found = {
        n.module_id
        for n in (node, *iter_ancestors(node, nodes_by_id))
        if n.stage_id == group_stage
    }
    if len(found) == 1:
        return found.pop()
    return None


def _stage_ancestry(node, nodes_by_id) -> dict:
    """Map stage_id → node ids over `node` and its full (parents-aware)
    ancestry. The per-stage sets are the join-compatibility fingerprint."""
    out: dict = {node.stage_id: {node.id}}
    for ancestor in iter_ancestors(node, nodes_by_id):
        out.setdefault(ancestor.stage_id, set()).add(ancestor.id)
    return out


def _lineages_consistent(a: dict, b: dict) -> bool:
    """Two ancestries may be joined iff they agree wherever they overlap: for
    every stage present in both, they share at least one node. Divergence at
    any common stage (different nodes of the same stage) is a cross-lineage
    pairing and must not be joined — sharing a root is not enough (#289)."""
    for stage_id, ids in a.items():
        other = b.get(stage_id)
        if other is not None and ids.isdisjoint(other):
            return False
    return True


def select_input_bundles(
    declared_input_ids: list,
    output_to_nodes: dict,
    resolved_nodes: list,
    stage_ids_in_order: list,
    previous_stage_nodes: list,
    nodes_by_id: dict,
) -> list:
    """Return input *bundles* — each a tuple of producer nodes to join into one
    downstream node (issue #289).

    Fast path (linear / single producing stage): every bundle is a 1-tuple and
    the result is identical to wrapping `select_input_nodes` — no behaviour
    change for existing benchmarks. Join path (a declared input is produced only
    on a divergent branch): the anchor is paired with every producer of each
    missing input whose ancestry is *consistent* with the anchor's (and with
    the other partners') — they agree at every shared stage, not merely at the
    root — and the cartesian product of those partners yields one bundle per
    emitted node.
    """
    anchors = select_input_nodes(
        declared_input_ids,
        output_to_nodes,
        resolved_nodes,
        stage_ids_in_order,
        previous_stage_nodes,
        nodes_by_id,
    )
    if not declared_input_ids:
        return [(a,) for a in anchors]

    producers: dict = {}
    for input_id in declared_input_ids:
        producers[input_id] = [
            nodes_by_id[node_id]
            for node_id, _path in output_to_nodes.get(input_id, [])
            if node_id in nodes_by_id
        ]

    # Shadowed producers must not become join partners. `select_input_nodes`
    # already drops them when picking anchors; without the same filter here the
    # cartesian product below emits one bundle per producer of a chain that
    # should have collapsed to its deepest (010 §3.1, "shadowed").
    maximal_stages: dict = {
        input_id: _maximal_producer_stages(
            {p.stage_id for p in nodes}, resolved_nodes, nodes_by_id
        )
        for input_id, nodes in producers.items()
    }

    ancestry_cache: dict = {}

    def _ancestry(n):
        cached = ancestry_cache.get(n.id)
        if cached is None:
            cached = _stage_ancestry(n, nodes_by_id)
            ancestry_cache[n.id] = cached
        return cached

    bundles: list = []
    for anchor in anchors:
        anchor_map = _ancestry(anchor)
        anchor_ids = {nid for ids in anchor_map.values() for nid in ids}
        missing = [
            input_id
            for input_id in declared_input_ids
            if not any(p.id in anchor_ids for p in producers.get(input_id, []))
        ]
        if not missing:
            bundles.append((anchor,))  # linear / covered — fast path
            continue
        # Join: for each uncovered input, take lineage-consistent partners.
        partner_lists = []
        for input_id in missing:
            partners = [
                p
                for p in producers[input_id]
                if p.stage_id in maximal_stages[input_id]
                and _lineages_consistent(anchor_map, _ancestry(p))
            ]
            partner_lists.append(partners)
        if not all(partner_lists):
            unmatched = [i for i, pl in zip(missing, partner_lists) if not pl]
            logger.warning(
                f"      No lineage-compatible producer for input(s) {unmatched} "
                f"of {anchor.id}; the input will not resolve."
            )
            bundles.append((anchor,))
            continue
        emitted = 0
        for combo in product(*partner_lists):
            # Partners must also be consistent with EACH OTHER (3-way joins).
            if all(
                _lineages_consistent(_ancestry(x), _ancestry(y))
                for i, x in enumerate(combo)
                for y in combo[i + 1 :]
            ):
                bundles.append((anchor,) + combo)
                emitted += 1
        if not emitted:
            logger.warning(
                f"      No mutually consistent partner combination for "
                f"{anchor.id}; the inputs will not resolve."
            )
            bundles.append((anchor,))
    return bundles


def expansion_segment(param_id: str, members) -> str:
    """The directory segment identifying one expansion of a (stage, module).

    A linear node needs only the parameter hash — its ancestry is already in
    the path prefix. A fan-in node's prefix carries just its deepest input, so
    two joins sharing that branch would collide (Snakemake: AmbiguousRule);
    append the parent-set digest, as the node id already does.
    """
    if len(members) > 1:
        return f"{param_id}-{join_hash(m.id for m in members)}"
    return param_id


def resolve_label_value(label: str, module_provides, module_id: str) -> str:
    """Resolve a single `Stage.provides` label's value for one node.

    Two-level chain, most specific first:
      1. ``module.provides[label]`` — explicit benchmark-local binding.
      2. ``module_id`` — default for the "label = module identity" pattern.

    The parameter-name fallback prototyped earlier is intentionally absent:
    sourcing a label from a same-named CLI parameter couples the module's CLI
    contract to the benchmark's routing across repositories. See
    docs/design/008-filtering.md §3.5.
    """
    if module_provides and label in module_provides:
        return str(module_provides[label])
    return module_id


def build_template_context(
    stage,
    module_id: str,
    module_name: Optional[str] = None,
    input_nodes=(),
    params=None,
    module_provides=None,
    extra_provides=None,
) -> TemplateContext:
    """Build a TemplateContext for a node during expansion.

    `input_nodes` is the producer bundle feeding this node: empty for a root or
    a gather (the cut has no parent), a 1-tuple for a linear node, several for
    a fan-in join. Labels are inherited from all of them; `{module.parent.*}`
    names the anchor, the first, which is also the node's id spine.

    `module_provides` is the optional `Module.provides` dict (label → value);
    it sources the values for the labels this stage advertises via
    `Stage.provides`, defaulting to the module id.

    `extra_provides` is an optional dict layered on top (authoritative) — gather
    nodes use it to bind their group-key label (the `group_by` stage → group
    value) so output templates can reference it (design 010 §3.3).
    """
    anchor = input_nodes[0] if input_nodes else None
    module_attrs: dict[str, str] = {
        "id": module_id,
        "stage": stage.id,
        "name": module_name or module_id,
    }

    stage_provides = getattr(stage, "provides", None)

    # Inherit the upstream lineage labels first, then layer this stage's own.
    provides: dict[str, str] = inherited_provides(input_nodes)

    if stage_provides:
        for label in stage_provides:
            provides[label] = resolve_label_value(label, module_provides, module_id)

    if anchor is not None:
        module_attrs["parent.id"] = anchor.module_id
        module_attrs["parent.stage"] = anchor.stage_id
    else:
        if extra_provides is None:
            # Builtin `dataset` label: root identity, propagated downstream. A
            # gather node (the only caller passing extra_provides) binds exactly
            # its group key; leaking `dataset=<gather module id>` would silently
            # shadow the real dataset downstream (design 010 §3.3).
            if params is not None and "dataset" in params:
                provides.setdefault("dataset", str(params["dataset"]))
            else:
                provides.setdefault("dataset", module_id)

    if extra_provides:
        provides.update(extra_provides)

    # {name} always resolves to the current module's own ID, never inherited
    provides["name"] = module_id

    return TemplateContext(provides=provides, module_attrs=module_attrs)
