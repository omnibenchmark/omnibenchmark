"""Pure lineage/selection helpers over resolved nodes.

These operate on the ``parent_id`` chain of resolved nodes (see
``model.resolved.ResolvedNode``) and back the stage-expansion loop in
``cli/run.py``. Kept out of the CLI layer because they are pure and unit-tested
in isolation.
"""

import hashlib
from collections import deque


def join_hash(node_ids) -> str:
    """Short, order-independent digest of a fan-in node's parent set.

    Keys on the *set*, so bundle order never leaks into a name or a path.
    Shared by the resolver (which has the member nodes) and the backend (which
    has `ResolvedNode.parents`), so both derive the same string from one rule.
    """
    return hashlib.sha1("|".join(sorted(node_ids)).encode()).hexdigest()[:8]


def iter_ancestors(node, nodes_by_id, through_gather: bool = True):
    """Yield every distinct ancestor of `node`, nearest-first, walking the
    `parent_id` chain plus explicit `parents` edges (fan-in branches), so a
    node downstream of a join sees every branch (design 010 §3.9).

    `through_gather=False` is the GATING walk: it does not expand a gather
    node's members — the cut deliberately forgets (design 010 §3.3), so e.g.
    one excluded member among hundreds gathered must not poison every
    downstream node. Resolution/provenance walks use the default and fan
    through the cut.
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

    Shadowed producers — those on the same chain as a deeper one — drop out
    here, which is what reproduces nearest-ancestor-wins (design 010 §3.1
    rule 1a).
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

    * same id — alternatives (010 §3.1 rule 1b): each is its own expansion
      base, so the consumer runs once per producer.
    * different ids — fan-in partners (010 §3.9): the consumer needs both at
      once. Left to ``_select_input_bundles``; returning them here would anchor
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

    ``nodes_by_id`` (id → node) is an optional prebuilt index for the O(1) node
    lookup; it is built from ``resolved_nodes`` when omitted, so callers that
    already keep the index (the run loop) avoid the O(N) rebuild.

    Producers on one chain collapse to the deepest — nearest-ancestor-wins
    (design 010 §3.1 rule 1a). Producers on *parallel* branches that declare
    the same input id are alternatives (rule 1b): each becomes its own
    expansion base, so the consumer runs once per producer. Producers of
    *different* ids on parallel branches are a fan-in and are left to
    ``_select_input_bundles``.

    Deliberately NOT gated on api_version, unlike `gather` and the fan-in join
    (both 0.7.0). It needs no gate — with one producer per id the maximal set
    is a singleton, so existing plans are byte-identical — and a gate would
    hurt: benchmarks pinned to 0.4.0 by the `--name`/`{dataset}` workaround
    (see backend/snakemake.py) would have to migrate their modules to get a
    resolver fix. Gating later is cheap if a reason appears: thread a flag from
    `_expand_scatter_stage`, which already has the `benchmark`.
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

    maximal = _maximal_producer_stages(
        set(providing_stage_id_to_depth), resolved_nodes, nodes_by_id
    )
    if len(maximal) < 2:
        return [n for n in resolved_nodes if n.stage_id == deepest_stage_id]

    selected = _alternative_producer_stages(
        declared_input_ids,
        output_to_nodes,
        nodes_by_id,
        maximal,
        deepest_stage_id,
    )
    return [n for n in resolved_nodes if n.stage_id in selected]


def satisfies_requires(requires: dict, input_node) -> bool:
    """Return True if the input_node's lineage satisfies all requires constraints."""
    if not input_node.template_context:
        return False
    for label, required_value in requires.items():
        actual_value = input_node.template_context.provides.get(label)
        if actual_value != required_value:
            return False
    return True


def lineage_module_ids(input_node, nodes_by_id) -> set:
    """Return the module_ids along input_node's full lineage, inclusive.

    Walks the explicit ``parent_id`` chain rather than matching node-ID string
    prefixes, so it is unaffected by module/stage ids that contain the id
    separator. ``nodes_by_id`` maps node id → node.
    """
    lineage = set()
    node = input_node
    while node is not None:
        lineage.add(node.module_id)
        node = nodes_by_id.get(node.parent_id) if node.parent_id else None
    return lineage


def lineage_ancestors(input_node, nodes_by_id) -> list:
    """Strict ancestors of ``input_node``, shallowest first.

    Node ids are built as ``{parent.id}-...``, so this is exactly the set of
    nodes whose id is a strict prefix of ``input_node.id`` -- i.e. its linear
    lineage. Walks the explicit ``parent_id`` chain (O(depth)), like
    ``lineage_module_ids``.

    The prior implementation found these by scanning all resolved nodes for
    id-prefix matches and RE-RECURSING on every match with no memo -- O(2^depth)
    -- which hung Snakefile generation on deep lineages. Shallowest-first so a
    caller doing ``dict.update`` per ancestor lets the deepest win on any
    output-id collision.
    """
    ancestors = []
    node = nodes_by_id.get(input_node.parent_id) if input_node.parent_id else None
    while node is not None:
        ancestors.append(node)
        node = nodes_by_id.get(node.parent_id) if node.parent_id else None
    ancestors.reverse()
    return ancestors
