"""Pure lineage/selection helpers over resolved nodes.

These operate on the ``parent_id`` chain of resolved nodes (see
``model.resolved.ResolvedNode``) and back the stage-expansion loop in
``cli/run.py``. Kept out of the CLI layer because they are pure and unit-tested
in isolation.
"""


def select_input_nodes(
    declared_input_ids: list[str],
    output_to_nodes: dict,
    resolved_nodes: list,
    nodes_by_id: dict,
    stage_ids_in_order: list[str],
    previous_stage_nodes: list,
) -> list:
    """Return the node list to use as the cartesian expansion base for a stage.

    NOTE(#289): this collapses a stage's inputs down to a *single* parent stage
    (the deepest already-expanded producer). It therefore assumes every declared
    input is reachable along one linear ancestry chain from that parent. When a
    stage joins two sibling branches (a "diamond"), the other branch's outputs are
    not in that lineage and later fail to resolve. This is the crux of the
    multi-stage input-collection limitation.

    TODO(#289): once arbitrary multi-stage input collection is supported, this
    single-parent selection (and the linear-ancestor walk ``lineage_ancestors``)
    must be generalised to gather inputs from *all* producing stages, not just the
    deepest one on a single chain.
    """
    if not declared_input_ids:
        return previous_stage_nodes

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
    return [n for n in resolved_nodes if n.stage_id == deepest_stage_id]


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
