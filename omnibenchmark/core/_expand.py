"""Stage expansion: turn declared stages into ResolvedNodes.

Two shapes, one per stage. ``expand_scatter_stage`` chains: it pairs each
module/parameter combination with the producer node(s) feeding it, so a stage
fans out over its upstream lineage (and fans *in* on a diamond join, #289).
``expand_gather_stage`` cuts the chain instead: it collects every node
producing a shared output id and partitions them by an ancestor stage, one node
per group (design 010).

Both are pure planning over the model and the nodes resolved so far — no
filesystem, no Snakemake, no CLI. ``cli/run.py`` picks between them per stage
and owns the progress reporting around it.
"""

from itertools import product
from typing import Optional

from omnibenchmark.core._lineage import (
    ancestor_module_at_stage,
    build_template_context,
    inherited_provides,
    expansion_segment,
    iter_ancestors,
    join_hash,
    lineage_module_ids,
    satisfies_requires,
    select_input_bundles,
)
from omnibenchmark.core._paths import is_lineage_excluded, truncate_path_filename
from omnibenchmark.core._prune import empty_stage_warning, select_capable_modules
from omnibenchmark.logging import logger
from omnibenchmark.model.params import Params
from omnibenchmark.model.resolved import ResolvedNode


def expand_gather_stage(
    stage,
    benchmark,
    resolved_modules_cache: dict,
    output_to_nodes: dict,
    nodes_by_id: dict,
    path_exclusions=None,
    module_filter=None,
    target_stage=None,
    available_capabilities=None,
) -> list:
    """Expand a gather stage into ResolvedNodes (design 010 MVP).

    For each `gather` entry, collect every node producing its `from` output id
    and partition them by the ancestor module of the entry's `group_by` stage —
    structural grouping over the lineage chain, no `provides` labels. One node
    is emitted per (group value × module × parameter). The chain is cut
    (`parent_id=None`); outputs are written under `stage.prefix` and registered
    in `output_to_nodes` so downstream stages consume them normally.

    Cross-cutting contracts honoured like the scatter sibling: `module_filter`
    prunes to a single execution path (`ob run -m`), and an exclusion rule
    pairing this stage's module with a member's lineage drops that member from
    that module's gather (transitive exclude at member level).
    """
    nodes: list = []

    # group value -> [(member id, from id, path)]. Insertion order is producer
    # expansion order — deterministic; entries are grouped by from id because
    # the outer loop runs per spec.
    grouped: dict = {}
    # The group-key label: the single group_by stage, or None for a global
    # gather (every producer into one node — the `metric_collector` shape).
    # Stage.validate_gather guarantees one axis, so [0] is authoritative.
    group_label = stage.gather[0].group_by

    for spec in stage.gather:
        producers = output_to_nodes.get(spec.from_, [])
        if not producers:
            raise ValueError(
                f"Stage '{stage.id}' gathers from output id '{spec.from_}', "
                f"which no stage produces (design 010 §3.2)."
            )
        for member_id, path in producers:
            if group_label is None:
                grouped.setdefault(None, []).append((member_id, spec.from_, path))
                continue
            gval = ancestor_module_at_stage(member_id, spec.group_by, nodes_by_id)
            if gval is None:
                logger.warning(
                    f"      Gather '{stage.id}': member {member_id} has no "
                    f"unambiguous ancestor in stage '{spec.group_by}'; dropped "
                    f"from grouping."
                )
                continue
            grouped.setdefault(gval, []).append((member_id, spec.from_, path))

    if not grouped:
        logger.warning(f"Stage '{stage.id}' gathered no members.")

    # Module-filter: same single-execution-path contract as the scatter path.
    if module_filter:
        if target_stage is not None and stage.id == target_stage.id:
            modules_to_expand = [m for m in stage.modules if m.id == module_filter]
        else:
            modules_to_expand = stage.modules[:1]
    else:
        modules_to_expand, _ = select_capable_modules(
            stage.modules, module_filter, available_capabilities
        )

    for module in modules_to_expand:
        module_id = module.id
        cache_key = (stage.id, module_id)
        if cache_key not in resolved_modules_cache:
            logger.warning(f"      Module {module_id} not in cache, skipping")
            continue
        resolved_module = resolved_modules_cache[cache_key]

        if module.parameters:
            params_list = []
            for param in module.parameters:
                params_list.extend(Params.expand_from_parameter(param))
        else:
            params_list = [None]

        combos = [(gval, params) for gval in grouped for params in params_list]
        if module_filter:
            combos = combos[:1]

        for gval, params in combos:
            members = grouped[gval]
            if path_exclusions:
                kept = []
                for member_id, from_id, path in members:
                    member = nodes_by_id.get(member_id)
                    lineage = {module_id}
                    if member is not None:
                        lineage |= lineage_module_ids(member, nodes_by_id)
                    if is_lineage_excluded(lineage, path_exclusions):
                        logger.debug(
                            f"      Gather '{stage.id}': dropping member "
                            f"{member_id} for module {module_id} "
                            f"(exclusion rule)"
                        )
                        continue
                    kept.append((member_id, from_id, path))
            else:
                kept = members
            if not kept:
                where = f"group '{gval}'" if gval is not None else "the gather"
                logger.warning(
                    f"      Gather '{stage.id}': {where} has no members left "
                    f"for module {module_id}; skipping node."
                )
                continue

            param_id = f".{params.hash_short()}" if params else ".default"
            # A global gather has no group value, so neither its id nor its
            # path carries a group segment.
            group_seg = f"-{gval}" if gval is not None else ""
            node_id = f"{stage.id}-{module_id}{group_seg}{param_id}"

            # Enumerated keys + name mapping: same shape metric collectors
            # use, so _write_gather_shell emits `--<from_id> p1 p2 …`.
            inputs: dict = {}
            input_name_mapping: dict = {}
            for idx, (_mid, from_id, path) in enumerate(kept):
                inputs[f"input_{idx}"] = path
                input_name_mapping[f"input_{idx}"] = from_id
            member_ids = list(dict.fromkeys(mid for mid, _f, _p in kept))

            ctx = build_template_context(
                stage=stage,
                module_id=module_id,
                module_name=getattr(module, "name", None),
                params=params,
                module_provides=module.provides,
                # Always a dict for a gather, empty when global: that is what
                # keeps the builtin `dataset` from leaking in (a gather binds
                # its group key and nothing else, design 010 §3.3).
                extra_provides={group_label: gval} if group_label else {},
            )

            outputs = []
            for output_spec in stage.outputs:
                # Only the group label is bound: a template referencing any
                # other lineage label raises in substitute() — the plan-time
                # error design 010 §3.3 wants, never a silent empty sub.
                tmpl = ctx.substitute(output_spec.path, params=params)
                group_dir = f"{gval}/" if gval is not None else ""
                output_path = truncate_path_filename(
                    f"{stage.prefix}/{group_dir}{stage.id}/{module_id}/{param_id}/{tmpl}"
                )
                outputs.append(output_path)
                output_to_nodes.setdefault(output_spec.id, []).append(
                    (node_id, output_path)
                )

            node_resources = None
            if getattr(module, "resources", None):
                node_resources = module.resources
            elif getattr(stage, "resources", None):
                node_resources = stage.resources

            nodes.append(
                ResolvedNode(
                    id=node_id,
                    stage_id=stage.id,
                    module_id=module_id,
                    param_id=param_id,
                    module=resolved_module,
                    parameters=params,
                    parent_id=None,
                    inputs=inputs,
                    outputs=outputs,
                    input_name_mapping=input_name_mapping,
                    benchmark_name=benchmark.model.get_name(),
                    benchmark_version=benchmark.model.get_version(),
                    benchmark_author=benchmark.model.get_author(),
                    resources=node_resources,
                    template_context=ctx,
                    is_gather=True,
                    gathered_from=member_ids,
                    parents=list(member_ids),
                )
            )

    return nodes


def _get_output_ids_for_node(node, benchmark) -> dict:
    """Map each output id declared by `node`'s stage to that node's resolved
    path (index-matched). Used to resolve a downstream stage's declared inputs."""
    stage = benchmark.model.get_stage(node.stage_id)
    if stage is None:
        return {}
    return {spec.id: path for spec, path in zip(stage.outputs, node.outputs)}


def expand_scatter_stage(
    stage,
    benchmark,
    resolved_modules_cache: dict,
    resolved_nodes: list,
    nodes_by_id: dict,
    output_to_nodes: dict,
    previous_stage_nodes: list,
    stages_to_expand: list,
    path_exclusions,
    nesting_strategy: str,
    module_filter,
    target_stage,
    dag_errors: list,
    prune_counts: dict,
    quiet: bool,
    available_capabilities: Optional[set] = None,
) -> list:
    """Expand one ordinary (scatter/chain) stage into ResolvedNodes.

    Appends created nodes to `resolved_nodes`/`nodes_by_id` as it builds (the
    ancestry walk reads them mid-expansion) and registers outputs in
    `output_to_nodes`; returns this stage's nodes. Extracted from the former
    inline expansion loop, plus one addition: inputs come as *bundles*
    (`select_input_bundles`), so a diamond join (#289, api ≥ 0.7.0) expands to
    one node per bundle instead of dropping a branch. The gather sibling is
    `expand_gather_stage`.
    """
    current_stage_nodes: list = []
    # How many (input, params) combinations this stage actually attempted. Zero
    # means nothing was generated to prune (an input id resolved to no upstream
    # node); >0 with an empty stage means a filter rejected them all. The
    # empty-stage warning distinguishes the two.
    combinations_seen = 0

    # Module-filter: which modules to expand per stage
    if module_filter:
        is_target_stage = stage.id == target_stage.id
        if is_target_stage:
            modules_to_expand = [m for m in stage.modules if m.id == module_filter]
        else:
            modules_to_expand = stage.modules[:1]
    else:
        modules_to_expand, _ = select_capable_modules(
            stage.modules, module_filter, available_capabilities
        )

    # Input bundles depend only on stage-level state — compute once, not per
    # module. None (vs []) distinguishes "no inputs declared" from "inputs
    # declared but nothing resolvable".
    input_bundles = None
    if stage.inputs and previous_stage_nodes:
        declared_input_ids = [
            entry
            for input_col in stage.inputs
            if hasattr(input_col, "entries")
            for entry in input_col.entries
        ]
        input_bundles = select_input_bundles(
            declared_input_ids=declared_input_ids,
            output_to_nodes=output_to_nodes,
            resolved_nodes=resolved_nodes,
            stage_ids_in_order=[s.id for s in stages_to_expand],
            previous_stage_nodes=previous_stage_nodes,
            nodes_by_id=nodes_by_id,
        )

    # Input resolution depends only on the bundle (never on module/params):
    # cache (inputs, name mapping, base_path) per member-id tuple.
    resolution_cache: dict = {}

    for module in modules_to_expand:
        module_id = module.id
        cache_key = (stage.id, module_id)

        try:
            if cache_key not in resolved_modules_cache:
                logger.warning(f"      Module {module_id} not in cache, skipping")
                continue

            resolved_module = resolved_modules_cache[cache_key]

            if module.parameters:
                params_list = []
                for param in module.parameters:
                    params_list.extend(Params.expand_from_parameter(param))
            else:
                params_list = [None]

            if input_bundles is not None:
                node_combinations = list(product(input_bundles, params_list))
            else:
                node_combinations = [(None, params) for params in params_list]

            if module_filter:
                node_combinations = node_combinations[:1]

            combinations_seen += len(node_combinations)

            for input_bundle, params in node_combinations:
                # A bundle is the producer node(s) feeding this node: a
                # 1-tuple for a linear node, several for a diamond join
                # (#289). `input_node` is the anchor — the primary (deepest)
                # branch — which gives the id spine and `{module.parent.*}`.
                # Labels come from the whole bundle, not just the anchor.
                members = input_bundle if input_bundle else ()
                input_node = members[0] if members else None
                is_join = len(members) > 1

                # Exclusions are transitive over the full lineage, not just
                # the immediate predecessor: prune if any exclusion rule has
                # both endpoints present along this node's lineage. For a join
                # the lineage is the union over every branch.
                if members:
                    lineage = {module_id}
                    for member in members:
                        lineage |= lineage_module_ids(member, nodes_by_id)
                    if is_lineage_excluded(lineage, path_exclusions):
                        prune_counts["exclude"] += 1
                        logger.debug(
                            f"      Excluding combination: lineage {lineage} "
                            f"violates an exclusion rule"
                        )
                        continue

                if input_node and module.requires:
                    # Gate against every branch's labels, the same lineage the
                    # exclusion check above unions over (010 §5.2).
                    upstream = inherited_provides(members)
                    if not satisfies_requires(module.requires, upstream):
                        prune_counts["requires"] += 1
                        logger.debug(
                            f"      Skipping combination: requires not satisfied for {module_id} "
                            f"(upstream context: {upstream})"
                        )
                        continue

                param_id = f".{params.hash_short()}" if params else ".default"

                if is_join:
                    # No single prefix chain: id is a readable stem plus a
                    # short hash of the (sorted) parents (design 010 §5.2).
                    node_id = f"{stage.id}-{module_id}-{join_hash(m.id for m in members)}{param_id}"
                elif input_node:
                    node_id = f"{input_node.id}-{stage.id}-{module_id}{param_id}"
                else:
                    node_id = f"{stage.id}-{module_id}{param_id}"

                inputs = {}
                input_name_mapping = {}
                base_path = None

                if input_node:
                    bundle_key = tuple(m.id for m in members)
                    cached = resolution_cache.get(bundle_key)
                    if cached is None:
                        output_id_to_path = {}

                        # Collect resolvable outputs from every branch's
                        # producer node and its ancestors (union across the
                        # join). Farthest ancestors first so the NEAREST
                        # producer of a shared output id wins (design 010
                        # §3.1).
                        for member in members:
                            for ancestor in reversed(
                                list(iter_ancestors(member, nodes_by_id))
                            ):
                                output_id_to_path.update(
                                    _get_output_ids_for_node(ancestor, benchmark)
                                )
                            output_id_to_path.update(
                                _get_output_ids_for_node(member, benchmark)
                            )

                        if stage.inputs:
                            for input_collection in stage.inputs:
                                if not hasattr(input_collection, "entries"):
                                    continue
                                for input_id in input_collection.entries:
                                    if input_id in output_id_to_path:
                                        sanitized_id = input_id.replace(".", "_")
                                        inputs[sanitized_id] = output_id_to_path[
                                            input_id
                                        ]
                                        input_name_mapping[sanitized_id] = input_id
                                    else:
                                        logger.warning(
                                            f"      Could not resolve input {input_id} for node {node_id}"
                                        )

                        if inputs:
                            from pathlib import Path as PathLib

                            deepest_input = max(
                                inputs.values(), key=lambda p: len(PathLib(p).parts)
                            )
                            base_path = str(PathLib(deepest_input).parent)
                        resolution_cache[bundle_key] = (
                            inputs,
                            input_name_mapping,
                            base_path,
                        )
                    else:
                        inputs, input_name_mapping, base_path = cached
                    # Fresh dicts per node — cached ones must stay pristine.
                    inputs = dict(inputs)
                    input_name_mapping = dict(input_name_mapping)

                ctx = build_template_context(
                    stage=stage,
                    module_id=module_id,
                    module_name=getattr(module, "name", None),
                    input_nodes=members,
                    params=params,
                    module_provides=module.provides,
                )

                outputs = []
                for output_spec in stage.outputs:
                    output_path_template = ctx.substitute(
                        output_spec.path, params=params
                    )

                    seg = expansion_segment(param_id, members)
                    if nesting_strategy == "nested":
                        if base_path:
                            output_path = f"{base_path}/{stage.id}/{module_id}/{seg}/{output_path_template}"
                        else:
                            output_path = (
                                f"{stage.id}/{module_id}/{seg}/{output_path_template}"
                            )
                    elif nesting_strategy == "flat":
                        output_path = (
                            f"{stage.id}/{module_id}/{seg}/{output_path_template}"
                        )
                    else:
                        raise ValueError(
                            f"Unknown nesting strategy: {nesting_strategy}"
                        )

                    output_path = truncate_path_filename(output_path)

                    outputs.append(output_path)
                    if output_spec.id not in output_to_nodes:
                        output_to_nodes[output_spec.id] = []
                    output_to_nodes[output_spec.id].append((node_id, output_path))

                node_resources = None
                if hasattr(module, "resources") and module.resources:
                    node_resources = module.resources
                elif hasattr(stage, "resources") and stage.resources:
                    node_resources = stage.resources

                node = ResolvedNode(
                    id=node_id,
                    stage_id=stage.id,
                    module_id=module_id,
                    param_id=param_id,
                    module=resolved_module,
                    parameters=params,
                    parent_id=input_node.id if input_node else None,
                    parents=[m.id for m in members] if is_join else [],
                    inputs=inputs,
                    outputs=outputs,
                    input_name_mapping=input_name_mapping,
                    benchmark_name=benchmark.model.get_name(),
                    benchmark_version=benchmark.model.get_version(),
                    benchmark_author=benchmark.model.get_author(),
                    resources=node_resources,
                    template_context=ctx,
                )

                resolved_nodes.append(node)
                nodes_by_id[node.id] = node
                current_stage_nodes.append(node)

            if not quiet:
                logger.info(
                    f"      Created {len([n for n in current_stage_nodes if n.module_id == module_id])} nodes for {module_id}"
                )

        except ValueError as e:
            msg = str(e)
            logger.error(f"      Failed to resolve {module_id}: {msg}")
            dag_errors.append((stage.id, module_id, msg))

        except Exception as e:
            logger.error(f"      Failed to resolve {module_id}: {e}")
            import traceback

            if logger.level <= 10:
                traceback.print_exc()

    # A stage that had modules to expand but produced no nodes should not
    # cascade an empty set downstream silently.
    if modules_to_expand and not current_stage_nodes:
        logger.warning(empty_stage_warning(stage.id, combinations_seen))

    return current_stage_nodes
