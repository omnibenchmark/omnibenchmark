"""Pure pruning helpers: which stages and modules survive into a plan.

Host-capability gating (``--with-capability``), the ``--until`` stage cut, and
the metric-collector filter that follows from it. All decide *what to keep*
from a declared benchmark before expansion, and none of them touch the
filesystem or a resolved node, so they live here rather than in the CLI.
"""


def module_capabilities_met(module, available_capabilities) -> bool:
    """True unless the module declares a capability the host did not provide.

    Modules with no `requires_capabilities` always pass. `available_capabilities`
    of None is treated as the empty set (no host facts declared).
    """
    required = getattr(module, "requires_capabilities", None)
    if not required:
        return True
    return set(required).issubset(available_capabilities or set())


def select_capable_modules(modules, module_filter, available_capabilities):
    """Partition modules into (kept, pruned) by the host-capability gate.

    A module is pruned when its `requires_capabilities` are not all provided.
    `-m/--module` (dev mode) bypasses the gate entirely — the user asked for a
    specific module explicitly, so nothing is silently dropped from under them.
    """
    if module_filter:
        return list(modules), []
    kept, pruned = [], []
    for module in modules:
        (
            kept if module_capabilities_met(module, available_capabilities) else pruned
        ).append(module)
    return kept, pruned


def capability_prune_summary(pruned_modules, available_capabilities) -> str:
    """One-line warning naming what the capability gate dropped and how to keep it."""
    ids = sorted({m.id for m in pruned_modules})
    missing = sorted(
        {
            cap
            for m in pruned_modules
            for cap in m.requires_capabilities
            if cap not in (available_capabilities or set())
        }
    )
    flags = " ".join(f"--with-capability {cap}" for cap in missing)
    return (
        f"{len(ids)} module(s) pruned by capability gate: {', '.join(ids)}. "
        f"Rerun with {flags} to include them."
    )


def empty_stage_warning(stage_id: str, combinations_seen: int) -> str:
    """Message for a stage that had modules to expand but produced no nodes.

    `combinations_seen` is how many (input, params) combinations the stage
    actually attempted. Zero means none was ever generated — a declared input
    id matched no upstream output, which is a wiring problem, not a filter one.
    Blaming requires/exclude there sends the user to the wrong config.
    """
    cause = (
        "every module/input combination was pruned by a filter (requires/exclude)"
        if combinations_seen
        else (
            "no upstream output matched its declared inputs (check for a "
            "misspelled output id or a skipped upstream stage)"
        )
    )
    return (
        f"Stage '{stage_id}' produced no nodes: {cause}. Downstream stages "
        f"depending on it will also be empty."
    )


def apply_until_filter(stages, until_stage, parents):
    """Restrict a stage list to `until_stage` plus its transitive ancestors.

    `parents` maps a stage id to the set of its *direct* upstream stage ids
    (built from the model's declared inputs — see `stage_adjacency`). Stages are
    kept by declared lineage, not by their order of appearance in the YAML: a
    benchmark may declare an ancestor of `until_stage` *after* it, and every
    stage downstream of (or unrelated to) `until_stage` is pruned regardless of
    declaration order. Output preserves the original declaration order.

    Returns the original sequence as a list when `until_stage` is None.
    Raises ValueError when the named stage is not present.
    """
    stages = list(stages)
    if until_stage is None:
        return stages
    if not any(s.id == until_stage for s in stages):
        available = ", ".join(s.id for s in stages)
        raise ValueError(
            f"--until: stage '{until_stage}' not found. Available stages: {available}"
        )
    keep = {until_stage}
    queue = [until_stage]
    while queue:
        for up in parents.get(queue.pop(), ()):
            if up not in keep:
                keep.add(up)
                queue.append(up)
    return [s for s in stages if s.id in keep]


def filter_collectors_by_stages(collectors, included_stage_ids, benchmark):
    """Drop metric collectors whose declared inputs reference pruned stages.

    Returns a tuple (kept, dropped_ids). A collector is dropped when one of its
    declared input ids has producing stages but none of them survived (an id
    with no producer at all is left to the regular collector resolver to warn
    about). An id may be a shared contract across stages (design 010 §3.1), and
    `_gather_collector_inputs` collects from every producer, so one surviving
    producer is enough to keep the collector.
    """
    kept = []
    dropped = []
    for c in collectors or []:
        keep = True
        for input_ref in c.inputs:
            input_id = input_ref if isinstance(input_ref, str) else input_ref.id
            producers = benchmark.get_stages_by_output(input_id)
            if producers and not any(s.id in included_stage_ids for s in producers):
                keep = False
                break
        if keep:
            kept.append(c)
        else:
            dropped.append(c.id)
    return kept, dropped
