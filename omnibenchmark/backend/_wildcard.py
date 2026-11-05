"""Wildcard rules for the nodes below a `foreach` producer (spike for 011).

Below a producer every node's directory is `P/{name}/seg_1/…/seg_k`, one
segment per lineage step; an `over` node drops the `/{name}/`. Nodes of one
(stage, module, params) that sit under the same chain of *stages* differ only
in the chain's segment values, so they lower to one rule with a wildcard per
chain position: `{s_<stage>}`, constrained to `<stage>\\.[^/]+`.

Inputs whose path is the same pattern across the group are written as that
pattern. When the file name depends on an ancestor module (a `{module.id}`
template) the paths differ: they go into a table keyed by the chain values
and the rule reads it through `_ob_in`.

Everything else (producers, joins, gathers, nodes outside a producer) keeps
one explicit rule per node, and `rule all` keeps every node's targets.
"""

from dataclasses import dataclass, field
import re
from typing import Dict, List, Optional, Tuple

from omnibenchmark.core._paths import sanitize_rule_name


@dataclass
class Group:
    rule: str
    kind: str  # "wild": carries {name}; "over": closes it
    producer_id: str
    name: str  # the foreach capture
    nodes: list
    chain: Tuple[str, ...]  # stage ids of the chain positions
    outputs: List[str]  # patterns, declaration order
    # key -> pattern (uniform across the group) or None (read the table)
    inputs: Dict[str, Optional[str]] = field(default_factory=dict)
    # key -> {chain values: that node's literal input path}
    table: Dict[str, Dict[tuple, str]] = field(default_factory=dict)

    @property
    def wildcards(self) -> List[str]:
        return [f"s_{s}" for s in self.chain]


def _chain(node, prod_dir: str, name: str) -> Tuple[Tuple[str, ...], str, str]:
    """(chain segments, own segment, kind) of a node under a producer."""
    parts = node.node_dir[len(prod_dir) + 1 :].split("/")
    if parts[0] == f"{{{name}}}":
        return tuple(parts[1:-1]), parts[-1], "wild"
    return tuple(parts[:-1]), parts[-1], "over"


def _patternize(path: str, seg_to_wc: Dict[str, str]) -> str:
    for seg, wc in seg_to_wc.items():
        path = path.replace(f"/{seg}/", f"/{{{wc}}}/")
    return path


def group_subtree(nodes, producers, foreach_dir) -> Tuple[List[Group], set]:
    """Group the groupable nodes below each producer. Returns the groups and
    the ids of every node they cover."""
    by_key: Dict[tuple, Group] = {}
    covered: set = set()
    for producer in producers:
        name, _ = producer.foreach
        prod_dir = foreach_dir(producer)
        for node in nodes:
            if (
                node is producer
                or node.foreach
                or node.is_gather
                or node.is_collector
                or node.parents
                or not node.node_dir.startswith(prod_dir + "/")
            ):
                continue
            chain, own, kind = _chain(node, prod_dir, name)
            stages = tuple(seg.split(".", 1)[0] for seg in chain)
            seg_to_wc = {seg: f"s_{st}" for seg, st in zip(chain, stages)}
            outputs = [_patternize(p, seg_to_wc) for p in node.outputs.values()]
            key = (kind, producer.id, stages, own, tuple(outputs))
            group = by_key.get(key)
            if group is None:
                stem = "-".join((producer.id, *stages, node.stage_id, node.module_id))
                group = Group(
                    rule=sanitize_rule_name(stem + node.param_id),
                    kind=kind,
                    producer_id=producer.id,
                    name=name,
                    nodes=[],
                    chain=stages,
                    outputs=outputs,
                )
                by_key[key] = group
            group.nodes.append(node)
            covered.add(node.id)

    rules: Dict[str, int] = {}
    for group in by_key.values():
        # Same stem, different file names: keep both rules apart.
        n = rules.get(group.rule, 0)
        rules[group.rule] = n + 1
        if n:
            group.rule = f"{group.rule}_{n}"
        rep = group.nodes[0]
        for key in rep.inputs:
            patterns = {}
            for node in group.nodes:
                chain, _own, _kind = _chain(
                    node, foreach_dir(_by_id(producers, group.producer_id)), group.name
                )
                seg_to_wc = dict(zip(chain, group.wildcards))
                patterns[chain] = (
                    _patternize(node.inputs[key], seg_to_wc),
                    node.inputs[key],
                )
            distinct = {p for p, _ in patterns.values()}
            if len(distinct) == 1:
                group.inputs[key] = distinct.pop()
            else:
                group.inputs[key] = None
                group.table[key] = {c: lit for c, (_p, lit) in patterns.items()}
    return list(by_key.values()), covered


def _by_id(producers, node_id):
    return next(p for p in producers if p.id == node_id)


def constraints(groups) -> Dict[str, str]:
    """Global `wildcard_constraints` for every chain wildcard the groups use."""
    out: Dict[str, str] = {}
    for g in groups:
        for st in g.chain:
            out[f"s_{st}"] = rf"{re.escape(st)}\.[^/]+"
    return dict(sorted(out.items()))


def input_expr(
    group: Group, key: str, producer_rule: str, producer_pattern: str
) -> str:
    """The Snakefile text for one input of a group rule."""
    pattern = group.inputs[key]
    segs = ", ".join(f'"{w}"' for w in group.wildcards) + (
        "," if len(group.chain) == 1 else ""
    )
    lookup = f'_ob_in(wc, "{group.rule}", "{key}", ({segs}))'
    name = group.name
    if group.kind == "wild" or (pattern is not None and f"{{{name}}}" not in pattern):
        if pattern is not None:
            return f'"{pattern}"'
        return f"lambda wc: {lookup}"
    # An over input: expand the captured values into the chain the rule
    # matched. The template is an f-string over the chain wildcards, or the
    # table entry, which already names this chain's files.
    if pattern is not None:
        template = pattern.replace(f"{{{name}}}", f"{{{{{name}}}}}")
        for w in group.wildcards:
            template = template.replace(f"{{{w}}}", f"{{wc.{w}}}")
        template = f'f"{template}"'
    else:
        template = lookup
    return (
        f"lambda wc: _ob_foreach(checkpoints.{producer_rule}, "
        f'"{producer_pattern}", "{name}", {template})'
    )


_IN_HELPER = """\
def _ob_in(wc, rule, key, segs):
    path = _OB_IN[rule][key][tuple(getattr(wc, s) for s in segs)]
    for n in _OB_NAMES:
        if hasattr(wc, n):
            path = path.replace("{%s}" % n, getattr(wc, n))
    return path

"""


def table_text(groups, names) -> str:
    """`_OB_NAMES`, `_OB_IN` and the `_ob_in` helper; empty when no group
    needs a table."""
    tables = {g.rule: g.table for g in groups if g.table}
    if not tables:
        return ""
    lines = [f"_OB_NAMES = {sorted(names)!r}", "_OB_IN = {"]
    for rule, by_key in tables.items():
        lines.append(f'    "{rule}": {{')
        for key, entries in by_key.items():
            lines.append(f'        "{key}": {{')
            for chain, path in sorted(entries.items()):
                lines.append(f'            {chain!r}: "{path}",')
            lines.append("        },")
        lines.append("    },")
    lines.append("}")
    return "\n".join(lines) + "\n\n" + _IN_HELPER
