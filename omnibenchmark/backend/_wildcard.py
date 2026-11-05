"""One rule per chain of stages below a `foreach` producer (design 013 §4.1).

Below a producer a node's directory is `P/{name}/seg_1/…/seg_k`, one segment
per lineage step (`stage.module.params`); an `over` node's has no `/{name}/`.
Nodes of one (stage, module, params) that sit under the same sequence of
*stages* differ only in the segment values, so they share one rule with a
wildcard per chain position, `{s_<stage>}`, constrained to `<stage>.…`.

An input's file name may depend on which module wrote it (a `{module.id}`
template). The rule then reads it at run time from `_OB_FILES[segment][id]`.
"""

from dataclasses import dataclass, field
import re
from typing import Dict, List, Optional, Tuple

from omnibenchmark.core._paths import sanitize_rule_name
from omnibenchmark.model.resolved import ResolvedNode


@dataclass
class Input:
    dir: str  # directory pattern: chain segments replaced by wildcards
    output_id: str
    seg: Optional[str]  # wildcard of the segment that wrote it, when in the chain
    files: Dict[str, str] = field(default_factory=dict)  # segment -> file name

    @property
    def file(self) -> Optional[str]:
        """The file name when every node agrees on it, else None."""
        names = set(self.files.values())
        return names.pop() if len(names) == 1 else None


@dataclass
class Group:
    rule: str
    producer: ResolvedNode  # the foreach producer
    stages: Tuple[str, ...]  # chain stages, outermost first
    outputs: List[str]  # output patterns
    nodes: list = field(default_factory=list)
    inputs: Dict[str, Input] = field(default_factory=dict)


def wildcard(stage: str) -> str:
    return f"s_{stage}"


def _chain(node, name: str, prod_dir: str) -> List[Tuple[str, str]]:
    """(segment, stage) per lineage step below the producer, the node's own last."""
    parts = node.node_dir[len(prod_dir) + 1 :].split("/")
    if parts[0] == f"{{{name}}}":
        parts = parts[1:]
    return [(seg, seg.split(".", 1)[0]) for seg in parts]


def group_nodes(nodes, producers, foreach_dir) -> Tuple[List[Group], set]:
    """Group the nodes below each producer; returns (groups, ids covered).

    Joins, gathers, collectors and nodes outside a producer are left alone."""
    groups: Dict[tuple, Group] = {}
    covered: set = set()
    for producer in producers:
        name = producer.foreach[0]
        prod_dir = foreach_dir(producer)
        for node in nodes:
            if (
                node.foreach
                or node.is_gather
                or node.is_collector
                or node.parents
                or not (node.node_dir or "").startswith(prod_dir + "/")
            ):
                continue
            *chain, (own, _) = _chain(node, name, prod_dir)
            by_seg = {seg: f"{{{wildcard(stage)}}}" for seg, stage in chain}

            def pattern(path: str) -> str:
                return "/".join(by_seg.get(part, part) for part in path.split("/"))

            stages = tuple(stage for _, stage in chain)
            outputs = [pattern(p) for p in node.outputs.values()]
            key = (producer.id, stages, own, tuple(outputs))
            if key not in groups:
                stem = "-".join((producer.id, *stages, node.stage_id, node.module_id))
                groups[key] = Group(
                    sanitize_rule_name(stem + node.param_id), producer, stages, outputs
                )
            group = groups[key]
            group.nodes.append(node)
            covered.add(node.id)
            for key_, path in node.inputs.items():
                d, file = path.rsplit("/", 1)
                written_by = d.rsplit("/", 1)[-1]
                inp = group.inputs.setdefault(
                    key_,
                    Input(
                        pattern(d),
                        node.input_name_mapping.get(key_, key_),
                        by_seg.get(written_by, "").strip("{}") or None,
                    ),
                )
                inp.files[written_by] = file
    return list(groups.values()), covered


def constraints(groups) -> Dict[str, str]:
    """`wildcard_constraints` for every chain wildcard the groups use."""
    pairs = {wildcard(s): rf"{re.escape(s)}\.[^/]+" for g in groups for s in g.stages}
    return dict(sorted(pairs.items()))


_HELPERS = '''\
def _ob_path(wc, template):
    """Fill the rule's wildcards into *template*; a foreach name is left in."""
    return expand(template, allow_missing=True, **dict(wc.items()))[0]


def _ob_file(wc, seg, output_id):
    """The file the module at chain position *seg* wrote for *output_id*."""
    return _OB_FILES[getattr(wc, seg)][output_id]

'''


def header(groups) -> str:
    """`_OB_FILES` and its helpers; empty when every input has one file name."""
    files: Dict[str, Dict[str, str]] = {}
    for group in groups:
        for inp in group.inputs.values():
            if inp.file is None:
                for seg, name in inp.files.items():
                    files.setdefault(seg, {})[inp.output_id] = name
    if not files:
        return ""
    lines = [
        "# File a module wrote for an output id, by its node segment.",
        "_OB_FILES = {",
    ]
    lines += [f"    {seg!r}: {by_id!r}," for seg, by_id in sorted(files.items())]
    lines.append("}")
    return "\n".join(lines) + "\n\n" + _HELPERS
