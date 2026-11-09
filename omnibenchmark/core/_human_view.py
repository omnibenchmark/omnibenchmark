"""The readable view of a flat (api ≥ 0.8.0) output tree: ``out/human/``.

Mirrors the result tree with each node segment's parameter hash replaced by
the readable parameter name (007 §3.1.2). Derived and disposable: rebuilt from
the node directories on disk after every run, skipped by archives and remote
storage. Reading the disk rather than the resolved nodes keeps the view whole
after a filtered run (``--module``, ``--until``, ``--filter``).
"""

import json
import os
import shutil
from pathlib import Path

from omnibenchmark.core._lineage import GATHERS_DIR, parse_node_segment
from omnibenchmark.core._paths import MAX_FILENAME_LEN, make_human_name
from omnibenchmark.model.params import Params

HUMAN_DIR = "human"


def _readable_segment(node_dir: Path, parsed: tuple) -> str:
    """``stage.module[.group].<name>[-join]``, the name read from parameters.json."""
    stage, module, group, param, join = parsed
    prefix = ".".join(f for f in (stage, module, group) if f) + "."
    suffix = f"-{join}" if join else ""
    if param == "default":
        return f"{prefix}default{suffix}"
    try:
        params = Params(json.loads((node_dir / "parameters.json").read_text()))
    except (OSError, ValueError):
        return node_dir.name  # not written yet: keep the hash
    # The whole segment must fit one path component.
    limit = MAX_FILENAME_LEN - len(prefix) - len(suffix)
    return f"{prefix}{make_human_name(params, limit=limit)}{suffix}"


def _holds_node(directory: Path) -> bool:
    return any(p.is_dir() and parse_node_segment(p.name) for p in directory.rglob("*"))


def _mirror(
    src: Path, dest: Path, top: bool = False, real: "Path | None" = None
) -> None:
    """Link *src*'s entries into *dest*, giving child nodes readable copies.

    A template subdirectory that holds child nodes is mirrored verbatim so its
    nodes get readable copies too; any other subdirectory is linked whole. At
    the top (*src* is the out dir) only node directories and ``gathers/`` are
    mirrored. *real* is where the nodes under ``gathers/<STAGE>/by-<S>/``
    really live (their parameters.json), up to the gather itself.
    """
    dest.mkdir(parents=True, exist_ok=True)
    for entry in sorted(src.iterdir()):
        parsed = parse_node_segment(entry.name) if entry.is_dir() else None
        if parsed:
            node = real / entry.name if real is not None else entry
            if not node.is_dir():  # the gather itself: real from here down
                node = entry
            target = dest / _readable_segment(node, parsed)
            if target.exists():  # two parameter sets sanitise to one name
                target = dest / entry.name
            _mirror(entry, target, real=node if node != entry else None)
        elif top and entry.name != GATHERS_DIR:
            continue
        elif entry.is_dir() and _holds_node(entry):
            by_group = entry.name.startswith("by-") and src.parent.name == GATHERS_DIR
            _mirror(
                entry, dest / entry.name, real=src.parent.parent if by_group else real
            )
        else:
            (dest / entry.name).symlink_to(os.path.relpath(entry, dest))


def build_human_view(out_dir: Path) -> None:
    """Rebuild ``out_dir/human`` from the node directories under *out_dir*."""
    out_dir = Path(out_dir)
    root = out_dir / HUMAN_DIR
    shutil.rmtree(root, ignore_errors=True)
    _mirror(out_dir, root, top=True)
