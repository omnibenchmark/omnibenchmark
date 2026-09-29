"""The readable view of a flat (api ≥ 0.8.0) output tree: ``out/human/``.

Mirrors the result tree with each node segment's parameter hash replaced by
the readable parameter name (007 §3.1.2). Derived and disposable: rebuilt from
the resolved nodes after every run, skipped by archives and remote storage.
"""

import os
import shutil
from pathlib import Path

from omnibenchmark.core._paths import MAX_FILENAME_LEN, make_human_name

HUMAN_DIR = "human"


def _readable_segment(segment: str, node) -> str:
    """``stage.module[.group].<hash>[-join]`` with the hash made readable."""
    *head, last = segment.split(".")
    _param, sep, join = last.partition("-")
    prefix, suffix = "".join(f"{h}." for h in head), f"{sep}{join}"
    name = "default"
    if node.parameters:  # the whole segment must fit one path component
        name = make_human_name(
            node.parameters, limit=MAX_FILENAME_LEN - len(prefix) - len(suffix)
        )
    return f"{prefix}{name}{suffix}"


def build_human_view(out_dir: Path, nodes) -> None:
    """Rebuild ``out_dir/human`` from the nodes whose directories exist.

    A node's directory holds its files and its children's directories. The
    readable copy links the files and nests the children's readable copies.
    """
    root = Path(out_dir) / HUMAN_DIR
    shutil.rmtree(root, ignore_errors=True)

    real_dirs = {n.param_dir_template: n for n in nodes if n.param_dir_template}
    # Every directory on the way to a node, so a template subdirectory that
    # holds child nodes is not linked as a whole (the children get readable
    # copies instead). Files beside child nodes in such a subdirectory are
    # therefore not linked; no plan puts files there today.
    on_path = {
        "/".join(p.split("/")[:i])
        for p in real_dirs
        for i in range(1, p.count("/") + 2)
    }
    human_of: dict = {}
    taken: set = set()

    for real in sorted(real_dirs, key=lambda p: p.count("/")):
        src = Path(out_dir) / real
        if not src.is_dir():
            continue  # not produced yet, so neither are its children

        # Nearest ancestor node with a readable copy; template subdirs between
        # it and this node are kept verbatim.
        parent, _, segment = real.rpartition("/")
        head, rest = parent, []
        while head and head not in human_of:
            head, _, tail = head.rpartition("/")
            rest.insert(0, tail)
        base = human_of.get(head, root).joinpath(*rest)

        dest = base / _readable_segment(segment, real_dirs[real])
        if dest in taken:  # two parameter sets sanitise to one name
            dest = base / segment
        taken.add(dest)
        human_of[real] = dest

        dest.mkdir(parents=True, exist_ok=True)
        for entry in src.iterdir():
            if f"{real}/{entry.name}" in on_path:
                continue
            (dest / entry.name).symlink_to(os.path.relpath(entry, dest))
