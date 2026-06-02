#!/usr/bin/env python3
"""Generate ARCHITECTURE.md from package docstrings + real import edges.

Single source of truth:
- *Semantics* (role, layer) come from each top-level package's ``__init__``
  docstring, which must carry the schema::

      <one-line role>
      Represents: <what this package embodies>
      Layer: <foundation|interface|engine|infrastructure|service|backend|cli>
      Depends on: <free text>

- *Structure* (the dependency edges) is computed from the actual imports via
  AST, collapsed to top-level packages — so the diagram can never drift from
  the code.

A Tarjan strongly-connected-components pass detects import cycles. With
``--check`` the script makes no changes and exits non-zero if the doc is stale
or any cycle exists, which makes it usable as a CI guard.

Usage::

    python scripts/gen_architecture.py            # (re)write ARCHITECTURE.md
    python scripts/gen_architecture.py --check     # verify, don't write
"""

from __future__ import annotations

import argparse
import ast
import re
import sys
from pathlib import Path

PKG = "omnibenchmark"
REPO = Path(__file__).resolve().parent.parent
PKG_DIR = REPO / PKG
OUT = REPO / "ARCHITECTURE.md"

# High (user-facing) → low (foundation). Matches the real dependency direction:
# cli → service/backend → interface(storage) → engine(core) → infrastructure(git) → foundation
LAYER_ORDER = [
    "cli",
    "backend",
    "service",
    "interface",
    "engine",
    "infrastructure",
    "foundation",
]


def _field(doc: str, key: str) -> str | None:
    """Capture a schema field value, joining wrapped lines.

    A value runs from ``Key:`` until the next blank line or the next
    ``Another-Key:`` line, so multi-line descriptions are captured whole.
    """
    lines = doc.splitlines()
    for i, line in enumerate(lines):
        m = re.match(rf"^{key}:\s*(.*)$", line)
        if not m:
            continue
        parts = [m.group(1).strip()]
        for nxt in lines[i + 1 :]:
            if not nxt.strip() or re.match(r"^[A-Z][\w ]*:", nxt):
                break
            parts.append(nxt.strip())
        return re.sub(r"\s+", " ", " ".join(parts)).strip()
    return None


def _first_sentence(text: str) -> str:
    """First sentence (for compact table cells); full text stays in the doc."""
    m = re.search(r"^(.+?[.])(\s|$)", text)
    return (m.group(1) if m else text).strip()


def discover():
    """Return {package: {role, represents, layer}} for documented packages."""
    pkgs = {}
    for init in sorted(PKG_DIR.glob("*/__init__.py")):
        name = init.parent.name
        doc = ast.get_docstring(ast.parse(init.read_text(encoding="utf-8"))) or ""
        layer = _field(doc, "Layer")
        represents = _field(doc, "Represents")
        if not layer or not represents:
            continue  # only packages that opt into the schema
        role = doc.strip().splitlines()[0].strip()
        # normalise "engine (central domain)" -> key "engine"
        layer_key = re.split(r"[ (]", layer)[0].strip().lower()
        pkgs[name] = {
            "role": role,
            "represents": represents,
            "layer": layer,
            "layer_key": layer_key,
        }
    return pkgs


def import_edges(known: set[str]) -> dict[str, set[str]]:
    """Actual import edges between top-level packages (AST, collapsed)."""
    edges: dict[str, set[str]] = {p: set() for p in known}
    for f in PKG_DIR.rglob("*.py"):
        src_pkg = f.relative_to(PKG_DIR).parts[0]
        if src_pkg not in known:
            continue
        try:
            tree = ast.parse(f.read_text(encoding="utf-8"))
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            mod = None
            if isinstance(node, ast.ImportFrom) and node.level == 0:
                mod = node.module
            elif isinstance(node, ast.Import):
                for a in node.names:
                    parts = a.name.split(".")
                    if len(parts) >= 2 and parts[0] == PKG and parts[1] in known:
                        if parts[1] != src_pkg:
                            edges[src_pkg].add(parts[1])
                continue
            if not mod:
                continue
            parts = mod.split(".")
            if len(parts) >= 2 and parts[0] == PKG and parts[1] in known:
                dst = parts[1]
                if dst != src_pkg:
                    edges[src_pkg].add(dst)
    return edges


def find_cycles(edges: dict[str, set[str]]) -> list[list[str]]:
    """Tarjan SCC; return components with more than one node (i.e. cycles)."""
    idx: dict[str, int] = {}
    low: dict[str, int] = {}
    on: dict[str, bool] = {}
    stack: list[str] = []
    counter = [0]
    out: list[list[str]] = []
    sys.setrecursionlimit(10000)

    def strong(v: str):
        idx[v] = low[v] = counter[0]
        counter[0] += 1
        stack.append(v)
        on[v] = True
        for w in edges.get(v, ()):
            if w not in idx:
                strong(w)
                low[v] = min(low[v], low[w])
            elif on.get(w):
                low[v] = min(low[v], idx[w])
        if low[v] == idx[v]:
            comp = []
            while True:
                w = stack.pop()
                on[w] = False
                comp.append(w)
                if w == v:
                    break
            if len(comp) > 1:
                out.append(sorted(comp))

    for v in edges:
        if v not in idx:
            strong(v)
    return out


def render(pkgs: dict, edges: dict[str, set[str]], cycles: list[list[str]]) -> str:
    by_layer: dict[str, list[str]] = {}
    for name, meta in pkgs.items():
        by_layer.setdefault(meta["layer_key"], []).append(name)

    lines: list[str] = []
    lines.append("# Architecture")
    lines.append("")
    lines.append(
        "> **Generated** by `scripts/gen_architecture.py` from package "
        "`__init__` docstrings and actual imports. Do not edit by hand; run "
        "`python scripts/gen_architecture.py` to refresh."
    )
    lines.append("")
    lines.append(
        "Packages are layered top (user-facing) to bottom (foundation). "
        "Imports only ever point downward — the graph is acyclic."
    )
    lines.append("")

    # --- Mermaid diagram (GitHub renders this natively) ---
    lines.append("```mermaid")
    lines.append("flowchart TB")
    for layer in LAYER_ORDER:
        members = by_layer.get(layer, [])
        if not members:
            continue
        lines.append(f"  subgraph layer_{layer}[{layer}]")
        for m in sorted(members):
            lines.append(f"    {m}")
        lines.append("  end")
    lines.append("")
    for src in sorted(edges):
        for dst in sorted(edges[src]):
            lines.append(f"  {src} --> {dst}")
    lines.append("```")
    lines.append("")

    # --- Cycle status ---
    if cycles:
        lines.append("## ⚠️ Import cycles detected")
        lines.append("")
        for c in cycles:
            lines.append(f"- `{' ↔ '.join(c)}`")
        lines.append("")
    else:
        lines.append(
            "**Import cycles:** none — every package is its own strongly-"
            "connected component (clean DAG)."
        )
        lines.append("")

    # --- Package table ---
    lines.append("## Packages")
    lines.append("")
    lines.append("| Package | Layer | Represents | Depends on |")
    lines.append("|---|---|---|---|")
    order = {layer: i for i, layer in enumerate(LAYER_ORDER)}
    for name in sorted(pkgs, key=lambda n: (order.get(pkgs[n]["layer_key"], 99), n)):
        deps = ", ".join(f"`{d}`" for d in sorted(edges.get(name, set()))) or "—"
        rep = _first_sentence(pkgs[name]["represents"]).rstrip(".")
        layer = _first_sentence(pkgs[name]["layer"])
        lines.append(f"| `{name}` | {layer} | {rep} | {deps} |")
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--check",
        action="store_true",
        help="verify ARCHITECTURE.md is current and cycle-free; do not write",
    )
    args = ap.parse_args()

    pkgs = discover()
    edges = import_edges(set(pkgs))
    cycles = find_cycles(edges)
    content = render(pkgs, edges, cycles) + "\n"

    if args.check:
        problems = []
        if cycles:
            problems.append(
                "import cycle(s): " + "; ".join(" <-> ".join(c) for c in cycles)
            )
        if not OUT.exists() or OUT.read_text(encoding="utf-8") != content:
            problems.append(
                "ARCHITECTURE.md is stale — run scripts/gen_architecture.py"
            )
        if problems:
            for p in problems:
                print(f"FAIL: {p}", file=sys.stderr)
            return 1
        print("OK: ARCHITECTURE.md current, no import cycles")
        return 0

    OUT.write_text(content, encoding="utf-8")
    status = "with CYCLES" if cycles else "acyclic"
    print(f"wrote {OUT.relative_to(REPO)} ({len(pkgs)} packages, {status})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
