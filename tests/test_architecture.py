"""Architecture fitness checks.

Guards the package layering established by the cycle-breaking refactor:
- the import graph must stay acyclic (no package-level import cycles)
- ARCHITECTURE.md must stay in sync with the code

Both are derived from the single source of truth (package docstrings + actual
imports) by ``scripts/gen_architecture.py``. A failure here breaks the build,
which is the point: it stops an import cycle (or a stale diagram) from landing.
"""

import importlib.util
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
GEN = REPO / "scripts" / "gen_architecture.py"


def _load_generator():
    spec = importlib.util.spec_from_file_location("gen_architecture", GEN)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.short
def test_no_package_import_cycles():
    ga = _load_generator()
    pkgs = ga.discover()
    edges = ga.import_edges(set(pkgs))
    cycles = ga.find_cycles(edges)
    assert not cycles, "package import cycle(s) detected: " + "; ".join(
        " <-> ".join(c) for c in cycles
    )


@pytest.mark.short
def test_architecture_doc_is_current():
    ga = _load_generator()
    pkgs = ga.discover()
    edges = ga.import_edges(set(pkgs))
    expected = ga.render(pkgs, edges, ga.find_cycles(edges)) + "\n"
    assert ga.OUT.exists(), "ARCHITECTURE.md missing — run scripts/gen_architecture.py"
    assert (
        ga.OUT.read_text(encoding="utf-8") == expected
    ), "ARCHITECTURE.md is stale — run `python scripts/gen_architecture.py`"
