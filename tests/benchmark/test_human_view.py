"""The `human/` readable view of a flat output tree (007 §3.1.2)."""

from types import SimpleNamespace

import pytest

from omnibenchmark.core._human_view import build_human_view
from omnibenchmark.model.params import Params


def _node(real, params=None):
    return SimpleNamespace(node_dir=real, parameters=Params(params) if params else None)


@pytest.mark.short
def test_view_mirrors_the_tree_with_readable_segments(tmp_path):
    data = "data.D1.default"
    method = f"{data}/method.M1.{Params({'k': 5}).hash_short()}"
    unrun = f"{data}/method.M1.{Params({'k': 9}).hash_short()}"
    for d, f in ((data, "d.txt"), (method, "m.txt")):
        (tmp_path / d).mkdir(parents=True)
        (tmp_path / d / f).write_text(d)
    (tmp_path / "human" / "stale").mkdir(parents=True)

    nodes = [_node(data), _node(method, {"k": 5}), _node(unrun, {"k": 9})]
    build_human_view(tmp_path, nodes)

    h = tmp_path / "human"
    assert sorted(p.name for p in h.iterdir()) == ["data.D1.default"]
    # The child node's hashed directory is not linked; its readable copy is.
    assert sorted(p.name for p in (h / data).iterdir()) == ["d.txt", "method.M1.k-5"]
    linked = h / data / "method.M1.k-5" / "m.txt"
    assert linked.is_symlink() and linked.read_text() == method

    build_human_view(tmp_path, nodes)  # idempotent
    assert linked.read_text() == method


@pytest.mark.short
def test_join_digest_survives_and_colliding_names_fall_back(tmp_path):
    j1 = "data.D1.default/join.J.default-aaaaaaaa"
    j2 = "data.D1.default/join.J.default-bbbbbbbb"
    # "a b" and "a_b" sanitise to the same readable name.
    p1 = "data.D1.default/method.M1.11111111"
    p2 = "data.D1.default/method.M1.22222222"
    for d in ("data.D1.default", j1, j2, p1, p2):
        (tmp_path / d).mkdir(parents=True, exist_ok=True)
    nodes = [
        _node("data.D1.default"),
        _node(j1),
        _node(j2),
        _node(p1, {"x": "a b"}),
        _node(p2, {"x": "a_b"}),
    ]
    build_human_view(tmp_path, nodes)
    names = sorted(p.name for p in (tmp_path / "human" / "data.D1.default").iterdir())
    assert "join.J.default-aaaaaaaa" in names and "join.J.default-bbbbbbbb" in names
    assert "method.M1.x-a_b" in names
    assert len([n for n in names if n.startswith("method.M1")]) == 2
