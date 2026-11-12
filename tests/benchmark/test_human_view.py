"""The `by-params/` readable view of a flat output tree (007 §3.1.2)."""

import json

import pytest

from omnibenchmark.core._human_view import build_human_view
from omnibenchmark.core._lineage import node_depth, parse_node_segment
from omnibenchmark.model.params import Params


def _mkdir(root, real, params=None, files=()):
    d = root / real
    d.mkdir(parents=True, exist_ok=True)
    if params:
        (d / "parameters.json").write_text(json.dumps(Params(params)._params))
    for f in files:
        (d / f).write_text(real)
    return real


@pytest.mark.short
def test_view_mirrors_the_tree_with_readable_segments(tmp_path):
    data = _mkdir(tmp_path, "data.D1.default", files=["d.txt"])
    method = _mkdir(
        tmp_path,
        f"{data}/method.M1.{Params({'k': 5}).hash_short()}",
        {"k": 5},
        ["m.txt"],
    )
    (tmp_path / "Snakefile").write_text("")  # not a node: not mirrored
    (tmp_path / "by-params" / "stale").mkdir(parents=True)

    build_human_view(tmp_path)

    h = tmp_path / "by-params"
    assert sorted(p.name for p in h.iterdir()) == ["data.D1.default"]
    # The child node's hashed directory is not linked; its readable copy is.
    assert sorted(p.name for p in (h / data).iterdir()) == ["d.txt", "method.M1.k-5"]
    linked = h / data / "method.M1.k-5" / "m.txt"
    assert linked.is_symlink() and linked.read_text() == method

    build_human_view(tmp_path)  # idempotent
    assert linked.read_text() == method


@pytest.mark.short
def test_join_digest_survives_and_colliding_names_fall_back(tmp_path):
    data = _mkdir(tmp_path, "data.D1.default")
    _mkdir(tmp_path, f"{data}/join.J.default-aaaaaaaa")
    _mkdir(tmp_path, f"{data}/join.J.default-bbbbbbbb")
    # "a b" and "a_b" sanitise to the same readable name.
    _mkdir(tmp_path, f"{data}/method.M1.11111111", {"x": "a b"})
    _mkdir(tmp_path, f"{data}/method.M1.22222222", {"x": "a_b"})

    build_human_view(tmp_path)

    names = sorted(p.name for p in (tmp_path / "by-params" / data).iterdir())
    assert "join.J.default-aaaaaaaa" in names and "join.J.default-bbbbbbbb" in names
    assert "method.M1.x-a_b" in names
    assert len([n for n in names if n.startswith("method.M1")]) == 2


@pytest.mark.short
def test_template_subdir_holding_a_node_is_mirrored(tmp_path):
    data = _mkdir(tmp_path, "data.D1.default")
    _mkdir(tmp_path, f"{data}/v1.0.0/method.M1.11111111", {"k": 1}, ["m.txt"])
    _mkdir(tmp_path, f"{data}/plots", files=["p.png"])

    build_human_view(tmp_path)

    h = tmp_path / "by-params" / data
    assert (h / "plots").is_symlink()
    assert (h / "v1.0.0" / "method.M1.k-1" / "m.txt").is_symlink()


@pytest.mark.short
@pytest.mark.parametrize(
    "part, parsed",
    [
        ("data.D1.default", ("data", "D1", None, "default", None)),
        ("m.M1.abcdef12", ("m", "M1", None, "abcdef12", None)),
        ("s.c.D1.default", None),  # no group field on 0.8.0
        ("data.D1", ("data", "D1", None, "default", None)),
        ("j.J-9f8e7d6c", ("j", "J", None, "default", "9f8e7d6c")),
        ("j.J.default-9f8e7d6c", ("j", "J", None, "default", "9f8e7d6c")),
        ("v1.0.0", None),  # template subdirectory
        ("a.b.c.d.e", None),
        (".logs", None),
    ],
)
def test_parse_node_segment(part, parsed):
    assert parse_node_segment(part) == parsed


@pytest.mark.short
def test_node_depth_ignores_template_subdirs():
    # 007 §3.1.1: a child nests under its deepest producer, counted in nodes.
    sub = "out/A.M.default/sub/dir/x/file"
    deeper = "out/A.M.default/B.N.default/file"
    assert node_depth(sub) == 1 and node_depth(deeper) == 2
