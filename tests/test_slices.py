import pytest

from omnibenchmark.slices import select, slice_hash, slice_label


def _meta():
    def n(id, stage, dataset, parents, key=None):
        return dict(
            id=id, stage=stage, dir=id, key=key or id, dataset=dataset, parents=parents
        )

    return {
        "stages": ["fetch", "data", "method", "summary"],
        "nodes": [
            n("f", "fetch", None, []),
            n("d1", "data", "D1", ["f"]),
            n("d2", "data", "D2", ["f"]),
            n("m1", "method", "D1", ["d1"]),
            n("m2", "method", "D2", ["d2"]),
            n("s", "summary", None, ["m1", "m2"]),
        ],
    }


@pytest.mark.short
def test_select_is_ancestor_closed():
    ids = lambda nodes: sorted(n["id"] for n in nodes)  # noqa: E731
    meta = _meta()
    assert ids(select(meta, "data", None)) == ["d1", "d2", "f"]
    # the trunk comes along; the cross-dataset gather does not
    assert ids(select(meta, None, "D1")) == ["d1", "f", "m1"]
    assert ids(select(meta, "data", "D1")) == ["d1", "f"]
    assert len(select(meta, None, None)) == 6
    assert slice_label("data", "D1") == "D1.until-data"


@pytest.mark.short
def test_slice_hash_follows_keys_not_order():
    meta = _meta()
    cut = slice_hash(select(meta, "data", None))
    meta["nodes"].reverse()
    assert slice_hash(select(meta, "data", None)) == cut
    meta["nodes"][0]["key"] = "changed"  # the summary node: outside the cut
    assert slice_hash(select(meta, "data", None)) == cut
    next(n for n in meta["nodes"] if n["id"] == "f")["key"] = "changed"
    assert slice_hash(select(meta, "data", None)) != cut
