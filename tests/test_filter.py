"""Unit tests for the obfilter spike (omnibenchmark/filter.py)."""

import pytest

from omnibenchmark import filter as f
from omnibenchmark.model.params import Params


# --------------------------------------------------------------------------- codec


@pytest.mark.short
class TestBlobCodec:
    def test_round_trip(self):
        picks = {"data": {"*": "all"}, "methods": {"M1": "all", "M2": ["abc12345"]}}
        parent = {"sha256": "de" * 32, "url": "https://x/b"}
        blob = f.unpack_blob(f.pack_blob(picks, parent))
        assert blob["v"] == 3
        assert blob["picks"] == picks
        assert blob["parent"] == parent

    def test_pack_is_deterministic(self):
        picks = {"s": {"m": "all"}}
        parent = {"sha256": "00" * 32}
        assert f.pack_blob(picks, parent) == f.pack_blob(picks, parent)

    def test_unpack_rejects_garbage(self):
        with pytest.raises(f.FilterError):
            f.unpack_blob("not-a-real-blob!!")

    def test_unpack_rejects_decompression_bomb(self):
        import base64
        import gzip

        import tracemalloc

        inflated = 256 * 1024 * 1024  # 256 MB once decompressed
        bomb = gzip.compress(b"\0" * inflated, mtime=0)
        packed = base64.urlsafe_b64encode(bomb).decode().rstrip("=")
        assert len(packed) < 1_000_000  # small on the wire, huge once inflated

        tracemalloc.start()
        try:
            with pytest.raises(f.FilterError, match="more than"):
                f.unpack_blob(packed)
            _, peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()
        # Bounded: the ceiling, not the 256 MB the blob asked for.
        assert peak < 3 * f.MAX_BLOB_BYTES, peak

    def test_unpack_accepts_a_blob_at_the_ceiling(self):
        # A legitimately large-but-bounded blob still round-trips.
        picks = {f"stage{i}": {"*": "all"} for i in range(2000)}
        assert f.unpack_blob(f.pack_blob(picks, {}))["picks"] == picks

    def test_unpack_rejects_wrong_version(self):
        import base64
        import gzip
        import json

        raw = json.dumps({"v": 2, "parent": {}, "picks": {}}).encode()
        packed = (
            base64.urlsafe_b64encode(gzip.compress(raw, mtime=0)).decode().rstrip("=")
        )
        with pytest.raises(f.FilterError, match="version"):
            f.unpack_blob(packed)


# --------------------------------------------------------------------------- schema


@pytest.mark.short
class TestLoadFilter:
    PICKS = {"data": {"*": "all"}, "methods": {"M1": "first"}}

    def test_yaml_needs_only_picks(self):
        text = "picks:\n  data: {'*': all}\n  methods:\n    M1: first\n"
        assert f.load_filter(text) == {"v": 3, "picks": self.PICKS}

    def test_json(self):
        import json

        text = json.dumps({"picks": self.PICKS})
        assert f.load_filter(text)["picks"] == self.PICKS

    def test_packed_blob(self):
        packed = f.pack_blob(self.PICKS, {"sha256": "00" * 32})
        assert f.load_filter(packed) == f.unpack_blob(packed)

    def test_yaml_wrong_version(self):
        with pytest.raises(f.FilterError, match="version"):
            f.load_filter("v: 2\npicks:\n  data: {'*': all}\n")

    def test_yaml_bad_picks(self):
        with pytest.raises(f.FilterError, match="spec must be"):
            f.load_filter("picks:\n  methods:\n    M1: some\n")


@pytest.mark.short
class TestValidatePicks:
    def test_accepts_specs(self):
        f.validate_picks({"s": {"m1": "all", "m2": "first", "m3": ["a", "b"]}})

    def test_wildcard_alone_ok(self):
        f.validate_picks({"s": {"*": "all"}})

    def test_wildcard_mixed_raises(self):
        with pytest.raises(f.FilterError, match="wildcard cannot be mixed"):
            f.validate_picks({"s": {"*": "all", "m1": "first"}})

    def test_bad_spec_raises(self):
        with pytest.raises(f.FilterError, match="spec must be"):
            f.validate_picks({"s": {"m": "second"}})

    def test_non_mapping_raises(self):
        with pytest.raises(f.FilterError):
            f.validate_picks({"s": ["m1", "m2"]})


# --------------------------------------------------------------------------- apply


@pytest.mark.short
class TestModuleSpec:
    PICKS = {"data": {"*": "first"}, "methods": {"M1": "all"}}

    def test_explicit_module(self):
        assert f.module_spec(self.PICKS, "methods", "M1") == "all"

    def test_wildcard_matches_any(self):
        assert f.module_spec(self.PICKS, "data", "whatever") == "first"

    def test_unpicked_module_drops(self):
        assert f.module_spec(self.PICKS, "methods", "M2") is f.DROP

    def test_unpicked_stage_drops(self):
        assert f.module_spec(self.PICKS, "ghost", "X") is f.DROP

    def test_keeps_module(self):
        assert f.keeps_module(self.PICKS, "methods", "M1") is True
        assert f.keeps_module(self.PICKS, "methods", "M2") is False


@pytest.mark.short
class TestFilterParams:
    def _params(self):
        return [Params({"n": str(i)}) for i in range(3)]

    def test_all_unchanged(self):
        ps = self._params()
        assert f.filter_params(ps, "all") == ps

    def test_no_params_unchanged(self):
        assert f.filter_params([None], "first") == [None]

    def test_first_keeps_one(self):
        ps = self._params()
        assert f.filter_params(ps, "first") == ps[:1]

    def test_hash_list_keeps_matching(self):
        ps = self._params()
        wanted = ps[2].hash_short()
        kept = f.filter_params(ps, [wanted])
        assert [p.hash_short() for p in kept] == [wanted]

    def test_hash_list_no_match_empty(self):
        assert f.filter_params(self._params(), ["deadbeef"]) == []


# --------------------------------------------------------------------------- drift


class _StubParam:
    """Minimal Parameter stand-in for Params.expand_from_parameter."""

    def __init__(self, values=None, params=None):
        self.values = values
        self.params = params


@pytest.mark.short
class TestFindOrphans:
    def _model(self):
        # Build a tiny stand-in model: stages with .id and .modules (.id, .parameters)
        from types import SimpleNamespace as NS

        m1 = NS(id="M1", parameters=None)
        m2 = NS(id="M2", parameters=None)
        data = NS(id="data", modules=[NS(id="D1", parameters=None)])
        methods = NS(id="methods", modules=[m1, m2])
        return NS(stages=[data, methods])

    def _params_model(self):
        """Same shape, but M1 and M2 each expand to one distinct combo."""
        from types import SimpleNamespace as NS

        m1 = NS(id="M1", parameters=[_StubParam(params={"k": "m1"})])
        m2 = NS(id="M2", parameters=[_StubParam(params={"k": "m2"})])
        data = NS(id="data", modules=[NS(id="D1", parameters=None)])
        return NS(stages=[data, NS(id="methods", modules=[m1, m2])])

    @staticmethod
    def _combo_hash(module):
        return Params.expand_from_parameter(module.parameters[0])[0].hash_short()

    def test_no_orphans_when_all_resolve(self):
        picks = {"data": {"*": "all"}, "methods": {"M1": "all"}}
        assert f.find_orphans(picks, self._model()) == []

    def test_missing_stage(self):
        orphans = f.find_orphans({"ghost": {"X": "all"}}, self._model())
        assert any("stage 'ghost'" in o for o in orphans)

    def test_missing_module(self):
        orphans = f.find_orphans({"methods": {"M9": "all"}}, self._model())
        assert any("methods/M9" in o for o in orphans)

    def test_wildcard_never_orphans(self):
        assert f.find_orphans({"methods": {"*": "all"}}, self._model()) == []

    def test_wildcard_hash_orphans_when_no_module_expands_it(self):
        orphans = f.find_orphans({"methods": {"*": ["deadbeef"]}}, self._params_model())
        assert any("methods/*/deadbeef" in o for o in orphans)

    def test_wildcard_hash_satisfied_by_any_module_in_stage(self):
        model = self._params_model()
        live = self._combo_hash(model.stages[1].modules[1])
        assert f.find_orphans({"methods": {"*": [live]}}, model) == []

    def test_explicit_module_hash_not_satisfied_by_sibling(self):
        model = self._params_model()
        sibling = self._combo_hash(model.stages[1].modules[1])
        orphans = f.find_orphans({"methods": {"M1": [sibling]}}, model)
        assert any(f"methods/M1/{sibling}" in o for o in orphans)
