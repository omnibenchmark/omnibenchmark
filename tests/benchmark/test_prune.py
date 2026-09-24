"""Unit tests for the pure pruning helpers in core/_prune.py."""

from dataclasses import dataclass
from unittest.mock import MagicMock

import pytest

from omnibenchmark.core._prune import (
    apply_until_filter,
    capability_prune_summary,
    empty_stage_warning,
    filter_collectors_by_stages,
    module_capabilities_met,
    select_capable_modules,
)

# ---------------------------------------------------------------------------
# module_capabilities_met (capability gating)
# ---------------------------------------------------------------------------


@dataclass
class _StubModule:
    id: str
    requires_capabilities: list = None


@pytest.mark.short
class TestModuleCapabilitiesMet:
    def test_no_requirements_always_passes(self):
        m = _StubModule("M1", requires_capabilities=None)
        assert module_capabilities_met(m, set()) is True
        assert module_capabilities_met(m, {"gpu"}) is True

    def test_empty_requirements_always_passes(self):
        m = _StubModule("M1", requires_capabilities=[])
        assert module_capabilities_met(m, set()) is True

    def test_required_capability_absent_is_pruned(self):
        m = _StubModule("gpu_pca", requires_capabilities=["gpu"])
        assert module_capabilities_met(m, set()) is False
        assert module_capabilities_met(m, {"large_mem"}) is False

    def test_required_capability_present_passes(self):
        m = _StubModule("gpu_pca", requires_capabilities=["gpu"])
        assert module_capabilities_met(m, {"gpu"}) is True

    def test_all_of_several_required(self):
        m = _StubModule("heavy", requires_capabilities=["gpu", "large_mem"])
        assert module_capabilities_met(m, {"gpu"}) is False
        assert module_capabilities_met(m, {"gpu", "large_mem"}) is True

    def test_none_available_treated_as_empty(self):
        m = _StubModule("gpu_pca", requires_capabilities=["gpu"])
        assert module_capabilities_met(m, None) is False


@pytest.mark.short
class TestSelectCapableModules:
    def _modules(self):
        return [
            _StubModule("M_cpu", requires_capabilities=None),
            _StubModule("M_gpu", requires_capabilities=["gpu"]),
        ]

    def test_prunes_module_missing_capability(self):
        kept, pruned = select_capable_modules(self._modules(), None, set())
        assert [m.id for m in kept] == ["M_cpu"]
        assert [m.id for m in pruned] == ["M_gpu"]

    def test_keeps_all_when_capability_provided(self):
        kept, pruned = select_capable_modules(self._modules(), None, {"gpu"})
        assert [m.id for m in kept] == ["M_cpu", "M_gpu"]
        assert pruned == []

    def test_module_filter_bypasses_gate(self):
        # -m/--module set: nothing is pruned, even the gpu module on a host
        # without gpu — the user asked for it explicitly.
        kept, pruned = select_capable_modules(self._modules(), "M_gpu", set())
        assert [m.id for m in kept] == ["M_cpu", "M_gpu"]
        assert pruned == []

    def test_none_available_prunes_requirers(self):
        kept, pruned = select_capable_modules(self._modules(), None, None)
        assert [m.id for m in kept] == ["M_cpu"]
        assert [m.id for m in pruned] == ["M_gpu"]


@pytest.mark.short
class TestCapabilityPruneSummary:
    def test_names_modules_and_missing_flags(self):
        pruned = [
            _StubModule("M_gpu", requires_capabilities=["gpu"]),
            _StubModule("M_heavy", requires_capabilities=["gpu", "large_mem"]),
        ]
        msg = capability_prune_summary(pruned, {"large_mem"})
        assert "2 module(s) pruned" in msg
        assert "M_gpu, M_heavy" in msg
        assert "--with-capability gpu" in msg
        # already-provided capabilities are not re-suggested
        assert "--with-capability large_mem" not in msg


# ---------------------------------------------------------------------------
# empty_stage_warning
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestEmptyStageWarning:
    def test_attempted_combinations_blame_the_filters(self):
        msg = empty_stage_warning("methods", combinations_seen=4)
        assert "'methods'" in msg
        assert "pruned by a filter (requires/exclude)" in msg

    def test_no_combinations_blames_the_input_wiring(self):
        """Zero combinations means nothing was generated to prune; pointing at
        requires/exclude would send the user to the wrong part of the YAML."""
        msg = empty_stage_warning("methods", combinations_seen=0)
        assert "no upstream output matched its declared inputs" in msg
        assert "filter" not in msg


# ---------------------------------------------------------------------------
# apply_until_filter
# ---------------------------------------------------------------------------


@dataclass
class _StubStage:
    id: str


@pytest.mark.short
class TestApplyUntilFilter:
    def _stages(self, *ids):
        return [_StubStage(i) for i in ids]

    # parents: stage_id -> set of direct upstream stage ids (declared lineage)
    def _linear(self, *ids):
        """Parents map for a simple chain ids[0] -> ids[1] -> ... ."""
        return {ids[i]: {ids[i - 1]} for i in range(1, len(ids))}

    def test_none_returns_all_stages_as_list(self):
        stages = self._stages("a", "b", "c")
        result = apply_until_filter(stages, None, {})
        assert [s.id for s in result] == ["a", "b", "c"]
        # must be a fresh list, not the same object
        assert result is not stages

    def test_initial_stage_keeps_only_first(self):
        result = apply_until_filter(
            self._stages("a", "b", "c"), "a", self._linear("a", "b", "c")
        )
        assert [s.id for s in result] == ["a"]

    def test_middle_stage_includes_named_and_ancestors(self):
        result = apply_until_filter(
            self._stages("a", "b", "c", "d"), "b", self._linear("a", "b", "c", "d")
        )
        assert [s.id for s in result] == ["a", "b"]

    def test_terminal_stage_returns_full_list(self):
        result = apply_until_filter(
            self._stages("a", "b", "c"), "c", self._linear("a", "b", "c")
        )
        assert [s.id for s in result] == ["a", "b", "c"]

    def test_ancestor_declared_after_target_is_kept(self):
        # Lineage a -> b -> c, but declared out of topological order [c, b, a].
        # Old declaration-order slicing kept only [c]; transitive closure keeps all.
        result = apply_until_filter(
            self._stages("c", "b", "a"), "c", self._linear("a", "b", "c")
        )
        assert [s.id for s in result] == ["c", "b", "a"]  # declaration order preserved

    def test_unrelated_branch_is_pruned(self):
        # root feeds both branch_a and branch_b; --until branch_a drops branch_b
        # even though it is declared before branch_a.
        parents = {"branch_a": {"root"}, "branch_b": {"root"}}
        result = apply_until_filter(
            self._stages("root", "branch_b", "branch_a"), "branch_a", parents
        )
        assert [s.id for s in result] == ["root", "branch_a"]

    def test_diamond_keeps_both_paths(self):
        # a -> {b, c} -> d ; --until d keeps everything, --until b drops c.
        #
        # Algorithm-only: `d` has two divergent direct parents, which
        # detect_diamond_input_joins rejects at validate/run time (#289), so this
        # topology never reaches apply_until_filter in a real benchmark. (Two
        # parents on ONE chain are legal and do reach it; two on divergent
        # branches, as here, do not.) Kept because the traversal must stay
        # correct for when #289 lifts and divergent joins become legal.
        parents = {"b": {"a"}, "c": {"a"}, "d": {"b", "c"}}
        stages = self._stages("a", "b", "c", "d")
        assert [s.id for s in apply_until_filter(stages, "d", parents)] == [
            "a",
            "b",
            "c",
            "d",
        ]
        assert [s.id for s in apply_until_filter(stages, "b", parents)] == ["a", "b"]

    def test_root_stage_keeps_only_itself(self):
        result = apply_until_filter(
            self._stages("a", "b", "c"), "a", {"b": {"a"}, "c": {"b"}}
        )
        assert [s.id for s in result] == ["a"]

    def test_unknown_stage_raises(self):
        with pytest.raises(ValueError, match="stage 'nope' not found"):
            apply_until_filter(self._stages("a", "b"), "nope", {})

    def test_unknown_stage_lists_available(self):
        with pytest.raises(ValueError, match="a, b"):
            apply_until_filter(self._stages("a", "b"), "missing", {})

    def test_empty_stages_with_until_raises(self):
        with pytest.raises(ValueError):
            apply_until_filter([], "anything", {})

    def test_empty_stages_without_until_returns_empty(self):
        assert apply_until_filter([], None, {}) == []


# ---------------------------------------------------------------------------
# filter_collectors_by_stages
# ---------------------------------------------------------------------------


@dataclass
class _StubCollector:
    id: str
    inputs: list


@pytest.mark.short
class TestFilterCollectorsByStages:
    def _benchmark(self, output_to_stage):
        """Build a stub benchmark whose get_stages_by_output returns mapped stages.

        Values may be a single stage or a list (an output id is a contract that
        several stages may produce — design 010 §3.1).
        """
        bench = MagicMock()

        def _producers(oid):
            found = output_to_stage.get(oid)
            if found is None:
                return []
            return found if isinstance(found, list) else [found]

        bench.get_stages_by_output.side_effect = _producers
        return bench

    def test_keeps_collector_when_all_inputs_in_included_stages(self):
        bench = self._benchmark({"methods.result": _StubStage("methods")})
        collectors = [_StubCollector("MC1", inputs=["methods.result"])]
        kept, dropped = filter_collectors_by_stages(
            collectors, included_stage_ids={"data", "methods"}, benchmark=bench
        )
        assert [c.id for c in kept] == ["MC1"]
        assert dropped == []

    def test_drops_collector_when_any_input_in_pruned_stage(self):
        bench = self._benchmark(
            {
                "methods.result": _StubStage("methods"),
                "data.raw": _StubStage("data"),
            }
        )
        collectors = [_StubCollector("MC1", inputs=["methods.result", "data.raw"])]
        kept, dropped = filter_collectors_by_stages(
            collectors, included_stage_ids={"data"}, benchmark=bench
        )
        assert kept == []
        assert dropped == ["MC1"]

    def test_unknown_input_does_not_drop(self):
        # Unknown outputs are left to the regular collector resolver to warn about.
        bench = self._benchmark({})
        collectors = [_StubCollector("MC1", inputs=["unknown.id"])]
        kept, dropped = filter_collectors_by_stages(
            collectors, included_stage_ids={"data"}, benchmark=bench
        )
        assert [c.id for c in kept] == ["MC1"]
        assert dropped == []

    def test_keeps_collector_when_one_producer_of_a_shared_id_survives(self):
        # A shared output id (design 010 §3.1) still has data to collect as long
        # as one of its producing stages was kept.
        bench = self._benchmark(
            {"methods.result": [_StubStage("methods_a"), _StubStage("methods_b")]}
        )
        collectors = [_StubCollector("MC1", inputs=["methods.result"])]
        kept, dropped = filter_collectors_by_stages(
            collectors, included_stage_ids={"methods_a"}, benchmark=bench
        )
        assert [c.id for c in kept] == ["MC1"]
        assert dropped == []

    def test_handles_iofile_inputs(self):
        # Inputs may be IOFile objects (with `.id`) rather than bare strings.
        from types import SimpleNamespace

        bench = self._benchmark({"methods.result": _StubStage("methods")})
        iofile = SimpleNamespace(id="methods.result")
        collectors = [_StubCollector("MC1", inputs=[iofile])]
        kept, dropped = filter_collectors_by_stages(
            collectors, included_stage_ids={"data"}, benchmark=bench
        )
        assert kept == []
        assert dropped == ["MC1"]

    def test_empty_collector_list(self):
        bench = self._benchmark({})
        kept, dropped = filter_collectors_by_stages(
            [], included_stage_ids={"data"}, benchmark=bench
        )
        assert kept == []
        assert dropped == []

    def test_none_collector_list(self):
        bench = self._benchmark({})
        kept, dropped = filter_collectors_by_stages(
            None, included_stage_ids={"data"}, benchmark=bench
        )
        assert kept == []
        assert dropped == []

    def test_partial_pruning_keeps_unrelated_collector(self):
        bench = self._benchmark(
            {
                "methods.result": _StubStage("methods"),
                "data.raw": _StubStage("data"),
            }
        )
        collectors = [
            _StubCollector("MC1", inputs=["methods.result"]),
            _StubCollector("MC2", inputs=["data.raw"]),
        ]
        kept, dropped = filter_collectors_by_stages(
            collectors, included_stage_ids={"data"}, benchmark=bench
        )
        assert [c.id for c in kept] == ["MC2"]
        assert dropped == ["MC1"]
