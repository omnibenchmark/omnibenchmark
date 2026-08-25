"""Unit tests for select_input_nodes (core/_lineage.py)."""

from dataclasses import dataclass

import pytest

from omnibenchmark.core._lineage import select_input_nodes


# ---------------------------------------------------------------------------
# Minimal stub that satisfies the .id / .stage_id duck-typing used by
# select_input_nodes.  Using a plain dataclass avoids importing heavy
# ResolvedNode infrastructure.
# ---------------------------------------------------------------------------


@dataclass
class _StubNode:
    id: str
    stage_id: str


# ---------------------------------------------------------------------------
# Helper: build output_to_nodes registry from a mapping
#   { output_id: [(node_id, path_str), ...] }
# ---------------------------------------------------------------------------


def _reg(*pairs):
    """Build an output_to_nodes dict.  Each pair is (output_id, node_id)."""
    reg = {}
    for output_id, node_id in pairs:
        reg.setdefault(output_id, []).append((node_id, f"/out/{node_id}/{output_id}"))
    return reg


# ===========================================================================
# TestSelectInputNodes
# ===========================================================================


@pytest.mark.short
class TestSelectInputNodes:
    """Tests for select_input_nodes pure function."""

    def test_no_declared_inputs_returns_previous(self):
        """When declared_input_ids is empty fall back to previous_stage_nodes."""
        prev = [_StubNode("n1", "stage_a")]
        result = select_input_nodes(
            declared_input_ids=[],
            output_to_nodes={},
            resolved_nodes=[],
            stage_ids_in_order=["stage_a"],
            previous_stage_nodes=prev,
        )
        assert result is prev

    def test_unresolvable_inputs_return_previous(self):
        """When no declared input exists in the registry fall back to previous."""
        prev = [_StubNode("n1", "stage_a")]
        result = select_input_nodes(
            declared_input_ids=["unknown_output"],
            output_to_nodes={},
            resolved_nodes=[],
            stage_ids_in_order=["stage_a"],
            previous_stage_nodes=prev,
        )
        assert result is prev

    def test_simple_linear_chain(self):
        """Inputs produced by previous stage → return that stage's nodes."""
        n1 = _StubNode("n1", "stage_a")
        reg = _reg(("pca.csv", "n1"))

        result = select_input_nodes(
            declared_input_ids=["pca.csv"],
            output_to_nodes=reg,
            resolved_nodes=[n1],
            stage_ids_in_order=["stage_a"],
            previous_stage_nodes=[],
        )
        assert result == [n1]

    def test_skip_dependency_picks_producing_stage_not_predecessor(self):
        """
        Skip-dependency scenario:
          data → pca → (leiden) → umap
        umap declares inputs produced by *pca*, not leiden.
        select_input_nodes must return pca nodes even though leiden is the
        most recent previous stage.
        """
        n_pca = _StubNode("pca-M1-default", "pca")
        n_leiden = _StubNode("leiden-M1-default", "leiden")
        # leiden is the previous_stage_nodes (most recent predecessor)
        # umap declares "neighbors.h5ad" which was produced by pca
        reg = _reg(("neighbors.h5ad", "pca-M1-default"))

        result = select_input_nodes(
            declared_input_ids=["neighbors.h5ad"],
            output_to_nodes=reg,
            resolved_nodes=[n_pca, n_leiden],
            stage_ids_in_order=["data", "pca", "leiden", "umap"],
            previous_stage_nodes=[n_leiden],
        )
        assert result == [n_pca]
        assert n_leiden not in result

    def test_inputs_from_two_stages_picks_deepest(self):
        """
        When declared inputs span two stages, pick the deepest one (later in
        document order) so its nodes carry the fullest TemplateContext.
        """
        n_data = _StubNode("data-D1-default", "data")
        n_pca = _StubNode("pca-M1-default", "pca")
        reg = _reg(
            ("data.h5ad", "data-D1-default"),
            ("pca.csv", "pca-M1-default"),
        )

        result = select_input_nodes(
            declared_input_ids=["data.h5ad", "pca.csv"],
            output_to_nodes=reg,
            resolved_nodes=[n_data, n_pca],
            stage_ids_in_order=["data", "pca", "harmony"],
            previous_stage_nodes=[],
        )
        # pca is deeper (index 1 > 0) → return pca nodes
        assert result == [n_pca]

    def test_multiple_nodes_in_deepest_stage_all_returned(self):
        """All nodes from the deepest stage are returned, not just one."""
        n1 = _StubNode("pca-M1-default", "pca")
        n2 = _StubNode("pca-M2-default", "pca")
        reg = _reg(("pca.csv", "pca-M1-default"), ("pca.csv", "pca-M2-default"))

        result = select_input_nodes(
            declared_input_ids=["pca.csv"],
            output_to_nodes=reg,
            resolved_nodes=[n1, n2],
            stage_ids_in_order=["data", "pca"],
            previous_stage_nodes=[],
        )
        assert set(n.id for n in result) == {"pca-M1-default", "pca-M2-default"}

    def test_node_not_in_resolved_nodes_is_ignored(self):
        """If the registry references a node_id not in resolved_nodes, skip it."""
        n1 = _StubNode("pca-M1-default", "pca")
        reg = _reg(
            ("pca.csv", "pca-M1-default"),
            ("pca.csv", "ghost-node"),  # ghost not in resolved_nodes
        )

        result = select_input_nodes(
            declared_input_ids=["pca.csv"],
            output_to_nodes=reg,
            resolved_nodes=[n1],
            stage_ids_in_order=["data", "pca"],
            previous_stage_nodes=[],
        )
        # Only n1 should be in result; ghost is silently skipped
        assert result == [n1]

    def test_stage_not_in_order_list_treated_as_depth_minus_one(self):
        """
        A node whose stage_id is absent from stage_ids_in_order gets depth -1
        and should lose to any real stage.
        """
        n_orphan = _StubNode("orphan-node", "orphan_stage")
        n_real = _StubNode("pca-M1-default", "pca")
        reg = _reg(
            ("output.csv", "orphan-node"),
            ("output.csv", "pca-M1-default"),
        )

        result = select_input_nodes(
            declared_input_ids=["output.csv"],
            output_to_nodes=reg,
            resolved_nodes=[n_orphan, n_real],
            stage_ids_in_order=["data", "pca"],
            previous_stage_nodes=[],
        )
        assert result == [n_real]

    def test_umap_harmony_neighbors_scenario(self):
        """
        Full scenario from the scrna benchmark:
          data → pca → harmony → leiden → umap
        umap inputs: harmony.csv, neighbors.h5ad
        harmony.csv produced by harmony stage
        neighbors.h5ad produced by harmony stage (same depth as harmony, deepest of the two)
        leiden is previous_stage_nodes.
        """
        n_pca = _StubNode("pca-S1-default", "pca")
        n_harmony = _StubNode("harmony-H1-default", "harmony")
        n_leiden_s = _StubNode("leiden-L1-scanpy", "leiden")
        n_leiden_r = _StubNode("leiden-L1-rapids", "leiden")
        reg = _reg(
            ("harmony.csv", "harmony-H1-default"),
            ("neighbors.h5ad", "harmony-H1-default"),
        )
        stage_order = ["data", "pca", "harmony", "leiden", "umap"]
        previous = [n_leiden_s, n_leiden_r]

        result = select_input_nodes(
            declared_input_ids=["harmony.csv", "neighbors.h5ad"],
            output_to_nodes=reg,
            resolved_nodes=[n_pca, n_harmony, n_leiden_s, n_leiden_r],
            stage_ids_in_order=stage_order,
            previous_stage_nodes=previous,
        )
        assert result == [n_harmony]
        assert n_leiden_s not in result
        assert n_leiden_r not in result


# ===========================================================================
# Shared output ids: several stages declaring the same output id (design 010
# §3.1). Nodes carry parent_id so ancestry — and therefore domination — is
# real; the stub above has none, which is the single-producer case.
# ===========================================================================


@dataclass
class _ChainNode:
    """Stub with a parent link, so ancestor walks see a real lineage."""

    id: str
    stage_id: str
    parent_id: str = None


@pytest.mark.short
class TestSharedOutputIdProducers:
    """`pcas_tsv` declared by both PCA and CNTFCT — interchangeable for any
    consumer referencing that id, so the consumer expands once per producer."""

    def _two_branch_setup(self):
        # feat -> pca and feat -> cntfct; both declare `pcas.tsv`.
        n_feat = _ChainNode("feat-F1", "feat")
        n_pca = _ChainNode("feat-F1-pca-P1", "pca", parent_id="feat-F1")
        n_cnt = _ChainNode("feat-F1-cntfct-C1", "cntfct", parent_id="feat-F1")
        reg = _reg(
            ("feat.h5", "feat-F1"),
            ("pcas.tsv", "feat-F1-pca-P1"),
            ("pcas.tsv", "feat-F1-cntfct-C1"),
        )
        return n_feat, n_pca, n_cnt, reg

    def test_alternatives_on_parallel_branches_all_selected(self):
        n_feat, n_pca, n_cnt, reg = self._two_branch_setup()

        result = select_input_nodes(
            declared_input_ids=["pcas.tsv"],
            output_to_nodes=reg,
            resolved_nodes=[n_feat, n_pca, n_cnt],
            stage_ids_in_order=["feat", "pca", "cntfct", "embed"],
            previous_stage_nodes=[],
        )
        assert sorted(n.id for n in result) == sorted([n_pca.id, n_cnt.id])

    def test_shared_ancestor_is_dominated_and_excluded(self):
        """`feat` also produces a declared input, but it is an ancestor of both
        alternatives, so it must not become an expansion base of its own."""
        n_feat, n_pca, n_cnt, reg = self._two_branch_setup()

        result = select_input_nodes(
            declared_input_ids=["pcas.tsv", "feat.h5"],
            output_to_nodes=reg,
            resolved_nodes=[n_feat, n_pca, n_cnt],
            stage_ids_in_order=["feat", "pca", "cntfct", "embed"],
            previous_stage_nodes=[],
        )
        assert sorted(n.id for n in result) == sorted([n_pca.id, n_cnt.id])
        assert n_feat not in result

    def test_shadowing_producers_on_one_chain_pick_nearest(self):
        """Same id declared twice along a single chain is shadowing, not
        alternatives: the nearest (deepest) producer wins, as before."""
        n_up = _ChainNode("up-U1", "up")
        n_down = _ChainNode("up-U1-down-D1", "down", parent_id="up-U1")
        reg = _reg(("shared.tsv", "up-U1"), ("shared.tsv", "up-U1-down-D1"))

        result = select_input_nodes(
            declared_input_ids=["shared.tsv"],
            output_to_nodes=reg,
            resolved_nodes=[n_up, n_down],
            stage_ids_in_order=["up", "down", "consumer"],
            previous_stage_nodes=[],
        )
        assert result == [n_down]

    def test_fan_in_partners_are_not_alternatives(self):
        """Two parallel producers of *different* ids are a join, not
        alternatives: only one anchor, so the bundle logic emits the join once
        rather than once per branch."""
        n_root = _ChainNode("root-R1", "root")
        n_left = _ChainNode("root-R1-left-L1", "left", parent_id="root-R1")
        n_right = _ChainNode("root-R1-right-G1", "right", parent_id="root-R1")
        reg = _reg(("left.tsv", "root-R1-left-L1"), ("right.tsv", "root-R1-right-G1"))

        result = select_input_nodes(
            declared_input_ids=["left.tsv", "right.tsv"],
            output_to_nodes=reg,
            resolved_nodes=[n_root, n_left, n_right],
            stage_ids_in_order=["root", "left", "right", "join"],
            previous_stage_nodes=[],
        )
        assert len(result) == 1
        assert result == [n_right]
