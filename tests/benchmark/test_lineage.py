"""Unit tests for the pure lineage helpers in core/_lineage.py."""

from unittest.mock import MagicMock

import pytest

from omnibenchmark.core._lineage import (
    build_template_context,
    inherited_provides,
    lineage_module_ids,
    resolve_label_value,
    satisfies_requires,
    select_input_nodes,
)
from omnibenchmark.model.params import Params
from omnibenchmark.model.resolved import TemplateContext

# ---------------------------------------------------------------------------
# resolve_label_value
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestResolveLabelValue:
    def test_module_binding_wins(self):
        assert resolve_label_value("size", {"size": "lg"}, "M1") == "lg"

    def test_defaults_to_module_id_when_no_binding(self):
        assert resolve_label_value("size", None, "M1") == "M1"

    def test_defaults_to_module_id_when_label_absent_from_binding(self):
        assert resolve_label_value("size", {"other": "x"}, "M1") == "M1"

    def test_non_string_binding_coerced(self):
        assert resolve_label_value("n", {"n": 3}, "M1") == "3"


# ---------------------------------------------------------------------------
# build_template_context
# ---------------------------------------------------------------------------


def _make_stage(stage_id, provides=None):
    s = MagicMock()
    s.id = stage_id
    s.provides = provides
    return s


def _make_input_node(module_id, stage_id, template_context=None):
    n = MagicMock()
    n.module_id = module_id
    n.stage_id = stage_id
    n.template_context = template_context
    return n


@pytest.mark.short
class TestBuildTemplateContext:
    def test_root_node_default_dataset(self):
        stage = _make_stage("data", provides=None)
        ctx = build_template_context(stage, "D1")
        assert ctx.provides["dataset"] == "D1"
        assert ctx.module_attrs["id"] == "D1"
        assert ctx.module_attrs["stage"] == "data"

    def test_root_node_dataset_from_params(self):
        stage = _make_stage("data", provides=None)
        p = Params({"dataset": "pbmc3k"})
        ctx = build_template_context(stage, "D1", params=p)
        assert ctx.provides["dataset"] == "pbmc3k"

    def test_root_node_provides_label_from_module_binding(self):
        stage = _make_stage("data", provides=["treatment"])
        ctx = build_template_context(stage, "D1", module_provides={"treatment": "ctrl"})
        assert ctx.provides["treatment"] == "ctrl"

    def test_root_node_provides_label_ignores_same_named_param(self):
        # The parameter-name fallback was cut: a same-named param must NOT
        # become the label value (see 008 §3.5). Falls through to module id.
        stage = _make_stage("data", provides=["treatment"])
        p = Params({"treatment": "ctrl"})
        ctx = build_template_context(stage, "D1", params=p)
        assert ctx.provides["treatment"] == "D1"

    def test_root_node_provides_label_defaults_to_module_id(self):
        stage = _make_stage("data", provides=["treatment"])
        ctx = build_template_context(stage, "D1")
        assert ctx.provides["treatment"] == "D1"

    def test_child_node_inherits_parent_context(self):
        parent_ctx = TemplateContext(
            provides={"dataset": "pbmc3k"},
            module_attrs={"id": "D1", "stage": "data"},
        )
        input_node = _make_input_node("D1", "data", template_context=parent_ctx)
        stage = _make_stage("methods", provides=None)
        ctx = build_template_context(stage, "M1", input_nodes=(input_node,))
        assert ctx.provides["dataset"] == "pbmc3k"
        assert ctx.module_attrs["parent.id"] == "D1"
        assert ctx.module_attrs["parent.stage"] == "data"

    def test_child_node_no_parent_context(self):
        input_node = _make_input_node("D1", "data", template_context=None)
        stage = _make_stage("methods", provides=None)
        ctx = build_template_context(stage, "M1", input_nodes=(input_node,))
        assert "dataset" not in ctx.provides
        assert ctx.module_attrs["parent.id"] == "D1"

    def test_child_node_adds_provides_label(self):
        parent_ctx = TemplateContext(
            provides={"dataset": "pbmc3k"},
            module_attrs={"id": "D1", "stage": "data"},
        )
        input_node = _make_input_node("D1", "data", template_context=parent_ctx)
        stage = _make_stage("methods", provides=["method"])
        ctx = build_template_context(stage, "M1", input_nodes=(input_node,))
        assert ctx.provides["method"] == "M1"

    def test_child_node_provides_label_from_module_binding(self):
        parent_ctx = TemplateContext(
            provides={"dataset": "pbmc3k"},
            module_attrs={"id": "D1", "stage": "data"},
        )
        input_node = _make_input_node("D1", "data", template_context=parent_ctx)
        stage = _make_stage("methods", provides=["method"])
        ctx = build_template_context(
            stage, "M1", input_nodes=(input_node,), module_provides={"method": "kmeans"}
        )
        assert ctx.provides["method"] == "kmeans"
        # parent lineage labels still inherited
        assert ctx.provides["dataset"] == "pbmc3k"

    # --- {name} template variable: always the current module's own ID ---

    def test_root_node_provides_name_is_module_id(self):
        stage = _make_stage("data", provides=None)
        ctx = build_template_context(stage, "D1")
        assert ctx.provides["name"] == "D1"

    def test_child_node_provides_name_is_own_module_id_not_inherited(self):
        # {name} must always be the *current* module's ID, never inherited from parent
        parent_ctx = TemplateContext(
            provides={"dataset": "D1", "name": "D1"},
            module_attrs={"id": "D1", "stage": "data"},
        )
        input_node = _make_input_node("D1", "data", template_context=parent_ctx)
        stage = _make_stage("methods", provides=None)
        ctx = build_template_context(stage, "M1", input_nodes=(input_node,))
        assert ctx.provides["name"] == "M1"

    def test_name_template_variable_substitution(self):
        stage = _make_stage("methods", provides=None)
        ctx = build_template_context(stage, "M1")
        assert ctx.substitute("{name}_result.txt") == "M1_result.txt"

    # --- {module.name} template variable: module's human-readable name ---

    def test_module_name_attr_uses_provided_name(self):
        stage = _make_stage("data", provides=None)
        ctx = build_template_context(stage, "D1", module_name="Dataset 1")
        assert ctx.module_attrs["name"] == "Dataset 1"

    def test_module_name_attr_falls_back_to_module_id(self):
        stage = _make_stage("data", provides=None)
        ctx = build_template_context(stage, "D1", module_name=None)
        assert ctx.module_attrs["name"] == "D1"

    def test_module_name_template_variable_substitution(self):
        stage = _make_stage("data", provides=None)
        ctx = build_template_context(stage, "D1", module_name="Dataset 1")
        assert ctx.substitute("{module.name}_output.txt") == "Dataset 1_output.txt"


# ---------------------------------------------------------------------------
# inherited_provides / satisfies_requires
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestInheritedProvides:
    def _node(self, module_id, **provides):
        return _make_input_node(
            module_id,
            "s",
            template_context=TemplateContext(
                provides=provides, module_attrs={"id": module_id, "stage": "s"}
            ),
        )

    def test_no_members_is_empty(self):
        assert inherited_provides(()) == {}

    def test_single_member_passes_its_labels_through(self):
        node = self._node("D1", dataset="pbmc3k")
        assert inherited_provides((node,)) == {"dataset": "pbmc3k"}

    def test_member_without_context_contributes_nothing(self):
        node = _make_input_node("D1", "s", template_context=None)
        assert inherited_provides((node,)) == {}

    def test_join_unions_every_branch(self):
        """A fan-in node sees every branch's labels, not just the anchor's —
        the same lineage `exclude` unions over (design 010 §5.2)."""
        anchor = self._node("log1p", dataset="pbmc3k", normalization="log")
        partner = self._node("umap", dataset="pbmc3k", embedding="umap")
        assert inherited_provides((anchor, partner)) == {
            "dataset": "pbmc3k",
            "normalization": "log",
            "embedding": "umap",
        }


@pytest.mark.short
class TestSatisfiesRequires:
    def test_empty_requires_is_satisfied_by_anything(self):
        assert satisfies_requires({}, {}) is True
        assert satisfies_requires({}, {"dataset": "pbmc3k"}) is True

    def test_matching_single_constraint(self):
        assert satisfies_requires({"dataset": "pbmc3k"}, {"dataset": "pbmc3k"}) is True

    def test_mismatched_value(self):
        assert satisfies_requires({"dataset": "other"}, {"dataset": "pbmc3k"}) is False

    def test_absent_label(self):
        assert satisfies_requires({"dataset": "pbmc3k"}, {}) is False

    def test_multiple_constraints_all_match(self):
        provides = {"dataset": "pbmc3k", "treatment": "ctrl"}
        assert satisfies_requires(provides, provides) is True

    def test_multiple_constraints_partial_mismatch(self):
        assert (
            satisfies_requires(
                {"dataset": "pbmc3k", "treatment": "stim"},
                {"dataset": "pbmc3k", "treatment": "ctrl"},
            )
            is False
        )

    def test_gate_matches_a_label_from_a_non_anchor_branch(self):
        """The point of the union: a gate may name a label carried by any
        branch of the join, not only the one that won the id spine."""
        anchor = TemplateContext(
            provides={"normalization": "log"}, module_attrs={"id": "log1p"}
        )
        partner = TemplateContext(
            provides={"embedding": "umap"}, module_attrs={"id": "umap"}
        )
        upstream = inherited_provides(
            (
                _make_input_node("log1p", "norm", template_context=anchor),
                _make_input_node("umap", "embed", template_context=partner),
            )
        )
        assert satisfies_requires({"embedding": "umap"}, upstream) is True


# ---------------------------------------------------------------------------
# lineage_module_ids
# ---------------------------------------------------------------------------


def _make_lineage_node(module_id, parent_id=None, node_id=None):
    n = MagicMock()
    n.id = node_id if node_id is not None else module_id
    n.module_id = module_id
    n.parent_id = parent_id
    n.parents = []
    n.is_gather = False
    return n


def _by_id(*nodes):
    return {n.id: n for n in nodes}


@pytest.mark.short
class TestLineageModuleIds:
    def test_root_node_lineage_is_self(self):
        root = _make_lineage_node("adamson")
        assert lineage_module_ids(root, _by_id(root)) == {"adamson"}

    def test_immediate_predecessor_in_lineage(self):
        root = _make_lineage_node("adamson")
        child = _make_lineage_node("prep", parent_id=root.id)
        assert lineage_module_ids(child, _by_id(root, child)) == {"adamson", "prep"}

    def test_transitive_lineage_across_non_adjacent_stages(self):
        """Regression: a download module several stages upstream must appear in
        the lineage of a leaf node, not just the immediate predecessor."""
        download = _make_lineage_node("adamson")
        preprocess = _make_lineage_node("prep", parent_id=download.id)
        split = _make_lineage_node("sim", parent_id=preprocess.id)
        # split is the immediate predecessor of a `methods` node; its lineage
        # must still include the non-adjacent `adamson` download module.
        assert lineage_module_ids(split, _by_id(download, preprocess, split)) == {
            "adamson",
            "prep",
            "sim",
        }

    def test_unrelated_branch_excluded_from_lineage(self):
        download = _make_lineage_node("adamson")
        other = _make_lineage_node("norman")  # sibling root, not an ancestor
        child = _make_lineage_node("prep", parent_id=download.id)
        assert lineage_module_ids(child, _by_id(download, other, child)) == {
            "adamson",
            "prep",
        }

    def test_dashed_module_ids_resolve_via_parent_links(self):
        """Lineage follows explicit parent_id links, so module ids containing the
        id separator ('-') cannot confuse ancestry (the old prefix-match risk)."""
        download = _make_lineage_node(
            "adamson-v2", node_id="download-adamson-v2.default"
        )
        method = _make_lineage_node(
            "my-method",
            parent_id=download.id,
            node_id="download-adamson-v2.default-methods-my-method.default",
        )
        assert lineage_module_ids(method, _by_id(download, method)) == {
            "adamson-v2",
            "my-method",
        }


# ---------------------------------------------------------------------------
# select_input_nodes
# ---------------------------------------------------------------------------


def _make_node(node_id, stage_id):
    n = MagicMock()
    n.id = node_id
    n.stage_id = stage_id
    return n


@pytest.mark.short
class TestSelectInputNodes:
    def test_empty_declared_inputs_returns_previous(self):
        prev = [_make_node("n1", "data")]
        result = select_input_nodes([], {}, [], [], prev)
        assert result is prev

    def test_no_matching_outputs_returns_previous(self):
        prev = [_make_node("n1", "data")]
        result = select_input_nodes(["data.raw"], {}, [], ["data"], prev)
        assert result is prev

    def test_single_input_selects_correct_stage(self):
        node_data = _make_node("data-D1", "data")
        node_methods = _make_node("methods-M1", "methods")
        resolved = [node_data, node_methods]
        output_to_nodes = {"data.raw": [("data-D1", "data/D1/out.json")]}
        stage_ids = ["data", "methods"]
        prev = []
        result = select_input_nodes(
            ["data.raw"], output_to_nodes, resolved, stage_ids, prev
        )
        assert all(n.stage_id == "data" for n in result)

    def test_selects_deepest_providing_stage(self):
        n_data = _make_node("data-D1", "data")
        n_prep = _make_node("prep-P1", "preprocessing")
        resolved = [n_data, n_prep]
        output_to_nodes = {
            "data.raw": [("data-D1", "p1.json")],
            "prep.out": [("prep-P1", "p2.json")],
        }
        stage_ids = ["data", "preprocessing", "methods"]
        prev = []
        result = select_input_nodes(
            ["data.raw", "prep.out"], output_to_nodes, resolved, stage_ids, prev
        )
        assert all(n.stage_id == "preprocessing" for n in result)

    def test_node_not_in_resolved_skipped(self):
        output_to_nodes = {"data.raw": [("ghost-node", "p.json")]}
        prev = [_make_node("fallback", "data")]
        result = select_input_nodes(["data.raw"], output_to_nodes, [], ["data"], prev)
        assert result is prev
