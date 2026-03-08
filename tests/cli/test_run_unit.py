"""Short unit tests for cli/run.py pure and near-pure helper functions."""

import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch
from pydantic import ValidationError as PydanticValidationError

from omnibenchmark.cli.run import (
    format_pydantic_errors,
    _read_rule_log,
    _substitute_params_in_path,
    _build_template_context,
    _select_input_nodes,
    _satisfies_requires,
)
from omnibenchmark.model.resolved import TemplateContext
from omnibenchmark.model.params import Params


# ---------------------------------------------------------------------------
# format_pydantic_errors
# ---------------------------------------------------------------------------


def _make_pydantic_error(errors: list[dict]):
    """Create a PydanticValidationError from a list of raw error dicts."""
    from pydantic import BaseModel

    class _M(BaseModel):
        x: int

    try:
        _M(x="bad")
    except PydanticValidationError:
        # Patch internal errors for custom test data
        pass

    mock_err = MagicMock(spec=PydanticValidationError)
    mock_err.errors.return_value = errors
    return mock_err


@pytest.mark.short
class TestFormatPydanticErrors:
    def test_missing_field(self):
        err = _make_pydantic_error(
            [{"loc": ("stages",), "msg": "Field required", "type": "missing"}]
        )
        result = format_pydantic_errors(err)
        assert "Missing required field: 'stages'" in result

    def test_non_missing_field(self):
        err = _make_pydantic_error(
            [
                {
                    "loc": ("version",),
                    "msg": "Input should be a string",
                    "type": "string_type",
                }
            ]
        )
        result = format_pydantic_errors(err)
        assert "Field 'version'" in result
        assert "Input should be a string" in result

    def test_multiple_errors(self):
        err = _make_pydantic_error(
            [
                {"loc": ("id",), "msg": "Field required", "type": "missing"},
                {"loc": ("version",), "msg": "bad value", "type": "value_error"},
            ]
        )
        result = format_pydantic_errors(err)
        assert "Validation failed:" in result
        assert "Missing required field: 'id'" in result
        assert "Field 'version'" in result

    def test_nested_loc(self):
        err = _make_pydantic_error(
            [{"loc": ("stages", 0, "id"), "msg": "bad", "type": "value_error"}]
        )
        result = format_pydantic_errors(err)
        assert "stages -> 0 -> id" in result


# ---------------------------------------------------------------------------
# _read_rule_log
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestReadRuleLog:
    def test_returns_none_when_missing(self, tmp_path):
        assert _read_rule_log(tmp_path, "some_rule") is None

    def test_reads_existing_log(self, tmp_path):
        log_dir = tmp_path / ".logs"
        log_dir.mkdir()
        (log_dir / "my_rule.log").write_text("some log content")
        result = _read_rule_log(tmp_path, "my_rule")
        assert result == "some log content"

    def test_returns_none_on_read_error(self, tmp_path):
        log_dir = tmp_path / ".logs"
        log_dir.mkdir()
        log_path = log_dir / "bad_rule.log"
        log_path.write_text("x")
        with patch.object(Path, "read_text", side_effect=OSError("perm")):
            result = _read_rule_log(tmp_path, "bad_rule")
        assert result is None


# ---------------------------------------------------------------------------
# _substitute_params_in_path
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestSubstituteParamsInPath:
    def test_no_placeholder_unchanged(self):
        assert (
            _substitute_params_in_path("data/D1/out.json", None) == "data/D1/out.json"
        )

    def test_none_params_unchanged(self):
        assert (
            _substitute_params_in_path("data/{params.key}/out.json", None)
            == "data/{params.key}/out.json"
        )

    def test_substitutes_key(self):
        p = Params({"key": "value123"})
        result = _substitute_params_in_path("data/{params.key}/out.json", p)
        assert result == "data/value123/out.json"

    def test_missing_key_left_as_is(self):
        p = Params({"other": "x"})
        result = _substitute_params_in_path("data/{params.missing}/out.json", p)
        assert result == "data/{params.missing}/out.json"

    def test_multiple_substitutions(self):
        p = Params({"a": "1", "b": "2"})
        result = _substitute_params_in_path("{params.a}_{params.b}.txt", p)
        assert result == "1_2.txt"

    def test_no_params_placeholder_skips_regex(self):
        p = Params({"k": "v"})
        assert _substitute_params_in_path("plain/path.txt", p) == "plain/path.txt"


# ---------------------------------------------------------------------------
# _build_template_context
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
        ctx = _build_template_context(stage, "D1")
        assert ctx.provides["dataset"] == "D1"
        assert ctx.module_attrs["id"] == "D1"
        assert ctx.module_attrs["stage"] == "data"

    def test_root_node_dataset_from_params(self):
        stage = _make_stage("data", provides=None)
        p = Params({"dataset": "pbmc3k"})
        ctx = _build_template_context(stage, "D1", params=p)
        assert ctx.provides["dataset"] == "pbmc3k"

    def test_root_node_provides_label_from_params(self):
        stage = _make_stage("data", provides=["treatment"])
        p = Params({"treatment": "ctrl"})
        ctx = _build_template_context(stage, "D1", params=p)
        assert ctx.provides["treatment"] == "ctrl"

    def test_root_node_provides_label_defaults_to_module_id(self):
        stage = _make_stage("data", provides=["treatment"])
        ctx = _build_template_context(stage, "D1")
        assert ctx.provides["treatment"] == "D1"

    def test_child_node_inherits_parent_context(self):
        parent_ctx = TemplateContext(
            provides={"dataset": "pbmc3k"},
            module_attrs={"id": "D1", "stage": "data"},
        )
        input_node = _make_input_node("D1", "data", template_context=parent_ctx)
        stage = _make_stage("methods", provides=None)
        ctx = _build_template_context(stage, "M1", input_node=input_node)
        assert ctx.provides["dataset"] == "pbmc3k"
        assert ctx.module_attrs["parent.id"] == "D1"
        assert ctx.module_attrs["parent.stage"] == "data"

    def test_child_node_no_parent_context(self):
        input_node = _make_input_node("D1", "data", template_context=None)
        stage = _make_stage("methods", provides=None)
        ctx = _build_template_context(stage, "M1", input_node=input_node)
        assert "dataset" not in ctx.provides
        assert ctx.module_attrs["parent.id"] == "D1"

    def test_child_node_adds_provides_label(self):
        parent_ctx = TemplateContext(
            provides={"dataset": "pbmc3k"},
            module_attrs={"id": "D1", "stage": "data"},
        )
        input_node = _make_input_node("D1", "data", template_context=parent_ctx)
        stage = _make_stage("methods", provides=["method"])
        ctx = _build_template_context(stage, "M1", input_node=input_node)
        assert ctx.provides["method"] == "M1"

    def test_child_node_provides_label_from_params(self):
        parent_ctx = TemplateContext(
            provides={"dataset": "pbmc3k"},
            module_attrs={"id": "D1", "stage": "data"},
        )
        input_node = _make_input_node("D1", "data", template_context=parent_ctx)
        stage = _make_stage("methods", provides=["method"])
        p = Params({"method": "kmeans"})
        ctx = _build_template_context(stage, "M1", input_node=input_node, params=p)
        assert ctx.provides["method"] == "kmeans"


# ---------------------------------------------------------------------------
# _satisfies_requires
# ---------------------------------------------------------------------------


@pytest.mark.short
class TestSatisfiesRequires:
    def test_no_template_context_returns_false(self):
        node = _make_input_node("D1", "data", template_context=None)
        assert _satisfies_requires({"dataset": "pbmc3k"}, node) is False

    def test_matching_single_constraint(self):
        ctx = TemplateContext(
            provides={"dataset": "pbmc3k"},
            module_attrs={"id": "D1", "stage": "data"},
        )
        node = _make_input_node("D1", "data", template_context=ctx)
        assert _satisfies_requires({"dataset": "pbmc3k"}, node) is True

    def test_mismatched_constraint(self):
        ctx = TemplateContext(
            provides={"dataset": "pbmc3k"},
            module_attrs={"id": "D1", "stage": "data"},
        )
        node = _make_input_node("D1", "data", template_context=ctx)
        assert _satisfies_requires({"dataset": "other"}, node) is False

    def test_missing_label_returns_false(self):
        ctx = TemplateContext(
            provides={},
            module_attrs={"id": "D1", "stage": "data"},
        )
        node = _make_input_node("D1", "data", template_context=ctx)
        assert _satisfies_requires({"dataset": "pbmc3k"}, node) is False

    def test_empty_requires_returns_true(self):
        ctx = TemplateContext(
            provides={"dataset": "pbmc3k"},
            module_attrs={"id": "D1", "stage": "data"},
        )
        node = _make_input_node("D1", "data", template_context=ctx)
        assert _satisfies_requires({}, node) is True

    def test_multiple_constraints_all_match(self):
        ctx = TemplateContext(
            provides={"dataset": "pbmc3k", "treatment": "ctrl"},
            module_attrs={"id": "D1", "stage": "data"},
        )
        node = _make_input_node("D1", "data", template_context=ctx)
        assert (
            _satisfies_requires({"dataset": "pbmc3k", "treatment": "ctrl"}, node)
            is True
        )

    def test_multiple_constraints_partial_mismatch(self):
        ctx = TemplateContext(
            provides={"dataset": "pbmc3k", "treatment": "ctrl"},
            module_attrs={"id": "D1", "stage": "data"},
        )
        node = _make_input_node("D1", "data", template_context=ctx)
        assert (
            _satisfies_requires({"dataset": "pbmc3k", "treatment": "stim"}, node)
            is False
        )


# ---------------------------------------------------------------------------
# _select_input_nodes
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
        result = _select_input_nodes([], {}, [], [], prev)
        assert result is prev

    def test_no_matching_outputs_returns_previous(self):
        prev = [_make_node("n1", "data")]
        result = _select_input_nodes(["data.raw"], {}, [], ["data"], prev)
        assert result is prev

    def test_single_input_selects_correct_stage(self):
        node_data = _make_node("data-D1", "data")
        node_methods = _make_node("methods-M1", "methods")
        resolved = [node_data, node_methods]
        output_to_nodes = {"data.raw": [("data-D1", "data/D1/out.json")]}
        stage_ids = ["data", "methods"]
        prev = []
        result = _select_input_nodes(
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
        result = _select_input_nodes(
            ["data.raw", "prep.out"], output_to_nodes, resolved, stage_ids, prev
        )
        assert all(n.stage_id == "preprocessing" for n in result)

    def test_node_not_in_resolved_skipped(self):
        output_to_nodes = {"data.raw": [("ghost-node", "p.json")]}
        prev = [_make_node("fallback", "data")]
        result = _select_input_nodes(["data.raw"], output_to_nodes, [], ["data"], prev)
        assert result is prev
