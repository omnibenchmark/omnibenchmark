"""Tests for the `requires` (explicit plugs) feature on Module."""

import pytest

from omnibenchmark.model.validation import ValidationError
from omnibenchmark.model.resolved import TemplateContext
from omnibenchmark.cli.run import _satisfies_requires
from tests.model.factories import make_benchmark, make_module


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_ctx(provides: dict) -> "object":
    """Return a minimal stub with a template_context carrying the given provides."""

    class _Stub:
        template_context = TemplateContext(provides=provides)

    return _Stub()


# ---------------------------------------------------------------------------
# Unit: _satisfies_requires
# ---------------------------------------------------------------------------


class TestSatisfiesRequires:
    @pytest.mark.short
    def test_returns_true_when_all_constraints_match(self):
        node = _make_ctx({"pca": "scanpy", "dataset": "D1"})
        assert _satisfies_requires({"pca": "scanpy"}, node) is True

    @pytest.mark.short
    def test_returns_false_when_label_value_mismatch(self):
        node = _make_ctx({"pca": "rapids"})
        assert _satisfies_requires({"pca": "scanpy"}, node) is False

    @pytest.mark.short
    def test_returns_false_when_label_missing(self):
        node = _make_ctx({"dataset": "D1"})
        assert _satisfies_requires({"pca": "scanpy"}, node) is False

    @pytest.mark.short
    def test_returns_true_for_multiple_constraints_all_match(self):
        node = _make_ctx({"pca": "scanpy", "harmony": "scanpy", "dataset": "D1"})
        assert _satisfies_requires({"pca": "scanpy", "harmony": "scanpy"}, node) is True

    @pytest.mark.short
    def test_returns_false_when_one_of_multiple_constraints_fails(self):
        node = _make_ctx({"pca": "scanpy", "harmony": "rapids"})
        assert (
            _satisfies_requires({"pca": "scanpy", "harmony": "scanpy"}, node) is False
        )

    @pytest.mark.short
    def test_returns_false_when_template_context_is_none(self):
        class _NoCtx:
            template_context = None

        assert _satisfies_requires({"pca": "scanpy"}, _NoCtx()) is False

    @pytest.mark.short
    def test_empty_requires_always_returns_true(self):
        node = _make_ctx({})
        assert _satisfies_requires({}, node) is True


# ---------------------------------------------------------------------------
# Unit: Model field
# ---------------------------------------------------------------------------


class TestRequiresField:
    @pytest.mark.short
    def test_module_accepts_requires_field(self):
        m = make_module(requires={"pca": "scanpy"})
        assert m.requires == {"pca": "scanpy"}

    @pytest.mark.short
    def test_module_requires_defaults_to_none(self):
        m = make_module()
        assert m.requires is None

    @pytest.mark.short
    def test_module_accepts_multiple_requires_keys(self):
        m = make_module(requires={"pca": "scanpy", "harmony": "scanpy"})
        assert m.requires == {"pca": "scanpy", "harmony": "scanpy"}


# ---------------------------------------------------------------------------
# Unit: Validation
# ---------------------------------------------------------------------------


def _two_stage_benchmark(rapids_requires=None):
    """Helper: pca stage (scanpy + rapids), then harmony stage (scanpy + rapids).

    rapids module in harmony can be given `requires` via `rapids_requires`.
    """
    return make_benchmark(
        stages=[
            {
                "id": "pca",
                "modules": [
                    {"id": "scanpy", "software_environment": "python_env"},
                    {"id": "rapids", "software_environment": "python_env"},
                ],
                "outputs": [{"id": "pca.csv", "path": "pca.csv.gz"}],
                "provides": {"pca": "pca.csv"},
            },
            {
                "id": "harmony",
                "inputs": [["pca.csv"]],
                "modules": [
                    {"id": "scanpy", "software_environment": "python_env"},
                    {
                        "id": "rapids",
                        "software_environment": "python_env",
                        **({"requires": rapids_requires} if rapids_requires else {}),
                    },
                ],
                "outputs": [{"id": "harmony.csv", "path": "harmony.csv.gz"}],
            },
        ]
    )


class TestRequiresValidation:
    @pytest.mark.short
    def test_valid_requires_passes_validation(self):
        """A valid requires declaration should not raise during construction."""
        b = _two_stage_benchmark(rapids_requires={"pca": "scanpy"})
        assert b is not None

    @pytest.mark.short
    def test_requires_unknown_label_raises(self):
        """Requires referencing a label not declared by any ancestor raises."""
        with pytest.raises(ValidationError) as exc_info:
            _two_stage_benchmark(rapids_requires={"nonexistent_label": "scanpy"})
        assert "nonexistent_label" in str(exc_info.value)
        assert "not declared by any ancestor stage" in str(exc_info.value)

    @pytest.mark.short
    def test_requires_unknown_module_id_raises(self):
        """Requires referencing a module that doesn't exist in the providing stage raises."""
        with pytest.raises(ValidationError) as exc_info:
            _two_stage_benchmark(rapids_requires={"pca": "no_such_module"})
        assert "no_such_module" in str(exc_info.value)
        assert "is not a module" in str(exc_info.value)

    @pytest.mark.short
    def test_requires_self_reference_raises(self):
        """Requires referencing a label provided by the current stage raises."""
        with pytest.raises(ValidationError) as exc_info:
            make_benchmark(
                stages=[
                    {
                        "id": "pca",
                        "modules": [
                            {"id": "scanpy", "software_environment": "python_env"},
                            {
                                "id": "rapids",
                                "software_environment": "python_env",
                                "requires": {"pca": "scanpy"},  # self-reference
                            },
                        ],
                        "outputs": [{"id": "pca.csv", "path": "pca.csv.gz"}],
                        "provides": {"pca": "pca.csv"},
                    }
                ]
            )
        assert "self-reference" in str(exc_info.value)
