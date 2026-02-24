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


# ---------------------------------------------------------------------------
# Unit: Ambiguous label detection (§9.5 uniqueness constraint)
# ---------------------------------------------------------------------------


def _ambiguous_label_benchmark(requires_label: str = "method"):
    """Three-stage benchmark where two stages both provide the same label.

    stages: methods_fast (provides 'method') → methods_accurate (provides 'method')
            → postprocess (rapids requires 'method': 'M1')

    The label 'method' is provided by two ancestor stages, so requires is
    ambiguous and must be rejected.
    """
    return make_benchmark(
        stages=[
            {
                "id": "methods_fast",
                "provides": {"method": "fast.result"},
                "modules": [{"id": "M1", "software_environment": "python_env"}],
                "outputs": [{"id": "fast.result", "path": "fast.csv"}],
            },
            {
                "id": "methods_accurate",
                "provides": {"method": "accurate.result"},
                "modules": [{"id": "M2", "software_environment": "python_env"}],
                "outputs": [{"id": "accurate.result", "path": "accurate.csv"}],
            },
            {
                "id": "postprocess",
                "inputs": [["fast.result"]],
                "modules": [
                    {"id": "scanpy", "software_environment": "python_env"},
                    {
                        "id": "rapids",
                        "software_environment": "python_env",
                        "requires": {requires_label: "M1"},
                    },
                ],
                "outputs": [{"id": "post.result", "path": "post.csv"}],
            },
        ]
    )


class TestRequiresAmbiguousLabel:
    @pytest.mark.short
    def test_ambiguous_label_raises(self):
        """requires referencing a label provided by more than one ancestor is an error."""
        with pytest.raises(ValidationError) as exc_info:
            _ambiguous_label_benchmark(requires_label="method")
        msg = str(exc_info.value)
        assert "ambiguous" in msg
        assert "method" in msg

    @pytest.mark.short
    def test_ambiguous_label_error_names_all_providers(self):
        """The error message must list every stage that provides the conflicting label."""
        with pytest.raises(ValidationError) as exc_info:
            _ambiguous_label_benchmark(requires_label="method")
        msg = str(exc_info.value)
        assert "methods_fast" in msg
        assert "methods_accurate" in msg

    @pytest.mark.short
    def test_ambiguous_label_error_names_affected_module(self):
        """The error message must identify the module that declared the bad requires."""
        with pytest.raises(ValidationError) as exc_info:
            _ambiguous_label_benchmark(requires_label="method")
        assert "rapids" in str(exc_info.value)

    @pytest.mark.short
    def test_unique_label_across_stages_passes(self):
        """When each stage provides a distinct label, requires is unambiguous and valid."""
        b = make_benchmark(
            stages=[
                {
                    "id": "methods_fast",
                    "provides": {"method_fast": "fast.result"},
                    "modules": [{"id": "M1", "software_environment": "python_env"}],
                    "outputs": [{"id": "fast.result", "path": "fast.csv"}],
                },
                {
                    "id": "methods_accurate",
                    "provides": {"method_accurate": "accurate.result"},
                    "modules": [{"id": "M2", "software_environment": "python_env"}],
                    "outputs": [{"id": "accurate.result", "path": "accurate.csv"}],
                },
                {
                    "id": "postprocess",
                    "inputs": [["fast.result"]],
                    "modules": [
                        {"id": "scanpy", "software_environment": "python_env"},
                        {
                            "id": "rapids",
                            "software_environment": "python_env",
                            "requires": {"method_fast": "M1"},  # unambiguous
                        },
                    ],
                    "outputs": [{"id": "post.result", "path": "post.csv"}],
                },
            ]
        )
        assert b is not None

    @pytest.mark.short
    def test_gather_with_shared_label_is_unaffected(self):
        """gather: is allowed on shared labels; the ambiguity rule applies only to requires."""
        from omnibenchmark.model.benchmark import GatherInput

        b = make_benchmark(
            stages=[
                {
                    "id": "methods_fast",
                    "provides": {"method": "fast.result"},
                    "modules": [{"id": "M1", "software_environment": "python_env"}],
                    "outputs": [{"id": "fast.result", "path": "fast.csv"}],
                },
                {
                    "id": "methods_accurate",
                    "provides": {"method": "accurate.result"},
                    "modules": [{"id": "M2", "software_environment": "python_env"}],
                    "outputs": [{"id": "accurate.result", "path": "accurate.csv"}],
                },
                {
                    "id": "summary",
                    "inputs": [GatherInput(gather="method")],
                    "modules": [
                        {"id": "reporter", "software_environment": "python_env"}
                    ],
                    "outputs": [{"id": "report", "path": "report.json"}],
                },
            ]
        )
        assert b is not None
