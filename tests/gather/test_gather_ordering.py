"""Tests for gather stage ordering validation (Phase 2).

Tests that _validate_gather_ordering correctly enforces that all
provider stages appear before their gather consumers in document order.
"""

import pytest

from omnibenchmark.model.benchmark import (
    GatherInput,
    InputCollection,
    Stage,
    Module,
    IOFile,
    Repository,
)
from omnibenchmark.cli.run import _validate_gather_ordering, GatherOrderingError


# ============================================================================
# Factories (minimal, same as test_gather_model.py)
# ============================================================================


def _repo(**kwargs):
    defaults = {"url": "https://example.com/repo.git", "commit": "abc123"}
    return Repository(**{**defaults, **kwargs})


def _module(**kwargs):
    defaults = {
        "id": "mod1",
        "software_environment": "host",
        "repository": _repo(),
    }
    return Module(**{**defaults, **kwargs})


def _iofile(**kwargs):
    defaults = {"id": "out1", "path": "output.txt"}
    return IOFile(**{**defaults, **kwargs})


def _stage(**kwargs):
    defaults = {
        "id": "stage1",
        "modules": [_module()],
        "outputs": [_iofile()],
    }
    return Stage(**{**defaults, **kwargs})


# ============================================================================
# Valid orderings (should not raise)
# ============================================================================


@pytest.mark.short
class TestValidGatherOrdering:
    def test_provider_before_gather(self):
        """Provider appears before gather — valid."""
        stages = [
            _stage(id="data", provides=["dataset"]),
            _stage(
                id="summary",
                inputs=[GatherInput(gather="dataset")],
                outputs=[_iofile(id="summary.report", path="report.json")],
            ),
        ]
        _validate_gather_ordering(stages)  # Should not raise

    def test_multiple_providers_before_gather(self):
        """Multiple providers all appear before the gather — valid."""
        stages = [
            _stage(id="methods_fast", provides=["method"]),
            _stage(id="methods_slow", provides=["method"]),
            _stage(
                id="summary",
                inputs=[GatherInput(gather="method")],
                outputs=[_iofile(id="summary.report", path="report.json")],
            ),
        ]
        _validate_gather_ordering(stages)

    def test_no_gather_stages(self):
        """Pipeline with no gather stages — trivially valid."""
        stages = [
            _stage(id="data"),
            _stage(
                id="methods",
                inputs=[InputCollection(entries=["data.raw"])],
            ),
        ]
        _validate_gather_ordering(stages)

    def test_provides_without_gather(self):
        """Provides with no gather consumer — valid (unused provides is ok)."""
        stages = [
            _stage(id="data", provides=["dataset"]),
            _stage(id="methods"),
        ]
        _validate_gather_ordering(stages)

    def test_empty_stages(self):
        """Empty stages list — trivially valid."""
        _validate_gather_ordering([])

    def test_provider_immediately_before_gather(self):
        """Provider is the stage immediately preceding gather — valid."""
        stages = [
            _stage(id="data", provides=["dataset"]),
            _stage(
                id="gather",
                inputs=[GatherInput(gather="dataset")],
                outputs=[_iofile(id="g.out", path="out.json")],
            ),
        ]
        _validate_gather_ordering(stages)

    def test_gather_with_gap_between_provider(self):
        """Non-provider stages between provider and gather — valid."""
        stages = [
            _stage(id="data", provides=["dataset"]),
            _stage(id="preprocessing"),
            _stage(id="methods"),
            _stage(
                id="summary",
                inputs=[GatherInput(gather="dataset")],
                outputs=[_iofile(id="s.out", path="out.json")],
            ),
        ]
        _validate_gather_ordering(stages)

    def test_chained_gathers(self):
        """Gather -> provides -> gather chain — valid if ordered correctly."""
        stages = [
            _stage(
                id="methods",
                provides=["method"],
                outputs=[_iofile(id="m.result", path="result.json")],
            ),
            _stage(
                id="first_gather",
                provides=["intermediate"],
                inputs=[GatherInput(gather="method")],
                outputs=[_iofile(id="fg.out", path="intermediate.json")],
            ),
            _stage(
                id="second_gather",
                inputs=[GatherInput(gather="intermediate")],
                outputs=[_iofile(id="sg.out", path="final.json")],
            ),
        ]
        _validate_gather_ordering(stages)


# ============================================================================
# Invalid orderings (should raise GatherOrderingError)
# ============================================================================


@pytest.mark.short
class TestInvalidGatherOrdering:
    def test_provider_after_gather(self):
        """Provider appears after gather — invalid."""
        stages = [
            _stage(
                id="summary",
                inputs=[GatherInput(gather="dataset")],
                outputs=[_iofile(id="s.out", path="report.json")],
            ),
            _stage(id="data", provides=["dataset"]),
        ]
        with pytest.raises(GatherOrderingError, match="appears after it"):
            _validate_gather_ordering(stages)

    def test_error_message_includes_stage_names(self):
        """Error message mentions both the gather stage and the provider."""
        stages = [
            _stage(
                id="my_summary",
                inputs=[GatherInput(gather="method")],
                outputs=[_iofile(id="s.out", path="report.json")],
            ),
            _stage(id="my_methods", provides=["method"]),
        ]
        with pytest.raises(GatherOrderingError) as exc_info:
            _validate_gather_ordering(stages)
        msg = str(exc_info.value)
        assert "my_summary" in msg
        assert "my_methods" in msg
        assert "method" in msg

    def test_self_referencing_gather_rejected(self):
        """A stage that both provides and gathers the same label — invalid."""
        stages = [
            _stage(
                id="self_ref",
                provides=["dataset"],
                inputs=[GatherInput(gather="dataset")],
                outputs=[_iofile(id="s.out", path="out.json")],
            ),
        ]
        with pytest.raises(GatherOrderingError, match="same stage"):
            _validate_gather_ordering(stages)

    def test_one_provider_before_one_after(self):
        """One provider before, one after — error only for the late one."""
        stages = [
            _stage(id="methods_fast", provides=["method"]),
            _stage(
                id="summary",
                inputs=[GatherInput(gather="method")],
                outputs=[_iofile(id="s.out", path="report.json")],
            ),
            _stage(id="methods_slow", provides=["method"]),
        ]
        with pytest.raises(GatherOrderingError) as exc_info:
            _validate_gather_ordering(stages)
        msg = str(exc_info.value)
        assert "methods_slow" in msg
        # methods_fast is fine, should not appear in error
        assert "methods_fast" not in msg

    def test_multiple_errors_collected(self):
        """Multiple ordering violations produce multiple error messages."""
        stages = [
            _stage(
                id="gather1",
                inputs=[GatherInput(gather="a")],
                outputs=[_iofile(id="g1.out", path="g1.json")],
            ),
            _stage(
                id="gather2",
                inputs=[GatherInput(gather="b")],
                outputs=[_iofile(id="g2.out", path="g2.json")],
            ),
            _stage(id="provider_a", provides=["a"]),
            _stage(id="provider_b", provides=["b"]),
        ]
        with pytest.raises(GatherOrderingError) as exc_info:
            _validate_gather_ordering(stages)
        assert len(exc_info.value.errors) == 2
