"""Short unit tests for the shared lineage-exclusion predicate.

``is_lineage_excluded`` is the single source of truth used by both
execution-path pruning (``_graph.exclude_paths``) and Snakefile generation
(``cli/run.py``); these tests pin its semantics directly.
"""

import pytest
from unittest.mock import MagicMock

from omnibenchmark.core._paths import is_lineage_excluded
from omnibenchmark.core._graph import exclude_paths


@pytest.mark.short
class TestIsLineageExcluded:
    def test_no_exclusions_keeps_everything(self):
        assert is_lineage_excluded({"adamson", "prep", "additive"}, {}) is False

    def test_excludes_adjacent_pair(self):
        # Minimal case: declaring and excluded module on the same path.
        assert (
            is_lineage_excluded({"adamson", "additive"}, {"adamson": ["additive"]})
            is True
        )

    def test_non_adjacent_endpoints_still_excluded(self):
        # Declaring and excluded modules are several stages apart, with
        # intervening modules — the exact case the immediate-predecessor
        # check missed.
        assert (
            is_lineage_excluded(
                {"adamson", "prep", "sim", "additive"}, {"adamson": ["additive"]}
            )
            is True
        )

    def test_excluded_module_absent_keeps_lineage(self):
        assert (
            is_lineage_excluded(
                {"adamson", "prep", "double"}, {"adamson": ["additive"]}
            )
            is False
        )

    def test_declaring_module_absent_keeps_lineage(self):
        # additive only appears with a different dataset that has no exclusions.
        assert (
            is_lineage_excluded(
                {"norman", "prep", "additive"}, {"adamson": ["additive"]}
            )
            is False
        )

    def test_symmetric_declaration(self):
        # Declaring the rule on either endpoint prunes the identical lineage:
        # the predicate is pure set co-occurrence, direction-agnostic.
        lineage = {"adamson", "prep", "additive"}
        forward = is_lineage_excluded(lineage, {"adamson": ["additive"]})
        reverse = is_lineage_excluded(lineage, {"additive": ["adamson"]})
        assert forward is reverse is True

    def test_unknown_id_in_list_is_a_noop(self):
        # A stage id or a typo is just a module id that never matches.
        assert (
            is_lineage_excluded(
                {"adamson", "prep", "additive"}, {"adamson": ["methods"]}
            )
            is False
        )

    def test_multi_exclude_is_or_not_and(self):
        # `adamson: [additive, multiplicative]` is two independent rules:
        # co-occurrence with *either* prunes; neither must be present together.
        rule = {"adamson": ["additive", "multiplicative"]}
        assert is_lineage_excluded({"adamson", "additive"}, rule) is True
        assert is_lineage_excluded({"adamson", "multiplicative"}, rule) is True
        assert (
            is_lineage_excluded({"adamson", "additive", "multiplicative"}, rule) is True
        )
        # Neither excluded module present → kept.
        assert is_lineage_excluded({"adamson", "double"}, rule) is False

    def test_multiple_rules_any_match_excludes(self):
        # Rules from different declaring modules are independent (OR across rules).
        rules = {"adamson": ["additive"], "norman": ["double"]}
        assert is_lineage_excluded({"norman", "double"}, rules) is True
        assert is_lineage_excluded({"adamson", "double"}, rules) is False


def _path(*module_ids):
    return [MagicMock(module_id=m) for m in module_ids]


@pytest.mark.short
class TestExcludePaths:
    def test_prunes_violating_path_across_stages(self):
        keep = _path("norman", "prep", "additive")
        drop = _path("adamson", "prep", "sim", "additive")
        result = exclude_paths([keep, drop], {"adamson": ["additive"]})
        assert result == [keep]

    def test_no_exclusions_returns_all(self):
        paths = [_path("adamson", "additive"), _path("norman", "double")]
        assert exclude_paths(paths, {}) == paths

    def test_symmetric_declaration_prunes_same_paths(self):
        keep = _path("norman", "prep", "additive")
        drop = _path("adamson", "prep", "sim", "additive")
        forward = exclude_paths([keep, drop], {"adamson": ["additive"]})
        reverse = exclude_paths([keep, drop], {"additive": ["adamson"]})
        assert forward == reverse == [keep]
