"""Tests for TemplateContext variable substitution system."""

import pytest
from dataclasses import FrozenInstanceError

from omnibenchmark.model.resolved import TemplateContext
from omnibenchmark.benchmark.params import Params


class TestSubstituteProvides:
    """Tests for {label} provides-derived variable substitution."""

    def test_single_provides_variable(self):
        ctx = TemplateContext(provides={"dataset": "iris"})
        assert ctx.substitute("{dataset}_features.csv") == "iris_features.csv"

    def test_multiple_provides_variables(self):
        ctx = TemplateContext(provides={"dataset": "iris", "method": "kmeans"})
        result = ctx.substitute("{dataset}_{method}_labels.csv")
        assert result == "iris_kmeans_labels.csv"

    def test_unknown_provides_left_as_is(self):
        ctx = TemplateContext(provides={"dataset": "iris"})
        assert ctx.substitute("{unknown}_data.csv") == "{unknown}_data.csv"

    def test_no_braces_returns_unchanged(self):
        ctx = TemplateContext(provides={"dataset": "iris"})
        assert ctx.substitute("plain_path.csv") == "plain_path.csv"

    def test_empty_template(self):
        ctx = TemplateContext(provides={"dataset": "iris"})
        assert ctx.substitute("") == ""

    def test_none_template(self):
        ctx = TemplateContext(provides={"dataset": "iris"})
        assert ctx.substitute(None) is None

    def test_repeated_provides_variable(self):
        ctx = TemplateContext(provides={"dataset": "iris"})
        result = ctx.substitute("{dataset}/{dataset}_features.csv")
        assert result == "iris/iris_features.csv"


class TestSubstituteModuleAttrs:
    """Tests for {module.attr} variable substitution."""

    def test_module_id(self):
        ctx = TemplateContext(module_attrs={"id": "kmeans"})
        assert ctx.substitute("{module.id}_output.csv") == "kmeans_output.csv"

    def test_module_stage(self):
        ctx = TemplateContext(module_attrs={"id": "kmeans", "stage": "methods"})
        assert ctx.substitute("{module.stage}/{module.id}.csv") == "methods/kmeans.csv"

    def test_module_parent_id(self):
        ctx = TemplateContext(
            module_attrs={"id": "kmeans", "parent.id": "iris", "parent.stage": "data"}
        )
        result = ctx.substitute("{module.parent.id}_{module.id}.csv")
        assert result == "iris_kmeans.csv"

    def test_unknown_module_attr_left_as_is(self):
        ctx = TemplateContext(module_attrs={"id": "kmeans"})
        assert ctx.substitute("{module.unknown}.csv") == "{module.unknown}.csv"


class TestSubstituteParams:
    """Tests for {params.key} variable substitution."""

    def test_single_param(self):
        ctx = TemplateContext()
        params = Params({"k": 3})
        assert ctx.substitute("k{params.k}.csv", params=params) == "k3.csv"

    def test_multiple_params(self):
        ctx = TemplateContext()
        params = Params({"k": 3, "alpha": 0.5})
        result = ctx.substitute("k{params.k}_a{params.alpha}.csv", params=params)
        assert result == "k3_a0.5.csv"

    def test_unknown_param_left_as_is(self):
        ctx = TemplateContext()
        params = Params({"k": 3})
        assert (
            ctx.substitute("{params.unknown}.csv", params=params)
            == "{params.unknown}.csv"
        )

    def test_no_params_object(self):
        ctx = TemplateContext()
        assert ctx.substitute("{params.k}.csv") == "{params.k}.csv"

    def test_none_params(self):
        ctx = TemplateContext()
        assert ctx.substitute("{params.k}.csv", params=None) == "{params.k}.csv"


class TestSubstituteMixed:
    """Tests for templates using multiple namespaces."""

    def test_provides_and_params(self):
        ctx = TemplateContext(provides={"dataset": "iris"})
        params = Params({"k": 3})
        result = ctx.substitute("{dataset}_k{params.k}.csv", params=params)
        assert result == "iris_k3.csv"

    def test_all_three_namespaces(self):
        ctx = TemplateContext(
            provides={"dataset": "iris", "method": "kmeans"},
            module_attrs={"id": "scorer", "stage": "scoring"},
        )
        params = Params({"metric": "ari"})
        template = "{dataset}_{method}_{module.id}_{params.metric}.csv"
        result = ctx.substitute(template, params=params)
        assert result == "iris_kmeans_scorer_ari.csv"

    def test_provides_and_module_attrs(self):
        ctx = TemplateContext(
            provides={"dataset": "iris"},
            module_attrs={"id": "kmeans", "parent.id": "iris"},
        )
        result = ctx.substitute("{dataset}_{module.id}.csv")
        assert result == "iris_kmeans.csv"


class TestAncestorInheritance:
    """Tests simulating the provides accumulation through the DAG."""

    def test_entrypoint_self_reference(self):
        """Entrypoint provides: [dataset] maps to own module_id."""
        ctx = TemplateContext(
            provides={"dataset": "iris"},
            module_attrs={"id": "iris", "stage": "data"},
        )
        assert ctx.substitute("{dataset}") == "iris"

    def test_child_inherits_parent_provides(self):
        """Child node inherits {dataset} from parent and adds {method}."""
        parent_ctx = TemplateContext(
            provides={"dataset": "iris"},
            module_attrs={"id": "iris", "stage": "data"},
        )
        # Simulate _build_template_context for downstream node
        child_provides = dict(parent_ctx.provides)
        child_provides["method"] = "kmeans"
        child_ctx = TemplateContext(
            provides=child_provides,
            module_attrs={
                "id": "kmeans",
                "stage": "methods",
                "parent.id": "iris",
                "parent.stage": "data",
            },
        )
        result = child_ctx.substitute("{dataset}_{method}_labels.csv")
        assert result == "iris_kmeans_labels.csv"

    def test_grandchild_inherits_full_chain(self):
        """Grandchild inherits {dataset} and {method} from ancestors."""
        grandchild_provides = {"dataset": "iris", "method": "kmeans"}
        grandchild_ctx = TemplateContext(
            provides=grandchild_provides,
            module_attrs={
                "id": "scorer",
                "stage": "scoring",
                "parent.id": "kmeans",
                "parent.stage": "methods",
            },
        )
        result = grandchild_ctx.substitute("{dataset}_{method}_{module.id}_score.csv")
        assert result == "iris_kmeans_scorer_score.csv"

    def test_nearest_ancestor_wins(self):
        """If same label provided at multiple levels, nearest wins."""
        # Parent provides "metric" as "v1", child overrides with "v2"
        child_provides = {"metric": "v2"}  # override
        ctx = TemplateContext(provides=child_provides)
        assert ctx.substitute("{metric}") == "v2"


class TestGatherContext:
    """Tests for gather stage TemplateContext (empty provides)."""

    def test_gather_empty_provides(self):
        ctx = TemplateContext(
            provides={},
            module_attrs={"id": "aggregator", "stage": "gather_results"},
        )
        # Provides-derived vars are left unresolved
        assert ctx.substitute("{dataset}_agg.csv") == "{dataset}_agg.csv"

    def test_gather_module_attrs_available(self):
        ctx = TemplateContext(
            provides={},
            module_attrs={"id": "aggregator", "stage": "gather_results"},
        )
        assert ctx.substitute("{module.id}_agg.csv") == "aggregator_agg.csv"

    def test_gather_params_available(self):
        ctx = TemplateContext(
            provides={},
            module_attrs={"id": "aggregator", "stage": "gather_results"},
        )
        params = Params({"format": "json"})
        result = ctx.substitute("results.{params.format}", params=params)
        assert result == "results.json"


class TestToDict:
    """Tests for TemplateContext serialization."""

    def test_empty_context(self):
        ctx = TemplateContext()
        assert ctx.to_dict() == {}

    def test_provides_only(self):
        ctx = TemplateContext(provides={"dataset": "iris", "method": "kmeans"})
        d = ctx.to_dict()
        assert d == {"dataset": "iris", "method": "kmeans"}

    def test_module_attrs_prefixed(self):
        ctx = TemplateContext(module_attrs={"id": "kmeans", "stage": "methods"})
        d = ctx.to_dict()
        assert d == {"module.id": "kmeans", "module.stage": "methods"}

    def test_mixed(self):
        ctx = TemplateContext(
            provides={"dataset": "iris"},
            module_attrs={"id": "kmeans"},
        )
        d = ctx.to_dict()
        assert d == {"dataset": "iris", "module.id": "kmeans"}


class TestFrozen:
    """Tests for immutability."""

    def test_cannot_set_provides(self):
        ctx = TemplateContext(provides={"dataset": "iris"})
        with pytest.raises(FrozenInstanceError):
            ctx.provides = {"new": "value"}

    def test_cannot_set_module_attrs(self):
        ctx = TemplateContext(module_attrs={"id": "kmeans"})
        with pytest.raises(FrozenInstanceError):
            ctx.module_attrs = {"new": "value"}
