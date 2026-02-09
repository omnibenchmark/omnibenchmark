"""Tests for {params.key} template substitution in output paths."""

from omnibenchmark.benchmark.params import Params
from omnibenchmark.cli.run import _substitute_params_in_path


class TestSubstituteParamsInPath:
    """Unit tests for _substitute_params_in_path."""

    def test_single_param(self):
        params = Params({"k": "10"})
        result = _substitute_params_in_path("{dataset}_k{params.k}.csv", params)
        assert result == "{dataset}_k10.csv"

    def test_multiple_params(self):
        params = Params({"alpha": "0.5", "k": "10"})
        result = _substitute_params_in_path(
            "result_a{params.alpha}_k{params.k}.csv", params
        )
        assert result == "result_a0.5_k10.csv"

    def test_no_params_placeholder(self):
        params = Params({"k": "10"})
        result = _substitute_params_in_path("{dataset}_result.csv", params)
        assert result == "{dataset}_result.csv"

    def test_none_params(self):
        result = _substitute_params_in_path("{dataset}_k{params.k}.csv", None)
        assert result == "{dataset}_k{params.k}.csv"

    def test_unknown_param_key_left_as_is(self):
        params = Params({"k": "10"})
        result = _substitute_params_in_path("{params.missing}.csv", params)
        assert result == "{params.missing}.csv"

    def test_mixed_known_and_unknown(self):
        params = Params({"k": "10"})
        result = _substitute_params_in_path("{params.k}_{params.missing}.csv", params)
        assert result == "10_{params.missing}.csv"

    def test_empty_template(self):
        params = Params({"k": "10"})
        result = _substitute_params_in_path("", params)
        assert result == ""

    def test_no_braces_at_all(self):
        params = Params({"k": "10"})
        result = _substitute_params_in_path("plain_filename.csv", params)
        assert result == "plain_filename.csv"

    def test_numeric_param_value(self):
        params = Params({"n_clusters": 5})
        result = _substitute_params_in_path("clusters_{params.n_clusters}.csv", params)
        assert result == "clusters_5.csv"

    def test_param_at_start(self):
        params = Params({"method": "pca"})
        result = _substitute_params_in_path("{params.method}_output.csv", params)
        assert result == "pca_output.csv"

    def test_param_at_end(self):
        params = Params({"method": "pca"})
        result = _substitute_params_in_path("output_{params.method}", params)
        assert result == "output_pca"

    def test_dataset_and_params_together(self):
        """Both {dataset} and {params.x} can coexist — params substitution
        handles only {params.*}, leaving {dataset} for the caller."""
        params = Params({"k": "10"})
        result = _substitute_params_in_path("{dataset}_k{params.k}_result.csv", params)
        assert result == "{dataset}_k10_result.csv"

    def test_repeated_same_param(self):
        params = Params({"k": "10"})
        result = _substitute_params_in_path("k{params.k}_again{params.k}.csv", params)
        assert result == "k10_again10.csv"

    def test_underscore_in_param_name(self):
        params = Params({"n_components": "3"})
        result = _substitute_params_in_path("nc{params.n_components}.csv", params)
        assert result == "nc3.csv"
