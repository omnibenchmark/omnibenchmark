from pathlib import Path
from typing import Dict
import json
import os

from omni.benchmark import Benchmark
from tests.workflow.Snakemake_setup import SnakemakeSetup


def test_run_workflow_001():
    benchmark_file = Path("..") / "data" / "Benchmark_001.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    with SnakemakeSetup(benchmark_file_path) as setup:
        benchmark = setup.benchmark
        assert benchmark.get_benchmark_name() == "Benchmark_001"

        # First run the whole workflow
        success = setup.workflow.run_workflow(benchmark)

        assert success


def test_run_workflow_backends_missing():
    benchmark_file = Path("..") / "data" / "benchmark_some_backends_missing.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    with SnakemakeSetup(benchmark_file_path) as setup:
        benchmark = setup.benchmark
        assert benchmark.get_benchmark_name() == "some_backends_missing"

        # First run the whole workflow
        success = setup.workflow.run_workflow(benchmark)

        assert success


def test_single_run_workflow_with_parameters():
    benchmark_file = Path("..") / "data" / "Clustering.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    expected_D1_param_dict = {
        "param7303dd4b11fc73cda647ecb3f7804e697c7263a96b17f866b3d33b3c61c7616a": {
            "measure": "chebyshev"
        },
        "param35923507f72ce06b1717bc24a2ed9be2b402465c4216083ea2649dacf4e7cfa3": {
            "measure": "cosine"
        },
        "param17e490053067e42124d830017cf74c8831676bc88fdb744769d12d71b0d91b51": {
            "measure": "manhattan"
        },
        "paramfbc351b86dafa7dc3e6f57558c75811033e0bfe4800932ad7760bbbb3c18853f": {
            "measure": "euclidean"
        },
    }

    with SnakemakeSetup(benchmark_file_path) as setup:
        benchmark = setup.benchmark

        # assert there are 2 datasets in the benchmark
        data_nodes = benchmark.get_nodes_by_stage_id("data")
        datasets = set([node.module_id for node in data_nodes])
        assert len(datasets) == 2

        # assert there are 4 parameters for module D1
        D1_nodes = benchmark.get_nodes_by_module_id("D1")
        D1_params_unique = set([node.param_id for node in D1_nodes])
        assert len(D1_params_unique) == 4

        # assert benchmark run is successful
        success = setup.workflow.run_workflow(benchmark)
        assert success

        # for each dataset, assert the parameter serialization is correct
        for dataset in datasets:
            D1_output_folder_path = (
                Path(os.getcwd())
                / "out"
                / "data"
                / dataset
                / "default"
                / "distances"
                / "D1"
            )
            _assert_parameters_output_for_module(
                D1_output_folder_path, expected_D1_param_dict
            )


def test_multi_run_workflow_with_parameters():
    benchmark_file = Path("..") / "data" / "Clustering.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    D1_output_folder_path = (
        Path(os.getcwd()) / "out" / "data" / "iris" / "default" / "distances" / "D1"
    )

    expected_D1_param_dict = {
        "param7303dd4b11fc73cda647ecb3f7804e697c7263a96b17f866b3d33b3c61c7616a": {
            "measure": "chebyshev"
        },
        "param35923507f72ce06b1717bc24a2ed9be2b402465c4216083ea2649dacf4e7cfa3": {
            "measure": "cosine"
        },
        "param17e490053067e42124d830017cf74c8831676bc88fdb744769d12d71b0d91b51": {
            "measure": "manhattan"
        },
        "paramfbc351b86dafa7dc3e6f57558c75811033e0bfe4800932ad7760bbbb3c18853f": {
            "measure": "euclidean"
        },
    }

    with SnakemakeSetup(benchmark_file_path) as setup:
        benchmark = setup.benchmark

        # assert benchmark 1st run is successful
        success = setup.workflow.run_workflow(benchmark)
        assert success

        # assert the parameter serialization is correct after 1st run
        _assert_parameters_output_for_module(
            D1_output_folder_path, expected_D1_param_dict
        )

        # assert benchmark 2nd run is successful
        success = setup.workflow.run_workflow(benchmark)
        assert success

        # assert the parameter serialization is correct after 2nd run
        _assert_parameters_output_for_module(
            D1_output_folder_path, expected_D1_param_dict
        )


def test_multi_run_workflow_with_parameter_removal():
    benchmark_file = Path("..") / "data" / "Clustering.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file
    D1_output_folder_path = (
        Path(os.getcwd()) / "out" / "data" / "iris" / "default" / "distances" / "D1"
    )

    with SnakemakeSetup(benchmark_file_path) as setup:
        benchmark = setup.benchmark

        # assert benchmark 1st run is successful
        success = setup.workflow.run_workflow(benchmark)
        assert success

        expected_param_dict_before_removal = {
            "param7303dd4b11fc73cda647ecb3f7804e697c7263a96b17f866b3d33b3c61c7616a": {
                "measure": "chebyshev"
            },
            "param35923507f72ce06b1717bc24a2ed9be2b402465c4216083ea2649dacf4e7cfa3": {
                "measure": "cosine"
            },
            "param17e490053067e42124d830017cf74c8831676bc88fdb744769d12d71b0d91b51": {
                "measure": "manhattan"
            },
            "paramfbc351b86dafa7dc3e6f57558c75811033e0bfe4800932ad7760bbbb3c18853f": {
                "measure": "euclidean"
            },
        }

        # assert the parameter serialization is correct after 1st run
        _assert_parameters_output_for_module(
            D1_output_folder_path, expected_param_dict_before_removal
        )

        # remove previous parameters dict
        os.remove(D1_output_folder_path / "parameters_dict.txt")

        # assert benchmark 2nd run is successful
        benchmark_file_trimmed = Path("..") / "data" / "Clustering_trimmed.yaml"
        benchmark_file_trimmed_path = Path(__file__).parent / benchmark_file_trimmed
        benchmark_without_param = Benchmark(benchmark_file_trimmed_path)
        success = setup.workflow.run_workflow(benchmark_without_param)
        assert success

        expected_param_dict_after_removal = {
            "param7303dd4b11fc73cda647ecb3f7804e697c7263a96b17f866b3d33b3c61c7616a": {
                "measure": "chebyshev"
            },
            "param17e490053067e42124d830017cf74c8831676bc88fdb744769d12d71b0d91b51": {
                "measure": "manhattan"
            },
            "param35923507f72ce06b1717bc24a2ed9be2b402465c4216083ea2649dacf4e7cfa3": {
                "measure": "cosine"
            },
        }

        # assert the parameter serialization is correct after 2nd run
        _assert_parameters_output_for_module(
            D1_output_folder_path, expected_param_dict_after_removal
        )


def _assert_parameters_output_for_module(
    module_output_path: Path, expected_param_dict: Dict[str, any]
):
    module_param_dict = module_output_path / "parameters_dict.txt"
    assert module_param_dict.exists(), f"File not found: {module_param_dict}"

    with module_param_dict.open("r") as f:
        param_dict_lines = [line.strip() for line in f if line.strip()]

    actual_param_dict = {}
    for line in param_dict_lines:
        key, value = line.split(" ", 1)
        actual_param_dict[key] = json.loads(value)

    assert (
        actual_param_dict == expected_param_dict
    ), f"Mismatch in parameters:\nExpected:\n{expected_param_dict}\n\nGot:\n{actual_param_dict}"

    # check that the content of the subfolder actually contains the parameters specified in the dict
    for param_folder in actual_param_dict.keys():
        module_param_file_path = module_output_path / param_folder / "parameters.txt"
        assert (
            module_param_file_path.exists()
        ), f"File not found: {module_param_file_path}"

        with module_param_file_path.open("r") as f:
            param_txt_content = json.loads(f.read())

        expected_content = actual_param_dict.get(param_folder)
        assert (
            expected_content == param_txt_content
        ), f"Mismatch in {module_param_file_path}:\nExpected: {expected_content}\nGot: {param_txt_content}"
