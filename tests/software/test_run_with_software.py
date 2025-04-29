import json
import os
import shutil
import functools

from pathlib import Path

from typing import Dict
import pandas as pd
import pytest

from omni_schema.datamodel.omni_schema import SoftwareBackendEnum

from omnibenchmark.benchmark import Benchmark
from omnibenchmark.workflow.snakemake import SnakemakeEngine


def skip_if_no_conda(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        if shutil.which("conda") is None:
            pytest.skip("conda not found in system path")
        return func(*args, **kwargs)

    return wrapper


def skip_if_no_envmodules(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        if shutil.which("module") is None:
            pytest.skip("envmodules not found in system path")
        return func(*args, **kwargs)

    return wrapper


benchmark_data = Path("..") / "data"
benchmark_data_path = Path(__file__).parent / benchmark_data


# TODO: use tmp_path, no need to cleanup
TO_CLEANUP = [".snakemake", "out", "Snakefile", "snakemake.log"]

# TODO(ben): reuse flags in tests/io to inject --current-dir and --keep-files


def _cleanup_snakemake():
    current_dir = os.getcwd()
    for file in TO_CLEANUP:
        file_path = os.path.join(current_dir, file)
        if os.path.exists(file_path):
            if os.path.isfile(file_path):
                os.remove(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path, onerror=_remove_readonly)


def _remove_readonly(func, path, _):
    """Clear the readonly bit and reattempt the removal"""
    os.chmod(path, os.stat.S_IWRITE)
    func(path)


@skip_if_no_envmodules
def test_run_benchmark_with_software_envmodules():
    benchmark_file = benchmark_data_path / "Clustering.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    if os.path.exists(benchmark_file):
        benchmark = Benchmark(benchmark_file_path)
        workflow = SnakemakeEngine()

        if os.environ.get("GITHUB_WORKSPACE", None):
            override_module_path = os.path.join(
                os.environ["GITHUB_WORKSPACE"],
                "tests",
                "data",
                "envs",
            )
        else:
            override_module_path = os.environ.get("MODULEPATH", None)

        success = workflow.run_workflow(
            benchmark,
            backend=SoftwareBackendEnum.envmodules,
            module_path=override_module_path,
        )

        assert success

        # TODO: use tmp_path
        _cleanup_snakemake()
    else:
        raise FileNotFoundError(f"Benchmark file {benchmark_file} does not exist.")


@skip_if_no_conda
def test_run_benchmark_with_software_conda():
    benchmark_file = benchmark_data_path / "Clustering.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    if os.path.exists(benchmark_file):
        benchmark = Benchmark(benchmark_file_path)
        workflow = SnakemakeEngine()

        success = workflow.run_workflow(
            benchmark,
            backend=SoftwareBackendEnum.conda,
        )

        assert success

        # TODO: use tmp_path
        _cleanup_snakemake()
    else:
        raise FileNotFoundError(f"Benchmark file {benchmark_file} does not exist.")


@pytest.mark.skip(reason="Apptainer image is not available yet")
def test_run_benchmark_with_software_apptainer():
    benchmark_file = benchmark_data_path / "Clustering.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    if os.path.exists(benchmark_file):
        benchmark = Benchmark(benchmark_file_path)
        workflow = SnakemakeEngine()

        success = workflow.run_workflow(
            benchmark,
            backend=SoftwareBackendEnum.apptainer,
        )

        assert success

        # TODO: use tmp_path
        _cleanup_snakemake()
    else:
        raise FileNotFoundError(f"Benchmark file {benchmark_file} does not exist.")


@skip_if_no_conda
def test_single_run_workflow_with_parameters():
    benchmark_file = Path("..") / "data" / "Clustering.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    expected_D1_param_dict = {
        "7303dd4b11fc73cda647ecb3f7804e697c7263a96b17f866b3d33b3c61c7616a": {
            "measure": "chebyshev"
        },
        "35923507f72ce06b1717bc24a2ed9be2b402465c4216083ea2649dacf4e7cfa3": {
            "measure": "cosine"
        },
        "17e490053067e42124d830017cf74c8831676bc88fdb744769d12d71b0d91b51": {
            "measure": "manhattan"
        },
        "fbc351b86dafa7dc3e6f57558c75811033e0bfe4800932ad7760bbbb3c18853f": {
            "measure": "euclidean"
        },
    }

    if os.path.exists(benchmark_file):
        benchmark = Benchmark(benchmark_file_path)
        workflow = SnakemakeEngine()

        # assert there are 2 datasets in the benchmark
        data_nodes = benchmark.get_nodes_by_stage_id("data")
        datasets = set([node.module_id for node in data_nodes])
        assert len(datasets) == 2

        # assert there are 4 parameters for module D1
        D1_nodes = benchmark.get_nodes_by_module_id("D1")
        D1_params_unique = set([node.param_id for node in D1_nodes])
        assert len(D1_params_unique) == 4

        # assert benchmark run is successful
        success = workflow.run_workflow(benchmark)
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

        _cleanup_snakemake()


@skip_if_no_conda
def test_multi_run_workflow_with_parameters():
    benchmark_file = Path("..") / "data" / "Clustering.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    D1_output_folder_path = (
        Path(os.getcwd()) / "out" / "data" / "iris" / "default" / "distances" / "D1"
    )

    expected_D1_param_dict = {
        "7303dd4b11fc73cda647ecb3f7804e697c7263a96b17f866b3d33b3c61c7616a": {
            "measure": "chebyshev"
        },
        "35923507f72ce06b1717bc24a2ed9be2b402465c4216083ea2649dacf4e7cfa3": {
            "measure": "cosine"
        },
        "17e490053067e42124d830017cf74c8831676bc88fdb744769d12d71b0d91b51": {
            "measure": "manhattan"
        },
        "fbc351b86dafa7dc3e6f57558c75811033e0bfe4800932ad7760bbbb3c18853f": {
            "measure": "euclidean"
        },
    }

    if os.path.exists(benchmark_file):
        benchmark = Benchmark(benchmark_file_path)
        workflow = SnakemakeEngine()

        # assert benchmark 1st run is successful
        success = workflow.run_workflow(benchmark, backend=SoftwareBackendEnum.conda)
        assert success

        # assert the parameter serialization is correct after 1st run
        _assert_parameters_output_for_module(
            D1_output_folder_path, expected_D1_param_dict
        )

        # assert benchmark 2nd run is successful
        success = workflow.run_workflow(benchmark, backend=SoftwareBackendEnum.conda)
        assert success

        # assert the parameter serialization is correct after 2nd run
        _assert_parameters_output_for_module(
            D1_output_folder_path, expected_D1_param_dict
        )

        _cleanup_snakemake()


@skip_if_no_conda
def test_multi_run_workflow_with_parameter_removal():
    benchmark_file = Path("..") / "data" / "Clustering.yaml"
    benchmark_file_path = Path(__file__).parent / benchmark_file

    D1_output_folder_path = (
        Path(os.getcwd()) / "out" / "data" / "iris" / "default" / "distances" / "D1"
    )

    if os.path.exists(benchmark_file):
        benchmark = Benchmark(benchmark_file_path)
        workflow = SnakemakeEngine()

        # assert benchmark 1st run is successful
        success = workflow.run_workflow(benchmark, backend=SoftwareBackendEnum.conda)
        assert success

        expected_param_dict_before_removal = {
            "7303dd4b11fc73cda647ecb3f7804e697c7263a96b17f866b3d33b3c61c7616a": {
                "measure": "chebyshev"
            },
            "35923507f72ce06b1717bc24a2ed9be2b402465c4216083ea2649dacf4e7cfa3": {
                "measure": "cosine"
            },
            "17e490053067e42124d830017cf74c8831676bc88fdb744769d12d71b0d91b51": {
                "measure": "manhattan"
            },
            "fbc351b86dafa7dc3e6f57558c75811033e0bfe4800932ad7760bbbb3c18853f": {
                "measure": "euclidean"
            },
        }

        # assert the parameter serialization is correct after 1st run
        _assert_parameters_output_for_module(
            D1_output_folder_path, expected_param_dict_before_removal
        )

        # remove previous parameters dict
        os.remove(D1_output_folder_path / "parameters_dict.tsv")

        # assert benchmark 2nd run is successful
        benchmark_file_trimmed = Path("..") / "data" / "Clustering_trimmed.yaml"
        benchmark_file_trimmed_path = Path(__file__).parent / benchmark_file_trimmed
        benchmark_without_param = Benchmark(benchmark_file_trimmed_path)
        success = workflow.run_workflow(
            benchmark_without_param, backend=SoftwareBackendEnum.conda
        )
        assert success

        expected_param_dict_after_removal = {
            "7303dd4b11fc73cda647ecb3f7804e697c7263a96b17f866b3d33b3c61c7616a": {
                "measure": "chebyshev"
            },
            "17e490053067e42124d830017cf74c8831676bc88fdb744769d12d71b0d91b51": {
                "measure": "manhattan"
            },
            "35923507f72ce06b1717bc24a2ed9be2b402465c4216083ea2649dacf4e7cfa3": {
                "measure": "cosine"
            },
        }

        # assert the parameter serialization is correct after 2nd run
        _assert_parameters_output_for_module(
            D1_output_folder_path, expected_param_dict_after_removal
        )

        _cleanup_snakemake()


def _assert_parameters_output_for_module(
    module_output_path: Path, expected_param_dict: Dict[str, any]
):
    module_param_tsv = module_output_path / "parameters_dict.tsv"
    assert module_param_tsv.exists(), f"File not found: {module_param_tsv}"

    module_param_df = pd.read_csv(module_param_tsv, sep="\t")

    # check that the content of the subfolder actually contains the parameters specified in the dict
    for _, module_param_row in module_param_df.iterrows():
        base_path = module_param_row["base_path"]
        base_param_file_path = Path(base_path) / "parameters.json"
        assert (
            base_param_file_path.exists()
        ), f"File not found through base path: {base_param_file_path}"

        alias_path = module_param_row["base_path"]
        alias_param_file_path = Path(alias_path) / "parameters.json"
        assert (
            alias_param_file_path.exists()
        ), f"File not found through alias: {alias_param_file_path}"

        param_hash = module_param_row["id"]
        expected_content = expected_param_dict.get(param_hash)
        with open(base_param_file_path, "r") as file:
            actual_content = json.load(file)

        assert (
            expected_content == actual_content
        ), f"Mismatch in {base_param_file_path}:\nExpected: {expected_content}\nGot: {actual_content}"
