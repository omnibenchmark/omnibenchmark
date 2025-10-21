import json
from pathlib import Path
import os
import shutil
from typing import Dict, Any

import csv


from omnibenchmark.benchmark import Benchmark
from tests.workflow.Snakemake_setup import SnakemakeSetup

from .path import data
from ..fixtures import bundled_repos  # noqa: F401 - pytest fixture


def pytest_adoption(parser):
    parser.addoption(
        "--keep-files",
        action="store_true",
        default=False,
        help="Keep temporary files after test execution",
    )
    parser.addoption(
        "--current-dir",
        action="store_true",
        default=False,
        help="Use current directory instead of temporary one",
    )


def test_run_workflow_001(snakemake_env, tmp_path):
    # Use current directory if specified, otherwise use tmp_path
    # tmp_path is already provided by pytest fixture if not using current dir
    if snakemake_env["current_dir"]:
        tmp_path = Path(os.getcwd())

    benchmark_file = data / "Benchmark_001.yaml"
    keep_files = snakemake_env["keep_files"]

    # Copy benchmark YAML to tmp_path so workflow can find it
    benchmark_file_in_tmp = tmp_path / "Benchmark_001.yaml"
    with open(benchmark_file, "rb") as src, open(benchmark_file_in_tmp, "wb") as dst:
        dst.write(src.read())

    # Copy envs directory to tmp_path so environment files are available
    envs_src = data / "envs"
    envs_dst = tmp_path / "envs"
    if envs_src.exists():
        shutil.copytree(envs_src, envs_dst, dirs_exist_ok=True)

    with SnakemakeSetup(
        benchmark_file_in_tmp,
        keep_files=keep_files,
        cwd=tmp_path.as_posix(),
        out_dir=(tmp_path / "out").absolute(),
    ) as setup:
        benchmark = setup.benchmark
        assert benchmark.get_benchmark_name() == "Benchmark_001"

        # Run the whole workflow
        success = setup.workflow.run_workflow(
            benchmark,
            work_dir=tmp_path,
        )

        assert success is True


def test_single_run_workflow_with_parameters(snakemake_env, tmp_path):  # noqa: F811
    # Use current directory if specified, otherwise use tmp_path
    # tmp_path is already provided by pytest fixture if not using current dir
    if snakemake_env["current_dir"]:
        tmp_path = Path(os.getcwd())

    benchmark_file = data / "Benchmark_003.yaml"
    keep_files = snakemake_env["keep_files"]

    # Copy benchmark YAML to tmp_path so workflow can find it
    benchmark_file_in_tmp = tmp_path / "Benchmark_003.yaml"
    with open(benchmark_file, "rb") as src, open(benchmark_file_in_tmp, "wb") as dst:
        dst.write(src.read())

    # Copy envs directory to tmp_path so environment files are available
    envs_src = data / "envs"
    envs_dst = tmp_path / "envs"
    if envs_src.exists():
        shutil.copytree(envs_src, envs_dst, dirs_exist_ok=True)

    # The bundled_repos fixture already creates a symlink to the bundles directory
    # at tmp_path/bundles, so we don't need to copy it manually

    expected_P1_param_dict = {
        "dca5d812e76298d5365cbc90ae29d54f867c0b23b0a3d4098826595b9a023054": {
            "a": "0",
            "b": "0",
        },
        "47fae0c4eff9d60d70f80ba26f299ece136010217b7393c672091bfbd5cecb1c": {
            "a": "1",
            "b": "0.1",
        },
        "2dc1a8f315206ba9aa73fd3af58b824faecc348f86c176c7325d2765151289ce": {
            "a": "2",
            "b": "0.5",
        },
    }

    with SnakemakeSetup(
        benchmark_file_in_tmp,
        keep_files=keep_files,
        cwd=tmp_path.as_posix(),
        out_dir=(tmp_path / "out").absolute(),
    ) as setup:
        benchmark = setup.benchmark

        # assert there are 2 datasets in the benchmark
        data_nodes = benchmark.get_nodes_by_stage_id("data")
        datasets = set([node.module_id for node in data_nodes])
        assert len(datasets) == 2

        # assert there are 4 parameters for module D1
        D1_nodes = benchmark.get_nodes_by_module_id("P1")
        D1_params_unique = set([node.param_id for node in D1_nodes])
        assert len(D1_params_unique) == 3

        success = setup.workflow.run_workflow(benchmark, work_dir=tmp_path)
        assert success is True

        # for each dataset, assert the parameter serialization is correct
        for dataset in datasets:
            P1_output_folder_path = (
                tmp_path / "out" / "data" / dataset / "default" / "process" / "P1"
            )

            assert_parameters_output_for_module(
                tmp_path, P1_output_folder_path, expected_P1_param_dict
            )


def test_multi_run_workflow_with_parameters(snakemake_env, tmp_path, bundled_repos):  # noqa: F811
    # Use current directory if specified, otherwise use tmp_path
    # tmp_path is already provided by pytest fixture if not using current dir
    if snakemake_env["current_dir"]:
        tmp_path = Path(os.getcwd())

    benchmark_file = data / "Benchmark_003.yaml"
    keep_files = snakemake_env["keep_files"]

    # Copy benchmark YAML to tmp_path so workflow can find it
    benchmark_file_in_tmp = tmp_path / "Benchmark_003.yaml"
    with open(benchmark_file, "rb") as src, open(benchmark_file_in_tmp, "wb") as dst:
        dst.write(src.read())

    # Copy envs directory to tmp_path so environment files are available
    envs_src = data / "envs"
    envs_dst = tmp_path / "envs"
    if envs_src.exists():
        shutil.copytree(envs_src, envs_dst, dirs_exist_ok=True)

    # The bundled_repos fixture already creates a symlink to the bundles directory
    # at tmp_path/bundles, so we don't need to copy it manually

    P1_output_folder_path = (
        tmp_path / "out" / "data" / "D1" / "default" / "process" / "P1"
    )

    expected_P1_param_dict = {
        "dca5d812e76298d5365cbc90ae29d54f867c0b23b0a3d4098826595b9a023054": {
            "a": "0",
            "b": "0",
        },
        "47fae0c4eff9d60d70f80ba26f299ece136010217b7393c672091bfbd5cecb1c": {
            "a": "1",
            "b": "0.1",
        },
        "2dc1a8f315206ba9aa73fd3af58b824faecc348f86c176c7325d2765151289ce": {
            "a": "2",
            "b": "0.5",
        },
    }

    with SnakemakeSetup(
        benchmark_file_in_tmp,
        keep_files=keep_files,
        cwd=tmp_path.as_posix(),
        out_dir=(tmp_path / "out").absolute(),
    ) as setup:
        benchmark = setup.benchmark

        # assert benchmark 1st run is successful
        success = setup.workflow.run_workflow(benchmark, work_dir=tmp_path)
        assert success is True

        # assert the parameter serialization is correct after 1st run
        assert_parameters_output_for_module(
            tmp_path, P1_output_folder_path, expected_P1_param_dict
        )

        # assert benchmark 2nd run is successful
        success = setup.workflow.run_workflow(benchmark, work_dir=tmp_path)
        assert success is True

        # assert the parameter serialization is correct after 2nd run
        assert_parameters_output_for_module(
            tmp_path, P1_output_folder_path, expected_P1_param_dict
        )


def test_multi_run_workflow_with_parameter_removal(
    snakemake_env,
    tmp_path,
    bundled_repos,  # noqa: F811
):
    # Use current directory if specified, otherwise use tmp_path
    # tmp_path is already provided by pytest fixture if not using current dir
    if snakemake_env["current_dir"]:
        tmp_path = Path(os.getcwd())

    benchmark_file = data / "Benchmark_003.yaml"
    keep_files = snakemake_env["keep_files"]

    # Copy benchmark YAML to tmp_path so workflow can find it
    benchmark_file_in_tmp = tmp_path / "Benchmark_003.yaml"
    with open(benchmark_file, "rb") as src, open(benchmark_file_in_tmp, "wb") as dst:
        dst.write(src.read())

    # Copy envs directory to tmp_path so environment files are available
    envs_src = data / "envs"
    envs_dst = tmp_path / "envs"
    if envs_src.exists():
        shutil.copytree(envs_src, envs_dst, dirs_exist_ok=True)

    # The bundled_repos fixture already creates a symlink to the bundles directory
    # at tmp_path/bundles, so we don't need to copy it manually

    P1_output_folder_path = (
        tmp_path / "out" / "data" / "D1" / "default" / "process" / "P1"
    )

    with SnakemakeSetup(
        benchmark_file_in_tmp,
        keep_files=keep_files,
        cwd=tmp_path.as_posix(),
        out_dir=(tmp_path / "out").absolute(),
    ) as setup:
        benchmark = setup.benchmark

        # assert benchmark 1st run is successful
        success = setup.workflow.run_workflow(benchmark, work_dir=tmp_path)
        assert success is True

        expected_param_dict_before_removal = {
            "dca5d812e76298d5365cbc90ae29d54f867c0b23b0a3d4098826595b9a023054": {
                "a": "0",
                "b": "0",
            },
            "47fae0c4eff9d60d70f80ba26f299ece136010217b7393c672091bfbd5cecb1c": {
                "a": "1",
                "b": "0.1",
            },
            "2dc1a8f315206ba9aa73fd3af58b824faecc348f86c176c7325d2765151289ce": {
                "a": "2",
                "b": "0.5",
            },
        }

        # assert the parameter serialization is correct after 1st run
        assert_parameters_output_for_module(
            tmp_path, P1_output_folder_path, expected_param_dict_before_removal
        )

        # remove previous parameters dict
        os.remove(P1_output_folder_path / "parameters_dict.tsv")

        # assert benchmark 2nd run is successful
        benchmark_file_trimmed = data / "Benchmark_003_trimmed.yaml"
        # Copy trimmed benchmark YAML to tmp_path
        benchmark_file_trimmed_in_tmp = tmp_path / "Benchmark_003_trimmed.yaml"
        with (
            open(benchmark_file_trimmed, "rb") as src,
            open(benchmark_file_trimmed_in_tmp, "wb") as dst,
        ):
            dst.write(src.read())
        benchmark_without_param = Benchmark(benchmark_file_trimmed_in_tmp)
        success = setup.workflow.run_workflow(
            benchmark_without_param, work_dir=tmp_path
        )
        assert success is True

        expected_param_dict_after_removal = {
            "dca5d812e76298d5365cbc90ae29d54f867c0b23b0a3d4098826595b9a023054": {
                "a": "0",
                "b": "0",
            },
            "2dc1a8f315206ba9aa73fd3af58b824faecc348f86c176c7325d2765151289ce": {
                "a": "2",
                "b": "0.5",
            },
        }

        # assert the parameter serialization is correct after 2nd run
        assert_parameters_output_for_module(
            tmp_path, P1_output_folder_path, expected_param_dict_after_removal
        )


def assert_parameters_output_for_module(
    cwd: Path, module_output_path: Path, expected_param_dict: Dict[str, Any]
):
    module_param_tsv = module_output_path / "parameters_dict.tsv"
    assert module_param_tsv.exists(), f"File not found: {module_param_tsv}"

    # Read the TSV file using the csv module
    with open(module_param_tsv, newline="") as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter="\t")
        for module_param_row in reader:
            base_path = module_param_row["base_path"]
            base_param_file_path = cwd / Path(base_path) / "parameters.json"
            assert base_param_file_path.exists(), (
                f"File not found through base path: {base_param_file_path}"
            )

            alias_path = module_param_row["base_path"]
            alias_param_file_path = cwd / Path(alias_path) / "parameters.json"
            assert alias_param_file_path.exists(), (
                f"File not found through alias: {alias_param_file_path}"
            )

            param_hash = module_param_row["id"]
            expected_content = expected_param_dict.get(param_hash)
            with open(base_param_file_path, "r") as file:
                actual_content = json.load(file)

            assert expected_content == actual_content, (
                f"Mismatch in {base_param_file_path}:\nExpected: {expected_content}\nGot: {actual_content}"
            )
