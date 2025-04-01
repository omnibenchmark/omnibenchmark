from pathlib import Path

from tests.workflow.Snakemake_setup import SnakemakeSetup
import textwrap


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

    expected_D1_param_dict = textwrap.dedent(
        """param_0 ['--measure', 'cosine']
        param_1 ['--measure', 'euclidean']
        param_2 ['--measure', 'manhattan']
        param_3 ['--measure', 'chebyshev']
    """
    )

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
                Path(__file__).parent
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
        Path(__file__).parent / "out" / "data" / "iris" / "default" / "distances" / "D1"
    )

    expected_D1_param_dict = textwrap.dedent(
        """param_0 ['--measure', 'cosine']
        param_1 ['--measure', 'euclidean']
        param_2 ['--measure', 'manhattan']
        param_3 ['--measure', 'chebyshev']
    """
    )

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
        Path(__file__).parent / "out" / "data" / "iris" / "default" / "distances" / "D1"
    )

    with SnakemakeSetup(benchmark_file_path) as setup:
        benchmark = setup.benchmark

        # assert benchmark 1st run is successful
        success = setup.workflow.run_workflow(benchmark)
        assert success

        expected_param_dict_before_removal = textwrap.dedent(
            """param_0 ['--measure', 'cosine']
            param_1 ['--measure', 'euclidean']
            param_2 ['--measure', 'manhattan']
            param_3 ['--measure', 'chebyshev']
        """
        )

        # assert the parameter serialization is correct after 1st run
        _assert_parameters_output_for_module(
            D1_output_folder_path, expected_param_dict_before_removal
        )

        node_id_to_remove = "distances-D1-param_1-after_data"
        benchmark = benchmark.remove_node_and_dependents(node_id_to_remove)

        # assert benchmark 2nd run is successful
        success = setup.workflow.run_workflow(benchmark)
        assert success

        expected_param_dict_after_removal = textwrap.dedent(
            """param_0 ['--measure', 'cosine']
            param_1 ['--measure', 'manhattan']
            param_2 ['--measure', 'chebyshev']
        """
        )

        # assert the parameter serialization is correct after 2nd run
        _assert_parameters_output_for_module(
            D1_output_folder_path, expected_param_dict_after_removal
        )


def _assert_parameters_output_for_module(
    module_output_path: Path, expected_module_param_dict: str
):
    module_param_dict = module_output_path / "parameters_dict.txt"
    assert module_param_dict.exists(), f"File not found: {module_param_dict}"

    with module_param_dict.open("r") as f:
        param_dict_lines = [line.strip() for line in f if line.strip()]
        actual_module_param_dict = f.read()

    assert (
        actual_module_param_dict == expected_module_param_dict
    ), f"Mismatch in parameters:\nExpected:\n{expected_module_param_dict}\n\nGot:\n{actual_module_param_dict}"

    param_dict = {}
    for line in param_dict_lines:
        key, value = line.split(" ", 1)
        param_dict[key] = value

    # check that the content of the subfolder actually contains the parameters specified in the dict
    for param_folder in param_dict.keys():
        module_param_file_path = module_output_path / param_folder / "parameters.txt"
        assert (
            module_param_file_path.exists()
        ), f"File not found: {module_param_file_path}"

        with module_param_file_path.open("r") as f:
            param_txt_content = f.read()

        expected_content = param_dict.get(param_folder)
        assert (
            expected_content == param_txt_content
        ), f"Mismatch in {module_param_file_path}:\nExpected: {expected_content}\nGot: {param_txt_content}"
