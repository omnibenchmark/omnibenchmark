from tests.cli.cli_setup import OmniCLISetup


def test_default():
    expected_output = """
    Error: At least one option must be specified. Use '--input-dir', '--example', or '--all'.
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                "../data/mock_benchmark.yaml",
                "--module",
                "D1",
            ]
        )
        assert result.exit_code == 1
        assert clean(result.output) == clean(expected_output)


def test_multiple_behaviours_set():
    expected_output = """
    Error: Only one of '--input-dir', '--example', or '--all' should be set. Please choose only one option.
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                "../data/mock_benchmark.yaml",
                "--module",
                "D1",
                "--all",
                "--example",
            ]
        )
        assert result.exit_code == 1
        assert clean(result.output) == clean(expected_output)


def test_behaviour_example():
    expected_output = """
    Running module on a predefined remote example dataset.
    Error: Remote execution is not supported yet. Workflows can only be run in local mode.
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                "../data/mock_benchmark.yaml",
                "--module",
                "D1",
                "--example",
            ]
        )
        assert result.exit_code == 1
        assert clean(result.output) == clean(expected_output)


def test_behaviour_all():
    expected_output = """
    Running module on all available remote datasets.
    Error: Remote execution is not supported yet. Workflows can only be run in local mode.
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                "../data/mock_benchmark.yaml",
                "--module",
                "D1",
                "--all",
            ]
        )
        assert result.exit_code == 1
        assert clean(result.output) == clean(expected_output)


def test_benchmark_not_found():
    expected_output = """
    Running module on a dataset provided in a custom directory. (input_dir: ../data/)
    Error: Benchmark YAML file not found.
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                "../data/does_not_exist.yaml",
                "--module",
                "D1",
                "--input",
                "../data/",
            ]
        )
        assert result.exit_code == 1
        assert clean(result.output) == clean(expected_output)


def test_benchmark_format_incorrect():
    expected_output = """
    Running module on a dataset provided in a custom directory. (input_dir: ../data/)
    Error: Failed to parse YAML as a valid OmniBenchmark: software_environments must be supplied.
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                "../data/benchmark_format_incorrect.yaml",
                "--module",
                "D1",
                "--input",
                "../data/",
            ]
        )
        assert result.exit_code == 1
        assert clean(result.output) == clean(expected_output)


def test_module_not_found():
    expected_output = """
    Running module on a dataset provided in a custom directory. (input_dir: ../data/)
    Benchmark YAML file integrity check passed.
    Error: Could not find module with id `not-existing` in benchmark definition
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                "../data/mock_benchmark.yaml",
                "--module",
                "not-existing",
                "--input",
                "../data/",
            ]
        )
        assert result.exit_code == 1
        assert clean(result.output) == clean(expected_output)


def test_behaviour_input_dir():
    expected_output = """
    Running module on a dataset provided in a custom directory. (input_dir: ../data/)
    Benchmark YAML file integrity check passed.
    Found 1 workflow nodes for module D1.
    Running module benchmark...
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                "../data/mock_benchmark.yaml",
                "--module",
                "D1",
                "--input",
                "../data/",
            ]
        )
        assert result.exit_code == 0
        assert clean(result.output).startswith(clean(expected_output))


def test_behaviour_input_dir_dry():
    expected_output = """
    Running module on a dataset provided in a custom directory. (input_dir: ../data/)
    Benchmark YAML file integrity check passed.
    Found 1 workflow nodes for module D1.
    Running module benchmark...
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                "../data/mock_benchmark.yaml",
                "--module",
                "D1",
                "--input",
                "../data/",
                "--dry",
            ]
        )
        assert result.exit_code == 0
        assert clean(result.output).startswith(clean(expected_output))


def test_behaviour_input_dir_update_true():
    expected_output = """
    Running module on a dataset provided in a custom directory. (input_dir: ../data/)
    Are you sure you want to re-run the entire workflow? [y/N]: y
    Benchmark YAML file integrity check passed.
    Found 1 workflow nodes for module D1.
    Running module benchmark...
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                "../data/mock_benchmark.yaml",
                "--module",
                "D1",
                "--input",
                "../data/",
                "--update",
            ],
            input="y",
        )
        assert result.exit_code == 0
        assert clean(result.output).startswith(clean(expected_output))


def test_behaviour_input_dir_update_false():
    expected_output = """
    Running module on a dataset provided in a custom directory. (input_dir: ../data/)
    Are you sure you want to re-run the entire workflow? [y/N]: N
    Aborted
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                "../data/mock_benchmark.yaml",
                "--module",
                "D1",
                "--input",
                "../data/",
                "--update",
            ],
            input="N",
        )
        assert result.exit_code == 1
        assert clean(result.output).startswith(clean(expected_output))


def test_behaviour_input_dir_update_dry():
    expected_output = """
    Running module on a dataset provided in a custom directory. (input_dir: ../data/)
    Are you sure you want to re-run the entire workflow? [y/N]: y
    Benchmark YAML file integrity check passed.
    Found 1 workflow nodes for module D1.
    Running module benchmark...
    """
    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                "module",
                "--benchmark",
                "../data/mock_benchmark.yaml",
                "--module",
                "D1",
                "--input",
                "../data/",
                "--update",
                "--dry",
            ],
            input="y",
        )
        assert result.exit_code == 0
        assert clean(result.output).startswith(clean(expected_output))


def clean(output: str) -> str:
    return output.strip().replace("    ", "").replace("\t", "").replace("\n", "")
