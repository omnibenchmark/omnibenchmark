"""
Test to verify that inputs_map is correctly passed to modules.

This test ensures that when inputs_map is a lambda function (as it is for
intermediate nodes), it gets properly evaluated and passed to modules with
the correct command-line format.
"""

import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch
from omnibenchmark.workflow.snakemake.scripts.execution import execution


@pytest.mark.short
def test_inputs_map_single_input():
    """Test that a single input is correctly formatted as CLI args."""

    module_dir = Path("/fake/module/dir")
    module_name = "test_module"
    output_dir = Path("/fake/output")
    dataset = "test_dataset"

    # Single input
    inputs_map = {"data.matrix": "/path/to/matrix.gz"}

    from omnibenchmark.benchmark.params import Params

    parameters = Params()

    # Mock the config file reading and execution
    with (
        patch(
            "omnibenchmark.workflow.snakemake.scripts.execution._read_entrypoint"
        ) as mock_read_entrypoint,
        patch(
            "omnibenchmark.workflow.snakemake.scripts.execution.os.path.exists",
            return_value=True,
        ),
        patch(
            "omnibenchmark.workflow.snakemake.scripts.execution.subprocess.Popen"
        ) as mock_popen,
        patch("builtins.open", MagicMock()),
    ):
        # Setup mock entrypoint
        mock_read_entrypoint.return_value = "run_module.py"

        # Setup successful subprocess Popen
        mock_process = MagicMock()
        mock_process.wait.return_value = 0
        mock_process.poll.return_value = 0
        mock_popen.return_value = mock_process

        # Execute
        exit_code = execution(
            module_dir=module_dir,
            module_name=module_name,
            output_dir=output_dir,
            dataset=dataset,
            inputs_map=inputs_map,
            parameters=parameters,
            keep_module_logs=False,
            timeout=100,
        )

        # Verify the command was called correctly
        assert mock_popen.called
        called_command = mock_popen.call_args[0][0]

        # Should contain: --data.matrix /path/to/matrix.gz
        assert "--data.matrix" in called_command
        matrix_idx = called_command.index("--data.matrix")
        assert called_command[matrix_idx + 1] == "/path/to/matrix.gz"

        assert exit_code == 0


@pytest.mark.short
def test_inputs_map_multiple_inputs():
    """Test that multiple inputs are correctly formatted as CLI args."""

    module_dir = Path("/fake/module/dir")
    module_name = "test_module"
    output_dir = Path("/fake/output")
    dataset = "test_dataset"

    # Multiple inputs
    inputs_map = {
        "data.matrix": "/path/to/matrix.gz",
        "data.true_labels": "/path/to/labels.gz",
    }

    from omnibenchmark.benchmark.params import Params

    parameters = Params()

    with (
        patch(
            "omnibenchmark.workflow.snakemake.scripts.execution._read_entrypoint"
        ) as mock_read_entrypoint,
        patch(
            "omnibenchmark.workflow.snakemake.scripts.execution.os.path.exists",
            return_value=True,
        ),
        patch(
            "omnibenchmark.workflow.snakemake.scripts.execution.subprocess.Popen"
        ) as mock_popen,
        patch("builtins.open", MagicMock()),
    ):
        # Setup mock entrypoint
        mock_read_entrypoint.return_value = "run_module.py"

        mock_process = MagicMock()
        mock_process.wait.return_value = 0
        mock_process.poll.return_value = 0
        mock_popen.return_value = mock_process

        exit_code = execution(
            module_dir=module_dir,
            module_name=module_name,
            output_dir=output_dir,
            dataset=dataset,
            inputs_map=inputs_map,
            parameters=parameters,
            keep_module_logs=False,
            timeout=100,
        )

        assert mock_popen.called
        called_command = mock_popen.call_args[0][0]

        # Should contain both inputs
        assert "--data.matrix" in called_command
        assert "--data.true_labels" in called_command

        matrix_idx = called_command.index("--data.matrix")
        assert called_command[matrix_idx + 1] == "/path/to/matrix.gz"

        labels_idx = called_command.index("--data.true_labels")
        assert called_command[labels_idx + 1] == "/path/to/labels.gz"

        assert exit_code == 0


@pytest.mark.short
def test_inputs_map_list_of_inputs():
    """Test that a list of inputs is correctly formatted as CLI args."""

    module_dir = Path("/fake/module/dir")
    module_name = "test_module"
    output_dir = Path("/fake/output")
    dataset = "test_dataset"

    # Input that is a list (e.g., multiple clustering results)
    inputs_map = {
        "clustering.results": [
            "/path/to/result1.gz",
            "/path/to/result2.gz",
            "/path/to/result3.gz",
        ]
    }

    from omnibenchmark.benchmark.params import Params

    parameters = Params()

    with (
        patch(
            "omnibenchmark.workflow.snakemake.scripts.execution._read_entrypoint"
        ) as mock_read_entrypoint,
        patch(
            "omnibenchmark.workflow.snakemake.scripts.execution.os.path.exists",
            return_value=True,
        ),
        patch(
            "omnibenchmark.workflow.snakemake.scripts.execution.subprocess.Popen"
        ) as mock_popen,
        patch("builtins.open", MagicMock()),
    ):
        # Setup mock entrypoint
        mock_read_entrypoint.return_value = "run_module.py"

        mock_process = MagicMock()
        mock_process.wait.return_value = 0
        mock_process.poll.return_value = 0
        mock_popen.return_value = mock_process

        exit_code = execution(
            module_dir=module_dir,
            module_name=module_name,
            output_dir=output_dir,
            dataset=dataset,
            inputs_map=inputs_map,
            parameters=parameters,
            keep_module_logs=False,
            timeout=100,
        )

        assert mock_popen.called
        called_command = mock_popen.call_args[0][0]

        # Should contain: --clustering.results /path/to/result1.gz /path/to/result2.gz /path/to/result3.gz
        assert "--clustering.results" in called_command
        results_idx = called_command.index("--clustering.results")

        # The three paths should follow the flag
        assert called_command[results_idx + 1] == "/path/to/result1.gz"
        assert called_command[results_idx + 2] == "/path/to/result2.gz"
        assert called_command[results_idx + 3] == "/path/to/result3.gz"

        assert exit_code == 0


@pytest.mark.short
def test_inputs_map_empty_dict():
    """Test that an empty inputs_map works correctly (for initial nodes)."""

    module_dir = Path("/fake/module/dir")
    module_name = "test_module"
    output_dir = Path("/fake/output")
    dataset = "test_dataset"

    # Empty inputs_map (initial nodes have no inputs)
    inputs_map = {}

    from omnibenchmark.benchmark.params import Params

    parameters = Params()

    with (
        patch(
            "omnibenchmark.workflow.snakemake.scripts.execution._read_entrypoint"
        ) as mock_read_entrypoint,
        patch(
            "omnibenchmark.workflow.snakemake.scripts.execution.os.path.exists",
            return_value=True,
        ),
        patch(
            "omnibenchmark.workflow.snakemake.scripts.execution.subprocess.Popen"
        ) as mock_popen,
        patch("builtins.open", MagicMock()),
    ):
        # Setup mock entrypoint
        mock_read_entrypoint.return_value = "run_module.py"

        mock_process = MagicMock()
        mock_process.wait.return_value = 0
        mock_process.poll.return_value = 0
        mock_popen.return_value = mock_process

        exit_code = execution(
            module_dir=module_dir,
            module_name=module_name,
            output_dir=output_dir,
            dataset=dataset,
            inputs_map=inputs_map,
            parameters=parameters,
            keep_module_logs=False,
            timeout=100,
        )

        assert mock_popen.called
        called_command = mock_popen.call_args[0][0]

        # Should still have output_dir and name
        assert "--output_dir" in called_command
        assert "--name" in called_command

        # But no input-related flags
        assert not any(arg.startswith("--data.") for arg in called_command)
        assert not any(arg.startswith("--clustering.") for arg in called_command)

        assert exit_code == 0
