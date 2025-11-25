"""Unit tests for Stage inputs field validator.

This tests the validate_inputs field validator that handles both:
1. NEW format: flat list of strings ['data.matrix', 'data.true_labels']
2. OLD format: list of dicts [{'entries': ['data.matrix', 'data.true_labels']}]
"""

import pytest
from omnibenchmark.model.benchmark import Stage, InputCollection


@pytest.mark.short
class TestStageInputsValidator:
    """Test the Stage.validate_inputs field validator."""

    def test_new_format_flat_list_creates_single_input_collection(self):
        """Test that a flat list of strings creates one InputCollection with multiple entries.

        This is the NEW recommended format:
        inputs: ['data.matrix', 'data.true_labels']

        Should create ONE InputCollection with entries=['data.matrix', 'data.true_labels']
        NOT two separate InputCollections.

        Validator should not treat each string as a separate InputCollection,
        causing duplicate nodes in the DAG.
        """
        stage_data = {
            "id": "clustering",
            "modules": [
                {
                    "id": "test_module",
                    "name": "Test Module",
                    "software_environment": "test_env",
                    "repository": {
                        "url": "https://github.com/example/repo",
                        "commit": "abc123",
                    },
                }
            ],
            "inputs": ["data.matrix", "data.true_labels"],
            "outputs": [{"id": "clustering.output", "path": "output.txt"}],
        }

        stage = Stage(**stage_data)

        # Should create exactly ONE InputCollection
        assert len(stage.inputs) == 1
        assert isinstance(stage.inputs[0], InputCollection)

        # That InputCollection should have BOTH entries
        assert len(stage.inputs[0].entries) == 2
        assert "data.matrix" in stage.inputs[0].entries
        assert "data.true_labels" in stage.inputs[0].entries

    def test_old_format_entries_dict_creates_single_input_collection(self):
        """Test that old format with entries dict creates ONE InputCollection.

        This is the OLD deprecated format:
        inputs:
          - entries:
              - data.matrix
              - data.true_labels

        Should create ONE InputCollection with entries=['data.matrix', 'data.true_labels']
        """
        stage_data = {
            "id": "clustering",
            "modules": [
                {
                    "id": "test_module",
                    "name": "Test Module",
                    "software_environment": "test_env",
                    "repository": {
                        "url": "https://github.com/example/repo",
                        "commit": "abc123",
                    },
                }
            ],
            "inputs": [{"entries": ["data.matrix", "data.true_labels"]}],  # OLD format
            "outputs": [{"id": "clustering.output", "path": "output.txt"}],
        }

        # Should emit deprecation warning
        with pytest.warns(
            FutureWarning, match="Using 'entries' field in inputs is deprecated"
        ):
            stage = Stage(**stage_data)

        # Should create exactly ONE InputCollection
        assert len(stage.inputs) == 1
        assert isinstance(stage.inputs[0], InputCollection)

        # That InputCollection should have BOTH entries
        assert len(stage.inputs[0].entries) == 2
        assert "data.matrix" in stage.inputs[0].entries
        assert "data.true_labels" in stage.inputs[0].entries

    def test_new_format_single_input(self):
        """Test new format with a single input string."""
        stage_data = {
            "id": "metrics",
            "modules": [
                {
                    "id": "test_module",
                    "name": "Test Module",
                    "software_environment": "test_env",
                    "repository": {
                        "url": "https://github.com/example/repo",
                        "commit": "abc123",
                    },
                }
            ],
            "inputs": ["clustering.output"],  # NEW format with single input
            "outputs": [{"id": "metrics.result", "path": "result.txt"}],
        }

        stage = Stage(**stage_data)

        assert len(stage.inputs) == 1
        assert isinstance(stage.inputs[0], InputCollection)

        assert len(stage.inputs[0].entries) == 1
        assert stage.inputs[0].entries[0] == "clustering.output"

    def test_new_format_empty_list(self):
        """Test new format with empty list."""
        stage_data = {
            "id": "data",
            "modules": [
                {
                    "id": "test_module",
                    "name": "Test Module",
                    "software_environment": "test_env",
                    "repository": {
                        "url": "https://github.com/example/repo",
                        "commit": "abc123",
                    },
                }
            ],
            "inputs": [],  # Empty inputs
            "outputs": [{"id": "data.output", "path": "output.txt"}],
        }

        stage = Stage(**stage_data)

        # Empty list should create ONE InputCollection with empty entries
        assert len(stage.inputs) == 1
        assert isinstance(stage.inputs[0], InputCollection)
        assert len(stage.inputs[0].entries) == 0

    def test_none_inputs(self):
        """Test that None inputs is handled correctly."""
        stage_data = {
            "id": "data",
            "modules": [
                {
                    "id": "test_module",
                    "name": "Test Module",
                    "software_environment": "test_env",
                    "repository": {
                        "url": "https://github.com/example/repo",
                        "commit": "abc123",
                    },
                }
            ],
            "inputs": None,  # No inputs
            "outputs": [{"id": "data.output", "path": "output.txt"}],
        }

        stage = Stage(**stage_data)

        # None should remain None
        assert stage.inputs is None

    def test_multiple_old_format_input_collections(self):
        """Test multiple InputCollections in old format (unusual but valid)."""
        stage_data = {
            "id": "complex_stage",
            "modules": [
                {
                    "id": "test_module",
                    "name": "Test Module",
                    "software_environment": "test_env",
                    "repository": {
                        "url": "https://github.com/example/repo",
                        "commit": "abc123",
                    },
                }
            ],
            "inputs": [
                {"entries": ["input1", "input2"]},
                {"entries": ["input3", "input4"]},
            ],
            "outputs": [{"id": "output", "path": "output.txt"}],
        }

        with pytest.warns(FutureWarning):
            stage = Stage(**stage_data)

        # Should create TWO InputCollections.
        # Not that we care too much about this format, but we should still test it.
        assert len(stage.inputs) == 2
        assert all(isinstance(ic, InputCollection) for ic in stage.inputs)

        # First InputCollection has input1 and input2
        assert len(stage.inputs[0].entries) == 2
        assert set(stage.inputs[0].entries) == {"input1", "input2"}

        # Second InputCollection has input3 and input4
        assert len(stage.inputs[1].entries) == 2
        assert set(stage.inputs[1].entries) == {"input3", "input4"}
