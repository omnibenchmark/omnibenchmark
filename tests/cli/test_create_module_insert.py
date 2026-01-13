"""Tests for automatic module insertion into benchmark YAML during 'ob create module'."""

import pytest
import yaml

from tests.cli.cli_setup import OmniCLISetup


@pytest.fixture
def cli_setup():
    """Fixture for CLI setup."""
    return OmniCLISetup()


@pytest.fixture
def temp_benchmark_yaml(tmp_path):
    """Create a minimal benchmark YAML for testing."""
    benchmark_content = {
        "id": "test_benchmark",
        "description": "Test benchmark for module insertion",
        "version": "1.0.0",
        "benchmarker": "Test User <test@example.com>",
        "api_version": "0.3.0",
        "software_backend": "host",
        "software_environments": [
            {
                "id": "host",
                "description": "Host environment",
            }
        ],
        "stages": [
            {
                "id": "data",
                "modules": [
                    {
                        "id": "existing_module",
                        "name": "Existing Module",
                        "software_environment": "host",
                        "repository": {
                            "url": "https://github.com/example/existing.git",
                            "commit": "abc123",
                        },
                    }
                ],
                "outputs": [{"id": "data.output", "path": "output.txt"}],
            },
            {
                "id": "methods",
                "modules": [],
                "inputs": [{"entries": ["data.output"]}],
                "outputs": [{"id": "methods.result", "path": "result.txt"}],
            },
        ],
    }

    benchmark_path = tmp_path / "test_benchmark.yaml"
    with open(benchmark_path, "w") as f:
        yaml.dump(benchmark_content, f, default_flow_style=False, sort_keys=False)

    return benchmark_path


class TestCreateModuleWithInsertion:
    """Test module creation with automatic insertion into benchmark YAML."""

    def test_create_module_inserts_into_local_benchmark(
        self, cli_setup, tmp_path, temp_benchmark_yaml
    ):
        """Test that creating a module with --benchmark inserts it into the YAML."""
        module_dir = tmp_path / "new_method_module"

        result = cli_setup.call(
            [
                "create",
                "module",
                str(module_dir),
                "--benchmark",
                str(temp_benchmark_yaml),
                "--for-stage",
                "methods",
                "--name",
                "TestMethod",
                "--author-name",
                "Test Author",
                "--author-email",
                "author@example.com",
                "--non-interactive",
            ],
            cwd=str(tmp_path),
        )

        # Check command succeeded
        assert result.returncode == 0

        # Check module directory was created
        assert module_dir.exists()
        assert (module_dir / "omnibenchmark.yaml").exists()

        # Check module was inserted into benchmark YAML
        with open(temp_benchmark_yaml, "r") as f:
            updated_benchmark = yaml.safe_load(f)

        # Find the methods stage
        methods_stage = None
        for stage in updated_benchmark["stages"]:
            if stage["id"] == "methods":
                methods_stage = stage
                break

        assert methods_stage is not None, "Methods stage not found"
        assert "modules" in methods_stage
        assert len(methods_stage["modules"]) == 1, "Module was not inserted"

        # Check the inserted module
        new_module = methods_stage["modules"][0]
        assert new_module["name"] == "TestMethod"
        assert new_module["software_environment"] == "host"
        assert "repository" in new_module
        assert new_module["repository"]["commit"] == ""  # Empty for local dev

    def test_create_module_preserves_existing_modules(
        self, cli_setup, tmp_path, temp_benchmark_yaml
    ):
        """Test that inserting a module preserves existing modules in the stage."""
        module_dir = tmp_path / "new_data_module"

        result = cli_setup.call(
            [
                "create",
                "module",
                str(module_dir),
                "--benchmark",
                str(temp_benchmark_yaml),
                "--for-stage",
                "data",
                "--name",
                "NewDataModule",
                "--author-name",
                "Test Author",
                "--author-email",
                "author@example.com",
                "--non-interactive",
            ],
            cwd=str(tmp_path),
        )

        assert result.returncode == 0

        # Check existing module is preserved
        with open(temp_benchmark_yaml, "r") as f:
            updated_benchmark = yaml.safe_load(f)

        data_stage = next(s for s in updated_benchmark["stages"] if s["id"] == "data")
        assert len(data_stage["modules"]) == 2, "Should have both modules"

        module_ids = {m["id"] for m in data_stage["modules"]}
        assert "existing_module" in module_ids, "Existing module should be preserved"
        assert "newdatamodule" in module_ids, "New module should be added"

    def test_create_module_with_url_benchmark_skips_insertion(
        self, cli_setup, tmp_path
    ):
        """Test that URL benchmarks don't trigger insertion."""
        module_dir = tmp_path / "url_test_module"

        result = cli_setup.call(
            [
                "create",
                "module",
                str(module_dir),
                "--benchmark",
                "https://github.com/example/benchmark.yaml",
                "--for-stage",
                "methods",
                "--name",
                "URLTest",
                "--author-name",
                "Test Author",
                "--author-email",
                "author@example.com",
                "--non-interactive",
            ],
            cwd=str(tmp_path),
        )

        # Should succeed but skip insertion
        assert result.returncode == 0
        assert (
            "skipping auto-insertion" in result.output.lower()
            or "url" in result.output.lower()
        )

    def test_create_module_without_benchmark_skips_insertion(self, cli_setup, tmp_path):
        """Test that modules created without --benchmark don't attempt insertion."""
        module_dir = tmp_path / "standalone_module"

        result = cli_setup.call(
            [
                "create",
                "module",
                str(module_dir),
                "--name",
                "Standalone",
                "--author-name",
                "Test Author",
                "--author-email",
                "author@example.com",
                "--non-interactive",
            ],
            cwd=str(tmp_path),
        )

        assert result.returncode == 0
        assert module_dir.exists()

    def test_create_module_with_nonexistent_stage_fails_gracefully(
        self, cli_setup, tmp_path, temp_benchmark_yaml
    ):
        """Test that specifying a non-existent stage logs error but creates module."""
        module_dir = tmp_path / "bad_stage_module"

        _ = cli_setup.call(
            [
                "create",
                "module",
                str(module_dir),
                "--benchmark",
                str(temp_benchmark_yaml),
                "--for-stage",
                "nonexistent_stage",
                "--name",
                "BadStage",
                "--author-name",
                "Test Author",
                "--author-email",
                "author@example.com",
                "--non-interactive",
            ],
            cwd=str(tmp_path),
        )

        # Module should still be created
        assert module_dir.exists()

        # But insertion should fail with error message
        # The actual behavior depends on whether extraction succeeds
        # If stage doesn't exist, _extract_stage_inputs should fail first

    def test_module_path_is_relative_when_possible(
        self, cli_setup, tmp_path, temp_benchmark_yaml
    ):
        """Test that module repository path is relative when in same directory tree."""
        # Create module in subdirectory of benchmark location
        benchmark_dir = temp_benchmark_yaml.parent
        module_dir = benchmark_dir / "modules" / "relative_module"
        module_dir.parent.mkdir(exist_ok=True)

        result = cli_setup.call(
            [
                "create",
                "module",
                str(module_dir),
                "--benchmark",
                str(temp_benchmark_yaml),
                "--for-stage",
                "methods",
                "--name",
                "RelativeModule",
                "--author-name",
                "Test Author",
                "--author-email",
                "author@example.com",
                "--non-interactive",
            ],
            cwd=str(tmp_path),
        )

        assert result.returncode == 0

        # Check the repository URL is relative
        with open(temp_benchmark_yaml, "r") as f:
            updated_benchmark = yaml.safe_load(f)

        methods_stage = next(
            s for s in updated_benchmark["stages"] if s["id"] == "methods"
        )
        new_module = methods_stage["modules"][0]

        # Should be a relative path
        repo_url = new_module["repository"]["url"]
        assert not repo_url.startswith("/"), "Should be relative path"
        assert "modules/relative_module" in repo_url
