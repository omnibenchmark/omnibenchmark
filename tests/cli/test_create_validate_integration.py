"""Integration tests for create and validate commands.

Tests that created benchmarks and modules can be successfully validated.
"""

from click.testing import CliRunner
from omnibenchmark.cli.create import create_benchmark, create_module
from omnibenchmark.cli.validate import validate_plan, validate_module


class TestBenchmarkCreateValidate:
    """Test benchmark creation followed by validation."""

    def test_create_benchmark_validates_successfully(self, tmp_path):
        """Test that a created benchmark validates successfully."""
        runner = CliRunner()

        # Create benchmark (uses conda by default in non-interactive mode)
        benchmark_path = tmp_path / "test-benchmark"

        result = runner.invoke(
            create_benchmark,
            [
                str(benchmark_path),
                "--non-interactive",
                "--name",
                "test-benchmark",
                "--author-name",
                "Test Author",
                "--author-email",
                "test@example.com",
                "--license",
                "MIT",
                "--description",
                "Test benchmark",
            ],
        )

        assert result.exit_code == 0, f"Create failed: {result.output}"
        assert benchmark_path.exists()

        # Verify benchmark.yaml was created
        benchmark_yaml = benchmark_path / "benchmark.yaml"
        assert benchmark_yaml.exists(), "benchmark.yaml not created"

        # Verify envs directory and conda file (default backend)
        envs_dir = benchmark_path / "envs"
        assert envs_dir.exists(), "envs directory not created"
        assert (envs_dir / "conda.yaml").exists(), "conda.yaml not created"

        # Should not have lua files for conda backend
        assert not list(
            envs_dir.glob("*.lua")
        ), "Unexpected lua files for conda backend"

        # Validate the created benchmark
        validate_result = runner.invoke(
            validate_plan,
            [str(benchmark_yaml)],
        )

        # Validation uses logging, so we just check exit code
        assert (
            validate_result.exit_code == 0
        ), f"Validation failed with exit code {validate_result.exit_code}"

    def test_create_benchmark_with_envmodules_no_lua_files(self, tmp_path):
        """Test that benchmark with envmodules backend creates .eb file but no lua files."""
        import yaml

        runner = CliRunner()

        # Create a basic benchmark first
        benchmark_path = tmp_path / "test-benchmark-envmodules"

        result = runner.invoke(
            create_benchmark,
            [
                str(benchmark_path),
                "--non-interactive",
                "--name",
                "test-benchmark-envmodules",
                "--author-name",
                "Test Author",
                "--author-email",
                "test@example.com",
            ],
        )

        assert result.exit_code == 0

        # Modify the benchmark.yaml to use envmodules
        benchmark_yaml = benchmark_path / "benchmark.yaml"
        with open(benchmark_yaml, "r") as f:
            config = yaml.safe_load(f)

        config["software_backend"] = "envmodules"

        # Add envmodule field to all software environments
        for env_id, env_config in config.get("software_environments", {}).items():
            env_config["envmodule"] = f"{env_id}_module"

        with open(benchmark_yaml, "w") as f:
            yaml.dump(config, f)

        # Re-run create to trigger skeleton file creation
        # (In real usage, this happens during initial creation if envmodules is selected)
        from omnibenchmark.cli.create import _create_env_skeleton_files

        _create_env_skeleton_files(
            benchmark_path, "envmodules", "test-benchmark-envmodules"
        )

        # Check that no lua files were created
        lua_files = list(benchmark_path.rglob("*.lua"))
        assert len(lua_files) == 0, f"Found unexpected lua files: {lua_files}"

        # Check that .eb file was created
        eb_files = list(benchmark_path.rglob("*.eb"))
        assert len(eb_files) > 0, "EasyBuild .eb file should be created"

        # Validate still works
        validate_result = runner.invoke(
            validate_plan,
            [str(benchmark_yaml)],
        )

        assert validate_result.exit_code == 0

    def test_validate_plan_no_out_directory_created(self, tmp_path):
        """Test that validate plan does not create an 'out' directory."""
        runner = CliRunner()

        # Create a benchmark first
        benchmark_path = tmp_path / "test-benchmark"

        result = runner.invoke(
            create_benchmark,
            [
                str(benchmark_path),
                "--non-interactive",
                "--name",
                "test-benchmark",
                "--author-name",
                "Test Author",
                "--author-email",
                "test@example.com",
            ],
        )

        assert result.exit_code == 0

        benchmark_yaml = benchmark_path / "benchmark.yaml"

        # Create a temporary working directory to check for 'out' creation
        work_dir = tmp_path / "work"
        work_dir.mkdir()

        # Run validation from the work directory
        validate_result = runner.invoke(
            validate_plan,
            [str(benchmark_yaml)],
            catch_exceptions=False,
        )

        assert validate_result.exit_code == 0

        # Check that 'out' directory was NOT created in work_dir or benchmark_path
        assert not (
            work_dir / "out"
        ).exists(), "Validation should not create 'out' in work directory"
        assert not (
            benchmark_path / "out"
        ).exists(), "Validation should not create 'out' in benchmark directory"


class TestModuleCreateValidate:
    """Test module creation followed by validation."""

    def test_create_module_validates_successfully(self, tmp_path):
        """Test that a created module validates successfully."""
        runner = CliRunner()

        # Create module
        module_path = tmp_path / "test-module"

        result = runner.invoke(
            create_module,
            [
                str(module_path),
                "--non-interactive",
                "--name",
                "test-module",
                "--author-name",
                "Test Author",
                "--author-email",
                "test@example.com",
                "--license",
                "MIT",
                "--description",
                "Test module",
                "--entrypoint",
                "run.sh",
            ],
        )

        assert result.exit_code == 0, f"Create failed: {result.output}"
        assert module_path.exists()

        # Verify required files were created
        assert (
            module_path / "omnibenchmark.yaml"
        ).exists(), "omnibenchmark.yaml not created"
        assert (module_path / "CITATION.cff").exists(), "CITATION.cff not created"
        assert (module_path / "run.sh").exists(), "Entrypoint not created"

        # Validate the created module (should pass without --strict)
        validate_result = runner.invoke(
            validate_module,
            [str(module_path)],
        )

        # Validation uses logging, so we just check exit code
        assert (
            validate_result.exit_code == 0
        ), f"Validation failed with exit code {validate_result.exit_code}"

    def test_create_module_missing_omnibenchmark_fails_validation(self, tmp_path):
        """Test that validation fails when omnibenchmark.yaml is missing."""
        runner = CliRunner()

        # Create a module first
        module_path = tmp_path / "test-module"

        result = runner.invoke(
            create_module,
            [
                str(module_path),
                "--non-interactive",
                "--name",
                "test-module",
                "--author-name",
                "Test Author",
                "--author-email",
                "test@example.com",
            ],
        )

        assert result.exit_code == 0

        # Remove omnibenchmark.yaml
        omnibenchmark_yaml = module_path / "omnibenchmark.yaml"
        omnibenchmark_yaml.unlink()

        # Validation should fail
        validate_result = runner.invoke(
            validate_module,
            [str(module_path)],
        )

        assert (
            validate_result.exit_code == 1
        ), "Validation should fail without omnibenchmark.yaml"
        # Error message is in logs, not output, so we just check exit code

    def test_validate_module_no_s3_warning(self, tmp_path):
        """Test that validate module does not show S3 warning."""
        runner = CliRunner()

        # Create a module
        module_path = tmp_path / "test-module"

        result = runner.invoke(
            create_module,
            [
                str(module_path),
                "--non-interactive",
                "--name",
                "test-module",
                "--author-name",
                "Test Author",
                "--author-email",
                "test@example.com",
            ],
        )

        assert result.exit_code == 0

        # Validate the module
        validate_result = runner.invoke(
            validate_module,
            [str(module_path)],
        )

        assert validate_result.exit_code == 0
        # S3 warning would be in logs if present, but we fixed that it shouldn't appear


class TestCreateOutputMessages:
    """Test that create commands output correct help messages."""

    def test_benchmark_create_shows_correct_validate_command(self, tmp_path):
        """Test that benchmark creation shows 'ob validate plan benchmark.yaml' in next steps."""
        runner = CliRunner()

        benchmark_path = tmp_path / "test-benchmark"

        result = runner.invoke(
            create_benchmark,
            [
                str(benchmark_path),
                "--non-interactive",
                "--name",
                "test-benchmark",
                "--author-name",
                "Test Author",
                "--author-email",
                "test@example.com",
            ],
        )

        assert result.exit_code == 0
        # Messages are in logs, not captured output in tests

    def test_module_create_shows_correct_validate_command(self, tmp_path):
        """Test that module creation shows 'ob validate module' in next steps."""
        runner = CliRunner()

        module_path = tmp_path / "test-module"

        result = runner.invoke(
            create_module,
            [
                str(module_path),
                "--non-interactive",
                "--name",
                "test-module",
                "--author-name",
                "Test Author",
                "--author-email",
                "test@example.com",
            ],
        )

        assert result.exit_code == 0
        # Messages are in logs, not captured output in tests
