"""
Tests for ModuleResolver.

These tests verify the module resolution process:
- Cloning to cache
- Checking out to work directory
- Dereferencing entrypoints from omnibenchmark.yaml or config.cfg
- Named entrypoint resolution
"""

import pytest
from pathlib import Path
import tempfile
import shutil

import yaml

from omnibenchmark.backend.resolver import ModuleResolver
from omnibenchmark.model import Benchmark


@pytest.fixture
def simple_benchmark():
    """Load the simple e2e benchmark."""
    benchmark_path = (
        Path(__file__).parent.parent / "e2e" / "configs" / "01_data_and_methods.yaml"
    )
    return Benchmark.from_yaml(benchmark_path)


@pytest.fixture
def temp_work_dir():
    """Create a temporary work directory."""
    temp_dir = Path(tempfile.mkdtemp())
    yield temp_dir
    shutil.rmtree(temp_dir, ignore_errors=True)


class TestModuleResolver:
    """Tests for ModuleResolver."""

    def test_resolver_creation(self, temp_work_dir):
        """Test creating a resolver."""
        resolver = ModuleResolver(work_base_dir=temp_work_dir)
        assert resolver.work_base_dir == temp_work_dir

    @pytest.mark.slow
    def test_populate_cache(self, simple_benchmark, temp_work_dir):
        """Test populating the git cache."""
        resolver = ModuleResolver(work_base_dir=temp_work_dir)

        # Get first module from data stage
        module = simple_benchmark.stages[0].modules[0]

        # Populate cache (should clone the repository)
        cache_path = resolver.populate_cache(module)

        assert cache_path.exists()
        assert (cache_path / ".git").exists()

    @pytest.mark.slow
    def test_resolve_module_with_entrypoint(self, simple_benchmark, temp_work_dir):
        """Test resolving a module and dereferencing entrypoint."""
        resolver = ModuleResolver(work_base_dir=temp_work_dir)

        # Get first module from data stage
        module = simple_benchmark.stages[0].modules[0]

        # Resolve the module
        resolved = resolver.resolve(
            module=module,
            module_id=module.id,
            software_environment_id=module.software_environment,
            dirty=False,
        )

        # Verify resolved module
        assert resolved.repository_url == module.repository.url
        assert resolved.commit == module.repository.commit
        assert resolved.module_dir.exists()
        assert resolved.entrypoint is not None
        assert resolved.software_environment_id == module.software_environment
        assert resolved.dirty is False

        # Verify entrypoint exists
        entrypoint_full = resolved.module_dir / resolved.entrypoint
        assert entrypoint_full.exists()

    def test_resolve_without_commit_fails_in_production(self, temp_work_dir):
        """Test that resolve fails without commit in production mode."""
        from omnibenchmark.model import Module, Repository

        # Create a module without commit
        module = Module(
            id="test_module",
            repository=Repository(url="https://github.com/test/repo.git", commit=""),
            software_environment="host",
        )

        resolver = ModuleResolver(work_base_dir=temp_work_dir)

        # Should fail in production mode
        with pytest.raises(RuntimeError, match="has no commit specified"):
            resolver.resolve(
                module=module,
                module_id=module.id,
                software_environment_id=module.software_environment,
                dirty=False,
            )

    @pytest.mark.slow
    def test_work_dir_is_relative_and_portable(self, simple_benchmark, temp_work_dir):
        """Test that work directories use relative paths for portability."""
        resolver = ModuleResolver(work_base_dir=Path(".snakemake/repos"))

        # Get first module
        module = simple_benchmark.stages[0].modules[0]

        # Resolve
        resolved = resolver.resolve(
            module=module,
            module_id=module.id,
            software_environment_id=module.software_environment,
            dirty=False,
        )

        # Work dir should be relative (or we can make it relative)
        # The important part is that the Snakefile can reference it relatively
        work_dir_str = str(resolved.module_dir)

        # Should start with .snakemake/repos
        assert ".snakemake/repos" in work_dir_str or resolved.module_dir.is_relative_to(
            Path(".snakemake/repos")
        )

    @pytest.mark.slow
    def test_resolve_conda_environment(self, temp_work_dir):
        """Test that conda environments are properly resolved and copied."""
        # Load benchmark with conda backend
        benchmark_path = (
            Path(__file__).parent.parent / "e2e" / "configs" / "conda_test.yaml"
        )
        benchmark = Benchmark.from_yaml(benchmark_path)

        # Create output directory structure
        output_dir = temp_work_dir / "out"
        output_dir.mkdir(parents=True, exist_ok=True)

        # Get benchmark directory (where the yaml and envs live)
        benchmark_dir = benchmark_path.parent

        # Create resolver with environment support
        resolver = ModuleResolver(
            work_base_dir=temp_work_dir / ".repos",
            output_dir=output_dir,
            software_backend=benchmark.get_software_backend(),
            software_environments=benchmark.get_software_environments(),
            benchmark_dir=benchmark_dir,
        )

        # Get first module (should reference python_env)
        module = benchmark.stages[0].modules[0]

        # Resolve the module
        resolved = resolver.resolve(
            module=module,
            module_id=module.id,
            software_environment_id=module.software_environment,
            dirty=False,
        )

        # Verify that resolved_environment is populated
        assert resolved.resolved_environment is not None
        assert resolved.resolved_environment.backend_type.value == "conda"

        # Verify the conda file was copied to .envs/
        assert resolved.resolved_environment.reference is not None
        env_file = output_dir / resolved.resolved_environment.reference
        assert env_file.exists(), f"Expected conda env file at {env_file}"

        # Verify the path is relative to output dir (.envs/ not .out/envs/)
        assert resolved.resolved_environment.reference.startswith(".envs/")

        # Verify it's a valid YAML file
        import yaml

        with open(env_file, "r") as f:
            env_data = yaml.safe_load(f)
            assert env_data is not None
            assert "dependencies" in env_data or "name" in env_data

    @pytest.mark.slow
    def test_resolve_envmodules_environment(self, temp_work_dir):
        """Test that envmodules are properly resolved."""
        # Load benchmark with envmodules backend
        benchmark_path = (
            Path(__file__).parent.parent / "e2e" / "configs" / "envmodules_test.yaml"
        )
        benchmark = Benchmark.from_yaml(benchmark_path)

        # Create output directory structure
        output_dir = temp_work_dir / "out"
        output_dir.mkdir(parents=True, exist_ok=True)

        # Get benchmark directory
        benchmark_dir = benchmark_path.parent

        # Create resolver with environment support
        resolver = ModuleResolver(
            work_base_dir=temp_work_dir / ".repos",
            output_dir=output_dir,
            software_backend=benchmark.get_software_backend(),
            software_environments=benchmark.get_software_environments(),
            benchmark_dir=benchmark_dir,
        )

        # Get first module (should reference gcc_env)
        module = benchmark.stages[0].modules[0]

        # Resolve the module
        resolved = resolver.resolve(
            module=module,
            module_id=module.id,
            software_environment_id=module.software_environment,
            dirty=False,
        )

        # Verify that resolved_environment is populated
        assert resolved.resolved_environment is not None
        assert resolved.resolved_environment.backend_type.value == "envmodules"

        # For envmodules, reference should be the module name
        assert resolved.resolved_environment.reference == "GCC/11.3.0"
        assert not resolved.resolved_environment.is_url()

        # No files should be created for envmodules (just module names)
        _envs_dir = output_dir / ".envs"
        # Directory might exist but should be empty or not contain module-related files
        # (envmodules don't need file copies)

    @pytest.mark.slow
    def test_resolve_apptainer_environment_url(self, temp_work_dir):
        """Test that apptainer URL-based images are properly resolved."""
        # Load benchmark with apptainer backend
        benchmark_path = (
            Path(__file__).parent.parent / "e2e" / "configs" / "apptainer_test.yaml"
        )
        benchmark = Benchmark.from_yaml(benchmark_path)

        # Create output directory structure
        output_dir = temp_work_dir / "out"
        output_dir.mkdir(parents=True, exist_ok=True)

        # Get benchmark directory
        benchmark_dir = benchmark_path.parent

        # Create resolver with environment support
        resolver = ModuleResolver(
            work_base_dir=temp_work_dir / ".repos",
            output_dir=output_dir,
            software_backend=benchmark.get_software_backend(),
            software_environments=benchmark.get_software_environments(),
            benchmark_dir=benchmark_dir,
        )

        # Get first module (should reference remote_image with docker:// URL)
        module = benchmark.stages[0].modules[0]

        # Resolve the module
        resolved = resolver.resolve(
            module=module,
            module_id=module.id,
            software_environment_id=module.software_environment,
            dirty=False,
        )

        # Verify that resolved_environment is populated
        assert resolved.resolved_environment is not None
        assert resolved.resolved_environment.backend_type.value == "apptainer"

        # For URL-based images, reference should be the URL itself
        assert resolved.resolved_environment.reference == "docker://ubuntu:22.04"
        assert resolved.resolved_environment.is_url()

    @pytest.mark.slow
    def test_resolve_apptainer_environment_local_sif(self, temp_work_dir):
        """Test that apptainer local SIF files are properly symlinked."""
        # Load benchmark with apptainer backend
        benchmark_path = (
            Path(__file__).parent.parent / "e2e" / "configs" / "apptainer_test.yaml"
        )
        benchmark = Benchmark.from_yaml(benchmark_path)

        # Create output directory structure
        output_dir = temp_work_dir / "out"
        output_dir.mkdir(parents=True, exist_ok=True)

        # Get benchmark directory
        benchmark_dir = benchmark_path.parent

        # Create resolver with environment support
        resolver = ModuleResolver(
            work_base_dir=temp_work_dir / ".repos",
            output_dir=output_dir,
            software_backend=benchmark.get_software_backend(),
            software_environments=benchmark.get_software_environments(),
            benchmark_dir=benchmark_dir,
        )

        # Get second stage module (should reference local_image with local SIF)
        module = benchmark.stages[1].modules[0]

        # Resolve the module
        resolved = resolver.resolve(
            module=module,
            module_id=module.id,
            software_environment_id=module.software_environment,
            dirty=False,
        )

        # Verify that resolved_environment is populated
        assert resolved.resolved_environment is not None
        assert resolved.resolved_environment.backend_type.value == "apptainer"

        # For local SIF, should be a symlink in .envs/
        assert resolved.resolved_environment.reference.startswith(".envs/")
        assert resolved.resolved_environment.reference.endswith(".sif")
        assert not resolved.resolved_environment.is_url()

        # Verify symlink was created
        symlink_path = output_dir / resolved.resolved_environment.reference
        assert symlink_path.exists(), f"Expected symlink at {symlink_path}"
        assert symlink_path.is_symlink(), f"Expected {symlink_path} to be a symlink"

        # Verify symlink points to the correct source file
        source_path = benchmark_dir / "containers" / "myimage.sif"
        assert symlink_path.resolve() == source_path.resolve()


class TestReadEntrypoint:
    """Tests for _read_entrypoint with named entrypoint keys."""

    @pytest.fixture
    def resolver(self, temp_work_dir):
        return ModuleResolver(work_base_dir=temp_work_dir)

    @pytest.fixture
    def module_dir_with_entrypoints(self, tmp_path):
        """Create a module directory with multiple entrypoints in omnibenchmark.yaml."""
        module_dir = tmp_path / "test_module"
        module_dir.mkdir()

        config = {
            "entrypoints": {
                "default": "run.py",
                "preprocess": "preprocess.py",
                "validate": "scripts/validate.R",
            }
        }
        (module_dir / "omnibenchmark.yaml").write_text(yaml.dump(config))
        return module_dir

    @pytest.fixture
    def module_dir_with_config_cfg(self, tmp_path):
        """Create a module directory with legacy config.cfg."""
        module_dir = tmp_path / "legacy_module"
        module_dir.mkdir()

        (module_dir / "config.cfg").write_text("[DEFAULT]\nSCRIPT = legacy_run.py\n")
        return module_dir

    def test_default_entrypoint(self, resolver, module_dir_with_entrypoints):
        """Default entrypoint key resolves to 'default' entry."""
        result = resolver._read_entrypoint(module_dir_with_entrypoints, "mod1")
        assert result == "run.py"

    def test_explicit_default_key(self, resolver, module_dir_with_entrypoints):
        """Explicitly passing 'default' works the same."""
        result = resolver._read_entrypoint(
            module_dir_with_entrypoints, "mod1", entrypoint_key="default"
        )
        assert result == "run.py"

    def test_named_entrypoint(self, resolver, module_dir_with_entrypoints):
        """A named entrypoint key resolves to the correct script."""
        result = resolver._read_entrypoint(
            module_dir_with_entrypoints, "mod1", entrypoint_key="preprocess"
        )
        assert result == "preprocess.py"

    def test_another_named_entrypoint(self, resolver, module_dir_with_entrypoints):
        """Another named entrypoint resolves correctly."""
        result = resolver._read_entrypoint(
            module_dir_with_entrypoints, "mod1", entrypoint_key="validate"
        )
        assert result == "scripts/validate.R"

    def test_missing_entrypoint_key_returns_none(
        self, resolver, module_dir_with_entrypoints
    ):
        """A missing entrypoint key returns None."""
        result = resolver._read_entrypoint(
            module_dir_with_entrypoints, "mod1", entrypoint_key="nonexistent"
        )
        assert result is None

    def test_missing_entrypoint_key_logs_available(
        self, resolver, module_dir_with_entrypoints, caplog
    ):
        """A missing entrypoint key logs the available keys."""
        import logging

        with caplog.at_level(logging.ERROR):
            resolver._read_entrypoint(
                module_dir_with_entrypoints, "mod1", entrypoint_key="nonexistent"
            )
        assert "nonexistent" in caplog.text
        assert "Available entrypoints:" in caplog.text
        assert "default" in caplog.text
        assert "preprocess" in caplog.text

    def test_config_cfg_default_key_works(self, resolver, module_dir_with_config_cfg):
        """config.cfg fallback works for default entrypoint."""
        with pytest.warns(FutureWarning, match="deprecated config.cfg"):
            result = resolver._read_entrypoint(
                module_dir_with_config_cfg, "mod1", entrypoint_key="default"
            )
        assert result == "legacy_run.py"

    def test_config_cfg_named_key_rejected(
        self, resolver, module_dir_with_config_cfg, caplog
    ):
        """config.cfg fallback rejects non-default entrypoint keys."""
        import logging

        with caplog.at_level(logging.ERROR):
            result = resolver._read_entrypoint(
                module_dir_with_config_cfg, "mod1", entrypoint_key="preprocess"
            )
        assert result is None
        assert "only supports" in caplog.text
        assert "preprocess" in caplog.text


class TestRepositoryEntrypointField:
    """Tests for the entrypoint field on Repository model."""

    def test_repository_default_entrypoint_is_none(self):
        """Repository without entrypoint field defaults to None."""
        from omnibenchmark.model.benchmark import Repository

        repo = Repository(url="https://example.com/repo.git", commit="abc123")
        assert repo.entrypoint is None

    def test_repository_with_entrypoint(self):
        """Repository with entrypoint field stores the value."""
        from omnibenchmark.model.benchmark import Repository

        repo = Repository(
            url="https://example.com/repo.git",
            commit="abc123",
            entrypoint="preprocess",
        )
        assert repo.entrypoint == "preprocess"

    def test_repository_empty_entrypoint_rejected(self):
        """Repository with empty string entrypoint is rejected."""
        from pydantic import ValidationError as PydanticValidationError
        from omnibenchmark.model.benchmark import Repository

        with pytest.raises(PydanticValidationError, match="entrypoint"):
            Repository(
                url="https://example.com/repo.git",
                commit="abc123",
                entrypoint="",
            )

    def test_repository_whitespace_entrypoint_rejected(self):
        """Repository with whitespace-only entrypoint is rejected."""
        from pydantic import ValidationError as PydanticValidationError
        from omnibenchmark.model.benchmark import Repository

        with pytest.raises(PydanticValidationError, match="entrypoint"):
            Repository(
                url="https://example.com/repo.git",
                commit="abc123",
                entrypoint="   ",
            )

    def test_repository_entrypoint_from_yaml(self):
        """Repository entrypoint parses correctly from YAML."""
        yaml_content = """\
id: test
benchmarker: tester
version: "1.0.0"
software_backend: host
software_environments:
  host:
    description: host
stages:
  - id: data
    modules:
      - id: D1
        software_environment: host
        repository:
          url: https://example.com/repo.git
          commit: abc123
          entrypoint: preprocess
    outputs:
      - id: data.raw
        path: output.json
"""
        benchmark = Benchmark.from_yaml(yaml_content)
        module = benchmark.stages[0].modules[0]
        assert module.repository.entrypoint == "preprocess"

    def test_repository_without_entrypoint_from_yaml(self):
        """Repository without entrypoint in YAML defaults to None."""
        yaml_content = """\
id: test
benchmarker: tester
version: "1.0.0"
software_backend: host
software_environments:
  host:
    description: host
stages:
  - id: data
    modules:
      - id: D1
        software_environment: host
        repository:
          url: https://example.com/repo.git
          commit: abc123
    outputs:
      - id: data.raw
        path: output.json
"""
        benchmark = Benchmark.from_yaml(yaml_content)
        module = benchmark.stages[0].modules[0]
        assert module.repository.entrypoint is None
