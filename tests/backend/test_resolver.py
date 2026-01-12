"""
Tests for ModuleResolver.

These tests verify the module resolution process:
- Cloning to cache
- Checking out to work directory
- Dereferencing entrypoints from omnibenchmark.yaml or config.cfg
"""

import pytest
from pathlib import Path
import tempfile
import shutil

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
