from pathlib import Path

import pytest

from omnibenchmark.benchmark import BenchmarkExecution
from omnibenchmark.model import SoftwareBackendEnum


class TestBenchmarkNodeEnvironmentPath:
    """Test BenchmarkNode.get_environment_path() method."""

    @pytest.fixture
    def benchmark_execution(self):
        """Load a real benchmark from test data."""
        benchmark_file_path = Path(__file__).parent / "../data/Benchmark_001.yaml"
        return BenchmarkExecution(benchmark_file_path)

    @pytest.mark.short
    def test_node_get_environment_path_with_conda_backend(self, benchmark_execution):
        """Test that node.get_environment_path() returns conda path for conda backend."""
        # Get a node that has a software environment
        nodes = benchmark_execution.get_nodes()
        node = next(n for n in nodes if n.get_software_environment() == "zlib_old")

        # Act
        result = node.get_environment_path(
            "zlib_old",
            SoftwareBackendEnum.conda,
            benchmark_execution.context.directory,
        )

        # Assert
        assert result is not None
        assert "envs/zlib_1.2.11.yaml" in result

    @pytest.mark.short
    def test_node_get_environment_path_with_envmodules_backend(
        self, benchmark_execution
    ):
        """Test that node.get_environment_path() returns envmodule for envmodules backend."""
        nodes = benchmark_execution.get_nodes()
        node = next(n for n in nodes if n.get_software_environment() == "zlib_new")

        # Act
        result = node.get_environment_path(
            "zlib_new",
            SoftwareBackendEnum.envmodules,
            benchmark_execution.context.directory,
        )

        # Assert
        assert result is not None
        assert result == "zlib/1.3.1"

    @pytest.mark.short
    def test_node_get_environment_path_with_apptainer_backend(
        self, benchmark_execution
    ):
        """Test that node.get_environment_path() returns apptainer path for apptainer backend."""
        nodes = benchmark_execution.get_nodes()
        node = next(n for n in nodes if n.get_software_environment() == "zlib_old")

        # Act
        result = node.get_environment_path(
            "zlib_old",
            SoftwareBackendEnum.apptainer,
            benchmark_execution.context.directory,
        )

        # Assert
        assert result is not None
        assert result == "http://registry.ch/notavailable.sif"

    @pytest.mark.short
    def test_node_get_environment_path_with_nonexistent_env(self, benchmark_execution):
        """Test that node.get_environment_path() returns None for non-existent environment."""
        nodes = benchmark_execution.get_nodes()
        node = nodes[0]

        # Act
        result = node.get_environment_path(
            "nonexistent_env",
            SoftwareBackendEnum.conda,
            benchmark_execution.context.directory,
        )

        # Assert
        assert result is None

    @pytest.mark.short
    def test_node_get_environment_path_with_host_backend(self, benchmark_execution):
        """Test that node.get_environment_path() returns None for host backend."""
        nodes = benchmark_execution.get_nodes()
        node = next(n for n in nodes if n.get_software_environment() == "zlib_old")

        # Act
        result = node.get_environment_path(
            "zlib_old", SoftwareBackendEnum.host, benchmark_execution.context.directory
        )

        # Assert - host backend should return None
        assert result is None


class TestBenchmarkExecutionEnvironmentPath:
    """Test BenchmarkExecution.get_environment_path() method."""

    @pytest.fixture
    def benchmark_execution(self):
        """Load a real benchmark from test data."""
        benchmark_file_path = Path(__file__).parent / "../data/Benchmark_001.yaml"
        return BenchmarkExecution(benchmark_file_path)

    @pytest.mark.short
    def test_execution_get_environment_path_with_conda_backend(
        self, benchmark_execution
    ):
        """Test that execution.get_environment_path() returns conda path for conda backend."""
        # Act
        result = benchmark_execution.get_environment_path(
            "zlib_old", SoftwareBackendEnum.conda
        )

        # Assert
        assert result is not None
        assert "envs/zlib_1.2.11.yaml" in result

    @pytest.mark.short
    def test_execution_get_environment_path_with_envmodules_backend(
        self, benchmark_execution
    ):
        """Test that execution.get_environment_path() returns envmodule for envmodules backend."""
        # Act
        result = benchmark_execution.get_environment_path(
            "zlib_new", SoftwareBackendEnum.envmodules
        )

        # Assert
        assert result is not None
        assert result == "zlib/1.3.1"

    @pytest.mark.short
    def test_execution_get_environment_path_with_docker_backend(
        self, benchmark_execution
    ):
        """Test that execution.get_environment_path() returns apptainer path for docker backend."""
        # Act
        result = benchmark_execution.get_environment_path(
            "zlib_new", SoftwareBackendEnum.docker
        )

        # Assert
        assert result is not None
        assert result == "http://registry.ch/notavailable.sif"

    @pytest.mark.short
    def test_execution_get_environment_path_with_nonexistent_env(
        self, benchmark_execution
    ):
        """Test that execution.get_environment_path() returns None for non-existent environment."""
        # Act
        result = benchmark_execution.get_environment_path(
            "nonexistent_env", SoftwareBackendEnum.conda
        )

        # Assert
        assert result is None

    @pytest.mark.short
    def test_execution_get_environment_path_with_host_backend(
        self, benchmark_execution
    ):
        """Test that execution.get_environment_path() returns None for host backend."""
        # Act
        result = benchmark_execution.get_environment_path(
            "zlib_old", SoftwareBackendEnum.host
        )

        # Assert - host backend should return None
        assert result is None
