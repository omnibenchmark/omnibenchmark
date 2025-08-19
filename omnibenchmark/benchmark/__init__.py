from .benchmark import BenchmarkExecution
from .benchmark import BenchmarkExecution as Benchmark  # Compatibility alias
from ._node import BenchmarkNode
from pathlib import Path
from omnibenchmark.model import SoftwareBackendEnum, SoftwareEnvironment


class Validator:
    """Utility class for validation and path resolution."""

    @staticmethod
    def get_environment_path(
        software_backend: SoftwareBackendEnum,
        environment: SoftwareEnvironment,
        benchmark_dir: Path,
    ) -> str:
        """Get the environment path based on software backend and environment configuration."""
        if software_backend == SoftwareBackendEnum.conda:
            env_path = environment.conda
        elif software_backend == SoftwareBackendEnum.docker:
            env_path = environment.docker
        elif software_backend == SoftwareBackendEnum.apptainer:
            env_path = environment.apptainer
        elif software_backend == SoftwareBackendEnum.envmodules:
            env_path = environment.envmodule
        else:
            env_path = None

        if not env_path:
            raise ValueError(
                f"No environment configuration found for backend {software_backend}"
            )

        # If it's a relative path, resolve it relative to benchmark directory
        if not Path(env_path).is_absolute() and "://" not in env_path:
            return str(benchmark_dir / env_path)
        else:
            return env_path


__all__ = ["BenchmarkExecution", "Benchmark", "BenchmarkNode", "Validator"]
