from pathlib import Path
from typing import Optional
from omnibenchmark.model import SoftwareBackendEnum, SoftwareEnvironment


class Validator:
    """Utility class for validation and path resolution."""

    @staticmethod
    def get_environment_path(
        software_backend: SoftwareBackendEnum,
        environment: SoftwareEnvironment,
        benchmark_dir: Path,
    ) -> Optional[str]:
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
            # Return None to allow Snakemake fallback (e.g., "conda_not_provided.yml")
            return None

        # For envmodules, return the module name as-is (it's not a file path)
        if software_backend == SoftwareBackendEnum.envmodules:
            return env_path

        # If it's a relative path, resolve it relative to benchmark directory
        if not Path(env_path).is_absolute() and "://" not in env_path:
            return str(benchmark_dir / env_path)
        else:
            return env_path
