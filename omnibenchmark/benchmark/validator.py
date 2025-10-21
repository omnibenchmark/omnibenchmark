"""Validator utilities - delegates to BenchmarkValidator in model.validation."""

from pathlib import Path
from typing import Optional
from omnibenchmark.model import SoftwareBackendEnum, SoftwareEnvironment
from omnibenchmark.model.validation import BenchmarkValidator


class Validator:
    """Utility class for validation and path resolution.

    This class now delegates to BenchmarkValidator in omnibenchmark.model.validation.
    Kept for backward compatibility.
    """

    @staticmethod
    def get_environment_path(
        software_backend: SoftwareBackendEnum,
        environment: SoftwareEnvironment,
        benchmark_dir: Path,
    ) -> Optional[str]:
        """Get the environment path based on software backend and environment configuration.

        Delegates to BenchmarkValidator.get_environment_path.
        """
        return BenchmarkValidator.get_environment_path(
            software_backend, environment, benchmark_dir
        )
