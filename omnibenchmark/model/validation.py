"""Validation logic for benchmark models."""

import os
import warnings
from typing import List, TYPE_CHECKING, Any

if TYPE_CHECKING:
    pass


class ValidationError(Exception):
    """Exception raised for validation errors."""

    def __init__(self, *errors: Any) -> None:
        self.errors = errors[0] if len(errors) == 1 else errors
        super().__init__(*errors)

    def __str__(self) -> str:
        if isinstance(self.errors, str):
            return self.errors
        else:
            return "\n".join(map(str, self.errors))


def _find_duplicates(items: List[str]) -> List[str]:
    """Find duplicate items in a list."""
    from collections import Counter

    counts = Counter(items)
    return [item for item, count in counts.items() if count > 1]


def _is_url(string: str) -> bool:
    """Check if the string is a valid URL using urlparse."""
    from urllib.parse import urlparse

    try:
        result = urlparse(string)
        return all([result.scheme, result.netloc])
    except ValueError:
        return False


class BenchmarkValidator:
    """Base validator class containing pure model validation for benchmarks.

    This is a mixin class that expects the implementing class to provide:
    - stages: List of stage objects
    - software_environments: List of software environment objects
    - metric_collectors: Optional list of metric collector objects
    - software_backend: Software backend configuration
    - get_modules(): Method returning dict of modules
    - get_outputs(): Method returning dict of outputs
    """

    def validate_model_structure(self) -> None:
        """Validate pure model structure without execution context."""
        errors: List[str] = []

        # 1. Validate unique IDs
        stage_ids = [stage.id for stage in self.stages]  # type: ignore
        duplicate_stage_ids = _find_duplicates(stage_ids)  # type: ignore
        if duplicate_stage_ids:
            errors.append(
                f"Found duplicate stage ids: {', '.join(duplicate_stage_ids)}"
            )

        all_modules = self.get_modules()  # type: ignore
        module_ids = list(all_modules.keys())  # type: ignore
        duplicate_module_ids = _find_duplicates(module_ids)  # type: ignore
        if duplicate_module_ids:
            errors.append(
                f"Found duplicate module ids: {', '.join(duplicate_module_ids)}"
            )

        all_outputs = self.get_outputs()  # type: ignore
        output_ids = list(all_outputs.keys())  # type: ignore
        duplicate_output_ids = _find_duplicates(output_ids)  # type: ignore
        if duplicate_output_ids:
            errors.append(
                f"Found duplicate output ids: {', '.join(duplicate_output_ids)}"
            )

        # 2. Validate output file paths are relative
        for output in all_outputs.values():  # type: ignore
            if output.path.strip() == "":  # type: ignore
                errors.append(f"Output path for file {output.id} is empty")  # type: ignore
            if os.path.isabs(output.path):  # type: ignore
                errors.append(
                    f"Output path for file {output.id} must be relative, not absolute: {output.path}"  # type: ignore
                )

        # 3. Validate stage inputs reference valid outputs
        for stage in self.stages:  # type: ignore
            if stage.inputs:  # type: ignore
                for input_collection in stage.inputs:  # type: ignore
                    for input_id in input_collection.entries:  # type: ignore
                        if input_id not in output_ids:
                            errors.append(
                                f"Input with id '{input_id}' in stage '{stage.id}' is not valid"  # type: ignore
                            )

        # 4. Validate that software environment references exist
        self._validate_environment_references(errors)

        # Raise error if any validation failed
        if errors:
            raise ValidationError(errors)

    def _validate_environment_references(self, errors: List[str]) -> None:
        """Validate software environment references without checking file paths."""
        env_ids = {env.id for env in self.software_environments}  # type: ignore

        # Check modules
        all_modules = self.get_modules()  # type: ignore
        for module in all_modules.values():  # type: ignore
            if module.software_environment not in env_ids:  # type: ignore
                errors.append(
                    f"Software environment with id '{module.software_environment}' is not defined."  # type: ignore
                )

        # Check metric collectors
        if self.metric_collectors:  # type: ignore
            for collector in self.metric_collectors:  # type: ignore
                if collector.software_environment not in env_ids:  # type: ignore
                    errors.append(
                        f"Software environment with id '{collector.software_environment}' for metric collector '{collector.id}' is not defined."  # type: ignore
                    )

                # Validate metric collector inputs
                all_outputs = self.get_outputs()  # type: ignore
                output_ids = list(all_outputs.keys())  # type: ignore
                for collector_input in collector.inputs:  # type: ignore
                    if collector_input.id not in output_ids:  # type: ignore
                        errors.append(
                            f"Input with id '{collector_input.id}' for metric collector '{collector.id}' is not valid."  # type: ignore
                        )

    def validate_software_environments(self) -> None:
        """
        Validate that all software environment references exist.

        Raises:
            ValueError: If any validation errors are found
        """
        # Import here to avoid circular imports
        from .benchmark import SoftwareBackendEnum

        errors: List[str] = []
        env_ids = {env.id for env in self.software_environments}  # type: ignore
        used_env_ids: set[str] = set()

        # Check modules
        for stage in self.stages:  # type: ignore
            for module in stage.modules:  # type: ignore
                if module.software_environment not in env_ids:  # type: ignore
                    errors.append(
                        f"Module '{module.id}' references undefined software environment: '{module.software_environment}'"  # type: ignore
                    )
                else:
                    used_env_ids.add(module.software_environment)  # type: ignore

        # Check metric collectors
        if self.metric_collectors:  # type: ignore
            for collector in self.metric_collectors:  # type: ignore
                if collector.software_environment not in env_ids:  # type: ignore
                    errors.append(
                        f"Metric collector '{collector.id}' references undefined software environment: '{collector.software_environment}'"  # type: ignore
                    )
                else:
                    used_env_ids.add(collector.software_environment)  # type: ignore

        # Validate backend-specific configurations
        for env in self.software_environments:  # type: ignore
            if self.software_backend == SoftwareBackendEnum.conda:  # type: ignore
                if not env.conda:  # type: ignore
                    errors.append(
                        f"Conda backend requires conda configuration for environment '{env.id}'"  # type: ignore
                    )
            elif self.software_backend == SoftwareBackendEnum.docker:  # type: ignore
                if not env.apptainer and not env.docker:  # type: ignore
                    errors.append(
                        f"Docker backend requires apptainer configuration for environment '{env.id}'"  # type: ignore
                    )
            elif self.software_backend == SoftwareBackendEnum.apptainer:  # type: ignore
                if not env.apptainer:  # type: ignore
                    errors.append(
                        f"Apptainer backend requires apptainer configuration for environment '{env.id}'"  # type: ignore
                    )
            elif self.software_backend == SoftwareBackendEnum.envmodules:  # type: ignore
                if not env.envmodule:  # type: ignore
                    errors.append(
                        f"Envmodules backend requires envmodule configuration for environment '{env.id}'"  # type: ignore
                    )

        # Check for unused environments
        unused_envs = env_ids - used_env_ids  # type: ignore
        for unused_env in unused_envs:  # type: ignore
            warnings.warn(
                f"Software environment '{unused_env}' is defined but not used",  # type: ignore
                UserWarning,
            )

        if errors:
            raise ValueError(
                f"Software environment validation failed: {'; '.join(errors)}. Environment not defined."
            )

    # Note: The old validate_structure method has been split into:
    # - validate_model_structure() for pure model validation (above)
    # - validate_execution_context() in the Benchmark class for path validation

    def is_initial(self, stage: Any) -> bool:
        """Check if a stage is an initial stage (has no inputs)."""
        return not stage.inputs or len(stage.inputs) == 0  # type: ignore
