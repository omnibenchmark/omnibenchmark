"""Validation logic for benchmark models."""

import os
from pathlib import Path
from typing import List, TYPE_CHECKING, Any, Optional
from urllib.parse import urlparse

if TYPE_CHECKING:
    from .benchmark import SoftwareEnvironment, SoftwareBackendEnum


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


class BenchmarkParseError(Exception):
    """Exception raised during benchmark parsing with location information.

    This exception carries additional context about where in the YAML file
    the error occurred, allowing the CLI layer to format it appropriately.
    """

    def __init__(
        self,
        message: str,
        yaml_file: Optional[Path] = None,
        line_number: Optional[int] = None,
        stage_id: Optional[str] = None,
        module_id: Optional[str] = None,
        parameter_index: Optional[int] = None,
        values: Optional[List[Any]] = None,
        original_error: Optional[Exception] = None,
    ):
        self.message = message
        self.yaml_file = yaml_file
        self.line_number = line_number
        self.stage_id = stage_id
        self.module_id = module_id
        self.parameter_index = parameter_index
        self.values = values
        self.original_error = original_error
        super().__init__(message)

    def __str__(self) -> str:
        """Format a compact, helpful error message with context."""
        parts = []

        # Add file and line number if available
        if self.yaml_file:
            location = str(self.yaml_file)
            if self.line_number:
                location += f":{self.line_number}"
            parts.append(f"Error in {location}")

        # Add stage/module context if available
        if self.stage_id or self.module_id:
            context_parts = []
            if self.stage_id:
                context_parts.append(f"stage '{self.stage_id}'")
            if self.module_id:
                context_parts.append(f"module '{self.module_id}'")
            if self.parameter_index is not None:
                context_parts.append(f"parameter #{self.parameter_index}")
            parts.append(f"  in {', '.join(context_parts)}")

        # Add the actual error message
        parts.append(f"  {self.message}")

        # Add parameter values if available (helps with debugging)
        if self.values:
            parts.append(f"  values: {self.values}")

        return "\n".join(parts)


def _find_duplicate_ids(items: List[str]) -> List[str]:
    """Find duplicate IDs in a list for validation purposes."""
    from collections import Counter

    counts = Counter(items)
    return [item for item, count in counts.items() if count > 1]


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
        duplicate_stage_ids = _find_duplicate_ids(stage_ids)  # type: ignore
        if duplicate_stage_ids:
            errors.append(
                f"Found duplicate stage ids: {', '.join(duplicate_stage_ids)}"
            )

        all_modules = self.get_modules()  # type: ignore
        module_ids = list(all_modules.keys())  # type: ignore
        duplicate_module_ids = _find_duplicate_ids(module_ids)  # type: ignore
        if duplicate_module_ids:
            errors.append(
                f"Found duplicate module ids: {', '.join(duplicate_module_ids)}"
            )

        all_outputs = self.get_outputs()  # type: ignore
        output_ids = list(all_outputs.keys())  # type: ignore
        duplicate_output_ids = _find_duplicate_ids(output_ids)  # type: ignore
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
        self._validate_software_environments(errors)

        # Raise error if any validation failed
        if errors:
            raise ValidationError(errors)

    def _validate_software_environments(self, errors: List[str]) -> None:
        """Validate software environment references without checking file paths."""
        env_ids = {env.id for env in self.software_environments}  # type: ignore

        # Check modules
        all_modules = self.get_modules()  # type: ignore
        for module in all_modules.values():  # type: ignore
            if module.software_environment not in env_ids:  # type: ignore
                errors.append(
                    f"Software environment with id '{module.software_environment}' is not declared. It should be listed as part of the stanza software_environments within the benchmarking YAML header."  # type: ignore
                )

        # Check metric collectors
        if self.metric_collectors:  # type: ignore
            for collector in self.metric_collectors:  # type: ignore
                if collector.software_environment not in env_ids:  # type: ignore
                    errors.append(
                        f"Software environment with id '{collector.software_environment}' for metric collector '{collector.id}' is not declared."  # type: ignore
                    )

                # Validate metric collector inputs
                all_outputs = self.get_outputs()  # type: ignore
                output_ids = list(all_outputs.keys())  # type: ignore
                for collector_input in collector.inputs:  # type: ignore
                    # Handle both string and IOFile inputs
                    if isinstance(collector_input, str):
                        input_id = collector_input
                    else:
                        input_id = collector_input.id  # type: ignore

                    if input_id not in output_ids:  # type: ignore
                        errors.append(
                            f"Input with id '{input_id}' for metric collector '{collector.id}' is not valid."  # type: ignore
                        )

    # Note: The old validate_structure method has been split into:
    # - validate_model_structure() for pure model validation (above)
    # - validate_execution_context() in the Benchmark class for path validation

    # Utility methods

    @staticmethod
    def is_initial(stage: Any) -> bool:
        """Check if a stage is an initial stage (has no inputs)."""
        return not stage.inputs or len(stage.inputs) == 0  # type: ignore

    @staticmethod
    def is_url(string: str) -> bool:
        """Check if the string is a valid URL using urlparse."""
        try:
            result = urlparse(string)
            return all(
                [result.scheme, result.netloc]
            )  # Valid if both scheme and netloc are present
        except ValueError:
            return False

    @staticmethod
    def is_absolute_path(string: str) -> bool:
        """Check if the string is an absolute path."""
        return Path(string).is_absolute()

    @staticmethod
    def get_environment_path(
        software_backend: "SoftwareBackendEnum",
        software: "SoftwareEnvironment",
        benchmark_dir: Path,
    ) -> Optional[str]:
        """Get the environment path based on software backend and environment configuration."""
        # Import here to avoid circular imports
        from .benchmark import SoftwareBackendEnum

        environment = None
        if (
            software_backend == SoftwareBackendEnum.apptainer
            or software_backend == SoftwareBackendEnum.docker
        ):
            environment = software.apptainer

        elif software_backend == SoftwareBackendEnum.conda:
            environment = software.conda

        elif software_backend == SoftwareBackendEnum.envmodules:
            environment = software.envmodule

        if not environment:
            return None

        if BenchmarkValidator.is_url(
            environment
        ) or BenchmarkValidator.is_absolute_path(environment):
            environment_path = environment
        elif software_backend == SoftwareBackendEnum.envmodules:
            environment_path = environment
        else:
            environment_path = os.path.join(benchmark_dir, environment)

        return environment_path
