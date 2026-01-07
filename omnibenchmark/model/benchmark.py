"""Pydantic models for Omnibenchmark."""

import os
import re
import warnings
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import yaml
from pydantic import (
    BaseModel,
    Field,
    field_validator,
    model_validator,
)
from pydantic import (
    ValidationError as PydanticValidationError,
)

from omnibenchmark.utils import merge_dict_list  # type: ignore[import]

from ._parsing import convert_pydantic_error_to_parse_error
from .validation import BenchmarkParseError, BenchmarkValidator, ValidationError


class LineNumberLoader(yaml.SafeLoader):
    """Custom YAML loader that tracks line numbers for better error reporting."""

    def compose_node(self, parent, index):
        # Get the line number where the node starts
        line = self.line + 1  # YAML lines are 0-indexed, we need 1-indexed for the user
        node = super().compose_node(parent, index)
        # I use setattr to avoid type checker errors with dynamic attributes (__line__ is our custom attribute)
        setattr(node, "__line__", line)
        return node


def load_yaml_with_lines(file_path: Path) -> tuple[Dict[str, Any], Dict[str, int]]:
    """
    Load YAML file and track line numbers for error reporting.

    Returns:
        tuple: (data dict, line_map dict mapping paths to line numbers)
    """
    with open(file_path, "r") as f:
        yaml_content = f.read()
        loader = LineNumberLoader(yaml_content)
        data = loader.get_single_data()

    # Build a map of paths to line numbers
    line_map = {}

    def traverse(node, path=""):
        """Traverse YAML nodes and build line map."""
        if hasattr(node, "__line__"):
            line_map[path] = node.__line__

        if node.id == "mapping":
            # Dictionary node
            for key_node, value_node in node.value:
                key = key_node.value
                new_path = f"{path}.{key}" if path else key
                if hasattr(value_node, "__line__"):
                    line_map[new_path] = value_node.__line__
                traverse(value_node, new_path)
        elif node.id == "sequence":
            # List node
            for i, item_node in enumerate(node.value):
                new_path = f"{path}[{i}]"
                if hasattr(item_node, "__line__"):
                    line_map[new_path] = item_node.__line__
                traverse(item_node, new_path)

    # Parse again with line tracking
    with open(file_path, "r") as f:
        loader = LineNumberLoader(f)
        root_node = loader.get_single_node()
        traverse(root_node)

    return data, line_map


def _is_environment_url(string: str) -> bool:
    """Check if an environment path string is a valid URL using urlparse."""
    from urllib.parse import urlparse

    try:
        result = urlparse(string)
        return all([result.scheme, result.netloc])
    except ValueError:
        return False


# Validators
def validate_non_empty_string(v: str) -> str:
    """Validate that a string is not empty."""
    if not v or not v.strip():
        raise ValueError("must be a non-empty string")
    return v


def validate_non_empty_commit(v: str) -> str:
    """Validate commit hash is not empty."""
    if not v or not v.strip():
        raise ValueError("Commit must be a non-empty string")
    return v


def validate_hex_string(v: str) -> str:
    """Validate that a string is a valid hex string."""
    if not v:
        raise ValueError("must be a valid hexadecimal string")
    try:
        int(v, 16)
    except ValueError:
        raise ValueError("must be a valid hexadecimal string")
    return v


# Enums
class APIVersion(str, Enum):
    """API version enum."""

    V0_1_0 = "0.1.0"
    V0_2_0 = "0.2.0"
    V0_3_0 = "0.3.0"

    @classmethod
    def latest(cls) -> "APIVersion":
        return cls.V0_3_0

    @classmethod
    def supported_versions(cls) -> set[str]:
        return {version.value for version in cls}


class SoftwareBackendEnum(str, Enum):
    """Software backend types."""

    apptainer = "apptainer"
    conda = "conda"
    docker = "docker"
    envmodules = "envmodules"
    host = "host"


class RepositoryType(str, Enum):
    """Repository types."""

    git = "git"


class StorageAPIEnum(str, Enum):
    """Storage API types. Currently only S3 is supported."""

    s3 = "S3"


# Base models
class IdentifiableEntity(BaseModel):
    """Base class for entities with an ID."""

    id: str = Field(..., description="Unique identifier")


class DescribableEntity(IdentifiableEntity):
    """Base class for entities with ID, name, and description."""

    name: Optional[str] = Field(None, description="Human-readable name")
    description: Optional[str] = Field(None, description="Description of the entity")


class Repository(BaseModel):
    """Repository information."""

    url: str = Field(..., description="Repository URL")
    commit: str = Field(..., description="Commit hash")

    @field_validator("url")
    @classmethod
    def validate_url(cls, v: str) -> str:
        return validate_non_empty_string(v)

    @field_validator("commit", mode="before")
    @classmethod
    def validate_commit(cls, v: Union[str, int, float]) -> str:
        # Convert numeric commit hashes to strings. Yes, they should not be numeric,
        # but we should handle it gracefully.
        v_str = str(v) if not isinstance(v, str) else v
        return validate_non_empty_commit(v_str)

    type: RepositoryType = Field(RepositoryType.git, description="Repository type")


class Storage(BaseModel):
    """
    Storage configuration for remote storage of benchmark artifacts.
    This is intended for the benchmarker role, since they will have to
    provide the needed credentials to write in the store.
    """

    api: StorageAPIEnum = Field(
        StorageAPIEnum.s3,
        description="Storage API type (currently only S3 is supported)",
    )
    endpoint: str = Field(..., description="Storage endpoint URL")
    bucket_name: Optional[str] = Field(
        None, description="Storage bucket name (only needed if using S3)"
    )


class Parameter(BaseModel):
    """Parameter definition - supports both CLI args format and key-value dict format.

    Two formats are supported:

    1. Legacy CLI args format (values):
       parameters:
         - values: ["--method", "cosine"]

    2. Key-value dict format (params):
       parameters:
         - method: genie
           threshold: 0.1
         - method: other
           k: [1, 2, 3]  # Expands to 3 runs with k=1, k=2, k=3

    IMPORTANT: In the dict format, any parameter value that is a list will trigger
    implicit Cartesian product expansion. Lists cannot be passed as literal values -
    they are always interpreted as grid expansion parameters.

    Example of Cartesian product:
       parameters:
         - alpha: [0.1, 0.5]  # 2 values
           beta: [1, 2, 3]    # 3 values
       # Expands to 6 parameter combinations: (0.1,1), (0.1,2), (0.1,3), (0.5,1), (0.5,2), (0.5,3)
    """

    id: str = Field(..., description="Parameter ID")
    values: Optional[List[str]] = Field(
        None, description="Parameter values as CLI args (legacy format)"
    )
    # NOTE: The 'params' field is an internal representation only. User should never
    # write 'params:' in YAML. When YAML like "method: genie" is parsed, the
    # @model_validator below automatically converts it to
    # params={"method": "genie"} internally.
    params: Optional[Dict[str, Any]] = Field(
        None,
        description="Parameter key-value pairs (supports lists for product expansion)",
    )

    @model_validator(mode="before")
    @classmethod
    def validate_parameter_format(cls, data: Any) -> Any:
        """Validate and normalize parameter format."""
        if isinstance(data, dict):
            has_values = "values" in data and data["values"] is not None
            has_params = "params" in data and data["params"] is not None

            # If neither is provided, check if data itself looks like a param dict
            if not has_values and not has_params:
                # Check if data has keys other than 'id'
                keys = set(data.keys()) - {"id"}
                if keys:
                    # Treat remaining keys as params
                    params_data = {k: data[k] for k in keys}
                    data["params"] = params_data
                    has_params = True

            if not has_values and not has_params:
                raise ValueError(
                    "Parameter must have either 'values' or parameter keys"
                )

            if has_values and has_params:
                raise ValueError(
                    "Parameter cannot have both 'values' and parameter keys"
                )

        return data

    @field_validator("values", mode="before")
    @classmethod
    def validate_values(cls, v: Optional[List[str]]) -> Optional[List[str]]:
        if v is None:
            return None
        if not v:
            raise ValueError("Parameter values cannot be empty")
        # Convert all values to strings to handle mixed types (float, int, str)
        return [str(val) for val in v]

    @field_validator("params")
    @classmethod
    def validate_params(cls, v: Optional[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
        if v is not None and not v:
            raise ValueError("Parameter params cannot be empty")
        return v


class IOFile(IdentifiableEntity):
    """Input/Output file definition."""

    path: str = Field(..., description="File path")

    @field_validator("path")
    @classmethod
    def validate_path(cls, v: str) -> str:
        return validate_non_empty_string(v)


class InputCollection(BaseModel):
    """Collection of input file IDs."""

    entries: List[str] = Field(..., description="List of input file IDs")


class SoftwareEnvironment(DescribableEntity):
    """Software environment configuration."""

    repository: Optional[Repository] = None
    conda: Optional[str] = Field(None, description="Conda environment file path")
    apptainer: Optional[str] = Field(
        None, description="Apptainer/Singularity image path"
    )
    docker: Optional[str] = Field(None, description="Docker image")
    envmodule: Optional[str] = Field(None, description="Environment module name")
    easyconfig: Optional[str] = Field(None, description="EasyBuild config path")

    @model_validator(mode="after")
    def validate_backend_config(self) -> "SoftwareEnvironment":
        """Ensure at least one backend configuration is provided."""
        # Temporarily disabled to allow benchmark-level validation
        # backends = [
        #     self.conda,
        #     self.apptainer,
        #     self.docker,
        #     self.envmodule,
        #     self.easyconfig,
        # ]
        # if not any(backends):
        #     raise ValueError(
        #         "At least one software backend must be configured (conda, apptainer, docker, envmodule, or easyconfig)"
        #     )
        return self


class SoftwareEnvironmentReference(BaseModel):
    """Reference to a software environment."""

    software_environment: str = Field(..., description="Software environment ID")


class Module(DescribableEntity, SoftwareEnvironmentReference):
    """Module definition."""

    repository: Repository = Field(..., description="Module repository")
    software_environment: str = Field(..., description="Software environment ID")
    parameters: Optional[List[Parameter]] = Field(None, description="Module parameters")
    exclude: Optional[List[str]] = Field(None, description="Paths to exclude")
    outputs: Optional[List[IOFile]] = Field(None, description="Module outputs")

    def has_environment_reference(self, env_id: Optional[str] = None) -> bool:
        """Check if module has a software environment reference."""
        if env_id is None:
            return bool(self.software_environment)
        return self.software_environment == env_id


class MetricCollector(DescribableEntity, SoftwareEnvironmentReference):
    """Metric collector definition."""

    repository: Repository = Field(..., description="Repository information")
    software_environment: str = Field(..., description="Software environment ID")
    inputs: List[Union[str, IOFile]] = Field(..., description="Input files")
    outputs: List[IOFile] = Field(..., description="Output files")

    def has_environment_reference(self, env_id: Optional[str] = None) -> bool:
        """Check if metric collector has a software environment reference."""
        if env_id is None:
            return bool(self.software_environment)
        return self.software_environment == env_id


class Stage(DescribableEntity):
    """Stage definition."""

    modules: List[Module] = Field(..., description="Modules in this stage")
    inputs: Optional[List[InputCollection]] = Field(None, description="Stage inputs")
    outputs: List[IOFile] = Field(..., description="Stage outputs")

    @field_validator("inputs", mode="before")
    @classmethod
    def validate_inputs(cls, v, info):
        """Convert legacy string format to InputCollection."""
        if v is None:
            return v

        # Check if all items are strings (non-legacy format: a simple list of input IDs)
        # This should be treated as a single InputCollection with all items
        if all(isinstance(item, str) for item in v):
            return [InputCollection(entries=v)]

        result = []
        for idx, item in enumerate(v):
            # If it's a dict with 'entries' (legacy format with explicit entries field)
            if isinstance(item, dict) and "entries" in item:
                import warnings

                # Try to get line number from context
                line_info = ""
                if info.context and "line_map" in info.context:
                    line_map = info.context["line_map"]
                    # Search for any key containing inputs[idx].entries
                    # We don't know the stage index, so we search through all keys
                    search_pattern = f"inputs[{idx}].entries"
                    for key in line_map:
                        if search_pattern in key:
                            line_info = f" (line {line_map[key]})"
                            break

                warnings.warn(
                    f"Using 'entries' field in inputs is deprecated{line_info}. "
                    "Please use a simple list of strings instead: "
                    "inputs: ['data.matrix', 'data.true_labels'] instead of "
                    "inputs: [{entries: ['data.matrix', 'data.true_labels']}]. "
                    "Support for 'entries' will be removed in a future release.",
                    FutureWarning,
                    stacklevel=2,
                )
                result.append(item)
            # If it's already an InputCollection or other dict, keep as-is
            else:
                result.append(item)

        return result


class Benchmark(DescribableEntity, BenchmarkValidator):
    """Main benchmark definition."""

    @model_validator(mode="after")
    def set_authors_if_missing(self):
        # If authors is missing or empty, set it to [benchmarker]
        if not self.authors or len(self.authors) == 0:
            self.authors = [self.benchmarker]
        return self

    benchmarker: str = Field(..., description="Benchmark author")
    authors: Optional[List[str]] = Field(None, description="List of benchmark authors")
    version: str = Field(..., description="Benchmark version")
    software_backend: SoftwareBackendEnum = Field(..., description="Software backend")
    software_environments: List[SoftwareEnvironment] = Field(
        ..., description="Available software environments"
    )
    stages: List[Stage] = Field(..., description="Benchmark stages")
    metric_collectors: Optional[List[MetricCollector]] = Field(
        None, description="Metric collectors"
    )
    storage: Optional[Storage] = Field(None, description="Remote storage configuration")
    # Legacy Compatibility Fields - Migration Strategy
    # These fields maintain backward compatibility during the LinkML → Pydantic transition.
    # They're handled via property methods that delegate to the new structured fields.
    # Future consideration: Remove after confirming no external YAML files use these fields.
    storage_api: Optional[StorageAPIEnum] = Field(
        None, description="Storage API type (deprecated, use storage.api)"
    )
    storage_bucket_name: Optional[str] = Field(
        None, description="Storage bucket name (deprecated, use storage.bucket_name)"
    )
    benchmark_yaml_spec: Optional[Union[str, float, int]] = Field(
        None,
        description="Benchmark YAML specification version (deprecated, use api_version)",
    )
    api_version: APIVersion = Field(APIVersion.V0_3_0, description="API version")

    @field_validator("api_version", mode="before")
    @classmethod
    def validate_api_version(cls, v: Union[str, APIVersion]) -> APIVersion:
        """Convert string API version to enum."""
        if isinstance(v, str):
            for version in APIVersion:
                if version.value == v:
                    return version
            # If no match found, raise ValueError
            raise ValueError(f"Invalid API version: {v}")
        return v

    @field_validator("version", mode="before")
    @classmethod
    def validate_version(cls, v: Union[str, int, float], info) -> str:
        """Validate that version follows strict semantic versioning format."""
        # Handle numeric types with deprecation warning
        if isinstance(v, (int, float)):
            # Try to get line number from context
            line_info = ""
            if info.context and "line_map" in info.context:
                line_map = info.context["line_map"]
                if "version" in line_map:
                    line_info = f" (line {line_map['version']})"

            warnings.warn(
                f"Field 'version' should be a string{line_info}. "
                f"Found {type(v).__name__} value: {v}. "
                f"This will not be valid in a future release. "
                f"Please update your YAML file to use a quoted string.",
                FutureWarning,
                stacklevel=4,
            )
            v = str(v)

        if not isinstance(v, str):
            raise ValueError("Version must be a string")

        if not v or not v.strip():
            raise ValueError("Version must be a non-empty string")

        v = v.strip()

        # Match semantic version pattern: x.y.z or x.y
        # No leading zeros allowed (except for "0" itself)
        pattern = r"^(0|[1-9]\d*)\.(0|[1-9]\d*)(?:\.(0|[1-9]\d*))?$"
        if not re.match(pattern, v):
            raise ValueError(
                f"Version '{v}' does not follow strict semantic versioning format. "
                "Expected x.y.z or x.y where x, y, z are non-negative integers without leading zeros."
            )

        return v

    @field_validator("benchmark_yaml_spec", mode="before")
    @classmethod
    def validate_benchmark_yaml_spec(
        cls, v: Optional[Union[str, float, int]], info
    ) -> Optional[str]:
        """Validate benchmark_yaml_spec format if provided."""
        if v is None:
            return v

        if not isinstance(v, str):
            # Allow float/int for backwards compatibility but warn
            if isinstance(v, (float, int)):
                # Try to get line number from context
                line_info = ""
                if info.context and "line_map" in info.context:
                    line_map = info.context["line_map"]
                    if "benchmark_yaml_spec" in line_map:
                        line_info = f" (line {line_map['benchmark_yaml_spec']})"

                warnings.warn(
                    f"Field 'benchmark_yaml_spec' should be a string{line_info}. "
                    f"Found {type(v).__name__} value: {v}. "
                    f"This will not be valid in a future release. "
                    f"Please update your YAML file to use a quoted string.",
                    FutureWarning,
                    stacklevel=4,
                )
                return str(v)
            else:
                raise ValueError("benchmark_yaml_spec must be a string")

        return v

    # _benchmark_dir: Optional[Path] = None

    @classmethod
    def from_yaml(cls, path_or_content: Union[str, Path]) -> "Benchmark":
        """Load benchmark from YAML file or string content."""
        line_map = {}
        yaml_file_path = None

        if isinstance(path_or_content, (str, Path)) and "\n" not in str(
            path_or_content
        ):
            # Treat as file path
            yaml_file_path = Path(path_or_content)
            data, line_map = load_yaml_with_lines(yaml_file_path)
        else:
            # Treat as YAML content string
            data = yaml.safe_load(str(path_or_content))

        # Convert dict-style software_environments to list format
        if "software_environments" in data and isinstance(
            data["software_environments"], dict
        ):
            envs: List[Dict[str, Any]] = []
            for env_id, env_config in data["software_environments"].items():
                env_dict = dict(env_config) if env_config else {}
                env_dict["id"] = env_id
                envs.append(env_dict)
            data["software_environments"] = envs

        # Generate parameter IDs for modules
        if "stages" in data:
            for stage_idx, stage in enumerate(data["stages"]):
                if "modules" in stage:
                    for module_idx, module in enumerate(stage["modules"]):
                        if "parameters" in module and module["parameters"]:
                            for param_idx, param in enumerate(module["parameters"]):
                                # Emit deprecation warning for old 'values' format
                                if "values" in param:
                                    # Convert all parameter values to strings to handle mixed types (float, int, str)
                                    param["values"] = [str(v) for v in param["values"]]

                                    # Emit structured warning with line context
                                    from omnibenchmark.cli.error_formatting import (
                                        format_yaml_warning,
                                    )
                                    from omnibenchmark.cli.utils.logging import logger

                                    stage_id = stage.get(
                                        "id", f"<stage index {stage_idx}>"
                                    )
                                    module_id = module.get(
                                        "id", f"<module index {module_idx}>"
                                    )

                                    # Find line number for this parameter
                                    param_line = None
                                    if line_map and yaml_file_path:
                                        possible_paths = [
                                            f"stages[{stage_idx}].modules[{module_idx}].parameters[{param_idx}].values",
                                            f"stages[{stage_idx}].modules[{module_idx}].parameters[{param_idx}]",
                                        ]
                                        for path in possible_paths:
                                            if path in line_map:
                                                param_line = line_map[path]
                                                break

                                    warning_msg = format_yaml_warning(
                                        message="Use dict format instead",
                                        yaml_file=yaml_file_path,
                                        line_number=param_line,
                                        stage_id=stage_id,
                                        module_id=module_id,
                                    )
                                    logger.warning(warning_msg)

                                if "id" not in param:
                                    # Generate a hash-based ID from the parameter
                                    import hashlib

                                    try:
                                        if "values" in param:
                                            # Legacy format: hash the values
                                            param_str = str(sorted(param["values"]))
                                        else:
                                            # New format: hash the entire param dict (excluding 'id')
                                            param_items = sorted(
                                                (k, v)
                                                for k, v in param.items()
                                                if k != "id"
                                            )
                                            param_str = str(param_items)

                                        param["id"] = hashlib.sha256(
                                            param_str.encode()
                                        ).hexdigest()[:8]
                                    except Exception as e:
                                        # Provide detailed error information
                                        stage_id = stage.get(
                                            "id", f"<stage index {stage_idx}>"
                                        )
                                        module_id = module.get(
                                            "id", f"<module index {module_idx}>"
                                        )

                                        # Try to find the line number in the YAML file
                                        line_num = None
                                        if line_map and yaml_file_path:
                                            # Try different path patterns to find the line
                                            possible_paths = [
                                                f"stages[{stage_idx}].modules[{module_idx}].parameters[{param_idx}]",
                                                f"stages[{stage_idx}].modules[{module_idx}].parameters[{param_idx}].values",
                                            ]
                                            for path in possible_paths:
                                                if path in line_map:
                                                    line_num = line_map[path]
                                                    break

                                        # Raise a BenchmarkParseError with all context information
                                        # The CLI layer will handle formatting
                                        raise BenchmarkParseError(
                                            message="Failed to process parameter values",
                                            yaml_file=yaml_file_path,
                                            line_number=line_num,
                                            stage_id=stage_id,
                                            module_id=module_id,
                                            parameter_index=param_idx,
                                            values=param.get("values"),
                                            original_error=e,
                                        ) from e

        # Convert string storage to Storage object and handle backward compatibility
        if "storage" in data:
            # Detect if using old format (top-level storage_api or storage_bucket_name)
            has_old_format_fields = (
                "storage_api" in data or "storage_bucket_name" in data
            )
            is_storage_dict = isinstance(data["storage"], dict)

            if isinstance(data["storage"], str):
                # Old format: storage: "https://endpoint.url"
                if not has_old_format_fields:
                    # Allow old format without storage_api/storage_bucket_name (use defaults)
                    pass
                storage_dict = {
                    "api": data.get("storage_api", "S3"),
                    "endpoint": data["storage"],
                }
                if "storage_bucket_name" in data:
                    storage_dict["bucket_name"] = data["storage_bucket_name"]
                data["storage"] = storage_dict
            elif is_storage_dict:
                # New format: storage: {endpoint, api, bucket_name}
                # Check for mixed format (both new dict and old top-level fields)
                if has_old_format_fields:
                    # Reject mixed format
                    line_num = None
                    if line_map and yaml_file_path:
                        # Try to find storage field line number
                        for path in ["storage", "storage_api", "storage_bucket_name"]:
                            if path in line_map:
                                line_num = line_map[path]
                                break

                    raise BenchmarkParseError(
                        message="Mixed storage format detected. Use either old format (storage: 'url', storage_api: 'S3', storage_bucket_name: 'name') OR new format (storage: {endpoint: 'url', api: 'S3', bucket_name: 'name'}), but not both.",
                        yaml_file=yaml_file_path,
                        line_number=line_num,
                        original_error=None,
                    )

        try:
            # Pass line_map through validation context for better error messages
            return cls.model_validate(data, context={"line_map": line_map})
        except PydanticValidationError as e:
            # Convert Pydantic validation error to BenchmarkParseError with context
            parse_error = convert_pydantic_error_to_parse_error(
                e, data, line_map, yaml_file_path
            )
            raise parse_error from e

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Benchmark":
        """Create benchmark from dictionary."""
        return cls(**data)

    def merge_with(self, other: "Benchmark") -> "Benchmark":
        """Merge this benchmark with another benchmark."""
        # Create a new benchmark with combined data
        merged_data = self.model_dump()
        other_data = other.model_dump()

        # Update with other's data (other takes precedence for scalar fields)
        for key, value in other_data.items():
            if key in ["stages", "metric_collectors", "software_environments"]:
                # For list fields, combine them
                if key == "stages":
                    existing_stages = {
                        stage["id"]: stage for stage in merged_data.get(key, [])
                    }
                    other_stages = {
                        stage["id"]: stage for stage in other_data.get(key, [])
                    }
                    existing_stages.update(other_stages)
                    merged_data[key] = list(existing_stages.values())
                elif key == "metric_collectors":
                    existing_collectors = {
                        mc["id"]: mc for mc in merged_data.get(key, [])
                    }
                    other_collectors = {mc["id"]: mc for mc in other_data.get(key, [])}
                    existing_collectors.update(other_collectors)
                    merged_data[key] = list(existing_collectors.values())
                elif key == "software_environments":
                    existing_envs = {env["id"]: env for env in merged_data.get(key, [])}
                    other_envs = {env["id"]: env for env in other_data.get(key, [])}
                    existing_envs.update(other_envs)
                    merged_data[key] = list(existing_envs.values())
            else:
                # For scalar fields, other takes precedence if not None
                if value is not None:
                    merged_data[key] = value

        return Benchmark(**merged_data)

    # Optional storage, with API compatibility
    # (Methods moved to StorageAccessors mixin)

    # API migration

    def upgrade_to_latest(self) -> "Benchmark":
        """Upgrade benchmark to latest API version."""
        # Check current version from either field
        current_version = getattr(self, "benchmark_yaml_spec", None) or str(
            self.api_version
        )
        latest_version = APIVersion.latest()

        if current_version == latest_version:
            return self

        # Create upgraded version
        data = self.model_dump()
        data["api_version"] = latest_version
        data["benchmark_yaml_spec"] = latest_version

        return Benchmark(**data)

    # Compatibility methods for LinkMLConverter interface
    # They can be removed when API migration is complete

    def get_storage_api(self) -> Optional[str]:
        """Get storage API with backward compatibility."""
        if self.storage and self.storage.api:
            return str(self.storage.api.value)
        return str(self.storage_api.value) if self.storage_api else None

    def get_storage_bucket_name(self) -> Optional[str]:
        """Get storage bucket name with backward compatibility."""
        if self.storage and self.storage.bucket_name:
            return self.storage.bucket_name
        return self.storage_bucket_name

    def get_storage_endpoint(self) -> Optional[str]:
        """Get storage endpoint."""
        if self.storage and self.storage.endpoint:
            return self.storage.endpoint
        return None

    def get_name(self) -> str:
        """Get name of the benchmark."""
        return self.name if self.name else self.id

    def get_version(self) -> str:
        """Get version of the benchmark."""
        return self.version

    def get_author(self) -> str:
        """Get author of the benchmark."""
        return self.benchmarker

    def get_software_backend(self) -> SoftwareBackendEnum:
        """Get software backend of the benchmark."""
        return self.software_backend

    def get_software_environments(self) -> Dict[str, SoftwareEnvironment]:
        """Get software environments as a dict."""
        return {env.id: env for env in self.software_environments}

    def get_stages(self) -> Dict[str, Stage]:
        """Get benchmark stages as a dict."""
        return {stage.id: stage for stage in self.stages}

    def get_stage(self, stage_id: str) -> Optional[Stage]:
        """Get stage by stage_id."""
        stages = self.get_stages()
        return stages.get(stage_id)

    def get_stage_by_output(self, output_id: str) -> Optional[Stage]:
        """Get stage that returns output with output_id."""
        for stage in self.stages:
            for output in stage.outputs:
                if output.id == output_id:
                    return stage
        return None

    def get_modules_by_stage(self, stage: Union[str, Stage]) -> Dict[str, Module]:
        """Get modules by stage/stage_id."""
        if isinstance(stage, str):
            stage_obj = self.get_stage(stage)
            if not stage_obj:
                return {}
            return {module.id: module for module in stage_obj.modules}
        return {module.id: module for module in stage.modules}

    def get_stage_implicit_inputs(self, stage: Union[str, Stage]) -> List[List[str]]:
        """Get implicit inputs of a stage by stage/stage_id."""
        if isinstance(stage, str):
            stage_obj = self.get_stage(stage)
            if not stage_obj or not stage_obj.inputs:
                return []
            return [input_col.entries for input_col in stage_obj.inputs]
        if not stage.inputs:
            return []
        return [input_col.entries for input_col in stage.inputs]

    def get_explicit_inputs(self, input_ids: List[str]) -> Dict[str, str]:
        """Get explicit inputs of a stage by input_id(s)."""
        all_stages_outputs: List[Dict[str, str]] = []
        for stage in self.stages:
            stage_outputs = self.get_stage_outputs(stage)
            # Substitute the actual stage_id into the template like the old LinkML code
            stage_outputs = {
                key: value.format(
                    input="{input}",
                    stage=stage.id,  # Substitute actual stage_id here
                    module="{module}",
                    params="{params}",
                    dataset="{dataset}",
                )
                for key, value in stage_outputs.items()
            }
            all_stages_outputs.append(stage_outputs)

        # Merge all stage outputs
        all_outputs = merge_dict_list(all_stages_outputs)

        result: Dict[str, str] = {}
        for input_id in input_ids:
            if input_id in all_outputs:
                result[input_id] = all_outputs[input_id]
            else:
                # Default to a placeholder
                result[input_id] = "{input_id}"

        return result

    def get_stage_outputs(self, stage: Union[str, Stage]) -> Dict[str, str]:
        """Get outputs of a stage by stage/stage_id."""
        if isinstance(stage, str):
            stage_obj = self.get_stage(stage)
            if not stage_obj:
                return {}
            stage = stage_obj

        return {output.id: expand_output_path(output) for output in stage.outputs}

    def get_output_stage(self, output_id: str) -> Optional[Stage]:
        """Get stage that returns output with output_id."""
        return self.get_stage_by_output(output_id)

    def _resolve_module_attr(self, module: Union[str, Module], attr: str):
        """Resolve a module attribute by module/module_id."""
        module_obj = (
            module if isinstance(module, Module) else self.get_modules().get(module)
        )
        return getattr(module_obj, attr, None) if module_obj else None

    def get_module_excludes(self, module: Union[str, Module]) -> Optional[List[str]]:
        """Get module excludes by module/module_id."""
        match module:
            case str():
                module_obj = self.get_modules().get(module)
                if not module_obj or not module_obj.exclude:
                    return None
                return module_obj.exclude
            case Module():
                if not module.exclude:
                    return None
                return module.exclude

    def get_module_parameters(self, module: Union[str, Module]) -> Optional[List[Any]]:
        """Get module parameters by module/module_id."""
        from omnibenchmark.benchmark.params import Params

        module_obj = (
            module if isinstance(module, Module) else self.get_modules().get(module)
        )
        if not module_obj or not module_obj.parameters:
            return None

        result: List[Any] = []
        for p in module_obj.parameters:
            result.extend(Params.expand_from_parameter(p))

        return result

    def get_module_repository(self, module: Union[str, Module]) -> Optional[Repository]:
        """Get module repository by module/module_id."""
        return self._resolve_module_attr(module, "repository")

    def get_module_environment(self, module: Union[str, Module]) -> Optional[str]:
        """Get module software environment by module/module_id."""
        return self._resolve_module_attr(module, "software_environment")

    def get_metric_collectors(self) -> List[MetricCollector]:
        """Get metric collectors."""
        return self.metric_collectors or []

    def is_initial(self, stage: Stage) -> bool:
        """Check if a stage is initial (has no inputs)."""
        return stage.inputs is None or len(stage.inputs) == 0

    def get_outputs(self) -> Dict[str, IOFile]:
        """Get all outputs from all stages."""
        outputs: Dict[str, IOFile] = {}
        for stage in self.stages:
            for output in stage.outputs:
                outputs[output.id] = output
        return outputs

    def get_modules(self) -> Dict[str, Module]:
        """Get all modules from all stages."""
        modules: Dict[str, Module] = {}
        for stage in self.stages:
            for module in stage.modules:
                modules[module.id] = module
        return modules

    def get_easyconfigs(self) -> List[Optional[str]]:
        """Get easyconfigs from software environments."""
        return [env.easyconfig for env in self.software_environments]

    def get_conda_envs(self) -> List[Optional[str]]:
        """Get conda environment files from software environments."""
        return [env.conda for env in self.software_environments]

    # Validations

    @model_validator(mode="after")
    def validate_model_structure_post_init(self) -> "Benchmark":
        """Validate pure model structure after initialization (no execution context)."""
        # Call the pure model validation from the validator base class
        self.validate_model_structure()
        return self

    def validate_execution_context(self, benchmark_dir: Path) -> None:
        """Validate execution context including file paths and environment availability.

        ARCHITECTURAL WARNING: This method performs system-specific validation
        that violates the principle of keeping models abstract and declarative.

        TODO: Move this entire method to BenchmarkExecution class or a separate
        validation layer. The model should only validate data structure and
        logical consistency, not system state or file existence.
        """
        errors: List[str] = []

        # Validate software environment paths (if benchmark_dir provided)
        if self.software_backend != SoftwareBackendEnum.host:
            for env in self.software_environments:
                env_errors = self._validate_environment_path(env, benchmark_dir)
                errors.extend(env_errors)

        # Raise error if any validation failed
        if errors:
            raise ValidationError(errors)

    def _validate_environment_path(
        self, env: SoftwareEnvironment, benchmark_dir: Path
    ) -> List[str]:
        """Validate software environment path based on backend.

        ARCHITECTURAL WARNING: This method performs system-specific validation
        that violates the principle of keeping models abstract and declarative.

        TODO: Move system-specific checks (file existence, envmodule availability)
        to BenchmarkExecution or a separate validation layer. The model should
        only validate data structure and logical consistency, not system state.
        """
        errors: List[str] = []

        # Get the appropriate environment configuration
        if self.software_backend == SoftwareBackendEnum.envmodules:
            # For envmodules, just check that the envmodule field is present
            # We don't try to load it during validation - that happens at runtime
            if not env.envmodule:
                errors.append(
                    f"Software environment with id '{env.id}' does not have a valid backend definition for: '{self.software_backend.value}'. The 'envmodule' field must be specified."
                )

        elif self.software_backend in [
            SoftwareBackendEnum.conda,
            SoftwareBackendEnum.docker,
            SoftwareBackendEnum.apptainer,
        ]:
            # Get the environment path based on backend
            env_path = None
            if self.software_backend == SoftwareBackendEnum.conda:
                env_path = env.conda
            elif self.software_backend in [
                SoftwareBackendEnum.docker,
                SoftwareBackendEnum.apptainer,
            ]:
                env_path = env.apptainer

            if not env_path:
                errors.append(
                    f"Software environment with id '{env.id}' does not have a valid backend definition for: '{self.software_backend.value}'."
                )
            elif not _is_environment_url(env_path) and not Path(env_path).is_absolute():
                # TODO: ARCHITECTURAL ISSUE - File system checks do not belong in abstract model
                # This should be moved to BenchmarkExecution or a separate validation layer
                # The model should only validate structure, not check file existence
                # Relative path - check relative to benchmark_dir
                full_path = benchmark_dir / env_path
                if not full_path.exists():
                    errors.append(
                        f"Software environment path for '{self.software_backend.value}' does not exist: '{full_path}'."
                    )
            elif not _is_environment_url(env_path) and Path(env_path).is_absolute():
                # TODO: ARCHITECTURAL ISSUE - File system checks do not belong in abstract model
                # This should be moved to BenchmarkExecution or a separate validation layer
                # The model should only validate structure, not check file existence
                # Absolute path - check directly
                if not Path(env_path).exists():
                    errors.append(
                        f"Software environment path for '{self.software_backend.value}' does not exist: '{env_path}'."
                    )

        return errors


def expand_output_path(file: IOFile) -> str:
    """
    Expands a relative output path into a standardized templated format.

    This function ensures output paths follow a consistent structure by:
    1. Prepending the standard OUTPUT_PATH_PREFIX if not already present
    2. Adding a {dataset} placeholder to the filename if not already included
    """
    # Import here to avoid circular imports
    OUTPUT_PATH_PREFIX = os.path.join("{input}", "{stage}", "{module}", "{params}")

    output_path = file.path
    if output_path.strip() == "":
        raise ValueError(f"Output path for file {file.id} is empty")

    if os.path.isabs(output_path):
        raise ValueError(
            f"Output path for file {file.id} must be relative, not absolute: {output_path}"
        )

    if not output_path.startswith(OUTPUT_PATH_PREFIX):
        output_path = os.path.join(OUTPUT_PATH_PREFIX, output_path)

    if "{dataset}" not in output_path:
        parts = output_path.rsplit(os.path.sep, 1)
        if len(parts) == 2:
            directory, filename = parts
            output_path = os.path.join(directory, f"{{dataset}}.{filename}")
        else:
            filename = parts[0]
            output_path = f"{{dataset}}.{filename}"

    return output_path


# For backwards compatibility
SoftwareEnvironmentId = str
