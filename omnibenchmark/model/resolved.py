"""
Resolved entity classes for benchmark execution.

This module defines ResolvedNode and related classes that represent fully resolved
execution units after module resolution. These contain all concrete paths and
information needed for execution, bridging the gap between DAG building
and backend-specific workflow generation.

Design Principles:
- Abstract: No backend-specific logic (e.g., no Snakemake-specific methods)
- Concrete: All paths and references are resolved (no lazy evaluation)
- Immutable: Once created, a resolved entity should not change
- Self-contained: Contains all information needed for execution
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Any

from omnibenchmark.benchmark.params import Params
from omnibenchmark.model.benchmark import SoftwareBackendEnum


@dataclass(frozen=True)
class ResolvedEnvironment:
    """
    Resolved software environment reference.

    This is a generic container for environment references that can handle
    different backend types uniformly.

    Attributes:
        backend_type: The software backend type (conda, apptainer, envmodules, etc.)
        reference: Generic reference string that can be:
                   - For conda: relative path to yaml file (e.g., ".out/envs/myenv.yaml")
                   - For apptainer: URL (oras://..., docker://...) or path to SIF
                   - For envmodules: module name string
                   - For host: None/empty
    """

    backend_type: SoftwareBackendEnum
    reference: str

    def is_url(self) -> bool:
        """Check if reference is a URL (for apptainer images)."""
        from urllib.parse import urlparse

        try:
            result = urlparse(self.reference)
            return all([result.scheme, result.netloc])
        except ValueError:
            return False

    def is_file_path(self) -> bool:
        """Check if reference is a file path."""
        return not self.is_url() and self.reference and self.reference.strip()


@dataclass(frozen=True)
class ResolvedModule:
    """
    A fully resolved module with concrete paths.

    This represents a module after cloning and entrypoint detection.

    Attributes:
        repository_url: Git repository URL
        commit: Commit hash (may be empty/None if dirty=True for working copy)
        module_dir: Path to cloned module (can be absolute or relative)
        entrypoint: Path to entrypoint script (can be absolute or relative to module_dir)
        software_environment_id: Reference to software environment (ID string)
        resolved_environment: Resolved environment with backend type and reference
                              - For conda: reference = path to yaml file in .out/envs/
                              - For apptainer: reference = URL or path to SIF
                              - For envmodules: reference = module name string
                              - For host: None
        dirty: If True, module is from working copy without specific commit
               If False (production), commit must be specified
    """

    # Repository information
    repository_url: str
    commit: str

    # Resolved paths
    module_dir: Path
    entrypoint: Path

    # Software environment reference
    software_environment_id: str
    resolved_environment: Optional[ResolvedEnvironment] = None

    # Entrypoint execution info
    has_shebang: bool = False  # True if entrypoint has valid shebang
    interpreter: Optional[str] = (
        None  # Interpreter to use if no shebang (python3, Rscript, bash)
    )

    # Development flag
    dirty: bool = False

    def __post_init__(self):
        """Validate module configuration."""
        # In production mode (dirty=False), commit must be specified
        if not self.dirty and (not self.commit or self.commit.strip() == ""):
            raise ValueError(
                "commit must be specified when dirty=False (production mode). "
                "Set dirty=True for development with working copies."
            )


@dataclass(frozen=True)
class ResolvedNode:
    """
    A fully resolved benchmark node ready for execution.

    This class represents a node after all resolution has completed. It contains
    everything needed to generate a workflow rule, with all paths and
    references fully resolved.

    Design:
    - All paths are concrete (absolute or template strings ready for wildcards)
    - Parameters are resolved and hashed
    - Module is cloned and entrypoint is known
    - No backend-specific logic (pure data)
    """

    # ========== Identity (required fields) ==========
    id: str
    stage_id: str
    module_id: str
    param_id: str
    module: ResolvedModule

    # ========== Optional fields with defaults ==========
    # Parameters
    parameters: Optional[Params] = None
    param_dir_template: str = ""  # e.g., "{input}/stage/module/.abc123"
    param_symlink_template: str = ""  # e.g., "{input}/stage/module/method-X_k-10"

    # DAG Structure
    parent_id: Optional[str] = None
    inputs: Dict[str, str] = field(default_factory=dict)
    outputs: List[str] = field(default_factory=list)

    # Input name mapping (sanitized -> original)
    # Maps Snakemake-safe input names (data_matrix) to original names (data.matrix)
    input_name_mapping: Dict[str, str] = field(default_factory=dict)

    # Dataset resolution (for entrypoint nodes)
    dataset: Optional[str] = (
        None  # Resolved dataset value (module_id for entrypoint nodes)
    )

    # Execution Config
    timeout: Optional[int] = None

    # Benchmark metadata
    benchmark_name: str = ""
    benchmark_version: str = ""
    benchmark_author: str = ""

    def is_entrypoint(self) -> bool:
        """Check if this is an entrypoint node (no inputs)."""
        return not self.inputs or len(self.inputs) == 0

    def get_input_list(self) -> List[str]:
        """Get list of input paths (template strings)."""
        return list(self.inputs.values())

    def get_input_dict(self) -> Dict[str, str]:
        """Get input mappings as dict."""
        return self.inputs

    def get_output_list(self) -> List[str]:
        """Get list of output paths (template strings)."""
        return self.outputs

    def get_parameter_cli_args(self, style: str = "gnu") -> List[str]:
        """
        Generate CLI arguments from parameters.

        Args:
            style: Either 'gnu' (--key value) or 'equals' (--key=value)

        Returns:
            List of command line argument strings
        """
        if self.parameters is None:
            return []
        return self.parameters.to_cli_args(style=style)

    def get_parameter_json(self) -> str:
        """Get parameters as JSON string."""
        if self.parameters is None:
            return "{}"
        return self.parameters.serialize()

    def get_parameter_hash(self) -> str:
        """Get short hash of parameters."""
        if self.parameters is None:
            return "default"
        return self.parameters.hash_short()

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert to dictionary representation for serialization.

        Useful for debugging, logging, or passing to backend systems.
        """
        return {
            "id": self.id,
            "stage_id": self.stage_id,
            "module_id": self.module_id,
            "param_id": self.param_id,
            "module": {
                "repository_url": self.module.repository_url,
                "commit": self.module.commit,
                "module_dir": str(self.module.module_dir),
                "entrypoint": str(self.module.entrypoint),
                "software_environment_id": self.module.software_environment_id,
                "resolved_environment": {
                    "backend_type": self.module.resolved_environment.backend_type.value,
                    "reference": self.module.resolved_environment.reference,
                }
                if self.module.resolved_environment
                else None,
            },
            "parameters": self.get_parameter_json(),
            "param_dir_template": self.param_dir_template,
            "param_symlink_template": self.param_symlink_template,
            "parent_id": self.parent_id,
            "inputs": self.inputs,
            "outputs": self.outputs,
            "input_name_mapping": self.input_name_mapping,
            "dataset": self.dataset,
            "timeout": self.timeout,
            "benchmark_name": self.benchmark_name,
            "benchmark_version": self.benchmark_version,
            "benchmark_author": self.benchmark_author,
        }

    @staticmethod
    def create_id(
        stage_id: str,
        module_id: str,
        param_id: str,
        after_stage_id: Optional[str] = None,
    ) -> str:
        """
        Create a node ID from components.

        Format: stage-module-param_id[-after_stage]
        """
        node_id = f"{stage_id}-{module_id}-{param_id}"
        if after_stage_id:
            node_id += f"-after_{after_stage_id}"
        return node_id

    def __str__(self) -> str:
        """String representation."""
        return self.id

    def __repr__(self) -> str:
        """Detailed representation."""
        return f"ResolvedNode({self.id})"


@dataclass(frozen=True)
class ResolvedMetricCollector:
    """
    A fully resolved metric collector node.

    Similar to ResolvedNode but for metric collection tasks.
    """

    # ========== Required fields ==========
    id: str
    module: ResolvedModule

    # ========== Optional fields with defaults ==========
    name: Optional[str] = None
    parameters: Optional[Params] = None
    param_id: str = "default"
    input_patterns: List[str] = field(default_factory=list)
    outputs: List[str] = field(default_factory=list)
    timeout: Optional[int] = None
    benchmark_name: str = ""
    benchmark_version: str = ""
    benchmark_author: str = ""

    def get_input_list(self) -> List[str]:
        """Get list of input patterns."""
        return self.input_patterns

    def get_output_list(self) -> List[str]:
        """Get list of output paths."""
        return self.outputs

    def get_parameter_cli_args(self, style: str = "gnu") -> List[str]:
        """
        Generate CLI arguments from parameters.

        Args:
            style: Either 'gnu' (--key value) or 'equals' (--key=value)

        Returns:
            List of command line argument strings
        """
        if self.parameters is None:
            return []
        return self.parameters.to_cli_args(style=style)

    def get_parameter_json(self) -> str:
        """Get parameters as JSON string."""
        if self.parameters is None:
            return "{}"
        return self.parameters.serialize()

    def get_parameter_hash(self) -> str:
        """Get short hash of parameters."""
        if self.parameters is None:
            return "default"
        return self.parameters.hash_short()

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        return {
            "id": self.id,
            "name": self.name,
            "module": {
                "repository_url": self.module.repository_url,
                "commit": self.module.commit,
                "module_dir": str(self.module.module_dir),
                "entrypoint": str(self.module.entrypoint),
                "software_environment_id": self.module.software_environment_id,
                "resolved_environment": {
                    "backend_type": self.module.resolved_environment.backend_type.value,
                    "reference": self.module.resolved_environment.reference,
                }
                if self.module.resolved_environment
                else None,
            },
            "input_patterns": self.input_patterns,
            "outputs": self.outputs,
            "timeout": self.timeout,
            "benchmark_name": self.benchmark_name,
            "benchmark_version": self.benchmark_version,
            "benchmark_author": self.benchmark_author,
        }

    def __str__(self) -> str:
        """String representation."""
        return self.id

    def __repr__(self) -> str:
        """Detailed representation."""
        return f"ResolvedMetricCollector({self.id})"
