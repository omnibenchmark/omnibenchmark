import os.path

from pathlib import Path
from typing import Dict, List, Set, Optional

from omnibenchmark.benchmark._mermaid import generate_mermaid_diagram
from omnibenchmark.benchmark._paths import (
    collect_output_paths,
    collect_path_exclusions,
)

from omnibenchmark.benchmark import _graph as graph
from omnibenchmark.model import Benchmark as BenchmarkModel, SoftwareBackendEnum
from omnibenchmark.utils import format_mc_output

from ._dag_builder import DAGBuilder
from ._dot import export_to_dot


class ExecutionContext:
    """Encapsulates execution-specific context like paths and directories.

    Created during LinkML → Pydantic migration to handle file system concerns
    that the pure Benchmark model should not know about. This allows the model
    to remain path-agnostic while providing execution operations with the
    directory context they need for relative path resolution.
    """

    def __init__(self, benchmark_yaml: Path, out_dir: Path = Path("out")):
        # base path is always the location of the benchmark YAML file
        # TODO: rename to definition_file
        self.path = benchmark_yaml

        # directory is used to lookup files and directories expressed in relative paths,
        # like the environment definition paths
        self.directory = benchmark_yaml.parent.absolute()

        # where are we going to write the output files, in the case of local execution
        self.out_dir = out_dir
        # Note: We don't create out_dir here anymore - it will be created lazily when needed

    def ensure_out_dir(self):
        """Create the output directory if it doesn't exist. Call this before writing outputs."""
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir, exist_ok=True)


class BenchmarkExecution:
    """
    BenchmarkExecution contains the execution context and the benchmark model.

    This class was created during the LinkML → Pydantic migration to separate
    pure data models from execution-specific operations that require file system
    paths and directory context. It serves as an adapter between the pure
    Benchmark model and the execution environment.

    Architecture Notes:
    - The Benchmark model (in omnibenchmark.model) is now path-agnostic and purely declarative
    - BenchmarkExecution adds execution context (paths, directories, DAG building)
    - This separation allows the model to be used independently for validation,
      serialization, and other operations that don't require file system access

    Architecture Note on DAG Usage:
    This class builds a computational DAG for pre-execution analysis (path planning,
    validation, visualization). While Snakemake builds its own DAG from the generated
    rules, our DAG serves different purposes:
    - Early validation of benchmark structure
    - Execution path enumeration for output collection
    - Visualization exports

    Future consideration: Evaluate whether Snakemake's DAG APIs could replace
    some of this functionality to reduce duplication.

    Future Refactoring Considerations:
    - This class currently has mixed responsibilities (model access + execution coordination)
    - Consider splitting into separate ExecutionContext and BenchmarkOrchestrator classes
    - The numerous getter methods suggest this may be serving as an anti-corruption layer
    """

    def __init__(self, benchmark_yaml: Path, out_dir: Path = Path("out")):
        self.context = ExecutionContext(benchmark_yaml, out_dir)

        self.model = BenchmarkModel.from_yaml(benchmark_yaml)

        # Validate execution context separately from pure model validation
        # Pure model validation happens automatically via Pydantic
        self.model.validate_execution_context(self.context.directory)

        # Use DAGBuilder to construct the computational graph
        # TODO: why does the DAG needs out_dir?
        dag_builder = DAGBuilder(self.model, self.context.out_dir)
        self.G = dag_builder.build()

        self.execution_paths = None

    def get_storage_api(self) -> Optional[str]:
        """Get storage API with backward compatibility."""
        return self.model.get_storage_api()

    def get_storage_bucket_name(self) -> Optional[str]:
        """Get storage bucket name with backward compatibility."""
        return self.model.get_storage_bucket_name()

    def get_storage_endpoint(self) -> Optional[str]:
        """Get storage endpoint."""
        return self.model.get_storage_endpoint()

    def get_model(self):
        """Get the underlying Pydantic model.

        Note: Direct model access - consider whether callers should use this
        class's methods instead to maintain proper abstraction boundaries.
        """
        return self.model

    def get_benchmark_name(self):
        return self.model.get_name()

    def get_benchmark_version(self):
        return self.model.get_version()

    def get_benchmark_author(self):
        return self.model.get_author()

    def get_benchmark_software_backend(self):
        return self.model.get_software_backend()

    def get_benchmark_software_environments(self):
        return self.model.get_software_environments()

    def get_definition(self):
        """Legacy method - returns the model for backward compatibility."""
        return self.model

    def get_definition_file(self) -> Path:
        """Get the path to the benchmark definition file.

        This bridges the gap between the path-agnostic model and execution context.
        """
        return self.context.path

    def get_easyconfigs(self):
        return self.model.get_easyconfigs()

    def get_conda_envs(self):
        return self.model.get_conda_envs()

    def get_stages(self):
        return self.model.get_stages()

    def get_nodes(self):
        return list(self.G.nodes)

    def get_node_by_id(self, node_id):
        return graph.find_node_by_id(self.G, node_id)

    def get_nodes_by_module_id(self, module_id: str) -> List:
        return graph.get_nodes_by_module_id(self.G, module_id)

    def get_nodes_by_stage_id(self, stage_id: str) -> List:
        return graph.get_nodes_by_stage_id(self.G, stage_id)

    def get_benchmark_datasets(self) -> List[str]:
        return graph.get_benchmark_datasets(self.model, self.get_stages())

    def get_execution_paths(self):
        if self.execution_paths is None:
            self.execution_paths = self._generate_execution_paths()

        return self.execution_paths

    def get_output_paths(self) -> Set[str]:
        execution_paths = self.get_execution_paths()
        return collect_output_paths(execution_paths, self.context.out_dir, self.model)

    def get_metric_collector_output_paths(self):
        collectors = self.get_metric_collectors()
        output_paths = []

        for collector in collectors:
            for output in collector.outputs:
                output_paths.append(
                    format_mc_output(output, self.context.out_dir, collector.id)
                )

        return set(output_paths)

    def get_explicit_inputs(self, stage_id: str):
        stage = self.model.get_stage(stage_id)
        if stage is None:
            return []
        implicit_inputs = self.model.get_stage_implicit_inputs(stage)
        explicit_inputs = [self.model.get_explicit_inputs(i) for i in implicit_inputs]
        return explicit_inputs

    def get_explicit_input(self, input_ids: List[str]) -> Dict[str, str]:
        return self.model.get_explicit_inputs(input_ids)

    def get_explicit_outputs(self, stage_id: str):
        stage = self.model.get_stage(stage_id)
        if stage is None:
            return {}
        return self.model.get_stage_outputs(stage)

    def get_available_parameter(self, module_id: str):
        node = graph.find_node_with_module_id(self.G, module_id)
        return node.get_parameters() if node else None

    def get_metric_collectors(self):
        return self.model.get_metric_collectors()

    def get_stage_ids(self):
        """Get all stage IDs."""
        return [stage.id for stage in self.model.stages]

    def _generate_execution_paths(self):
        path_exclusions = collect_path_exclusions(self.model)
        initial_nodes, terminal_nodes = graph.find_initial_and_terminal_nodes(self.G)

        execution_paths = []
        for initial_node in initial_nodes:
            for terminal_node in terminal_nodes:
                paths = graph.list_all_paths(self.G, initial_node, terminal_node)
                paths_after_exclusion = graph.exclude_paths(paths, path_exclusions)

                execution_paths.extend(paths_after_exclusion)

        return execution_paths

    # Visualization functions

    def export_to_dot(self):
        return export_to_dot(self.G, title=self.get_benchmark_name())

    def get_environment_path(
        self, env_key: str, software_backend: SoftwareBackendEnum
    ) -> Optional[str]:
        from omnibenchmark.benchmark import Validator

        benchmark_dir = self.context.directory
        environment = self.get_benchmark_software_environments().get(env_key, None)
        environment_path = (
            Validator.get_environment_path(software_backend, environment, benchmark_dir)
            if environment
            else None
        )

        return environment_path

    def export_to_mermaid(self, show_params: bool = True) -> str:
        """Export the benchmark workflow as a Mermaid diagram.

        Args:
            show_params: Whether to include parameter subgraphs in the diagram

        Returns:
            A string containing the complete Mermaid diagram syntax
        """
        return generate_mermaid_diagram(self.model, self.G, show_params)

    def __str__(self):
        return f"Benchmark({self.get_definition})"
