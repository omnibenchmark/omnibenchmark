"""
Metric collector resolution module.

This module handles the conversion of MetricCollector entities into regular
ResolvedNode instances, allowing them to be processed through the standard
pipeline instead of requiring special treatment.

Design Philosophy:
- Metric collectors are conceptually just modules that aggregate outputs
- They should be treated as regular nodes in the DAG
- The only difference is their input gathering behavior (collect ALL outputs from a stage)

FUTURE DESIGN NOTE:
- Ideally, collectors should be regular stages, not requiring special "metric_collectors" stanza
- Could use a "gather" input modifier in regular stage inputs to collect all matching outputs
- This would eliminate the need for special collector resolution logic
- See METRIC_COLLECTOR_ISSUES.md "Migration Path" section for proposed syntax

Resolution Process:
1. Parse collector inputs to identify referenced stage outputs
2. Gather ALL actual output files from resolved nodes of those stages
3. Resolve the collector's repository/module using ModuleResolver
4. Create ResolvedNode instances (one per parameter combination)
5. Handle path templates (strip {input}, replace {name} with collector id)
"""

import logging
import warnings
from typing import List, Optional

from omnibenchmark.backend.resolver import ModuleResolver
from omnibenchmark.benchmark.params import Params
from omnibenchmark.model.benchmark import MetricCollector, Benchmark
from omnibenchmark.model.resolved import ResolvedNode


logger = logging.getLogger(__name__)


def resolve_metric_collectors(
    metric_collectors: List[MetricCollector],
    resolved_nodes: List[ResolvedNode],
    benchmark: Benchmark,
    resolver: ModuleResolver,
    quiet: bool = False,
    dirty: bool = False,
    dev: bool = False,
) -> List[ResolvedNode]:
    """
    Convert metric collectors to regular ResolvedNode instances.

    This function transforms MetricCollector definitions into ResolvedNodes
    that can be processed through the standard execution pipeline.

    Args:
        metric_collectors: List of MetricCollector from benchmark definition
        resolved_nodes: Already-resolved nodes from regular stages
        benchmark: Benchmark model (for looking up stages/outputs)
        resolver: ModuleResolver for resolving collector repositories
        quiet: If True, suppress info logging
        dirty: If True, allow local paths with uncommitted changes
        dev: If True, allow unpinned branch refs on remote repos

    Returns:
        List of ResolvedNode instances representing metric collectors

    Note:
        Each collector may produce multiple nodes if it has parameter expansion.
        Collector nodes depend on ALL outputs from the referenced stage(s).
    """
    collector_nodes = []

    if not metric_collectors:
        return collector_nodes

    if not quiet:
        logger.info(f"Resolving {len(metric_collectors)} metric collector(s)...")

    for collector in metric_collectors:
        if not quiet:
            logger.info(f"  Processing collector: {collector.id}")

        try:
            # Resolve the collector's module/repository
            resolved_module = _resolve_collector_module(
                collector=collector,
                resolver=resolver,
                quiet=quiet,
                dirty=dirty,
                dev=dev,
            )

            # Expand parameters (same as regular modules)
            parameters_list = _expand_collector_parameters(collector)

            # For each parameter combination, create a resolved node
            for params in parameters_list:
                param_id = params.hash_short() if params else "default"

                # Gather inputs: ALL outputs from referenced stages
                # Returns dict mapping input names (e.g., "metrics.scores") to list of file paths
                inputs_by_name = _gather_collector_inputs(
                    collector=collector,
                    resolved_nodes=resolved_nodes,
                    benchmark=benchmark,
                )

                # Resolve output paths - need input paths to determine {input} base dir
                outputs = _resolve_collector_outputs(
                    collector=collector,
                    param_id=param_id,
                    inputs_by_name=inputs_by_name,
                )

                # Create node ID: collector_{id}_{param_id}
                node_id = f"collector_{collector.id}_{param_id}"

                # Create ResolvedNode (treating collector as a regular node)
                # For collectors, enumerate inputs with their proper names
                # Each input spec gets multiple file paths that need to be passed together
                inputs_dict = {}
                input_name_mapping = {}
                idx = 0
                for input_name, file_paths in inputs_by_name.items():
                    # Store each file under an enumerated key for Snakemake
                    # but track which ones belong to which input name
                    for file_path in file_paths:
                        key = f"input_{idx}"
                        inputs_dict[key] = file_path
                        # Map sanitized key back to original input name
                        input_name_mapping[key] = input_name
                        idx += 1

                # Get resource requirements from collector (if specified)
                node_resources = None
                if hasattr(collector, "resources") and collector.resources:
                    node_resources = collector.resources

                node = ResolvedNode(
                    id=node_id,
                    stage_id=f"_collector_{collector.id}",  # Synthetic stage ID
                    module_id=collector.id,
                    param_id=param_id,
                    module=resolved_module,
                    parameters=params,
                    inputs=inputs_dict,
                    outputs=outputs,
                    input_name_mapping=input_name_mapping,
                    benchmark_name=benchmark.get_name(),
                    benchmark_version=benchmark.get_version(),
                    benchmark_author=benchmark.get_author(),
                    resources=node_resources,
                )

                collector_nodes.append(node)

                if not quiet:
                    total_files = sum(len(files) for files in inputs_by_name.values())
                    logger.info(
                        f"    Created collector node {node_id} with {total_files} inputs"
                    )

        except Exception as e:
            import traceback

            if logger.isEnabledFor(logging.DEBUG):
                traceback.print_exc()
            raise RuntimeError(
                f"Failed to resolve collector '{collector.id}': {e}"
            ) from e

    if not quiet:
        logger.info(f"Created {len(collector_nodes)} collector node(s)")

    return collector_nodes


def _resolve_collector_module(
    collector: MetricCollector,
    resolver: ModuleResolver,
    quiet: bool = False,
    dirty: bool = False,
    dev: bool = False,
):
    """
    Resolve a metric collector's repository/module.

    Args:
        collector: MetricCollector to resolve
        resolver: ModuleResolver instance
        quiet: If True, suppress logging output

    Returns:
        ResolvedModule with concrete paths

    Note:
        This creates a synthetic Module-like object for the resolver.
        TODO: This could be cleaner with a proper adapter pattern.
    """
    # Create a minimal module-like object for the resolver
    # The resolver expects a Module with repository and software_environment attributes
    from omnibenchmark.model.benchmark import Module

    synthetic_module = Module(
        id=collector.id,
        name=collector.name or collector.id,
        repository=collector.repository,
        software_environment=collector.software_environment,
        outputs=[],  # Collectors define outputs directly, not via module
    )

    # Suppress resolver logging when quiet=True (same as regular module resolution)
    if quiet:
        # Save current log levels
        root_logger = logging.getLogger()
        omnibenchmark_logger = logging.getLogger("omnibenchmark")
        git_logger = logging.getLogger("omnibenchmark.git")
        resolver_logger = logging.getLogger("omnibenchmark.backend.resolver")

        old_level = root_logger.level
        old_handlers = root_logger.handlers[:]
        old_omni_level = omnibenchmark_logger.level
        old_omni_propagate = omnibenchmark_logger.propagate
        old_omni_handlers = omnibenchmark_logger.handlers[:]

        # Suppress all logging (root, omnibenchmark, and child loggers)
        root_logger.handlers = []
        root_logger.setLevel(logging.CRITICAL + 1)
        omnibenchmark_logger.handlers = []
        omnibenchmark_logger.setLevel(logging.CRITICAL + 1)
        omnibenchmark_logger.propagate = False
        git_logger.setLevel(logging.CRITICAL + 1)
        resolver_logger.setLevel(logging.CRITICAL + 1)

    # Capture warnings instead of displaying them when quiet
    with warnings.catch_warnings(record=True):
        if quiet:
            warnings.simplefilter("ignore")

        try:
            resolved_module = resolver.resolve(
                module=synthetic_module,
                module_id=collector.id,
                software_environment_id=collector.software_environment,
                dirty=dirty,
                dev=dev,
            )
        finally:
            # Restore logging if it was suppressed
            if quiet:
                root_logger.handlers = old_handlers
                root_logger.setLevel(old_level)
                omnibenchmark_logger.handlers = old_omni_handlers
                omnibenchmark_logger.setLevel(old_omni_level)
                omnibenchmark_logger.propagate = old_omni_propagate

    return resolved_module


def _expand_collector_parameters(collector: MetricCollector) -> List[Optional[Params]]:
    """
    Expand collector parameters into list of Params instances.

    Args:
        collector: MetricCollector with optional parameters

    Returns:
        List of Params instances (or [None] if no parameters)
    """
    parameters_list = []

    if collector.parameters:
        for param in collector.parameters:
            parameters_list.extend(Params.expand_from_parameter(param))
    else:
        parameters_list = [None]

    return parameters_list


def _gather_collector_inputs(
    collector: MetricCollector,
    resolved_nodes: List[ResolvedNode],
    benchmark: Benchmark,
) -> dict[str, List[str]]:
    """
    Gather ALL output files from stages referenced in collector inputs.

    A collector input like "metrics.scores" means:
    - Find the stage that produces output "scores"
    - Collect ALL output files from ALL resolved nodes in that stage

    Args:
        collector: MetricCollector definition
        resolved_nodes: All resolved nodes from regular stages
        benchmark: Benchmark model (for stage lookups)

    Returns:
        Dictionary mapping input names (e.g., "metrics.scores") to lists of file paths
        This preserves the grouping needed for passing as named arguments to the collector

    Complexity Analysis:
        - For each collector input (usually 1-3): O(i)
        - Find stage by output: O(s) where s = number of stages
        - Filter nodes by stage: O(n) where n = number of nodes
        - Gather outputs: O(n * o) where o = outputs per node
        Total: O(i * s * n * o) - acceptable for typical benchmark sizes
    """
    inputs_by_name = {}

    for input_ref in collector.inputs:
        # Handle both string and IOFile inputs
        input_id = input_ref if isinstance(input_ref, str) else input_ref.id

        # Find the stage that produces this output
        stage = benchmark.get_stage_by_output(input_id)

        if stage:
            # Gather ALL outputs from this stage across all resolved nodes
            # Output paths are already fully resolved (no wildcards)
            stage_outputs = []
            for node in resolved_nodes:
                if node.stage_id == stage.id:
                    for output in node.outputs:
                        stage_outputs.append(output)

            if stage_outputs:
                inputs_by_name[input_id] = stage_outputs
            else:
                logger.warning(
                    f"Collector {collector.id} references stage {stage.id} but found no outputs"
                )
        else:
            logger.warning(
                f"Collector {collector.id} references unknown output: {input_id}"
            )

    return inputs_by_name


def _resolve_collector_outputs(
    collector: MetricCollector,
    param_id: str,
    inputs_by_name: dict = None,
) -> List[str]:
    """
    Resolve collector output paths by handling template variables.

    Template variables in collector outputs:
    - {name}: Replaced with collector.id
    - {input}: Should map to the base output directory (typically "out")
                For now, we strip it since outputs are relative to Snakefile location

    Path structure:
    - Without parameters: {collector.id}/output_file.html
    - With parameters: {collector.id}/{param_id}/output_file.html

    Args:
        collector: MetricCollector definition
        param_id: Parameter hash ID (or "default")
        inputs_by_name: Optional dict of input paths (not currently used)

    Returns:
        List of resolved output paths

    Corner Cases:
    - Empty path after stripping {input}
    - Multiple consecutive slashes
    - Path that is just {input} or {name}
    """
    outputs = []

    for output_spec in collector.outputs:
        # Start with the path template
        output_path = output_spec.path

        # Replace {name} with collector id
        output_path = output_path.replace("{name}", collector.id)

        # Strip {input}/ prefix - outputs are relative to Snakefile (in out/ dir)
        # So {input}/plotting/report.html becomes plotting/report.html
        output_path = output_path.replace("{input}/", "")
        output_path = output_path.replace("{input}", "")

        # Clean up any double slashes that might result
        while "//" in output_path:
            output_path = output_path.replace("//", "/")

        # Strip leading slash if present
        output_path = output_path.lstrip("/")

        # If there are parameters, insert param_id into the path
        if param_id != "default":
            # Insert param_id before the filename
            # e.g., "plotting/report.html" -> "plotting/abc123/report.html"
            parts = output_path.rsplit("/", 1)
            if len(parts) == 2:
                output_path = f"{parts[0]}/{param_id}/{parts[1]}"
            else:
                # No directory structure, just filename
                output_path = f"{param_id}/{output_path}"

        outputs.append(output_path)

    return outputs


def format_collector_inputs_for_snakemake(inputs: List[str]) -> str:
    """
    Format collector inputs for Snakemake rule input section.

    Since collectors gather many files, we format them as a list
    that Snakemake can expand.

    Args:
        inputs: List of input file paths

    Returns:
        Formatted input string for Snakemake

    Example:
        ["file1.txt", "file2.txt"] -> '["file1.txt", "file2.txt"]'
    """
    if not inputs:
        return "[]"

    # Format as Python list literal
    formatted = "[\n"
    for inp in inputs:
        formatted += f'        "{inp}",\n'
    formatted += "    ]"

    return formatted
