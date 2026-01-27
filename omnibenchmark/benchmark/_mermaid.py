"""Mermaid diagram generation utilities for omnibenchmark."""

from typing import Dict, Any, List
from omnibenchmark.model import Benchmark as BenchmarkModel
from omnibenchmark.benchmark._graph import get_nodes_by_module_id


class MermaidDiagramGenerator:
    """Generates Mermaid diagrams for benchmark workflows."""

    def __init__(self, model: BenchmarkModel, graph):
        self.model = model
        self.graph = graph

    def generate_diagram(self, show_params: bool = True) -> str:
        """Generate a complete Mermaid flowchart diagram.

        Args:
            show_params: Whether to include parameter subgraphs in the diagram

        Returns:
            A string containing the complete Mermaid diagram syntax
        """
        diagram_parts = [
            self._generate_header(),
            self._generate_flowchart_declaration(),
            self._generate_css_classes(),
            self._generate_stage_subgraphs(),
        ]

        if show_params:
            diagram_parts.append(self._generate_parameter_subgraphs())

        return "\n".join(diagram_parts)

    def _generate_header(self) -> str:
        """Generate the diagram header with title."""
        return "\n".join(["---", f"title: {self.model.get_name()}", "---"])

    def _generate_flowchart_declaration(self) -> str:
        """Generate the flowchart declaration."""
        return "flowchart LR"

    def _generate_css_classes(self) -> str:
        """Generate CSS class definitions for styling."""
        return "\tclassDef param fill:#f96"

    def _generate_stage_subgraphs(self) -> str:
        """Generate subgraphs for each stage and their modules."""
        subgraphs = []

        for stage_id, stage in self.model.get_stages().items():
            subgraph_lines = [f"\tsubgraph {stage_id}"]

            # Add modules to the stage subgraph
            for module_id, module in self.model.get_modules_by_stage(stage).items():
                subgraph_lines.append(f"\t\t{module_id}")

                # Add connections from predecessor modules
                module_connections = self._get_module_connections(module_id)
                subgraph_lines.extend(module_connections)

            subgraph_lines.append("\tend")
            subgraphs.append("\n".join(subgraph_lines))

        return "\n".join(subgraphs)

    def _get_module_connections(self, module_id: str) -> List[str]:
        """Get Mermaid connection strings for a module's dependencies."""
        connections = []
        module_nodes = self._get_nodes_by_module_id(module_id)

        from_nodes = [list(self.graph.predecessors[node]) for node in module_nodes]

        from_node_module_ids = sorted(
            list(set([node.module_id for flatten in from_nodes for node in flatten]))
        )

        for from_module_id in from_node_module_ids:
            connections.append(f"\t\t{from_module_id} --> {module_id}")

        return connections

    def _generate_parameter_subgraphs(self) -> str:
        """Generate parameter subgraphs and their connections to modules."""
        param_subgraphs = []
        module_id_params_dict = self._collect_module_parameters()

        for module_id, params in module_id_params_dict.items():
            # Create parameter subgraph
            subgraph_lines = [f"\tsubgraph params_{module_id}"]

            for param in params:
                param_node = f"{param.hash()}{param.to_cli_args()}"
                subgraph_lines.append(f"\t\t{param_node}")

            subgraph_lines.extend(
                ["\tend", f"\tparams_{module_id}:::param --o {module_id}"]
            )

            param_subgraphs.append("\n".join(subgraph_lines))

        return "\n".join(param_subgraphs)

    def _collect_module_parameters(self) -> Dict[str, Any]:
        """Collect parameters for all modules that have them."""
        module_id_params_dict = {}

        for stage_id, stage in self.model.get_stages().items():
            for module_id, module in self.model.get_modules_by_stage(stage).items():
                module_parameters = self.model.get_module_parameters(module)
                if module_parameters:
                    module_id_params_dict[module_id] = module_parameters

        return module_id_params_dict

    def _get_nodes_by_module_id(self, module_id: str) -> List:
        """Get all graph nodes for a given module ID."""
        return get_nodes_by_module_id(self.graph, module_id)


def generate_mermaid_diagram(
    model: BenchmarkModel, graph, show_params: bool = True
) -> str:
    """Generate a Mermaid diagram for a benchmark workflow.

    Args:
        model: The benchmark model containing workflow definition
        graph: The benchmark DAG graph
        show_params: Whether to include parameter subgraphs

    Returns:
        A string containing the complete Mermaid diagram syntax
    """
    generator = MermaidDiagramGenerator(model, graph)
    return generator.generate_diagram(show_params)
