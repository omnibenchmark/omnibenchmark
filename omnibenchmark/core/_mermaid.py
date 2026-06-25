"""Mermaid diagram generation utilities for omnibenchmark."""

from typing import Dict, List
from omnibenchmark.model import Benchmark as BenchmarkModel
from omnibenchmark.model.benchmark import Parameter
from omnibenchmark.core._graph import upstream_stage_ids
from omnibenchmark.model.params import Params


def _quote_label(text: str) -> str:
    """Wrap a node label so Mermaid treats it as literal text.

    Mermaid only tolerates quotes, commas, brackets and other punctuation
    inside a node label when the label is wrapped in double quotes; embedded
    double quotes themselves must be entity-escaped. Without this, labels such
    as ``--alpha 0.1`` or ``'--beta', '2'`` break the parser.
    """
    return '"' + text.replace('"', "&quot;") + '"'


class MermaidDiagramGenerator:
    """Generates Mermaid diagrams for benchmark workflows."""

    def __init__(self, model: BenchmarkModel, graph=None):
        # ``graph`` (the execution DAG) is accepted for backward compatibility
        # but no longer used to derive edges. The execution DAG connects each
        # node only to its single "latest" upstream stage (so that output-path
        # construction stays a linear ancestry); that drops edges whenever an
        # input group spans multiple upstream stages. The topology view instead
        # reads the model's declared inputs directly, so every declared
        # stage-to-stage dependency is rendered.
        self.model = model
        self.graph = graph

    def generate_diagram(
        self, show_params: bool = True, compact_params: bool = False
    ) -> str:
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
            diagram_parts.append(
                self._generate_parameter_subgraphs(compact_params=compact_params)
            )

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

            # All modules in a stage share the stage's declared inputs, so the
            # set of upstream modules is the same for every module in the stage.
            upstream_module_ids = self._get_upstream_module_ids(stage)

            # Add modules to the stage subgraph
            for module_id, module in self.model.get_modules_by_stage(stage).items():
                subgraph_lines.append(f"\t\t{module_id}")

                # Add connections from every upstream module the stage depends on
                for from_module_id in upstream_module_ids:
                    subgraph_lines.append(f"\t\t{from_module_id} --> {module_id}")

            subgraph_lines.append("\tend")
            subgraphs.append("\n".join(subgraph_lines))

        return "\n".join(subgraphs)

    def _get_upstream_module_ids(self, stage) -> List[str]:
        """Module IDs of every stage feeding *stage* via its declared inputs.

        Collects the modules of every stage feeding *stage* (see
        :func:`upstream_stage_ids`). Unlike the execution DAG, a single input
        group spanning multiple upstream stages contributes edges from *all* of
        them, not just the topologically-latest one.
        """
        module_ids = set()
        for stage_id in upstream_stage_ids(self.model, stage):
            module_ids.update(self.model.get_modules_by_stage(stage_id).keys())

        return sorted(module_ids)

    def _generate_parameter_subgraphs(self, compact_params: bool) -> str:
        """Generate parameter subgraphs and their connections to modules."""
        param_subgraphs = []
        module_id_params_dict = self._collect_module_parameters()

        for module_id, parameters in module_id_params_dict.items():
            subgraph_lines = [f"\tsubgraph params_{module_id}"]

            if compact_params:
                # One node per Parameter before expansion (dict-format only)
                for paramset, parameter in enumerate(parameters):
                    if parameter.params is None:
                        continue
                    comp_str = " ".join(f"{k}:{v}" for k, v in parameter.params.items())
                    paramset_node = f"{module_id}_paramset{paramset}"
                    subgraph_lines.append(
                        f"\t\t{paramset_node}[{_quote_label(comp_str)}]"
                    )
            else:
                # Expand each Parameter and render one node per combination
                for parameter in parameters:
                    for param in Params.expand_from_parameter(parameter):
                        label = " ".join(param.to_cli_args())
                        param_node = f"{param.hash()}[{_quote_label(label)}]"
                        subgraph_lines.append(f"\t\t{param_node}")

            subgraph_lines.extend(
                ["\tend", f"\tparams_{module_id}:::param --o {module_id}"]
            )
            param_subgraphs.append("\n".join(subgraph_lines))

        return "\n".join(param_subgraphs)

    def _collect_module_parameters(self) -> Dict[str, List[Parameter]]:
        """Collect raw Parameter objects for all modules that have them."""
        module_id_params_dict = {}

        for stage_id, stage in self.model.get_stages().items():
            for module_id, module in self.model.get_modules_by_stage(stage).items():
                if module.parameters:
                    module_id_params_dict[module_id] = module.parameters

        return module_id_params_dict


def generate_mermaid_diagram(
    model: BenchmarkModel,
    graph=None,
    show_params: bool = True,
    compact_params: bool = False,
) -> str:
    """Generate a Mermaid diagram for a benchmark workflow.

    Args:
        model: The benchmark model containing workflow definition
        graph: Deprecated/unused — the execution DAG. Topology edges are now
            derived from the model's declared inputs.
        show_params: Whether to include parameter subgraphs

    Returns:
        A string containing the complete Mermaid diagram syntax
    """
    generator = MermaidDiagramGenerator(model, graph)
    return generator.generate_diagram(show_params, compact_params)
