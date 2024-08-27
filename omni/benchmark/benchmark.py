from omni.benchmark import dag
from omni.benchmark.converter import LinkMLConverter
from omni.benchmark.validation import Validator
from omni.utils import *


class Benchmark:
    def __init__(self, benchmark_yaml: Path, out_dir: str = "out"):
        converter = LinkMLConverter(benchmark_yaml)
        validator = Validator()
        converter = validator.validate(converter)

        self.converter = converter
        self.out_dir = out_dir
        self.G = dag.build_benchmark_dag(converter, self.out_dir)

        self.execution_paths = None

    def get_converter(self):
        return self.converter

    def get_benchmark_name(self):
        return self.converter.get_name()

    def get_benchmark_version(self):
        return self.converter.get_version()

    def get_benchmark_author(self):
        return self.converter.get_author()

    def get_definition(self):
        return self.converter.get_definition()

    def get_definition_file(self) -> Path:
        return self.converter.benchmark_file

    def get_nodes(self):
        return list(self.G.nodes)

    def get_stage_ids(self):
        return self.converter.get_stages().keys()

    def get_node_by_id(self, node_id):
        for node in self.G.nodes:
            if node.get_id() == node_id:
                return node

        return None

    def get_nodes_by_module_id(self, module_id: str) -> List:
        nodes = []
        for node in self.G.nodes:
            if node.module_id == module_id:
                nodes.append(node)

        return nodes

    def get_benchmark_datasets(self) -> List[str]:
        datasets = []
        for _, stage in self.converter.get_stages().items():
            # There should be only one initial stage
            if self.converter.is_initial(stage):
                for module in stage.modules:
                    datasets.append(module.id)

                break

        return datasets

    def get_execution_paths(self):
        if self.execution_paths is None:
            self.execution_paths = self._generate_execution_paths()

        return self.execution_paths

    def get_output_paths(self):
        execution_paths = self.get_execution_paths()

        output_paths = [
            format_name(output, self.out_dir)
            for path in execution_paths
            for output in self._construct_output_paths(prefix=self.out_dir, nodes=path)
        ]

        return set(output_paths)

    def get_explicit_inputs(self, stage_id: str):
        stage = self.converter.get_stage(stage_id)
        implicit_inputs = self.converter.get_stage_implicit_inputs(stage)
        explicit_inputs = [
            self.converter.get_explicit_inputs(i) for i in implicit_inputs
        ]
        return explicit_inputs

    def get_explicit_outputs(self, stage_id: str):
        stage = self.converter.get_stage(stage_id)
        return self.converter.get_stage_outputs(stage)

    def get_available_parameter(self, module_id: str):
        node = next(node for node in self.G.nodes if node.module_id == module_id)
        return node.get_parameters()

    def plot_benchmark_graph(self):
        dag.plot_graph(
            self.G, output_file="output_dag.png", scale_factor=1.5, node_spacing=0.2
        )

    def __str__(self):
        return f"Benchmark({self.get_definition})"

    def _generate_execution_paths(self):
        path_exclusions = self._get_path_exclusions()
        initial_nodes, terminal_nodes = dag.find_initial_and_terminal_nodes(self.G)

        execution_paths = []
        for initial_node in initial_nodes:
            for terminal_node in terminal_nodes:
                paths = dag.list_all_paths(self.G, initial_node, terminal_node)
                paths_after_exclusion = dag.exclude_paths(paths, path_exclusions)

                execution_paths.extend(paths_after_exclusion)

        return execution_paths

    def _get_path_exclusions(self):
        path_exclusions = {}
        stages = self.converter.get_stages()
        for stage_id in stages:
            stage = stages[stage_id]

            modules_in_stage = self.converter.get_modules_by_stage(stage)
            for module_id in modules_in_stage:
                module_excludes = self.converter.get_module_excludes(module_id)
                if module_excludes:
                    path_exclusions[module_id] = module_excludes

        return path_exclusions

    def _construct_output_paths(self, prefix, nodes):
        if nodes is None or len(nodes) == 0:
            return []
        else:
            head = nodes[0]
            tail = nodes[1:]
            stage_outputs = self.converter.get_stage_outputs(head.stage_id).values()

            current_path = f"{head.stage_id}/{head.module_id}"
            if any(["{params}" in o for o in stage_outputs]):
                current_path += f"/{head.param_id}"

            new_prefix = f"{prefix}/{current_path}"
            paths = [
                x.format(
                    input=prefix,
                    stage=head.stage_id,
                    module=head.module_id,
                    params=head.param_id,
                    dataset="{dataset}",
                )
                for x in stage_outputs
            ]

            return paths + self._construct_output_paths(new_prefix, tail)
