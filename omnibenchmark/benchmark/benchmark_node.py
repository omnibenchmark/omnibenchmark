import os.path
from pathlib import Path


class BenchmarkNode:
    def __init__(
        self,
        converter,
        stage,
        module,
        parameters,
        inputs,
        outputs,
        after=None,
    ):
        self.converter = converter
        self.stage = stage
        self.module = module
        self.parameters = parameters
        self.inputs = inputs
        self.outputs = outputs
        self.after = after

        self.stage_id = stage.id
        self.module_id = module.id
        self.param_id = (
            "default" if not self.parameters else f".{self.parameters.hash()}"
        )

    def is_entrypoint(self):
        return not self.inputs or len(self.inputs) == 0

    def get_id(self):
        return BenchmarkNode.to_id(
            self.stage_id, self.module_id, self.param_id, self.after
        )

    def get_benchmark_name(self):
        return self.converter.get_name()

    def get_benchmark_version(self):
        return self.converter.get_version()

    def get_benchmark_author(self):
        return self.converter.get_author()

    def get_benchmark_software_backend(self):
        return self.converter.get_software_backend()

    def get_benchmark_software_environments(self):
        return self.converter.get_software_environments()

    def get_definition(self):
        return self.converter.get_definition()

    def get_definition_file(self) -> Path:
        return self.converter.benchmark_file

    def get_inputs(self):
        return self.inputs.values() if self.inputs else []

    def get_inputs_dict(self):
        return self.inputs if self.inputs else {}

    def get_explicit_inputs(self):
        explicit_inputs = [
            self.converter.get_explicit_inputs(i)
            for i in self.converter.get_stage_implicit_inputs(self.stage)
        ]
        return explicit_inputs

    def get_input_paths(self, config, return_as_dict=False):
        input_paths = {}
        inputs = self.inputs.items() if self.inputs else {}
        for key, input in inputs:
            input = os.path.basename(input)
            input = input.replace("{dataset}", config["dataset"])
            input = os.path.join(config["input"], input)
            input_paths[key] = input

        if return_as_dict:
            return input_paths
        else:
            return list(input_paths.values())

    def get_outputs(self):
        return self.outputs if self.outputs else []

    def get_output_paths(self, config):
        output_paths = []

        pre = config.get("input", config.get("output"))
        dataset = config["dataset"]
        for output in self.get_outputs():
            output = output.format(
                pre=pre,
                dataset=dataset,
                stage=self.stage_id,
                module=self.module_id,
                params=self.param_id,
            )

            output_paths.append(output)

        return output_paths

    def get_benchmark_path(self, config):
        pre = config.get("input", config.get("output"))
        dataset = config["dataset"]

        output_paths = self.get_outputs()
        output_dir = os.path.commonpath(output_paths)
        if len(output_paths) == 1:
            output_dir = Path(os.path.dirname(output_dir))

        benchmark_path = Path(output_dir) / "{dataset}_performance.txt"

        benchmark_path = str(benchmark_path).format(
            pre=pre,
            dataset=dataset,
            stage=self.stage_id,
            module=self.module_id,
            params=self.param_id,
        )

        return benchmark_path

    def get_parameters(self):
        return self.parameters

    def get_repository(self):
        return self.converter.get_module_repository(module=self.module)

    def get_software_environment(self):
        return self.converter.get_module_environment(module=self.module)

    def is_initial(self):
        return self.converter.is_initial(self.stage)

    def get_stage(self):
        return self.stage

    def display_name(self):
        return BenchmarkNode.to_id(
            self.stage_id, self.module_id, self.param_id, self.after
        )

    def __str__(self):
        return self.display_name()

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        if isinstance(other, BenchmarkNode):
            return (self.stage_id, self.module_id, self.parameters, self.inputs) == (
                other.stage_id,
                other.module_id,
                other.parameters,
                other.inputs,
            )
        return False

    def __hash__(self):
        return hash(self.get_id())

    @staticmethod
    def to_id(stage_id, module_id, param_id, after_stage_id=None):
        node_id = f"{stage_id}-{module_id}-{param_id}"
        node_id += f"-after_{after_stage_id}" if after_stage_id else ""

        return node_id
