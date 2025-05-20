from pathlib import Path
from typing import Union, Optional, Dict, List

from omni_schema.datamodel import omni_schema

from omnibenchmark.benchmark import params
from omnibenchmark.utils import merge_dict_list, parse_instance


class LinkMLConverter:
    def __init__(self, benchmark_file: Path):
        """
        The constructor tries to parse the benchmark file and initialize the model.
        It also validates the model against the schema, and raises a ValidationError if
        """
        self.benchmark_file = benchmark_file
        self.model: omni_schema.Benchmark = parse_instance(
            benchmark_file, omni_schema.Benchmark
        )

    def get_name(self) -> str:
        """Get name of the benchmark"""

        return self.model.name if self.model.name else self.model.id

    def get_version(self) -> str:
        """Get version of the benchmark"""

        return self.model.version

    def get_author(self) -> str:
        """Get author of the benchmark"""

        return self.model.benchmarker

    def get_software_backend(self) -> omni_schema.SoftwareBackendEnum:
        """Get software backend of the benchmark"""

        return self.model.software_backend.code

    def get_software_environments(
        self,
    ) -> Dict[omni_schema.SoftwareEnvironmentId, omni_schema.SoftwareEnvironment]:
        """Get software environments"""

        environments = {}
        for env in self.model.software_environments:
            environments[env.id] = env

        return environments

    def get_definition(self) -> omni_schema.Benchmark:
        """Get underlying benchmark"""

        return self.model

    def get_easyconfigs(self):
        """Get easyconfigs"""
        easyconfigs = []
        for se in self.model.software_environments:
            easyconfigs.append(se["easyconfig"])
        return easyconfigs

    def get_conda_envs(self):
        """Get conda envs"""
        envs = []
        for se in self.model.software_environments:
            envs.append(se["conda"])
        return envs

    def get_stages(self) -> Dict[str, omni_schema.Stage]:
        """Get benchmark stages"""

        return dict([(x.id, x) for x in self.model.stages])

    def get_stage(self, stage_id: str) -> Optional[omni_schema.Stage]:
        """Get stage by stage_id"""

        return self.get_stages()[stage_id]

    def get_stage_by_output(self, output_id: str) -> Optional[omni_schema.Stage]:
        """Get stage that returns output with output_id"""

        stage_by_output: dict = {}
        for stage_id, stage in self.get_stages().items():
            stage_by_output.update({output.id: stage for output in stage.outputs})

        return stage_by_output.get(output_id)

    def get_modules_by_stage(
        self, stage: Union[str, omni_schema.Stage]
    ) -> Dict[str, omni_schema.Module]:
        """Get modules by stage/stage_id"""

        if isinstance(stage, str):
            stage = self.get_stages()[stage]

        return dict([(x.id, x) for x in stage.modules])

    def get_stage_implicit_inputs(
        self, stage: Union[str, omni_schema.Stage]
    ) -> List[str]:
        """Get implicit inputs of a stage by stage/stage_id"""

        if isinstance(stage, str):
            stage = self.get_stages()[stage]

        return [input.entries for input in stage.inputs]

    def get_explicit_inputs(self, input_ids: List[str]) -> Dict[str, str]:
        """Get explicit inputs of a stage by input_id(s)"""
        all_stages_outputs = []
        for stage_id in self.get_stages():
            outputs = self.get_stage_outputs(stage=stage_id)
            outputs = {
                key: value.format(
                    input="{input}",
                    stage=stage_id,
                    module="{module}",
                    params="{params}",
                    dataset="{dataset}",
                )
                for key, value in outputs.items()
            }
            all_stages_outputs.append(outputs)

        all_stages_outputs = merge_dict_list(all_stages_outputs)

        explicit = {key: None for key in input_ids}
        for in_deliverable in input_ids:
            # beware stage needs to be substituted
            curr_output = all_stages_outputs[in_deliverable]

            explicit[in_deliverable] = curr_output

        return explicit

    def get_stage_outputs(self, stage: Union[str, omni_schema.Stage]) -> Dict[str, str]:
        """Get outputs of a stage by stage/stage_id"""

        if isinstance(stage, str):
            stage = self.get_stages()[stage]

        return dict([(output.id, output.path) for output in stage.outputs])

    def get_output_stage(self, output_id: str) -> omni_schema.Stage:
        """Get stage that returns output with out_id"""

        stage_by_output: dict = {}
        for stage in self.model.stages:
            stage_by_output.update({out.id: stage for out in stage.outputs})

        return stage_by_output.get(output_id)

    def get_module_excludes(self, module: Union[str, omni_schema.Module]) -> List[str]:
        """Get module excludes by module/module_id"""

        if isinstance(module, str):
            module = self.get_modules()[module]

        return module.exclude

    def get_module_parameters(
        self, module: Union[str, omni_schema.Module]
    ) -> List[params.Params]:
        """Get module parameters by module/module_id"""

        if isinstance(module, str):
            module = self.get_modules()[module]

        parameters = None
        if module.parameters is not None:
            parameters = [
                params.Params.from_cli_args(x.values) for x in module.parameters
            ]

        return parameters

    def get_module_repository(
        self, module: Union[str, omni_schema.Module]
    ) -> omni_schema.Repository:
        """Get module repository by module/module_id"""

        if isinstance(module, str):
            module = self.get_modules()[module]

        return module.repository

    def get_module_environment(
        self, module: Union[str, omni_schema.Module]
    ) -> omni_schema.SoftwareEnvironmentId:
        """Get module software environment by module/module_id"""

        if isinstance(module, str):
            module = self.get_modules()[module]

        return module.software_environment

    def get_metric_collectors(self) -> List[omni_schema.MetricCollector]:
        return self.model.metric_collectors

    def is_initial(self, stage: omni_schema.Stage) -> bool:
        """Check if stage is initial"""

        if stage.inputs is None or len(stage.inputs) == 0:
            return True
        else:
            return False

    def get_outputs(self) -> Dict[str, str]:
        """Get outputs"""

        outputs = {}
        for stage_id, stage in self.get_stages().items():
            for output in stage.outputs:
                outputs[output.id] = output

        return outputs

    def get_modules(self) -> Dict[str, omni_schema.Module]:
        """Get modules"""

        modules = {}

        for stage_id, stage in self.get_stages().items():
            modules_in_stage = self.get_modules_by_stage(stage)
            modules.update(modules_in_stage)

        return modules
