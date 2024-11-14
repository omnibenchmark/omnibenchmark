import logging
import os.path
from collections import Counter
from pathlib import Path
from typing import Union, Optional
from urllib.parse import urlparse

from omni_schema.datamodel.omni_schema import SoftwareBackendEnum, SoftwareEnvironment

from omni.benchmark.converter import LinkMLConverter
from omni.benchmark.validation.error import ValidationError
from omni.utils import try_load_envmodule


class Validator:
    """Simple validator class for Benchmark."""

    def __init__(self):
        self.errors = []

    def validate(
        self, benchmark_dir: Path, converter: LinkMLConverter
    ) -> Union[ValidationError, LinkMLConverter]:
        # Validate ids are unique
        stage_ids = converter.get_stages().keys()
        duplicate_stage_ids = Validator.find_duplicate(stage_ids)
        if duplicate_stage_ids:
            self.errors.append(
                ValidationError(
                    f"Found duplicate stage ids: {', '.join(duplicate_stage_ids)}"
                )
            )

        module_ids = converter.get_modules().keys()
        duplicate_module_ids = Validator.find_duplicate(module_ids)
        if duplicate_module_ids:
            self.errors.append(
                ValidationError(
                    f"Found duplicate module ids: {', '.join(duplicate_module_ids)}"
                )
            )

        output_ids = converter.get_outputs().keys()
        duplicate_output_ids = Validator.find_duplicate(output_ids)
        if duplicate_output_ids:
            self.errors.append(
                ValidationError(
                    f"Found duplicate output ids: {', '.join(duplicate_output_ids)}"
                )
            )

        for stage_id in converter.get_stages():
            stage_inputs_set = converter.get_stage_implicit_inputs(stage_id)
            for stage_inputs in stage_inputs_set:
                for stage_input_id in stage_inputs:
                    if stage_input_id not in output_ids:
                        self.errors.append(
                            ValidationError(
                                f"Input with id '{stage_input_id}' in stage '{stage_id}' is not valid"
                            )
                        )

        # Validate software environments
        software_environments = converter.get_software_environments()
        for module_id, module in converter.get_modules().items():
            environment_id = module.software_environment
            if software_environments.get(environment_id) is None:
                self.errors.append(
                    ValidationError(
                        f"Software environment with id '{environment_id}' is not defined."
                    )
                )

        software_backend = converter.get_software_backend()
        for environment in software_environments.values():
            if software_backend == SoftwareBackendEnum.host:
                continue

            elif software_backend == SoftwareBackendEnum.envmodules:
                # available_modules = get_available_modules(environment.envmodule)
                # if len(available_modules) > 1:
                #     logging.warning(
                #         f"WARNING: Ambiguous envmodule name. Found the following modules matching the name: {available_modules}."
                #     )
                #     self.errors.append(
                #         ValidationError(
                #             f"Ambiguous envmodule name. Found the following modules for software environment with id '{environment.id}' matching the name: {available_modules}."
                #         )
                #     )
                # elif len(available_modules) == 0:
                #     self.errors.append(
                #         ValidationError(
                #             f"Modules with name `{environment.envmodule}` were not found for software environment with id '{environment.id}'."
                #         )
                #     )
                # else:
                load_result = try_load_envmodule(environment.envmodule)
                if not load_result:
                    self.errors.append(
                        ValidationError(
                            f"Software environment with id '{environment.id}' could not be loaded as a valid `envmodule`."
                        )
                    )

            elif (
                software_backend == SoftwareBackendEnum.conda
                or software_backend == SoftwareBackendEnum.docker
                or software_backend == SoftwareBackendEnum.apptainer
            ):
                environment_path = Validator.get_environment_path(
                    software_backend, environment, benchmark_dir
                )

                if not environment_path:
                    self.errors.append(
                        ValidationError(
                            f"Software environment with id '{environment.id}' does not have a valid backend definition for: '{software_backend.text}'."
                        )
                    )

                # else for other backends, attempt to find path on disk
                elif not Validator.is_url(environment_path) and not os.path.exists(
                    environment_path
                ):
                    self.errors.append(
                        ValidationError(
                            f"Software environment path for '{software_backend.text}' does not exist: '{environment_path}'."
                        )
                    )

        # Raise ValidationError if there are errors
        if self.errors:
            raise ValidationError(self.errors)
        else:
            return converter

    @staticmethod
    def find_duplicate(ids):
        id_counts = Counter(ids)
        duplicate_ids = [id for id, count in id_counts.items() if count > 1]

        return duplicate_ids

    @staticmethod
    def get_environment_path(
        software_backend: SoftwareBackendEnum,
        software: SoftwareEnvironment,
        benchmark_dir: Path,
    ) -> Optional[str]:
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

        if Validator.is_url(environment) or Validator.is_absolute_path(environment):
            environment_path = environment
        else:
            environment_path = os.path.join(benchmark_dir, environment)

        return environment_path

    @staticmethod
    def is_url(string: str) -> bool:
        # Check if the string is a valid URL using urlparse
        try:
            result = urlparse(string)
            return all(
                [result.scheme, result.netloc]
            )  # Valid if both scheme and netloc are present
        except ValueError:
            return False

    @staticmethod
    def is_absolute_path(string: str) -> bool:
        # Check if the string is an absolute path
        return Path(string).is_absolute()
