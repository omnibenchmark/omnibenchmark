import os
import re
from pathlib import Path
from typing import List

from omni_schema.datamodel.omni_schema import MetricCollector, SoftwareBackendEnum

from omnibenchmark.benchmark import Benchmark, Validator
from omnibenchmark.workflow.snakemake import scripts
from omnibenchmark.workflow.snakemake.format import formatter
from omnibenchmark.workflow.snakemake.scripts.parse_performance import write_combined_performance_file

import logging
logging.getLogger('snakemake').setLevel(logging.DEBUG)


def create_all_rule(paths: List[str], aggregate_performance: bool = False):
    if not aggregate_performance:
        rule all:
            input: paths,
            # "out/data/D2/default/D2.data.ext",
            # "out/data/D2/default/process/P1/default/D2.txt.gz",
            # "out/data/D1/default/process/P2/default/methods/M2/default/D1.model.out.gz"
            # "out/data/D1/default/process/P2/default/methods/M2/default/m1/default/D1.results.txt"
    else:
        rule all:
            input: paths
            output:
                f"{benchmark.out_dir}/performances.tsv"
            run:
                result = glob_wildcards(str(benchmark.out_dir) + "/{path}/{dataset}_performance.txt")
                performances = expand(str(benchmark.out_dir) + "/{path}/{dataset}_performance.txt", path=result.path, dataset=result.dataset)
                performances = sorted(list(set(performances)))

                output_dir = Path(str(os.path.commonpath(output)))
                if len(output) == 1:
                    output_dir = Path(os.path.dirname(output_dir))

                write_combined_performance_file(output_dir, performances)


def create_metric_collector_rule(benchmark: Benchmark, collector: MetricCollector, node_output_paths: List[str]):
    repository = collector.repository
    repository_url = repository.url if repository else None
    commit_hash = repository.commit if repository else None

    inputs_map = formatter.format_metric_collector_input(benchmark, collector, return_as_dict=True)

    updated_inputs_map = {}
    updated_inputs = []
    for key, pattern_str in inputs_map.items():
        pattern = _compile_regex_pattern_for_collectors_input(pattern_str)
        filtered_input_paths = [path for path in node_output_paths if pattern.match(path)]
        updated_inputs.extend(filtered_input_paths)
        updated_inputs_map[key] = filtered_input_paths

    software_backend = benchmark.get_benchmark_software_backend()

    # Only set environment directive for the backend that's actually being used
    if software_backend == SoftwareBackendEnum.conda:
        conda_env = _get_environment_paths(benchmark, collector, SoftwareBackendEnum.conda)
        rule:
            name: f"metric_collector_{{name}}".format(name=collector.id)
            input:
                expand(updated_inputs, allow_missing=True)
            output:
                formatter.format_metric_collector_output(benchmark.out_dir, collector)
            conda:
                conda_env
            params:
                inputs_map=updated_inputs_map,
                repository_url=repository_url,
                commit_hash=commit_hash,
                keep_module_logs=config['keep_module_logs']
            script: os.path.join(os.path.dirname(os.path.realpath(scripts.__file__)),'run_module.py')

    elif software_backend == SoftwareBackendEnum.envmodules:
        envmodules_env = _get_environment_paths(benchmark, collector, SoftwareBackendEnum.envmodules)
        rule:
            name: f"metric_collector_{{name}}".format(name=collector.id)
            input:
                expand(updated_inputs, allow_missing=True)
            output:
                formatter.format_metric_collector_output(benchmark.out_dir, collector)
            envmodules:
                envmodules_env
            params:
                inputs_map=updated_inputs_map,
                repository_url=repository_url,
                commit_hash=commit_hash,
                keep_module_logs=config['keep_module_logs']
            script: os.path.join(os.path.dirname(os.path.realpath(scripts.__file__)),'run_module.py')

    elif software_backend == SoftwareBackendEnum.apptainer or software_backend == SoftwareBackendEnum.docker:
        container_env = _get_environment_paths(benchmark, collector, SoftwareBackendEnum.apptainer)
        rule:
            name: f"metric_collector_{{name}}".format(name=collector.id)
            input:
                expand(updated_inputs, allow_missing=True)
            output:
                formatter.format_metric_collector_output(benchmark.out_dir, collector)
            container:
                container_env
            params:
                inputs_map=updated_inputs_map,
                repository_url=repository_url,
                commit_hash=commit_hash,
                keep_module_logs=config['keep_module_logs']
            script: os.path.join(os.path.dirname(os.path.realpath(scripts.__file__)),'run_module.py')

    else:
        # host or other backend - no environment directive needed
        rule:
            name: f"metric_collector_{{name}}".format(name=collector.id)
            input:
                expand(updated_inputs, allow_missing=True)
            output:
                formatter.format_metric_collector_output(benchmark.out_dir, collector)
            params:
                inputs_map=updated_inputs_map,
                repository_url=repository_url,
                commit_hash=commit_hash,
                keep_module_logs=config['keep_module_logs']
            script: os.path.join(os.path.dirname(os.path.realpath(scripts.__file__)),'run_module.py')


def _compile_regex_pattern_for_collectors_input(pattern: str) -> re.Pattern[str]:
    pattern_regex = (
        pattern
        .replace('{input}',r'.+')  # {input} can be anything before "/metrics/"
        .replace('{stage}',r'[^/]+') # {stage} is a single folder
        .replace('{module}',r'[^/]+')  # {module} is a single folder
        .replace('{params}',r'[^/]+')  # {params} is a single folder
        .replace('{dataset}',r'[^/]+')  # {dataset} is the filename prefix
    )

    pattern_regex = re.compile(r'^' + pattern_regex + r'$')
    return pattern_regex


def _get_environment_paths(benchmark: Benchmark, collector: MetricCollector, software_backend: SoftwareBackendEnum) -> str:
    benchmark_dir = benchmark.directory
    environment = benchmark.get_benchmark_software_environments()[collector.software_environment]
    environment_path = Validator.get_environment_path(software_backend, environment, benchmark_dir)

    return environment_path
