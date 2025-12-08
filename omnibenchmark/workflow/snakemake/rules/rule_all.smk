import os
import re
from importlib import resources
from pathlib import Path
from typing import List, Any

from omnibenchmark.model import MetricCollector, SoftwareBackendEnum

from omnibenchmark.benchmark import Benchmark, Validator
from omnibenchmark.workflow.snakemake import scripts
from omnibenchmark.workflow.snakemake.format import formatter

import logging
logging.getLogger('snakemake').setLevel(logging.DEBUG)

RUN_MODULE = "run_module.py"


def get_script_path(script_name: str) -> str:
    """Get the filesystem path to a script in the scripts package.
    """
    path = Path(resources.files(scripts) / script_name)
    return str(path)

def create_all_rule(config, paths: List[str], aggregate_performance: bool = False):
    out_dir = config['out_dir']

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
                f"{benchmark.context.out_dir}/performances.tsv"
            script:
                get_script_path("aggregate_performance.py")


def create_metric_collector_rule(benchmark: Benchmark, collector: MetricCollector, config: dict[str, Any], node_output_paths: List[str]):
    repository = collector.repository
    repository_url = repository.url if repository else None
    commit_hash = repository.commit if repository else None

    out_dir = config['out_dir']
    software_backend = config['backend']
    keep_module_logs = config.get('keep_module_logs', False)

    inputs_map = formatter.format_metric_collector_input(benchmark, collector, return_as_dict=True)

    updated_inputs_map = {}
    updated_inputs = []
    for key, pattern_str in inputs_map.items():
        pattern = _compile_regex_pattern_for_collectors_input(pattern_str)
        filtered_input_paths = [path for path in node_output_paths if pattern.match(path)]
        updated_inputs.extend(filtered_input_paths)
        updated_inputs_map[key] = filtered_input_paths

    # Only set environment directive for the backend that's actually being used
    # TODO https://github.com/omnibenchmark/omnibenchmark/issues/201
    # TODO Factor out the conditional rule generation when working on the above issue
    if software_backend == SoftwareBackendEnum.conda:
        conda_env = benchmark.get_environment_path(collector.software_environment, SoftwareBackendEnum.conda)
        rule:
            name: f"metric_collector_{{name}}".format(name=collector.id)
            input:
                expand(updated_inputs, allow_missing=True)
            output:
                formatter.format_metric_collector_output(out_dir, collector)
            conda:
                conda_env
            params:
                inputs_map=updated_inputs_map,
                repository_url=repository_url,
                commit_hash=commit_hash,
                keep_module_logs=keep_module_logs
            script: get_script_path(RUN_MODULE)

    elif software_backend == SoftwareBackendEnum.envmodules:
        envmodules_env = benchmark.get_environment_path(collector.software_environment, SoftwareBackendEnum.envmodules)
        rule:
            name: f"metric_collector_{{name}}".format(name=collector.id)
            input:
                expand(updated_inputs, allow_missing=True)
            output:
                formatter.format_metric_collector_output(out_dir, collector)
            envmodules:
                envmodules_env
            params:
                inputs_map=updated_inputs_map,
                repository_url=repository_url,
                commit_hash=commit_hash,
                keep_module_logs=keep_module_logs
            script: get_script_path(RUN_MODULE)

    elif software_backend == SoftwareBackendEnum.apptainer or software_backend == SoftwareBackendEnum.docker:
        container_env = benchmark.get_environment_path(collector.software_environment, SoftwareBackendEnum.apptainer)
        rule:
            name: f"metric_collector_{{name}}".format(name=collector.id)
            input:
                expand(updated_inputs, allow_missing=True)
            output:
                formatter.format_metric_collector_output(out_dir, collector)
            container:
                container_env
            params:
                inputs_map=updated_inputs_map,
                repository_url=repository_url,
                commit_hash=commit_hash,
                keep_module_logs=keep_module_logs
            script: get_script_path(RUN_MODULE)

    else:
        # host or other backend - no environment directive needed
        rule:
            name: f"metric_collector_{{name}}".format(name=collector.id)
            input:
                expand(updated_inputs, allow_missing=True)
            output:
                formatter.format_metric_collector_output(out_dir, collector)
            params:
                inputs_map=updated_inputs_map,
                repository_url=repository_url,
                commit_hash=commit_hash,
                keep_module_logs=keep_module_logs
            script: get_script_path(RUN_MODULE)


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
