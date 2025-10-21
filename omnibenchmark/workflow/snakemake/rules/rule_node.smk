from importlib import resources
from pathlib import Path

from omni_schema.datamodel.omni_schema import SoftwareBackendEnum, Benchmark

from omnibenchmark.benchmark import Validator, BenchmarkNode
from omnibenchmark.workflow.snakemake import scripts
from omnibenchmark.workflow.snakemake.format import formatter

RUN_MODULE = "run_module.py"


def get_script_path(script_name: str) -> str:
    """Get the filesystem path to a script in the scripts package.
    """
    path = Path(resources.files(scripts) / script_name)
    return str(path)

def create_node_rule(node, benchmark, config, local_timeout):
    if node.is_initial():
        return _create_initial_node(benchmark, node, config, local_timeout)
    else:
        return _create_intermediate_node(benchmark, node, config, local_timeout)

def _create_initial_node(benchmark, node, config, local_timeout):
    stage_id = node.stage_id
    module_id = node.module_id
    param_id = node.param_id

    out_dir = config['out_dir']

    repository = node.get_repository()
    repository_url = repository.url if repository else None
    commit_hash = repository.commit if repository else None

    software_backend = benchmark.get_benchmark_software_backend()

    # Only set environment directive for the backend that's actually being used
    if software_backend == SoftwareBackendEnum.conda:
        conda_env = _get_environment_path(benchmark, node, SoftwareBackendEnum.conda)
        rule:
            name: f"{{stage}}_{{module}}_{{param}}".format(stage=stage_id,module=module_id,param=param_id)
            wildcard_constraints:
                stage=stage_id,
                module=module_id,
                params=param_id,
                dataset=module_id
            benchmark:
                formatter.format_performance_file(node, out_dir)
            output:
                formatter.format_output_templates_to_be_expanded(node)
            conda:
                conda_env
            params:
                repository_url = repository_url,
                commit_hash = commit_hash,
                parameters = node.get_parameters(),
                dataset = module_id,
                keep_module_logs=config['keep_module_logs'],
                local_task_timeout=local_timeout
            script: get_script_path(RUN_MODULE)

    elif software_backend == SoftwareBackendEnum.envmodules:
        envmodules_env = _get_environment_path(benchmark, node, SoftwareBackendEnum.envmodules)
        rule:
            name: f"{{stage}}_{{module}}_{{param}}".format(stage=stage_id,module=module_id,param=param_id)
            wildcard_constraints:
                stage=stage_id,
                module=module_id,
                params=param_id,
                dataset=module_id
            benchmark:
                formatter.format_performance_file(node, out_dir)
            output:
                formatter.format_output_templates_to_be_expanded(node)
            envmodules:
                envmodules_env
            params:
                repository_url = repository_url,
                commit_hash = commit_hash,
                parameters = node.get_parameters(),
                dataset = module_id,
                keep_module_logs=config['keep_module_logs'],
                local_task_timeout=local_timeout
            script: get_script_path(RUN_MODULE)

    elif software_backend == SoftwareBackendEnum.apptainer or software_backend == SoftwareBackendEnum.docker:
        container_env = _get_environment_path(benchmark, node, SoftwareBackendEnum.apptainer)
        rule:
            name: f"{{stage}}_{{module}}_{{param}}".format(stage=stage_id,module=module_id,param=param_id)
            wildcard_constraints:
                stage=stage_id,
                module=module_id,
                params=param_id,
                dataset=module_id
            benchmark:
                formatter.format_performance_file(node, out_dir)
            output:
                formatter.format_output_templates_to_be_expanded(node)
            container:
                container_env
            params:
                repository_url = repository_url,
                commit_hash = commit_hash,
                parameters = node.get_parameters(),
                dataset = module_id,
                keep_module_logs=config['keep_module_logs'],
                local_task_timeout=local_timeout
            script: get_script_path(RUN_MODULE)

    else:
        # host or other backend - no environment directive needed
        rule:
            name: f"{{stage}}_{{module}}_{{param}}".format(stage=stage_id,module=module_id,param=param_id)
            wildcard_constraints:
                stage=stage_id,
                module=module_id,
                params=param_id,
                dataset=module_id
            benchmark:
                formatter.format_performance_file(node, out_dir)
            output:
                formatter.format_output_templates_to_be_expanded(node)
            params:
                repository_url = repository_url,
                commit_hash = commit_hash,
                parameters = node.get_parameters(),
                dataset = module_id,
                keep_module_logs=config['keep_module_logs'],
                local_task_timeout=local_timeout
            script: get_script_path(RUN_MODULE)


def _create_intermediate_node(benchmark, node, config, local_timeout):
    stage_id = node.stage_id
    module_id = node.module_id
    param_id = node.param_id

    outputs = node.get_outputs()

    out_dir = config['out_dir']

    post = stage_id + '/' + module_id
    if any(['{params}' in o for o in outputs]):
        post += '/' + param_id

    repository = node.get_repository()
    repository_url = repository.url if repository else None
    commit_hash = repository.commit if repository else None

    inputs_map = lambda wildcards: formatter.format_input_templates_to_be_expanded(benchmark, wildcards, return_as_dict=True)

    keep_going = config.get("keep_going", False)
    # just for intermediate notes, we allow the rules to touch their expected
    # outputs. the logic is that for initial nodes we want to fail hard,
    # but we want the collectors to meet the preconditions for execution.

    fmt_fn = formatter.format_output_templates_to_be_expanded

    # Get the software backend and corresponding environment path
    software_backend = benchmark.get_benchmark_software_backend()

    # Only set environment directive for the backend that's actually being used
    if software_backend == SoftwareBackendEnum.conda:
        conda_env = _get_environment_path(benchmark, node, SoftwareBackendEnum.conda)
        rule:
            name: f"{{stage}}_{{module}}_{{param}}".format(stage=stage_id,module=module_id,param=param_id)
            wildcard_constraints:
                post=post,
                stage=stage_id,
                module=module_id
            input:
                lambda wildcards: formatter.format_input_templates_to_be_expanded(benchmark, wildcards)
            benchmark:
                formatter.format_performance_file(node, out_dir)
            output: touch(fmt_fn(node)) if keep_going else fmt_fn(node)
            conda:
                conda_env
            params:
                inputs_map = inputs_map,
                repository_url = repository_url,
                commit_hash = commit_hash,
                parameters = node.get_parameters(),
                keep_module_logs=config['keep_module_logs'],
                local_task_timeout=local_timeout
            script: get_script_path(RUN_MODULE)
    elif software_backend == SoftwareBackendEnum.envmodules:
        envmodules_env = _get_environment_path(benchmark, node, SoftwareBackendEnum.envmodules)
        rule:
            name: f"{{stage}}_{{module}}_{{param}}".format(stage=stage_id,module=module_id,param=param_id)
            wildcard_constraints:
                post=post,
                stage=stage_id,
                module=module_id
            input:
                lambda wildcards: formatter.format_input_templates_to_be_expanded(benchmark, wildcards)
            benchmark:
                formatter.format_performance_file(node, out_dir)
            output: touch(fmt_fn(node)) if keep_going else fmt_fn(node)
            envmodules:
                envmodules_env
            params:
                inputs_map = inputs_map,
                repository_url = repository_url,
                commit_hash = commit_hash,
                parameters = node.get_parameters(),
                keep_module_logs=config['keep_module_logs'],
                local_task_timeout=local_timeout
            script: get_script_path(RUN_MODULE)
    elif software_backend == SoftwareBackendEnum.apptainer or software_backend == SoftwareBackendEnum.docker:
        container_env = _get_environment_path(benchmark, node, SoftwareBackendEnum.apptainer)
        rule:
            name: f"{{stage}}_{{module}}_{{param}}".format(stage=stage_id,module=module_id,param=param_id)
            wildcard_constraints:
                post=post,
                stage=stage_id,
                module=module_id
            input:
                lambda wildcards: formatter.format_input_templates_to_be_expanded(benchmark, wildcards)
            benchmark:
                formatter.format_performance_file(node, out_dir)
            output: touch(fmt_fn(node)) if keep_going else fmt_fn(node)
            container:
                container_env
            params:
                inputs_map = inputs_map,
                repository_url = repository_url,
                commit_hash = commit_hash,
                parameters = node.get_parameters(),
                keep_module_logs=config['keep_module_logs'],
                local_task_timeout=local_timeout
            script: get_script_path(RUN_MODULE)
    else:
        # host or other backend - no environment directive needed
        rule:
            name: f"{{stage}}_{{module}}_{{param}}".format(stage=stage_id,module=module_id,param=param_id)
            wildcard_constraints:
                post=post,
                stage=stage_id,
                module=module_id
            input:
                lambda wildcards: formatter.format_input_templates_to_be_expanded(benchmark, wildcards)
            benchmark:
                formatter.format_performance_file(node, out_dir)
            output: touch(fmt_fn(node)) if keep_going else fmt_fn(node)
            params:
                inputs_map = inputs_map,
                repository_url = repository_url,
                commit_hash = commit_hash,
                parameters = node.get_parameters(),
                keep_module_logs=config['keep_module_logs'],
                local_task_timeout=local_timeout
            script: get_script_path(RUN_MODULE)


def create_standalone_node_rule(node, config):
    stage_id = node.stage_id
    module_id = node.module_id
    param_id = node.param_id

    repository = node.get_repository()
    repository_url = repository.url if repository else None
    commit_hash = repository.commit if repository else None

    # TODO(ben): can factor out common parts?

    if node.is_initial():
        rule:
            name: f"{{stage}}_{{module}}_{{param}}".format(stage=stage_id,module=module_id,param=param_id)
            wildcard_constraints:
                dataset=config['dataset']
            output:
                node.get_output_paths(config)
            params:
                repository_url = repository_url,
                commit_hash = commit_hash,
                parameters=node.get_parameters(),
                dataset=config['dataset'],
                keep_module_logs=config['keep_module_logs'],
            benchmark:
                node.get_benchmark_path(config)
            script: get_script_path(RUN_MODULE)
    else:
        rule:
            name: f"{{stage}}_{{module}}_{{param}}".format(stage=stage_id,module=module_id,param=param_id)
            input:
                node.get_input_paths(config)
            output:
                node.get_output_paths(config)
            params:
                inputs_map = node.get_input_paths(config, return_as_dict=True),
                repository_url = repository_url,
                commit_hash = commit_hash,
                parameters=node.get_parameters(),
                dataset=config['dataset'],
                keep_module_logs=config['keep_module_logs'],
            benchmark:
                node.get_benchmark_path(config)
            script: get_script_path(RUN_MODULE)


def _get_environment_path(benchmark: Benchmark, node: BenchmarkNode, software_backend: SoftwareBackendEnum):
    benchmark_dir = benchmark.directory
    environment = benchmark.get_benchmark_software_environments()[node.get_software_environment()]
    environment_path = Validator.get_environment_path(software_backend, environment, benchmark_dir)

    return environment_path
