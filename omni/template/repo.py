from pathlib import Path
from typing import List
import os

from omni.benchmark import Benchmark
from omni.io.code import get_git_archive, check_remote_repo_existance

# parse benchmark

# get all repositories

# check if repository (and commit hash) exists (e.g. is reachable), and if yes if actually is omnibenchmark repo, i.e. contains config.cfg file and the specified file in the config.cfg file

# if not, create a new repository:
# - figure out which language: python, R (by description/name of software env used?)
# - based on template in /data/templates/main.(py|R) create a updated file that includes the expected argparse arguments (derived from input/output in yaml)
# - also add config.cfg file
# - create repo and add files


# template adaption:
from enum import Enum


class filetype(Enum):
    python = "py"
    R = "R"
    sh = "sh"


class languageregex(Enum):
    python = "(P|p)ython"
    R = "(R|r)(_|-| )"
    sh = "sh"


# offset of each line in the script
script_offsets = {
    "argparse": {
        "python": "    ",
        "R": "",
        "sh": "    ",
    },
    "variable": {
        "python": "    ",
        "R": "",
        "sh": "",
    },
}

# replacement line for argparse
argparse_dict_base = {
    "python": "parser.add_argument('--{argument}', type=str, help='input file {argument}')",
    "R": "parser$add_argument('--{argument}', type='character', help='input file {argument}')",
    "sh": '--{argument})\n      {argument}="$2"\n      shift\n      shift\n      ;;',
}

# replacement line for variable
replace_dict_base = {
    "implicit": {
        "input": {
            "python": "input_file_{iter} = snakemake.input['{argument}']",
            "R": "input_file_{iter} <- snakemake@input[['{argument}']]",
            "sh": "input_file_{iter}=${{snakemake_input[{argument}]}}\n",
        },
        "params": {
            "python": "params_file_{iter} = snakemake.params['{argument}']",
            "R": "params_file_{iter} <- snakemake@params[['{argument}']]",
            "sh": "params_file_{iter}=${{snakemake_params[{argument}]}}\n",
        },
        "output": {
            "python": "output_file_{iter} = os.path.join(snakemake.output['output_dir'], f'{{snakemake.output['name']}}{argument}')",
            "R": "output_file_{iter} <- file.path(snakemake@output[['output_dir']], paste0(snakemake@output[['name']], '{argument}'))",
            "sh": "output_file_{iter}=${{snakemake_output[output_dir]}}/${{snakemake_output[name]}}{argument}\n",
        },
    },
    "explicit": {
        "input": {
            "python": "input_file_{iter} = getattr(args, '{argument}')",
            "R": "input_file_{iter} = args[['{argument}']]",
            "sh": "input_file_{iter}=${{{argument}}}\n",
        },
        "params": {
            "python": "params_file_{iter} = getattr(args, '{argument}')",
            "R": "params_file_{iter} = args[['{argument}']]",
            "sh": "params_file_{iter}=${{{argument}}}\n",
        },
        "output": {
            "python": "output_file_{iter} = os.path.join(getattr(args, 'output_dir'), f'{{getattr(args, 'name')}}{argument}')",
            "R": "output_file_{iter} = file.path(args[['output_dir']], paste0(args[['name']], 'args[[{argument}]]'))",
            "sh": "output_file_{iter}=${{output_dir}}/${{name}}{argument}\n",
        },
    },
}


def create_argparse_expression(
    argument: str, description: str = "", language: str = "python"
) -> str:
    assert language in [l.name for l in list(filetype)]
    return argparse_dict_base[language].format(
        argument=argument, description=description
    )


def create_variable_expression(
    argument: str,
    iter: int = 1,
    vartype="input",
    language: str = "python",
    script_type: str = "explicit",
) -> str:
    assert vartype in ["input", "params", "output"]
    assert script_type in ["explicit", "implicit"]
    assert language in [l.name for l in list(filetype)]
    return replace_dict_base[script_type][vartype][language].format(
        argument=argument, iter=iter
    )


def prepare_template(
    node, language: str = "python", script_type: str = "explicit"
) -> str:
    assert language in [l.name for l in list(filetype)]
    input_argparse_str, input_variable_str = "", ""
    try:
        inputs = node.get_explicit_inputs()[0]
        for i, input in enumerate(inputs.keys()):
            input_argparse_str += f"{script_offsets['argparse'][language]}{create_argparse_expression(input, language)}\n"
            input_variable_str += f"{script_offsets['variable'][language]}{create_variable_expression(input, i + 1, 'input', language, script_type)}\n"
    except Exception:
        pass
    parameters_argparse_str, parameters_variable_str = "", ""
    try:
        parameters = node.get_parameters()[::2]
        for i, params in enumerate(parameters):
            parameters_argparse_str += f"{script_offsets['argparse'][language]}{create_argparse_expression(params.replace("--", ""), language)}\n"
            parameters_variable_str += f"{script_offsets['variable'][language]}{create_variable_expression(params, i+1, "param", language, script_type)}\n"
    except Exception:
        pass

    outputs_variable_str = ""
    try:
        outputs = node.get_outputs()
        for i, output in enumerate(outputs):
            output_basename = os.path.basename(output)
            ti = output_basename.find("}")
            file_ending = output_basename[(ti + 1) :]
            outputs_variable_str += f"{script_offsets['variable'][language]}{create_variable_expression(file_ending, i + 1, 'output', language, script_type)}\n"
    except Exception:
        pass

    with open(
        f"data/templates/main_{script_type}.{filetype[language].value}", "r"
    ) as f:
        template = f.read()
        replace_dict = {
            "INPUTS_ARGPARSE": input_argparse_str,
            "PARAMETERS_ARGPARSE": parameters_argparse_str,
            "VARIABLES_INPUTS": input_variable_str,
            "VARIABLES_PARAMETERS": parameters_variable_str,
            "VARIABLES_OUTPUTS": outputs_variable_str,
        }
        template = template.format(**replace_dict)
    return template


def get_missing_repos(benchmark: Benchmark) -> List:
    modules = {}
    for i, node in enumerate(benchmark.get_nodes()):
        modules[node.get_repository()["url"]] = i

    unique_module_ids = list(modules.values())
    nodes = [no for i, no in enumerate(benchmark.get_nodes()) if i in unique_module_ids]

    inexistent_nodes = []
    for node in nodes:
        repo = node.get_repository()
        if not get_git_archive("", repo["url"], repo["commit"], True):
            inexistent_nodes.append(node)

    repo_urls = [node.get_repository()["url"] for node in inexistent_nodes]
    return repo_urls, inexistent_nodes


def guess_language(node):
    search_list = [node.module.software_environment]
    import re

    for language in languageregex:
        if any(re.search(language.value, x) for x in search_list):
            return language.name
    return None


def create_template(
    node, language: str = None, outfile: str = None, script_type="explicit"
):
    if language is None:
        language = guess_language(node)
    if outfile is None:
        return prepare_template(node, language=language, script_type=script_type)
    else:
        with open(outfile, "w") as f:
            f.write(prepare_template(node, language=language, script_type=script_type))


def create_repo_files(
    node,
    language: str = None,
    output_dir: Path = Path("output"),
    script_type="explicit",
):
    if language is None:
        language = guess_language(node)
    assert language in [l.name for l in list(filetype)]

    create_template(
        node,
        language=language,
        outfile=output_dir / Path(f"main.{filetype[language].value}"),
        script_type=script_type,
    )

    config_cfg_str = f"[DEFAULT]\nSCRIPT=main.{filetype[language].value}"
    with open(output_dir / Path("config.cfg"), "w") as f:
        f.write(config_cfg_str)
    return output_dir


if __name__ == "__main__":
    benchmark = Benchmark(Path("tests/data/Clustering.yaml"))
    benchmark
    urls, nodes = get_missing_repos(benchmark)

    # print(prepare_template(nodes[0], language="R"))
    # print(prepare_template(nodes[0], language="python"))
    create_template(nodes[0])
    create_repo_files(
        nodes[0], output_dir=Path("tmp"), language="python", script_type="implicit"
    )
    create_repo_files(nodes[2], output_dir=Path("tmp"))
    # prepare_template(nodes[2])

    repository_url = urls[0]
    check_remote_repo_existance(repository_url)
