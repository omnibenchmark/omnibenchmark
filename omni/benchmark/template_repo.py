from pathlib import Path
from typing import List
import os

from omni.benchmark import Benchmark
from omni.io.code import get_git_archive

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


class languageregex(Enum):
    python = "(P|p)ython"
    R = "(R|r)(_|-| )"


def create_argparse_expression(
    argument: str, description: str = "", language: str = "python"
) -> str:
    assert language in [l.name for l in list(filetype)]
    if language == "python":
        return f"parser.add_argument('--{argument}', type=str, help='{description}', required = True)"
    elif language == "R":
        return f"parser$add_argument('--{argument}', type='character', help='{description}')"
    else:
        raise ValueError(f"Language {language} not supported.")


def create_variable_expression(
    argument: str, iter: int = 1, vartype="input", language: str = "python"
) -> str:
    assert vartype in ["input", "parameter", "output"]
    assert language in [l.name for l in list(filetype)]
    if language == "python":
        if vartype == "output":
            return f"{vartype}_file_{iter} = os.path.join(getattr(args, 'output_dir'), f'{{getattr(args, 'name')}}{argument}')"
        else:
            return f"{vartype}_file_{iter} = getattr(args, '{argument}')"
    elif language == "R":
        if vartype == "output":
            return f"{vartype}_file_{iter} = file.path(args[['output_dir']], paste0(args[['name']], 'args[[{argument}]]'))"
        else:
            return f"{vartype}_file_{iter} = args[['{argument}']]"
    else:
        raise ValueError(f"Language {language} not supported.")


def prepare_template(node, language: str = "python") -> str:
    assert language in [l.name for l in list(filetype)]
    input_argparse_str = ""
    input_variable_str = ""
    try:
        inputs = node.get_explicit_inputs()[0]
        for i, input in enumerate(inputs.keys()):
            input_argparse_str += (
                "    " + create_argparse_expression(input, language=language) + "\n"
            )
            input_variable_str += (
                "    "
                + create_variable_expression(input, i + 1, "input", language=language)
                + "\n"
            )
    except:
        pass

    parameters_argparse_str = ""
    parameters_variable_str = ""
    try:
        parameters = node.get_parameters()[::2]
        for i, parameter in enumerate(parameters):
            parameters_argparse_str += (
                "    "
                + create_argparse_expression(
                    parameter.replace("--", ""), language=language
                )
                + "\n"
            )
            parameters_variable_str += (
                "    "
                + create_variable_expression(
                    parameter, i + 1, "parameter", language=language
                )
                + "\n"
            )
    except:
        pass

    outputs_variable_str = ""
    outputs = node.get_outputs()
    try:
        for i, output in enumerate(outputs):
            output_basename = os.path.basename(output)
            ti = output_basename.find("}")
            file_ending = output_basename[(ti + 1) :]
            outputs_variable_str += (
                "    "
                + create_variable_expression(
                    file_ending, i + 1, "output", language=language
                )
                + "\n"
            )
    except:
        pass

    with open(f"data/templates/main.{filetype[language].value}", "r") as f:
        template = f.read()
        template = template.replace("    INPUTS_ARGPARSE", input_argparse_str)
        template = template.replace("    PARAMETERS_ARGPARSE", parameters_argparse_str)
        template = template.replace("    VARIABLES_INPUTS", input_variable_str)
        template = template.replace("    VARIABLES_PARAMETERS", parameters_variable_str)
        template = template.replace("    VARIABLES_OUTPUTS", outputs_variable_str)

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
    return inexistent_nodes


def guess_language(node):
    search_list = [node.module.software_environment]
    import re

    for language in languageregex:
        if any(re.search(language.value, x) for x in search_list):
            return language.name
    return None


def create_template(node, language: str = None, outfile: str = None):
    if language is None:
        language = guess_language(node)
    if outfile is None:
        return prepare_template(node, language=language)
    else:
        with open(outfile, "w") as f:
            f.write(prepare_template(node, language=language))


def create_repo_files(node, language: str = None, output_dir: Path = Path("output")):
    if language is None:
        language = guess_language(node)
    assert language in [l.name for l in list(filetype)]

    create_template(
        node,
        language=language,
        outfile=output_dir / Path(f"main.{filetype[language].value}"),
    )

    config_cfg_str = f"[DEFAULT]\nSCRIPT=main.{filetype[language].value}"
    with open(output_dir / Path("config.cfg"), "w") as f:
        f.write(config_cfg_str)
    return output_dir


if __name__ == "__main__":
    benchmark = Benchmark(Path("tests/data/Clustering.yaml"))
    benchmark
    nodes = get_missing_repos(benchmark)

    # print(prepare_template(nodes[0], language="R"))
    # print(prepare_template(nodes[0], language="python"))
    create_template(nodes[0])
    create_repo_files(nodes[0], output_dir=Path("tmp"))
    # prepare_template(nodes[2])
