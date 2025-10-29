from datetime import datetime
import os
from pathlib import Path
from typing import Union
from omnibenchmark.benchmark import Benchmark
from omnibenchmark.cli.utils.validation import validate_benchmark
from omnibenchmark.workflow.snakemake.scripts.utils import (
    generate_unique_repo_folder_name,
)

# ndstdf = pd.DataFrame(
#     {
#         "stage": [st for st in stages for nd in b.get_nodes_by_stage_id(st)],
#         "module": [
#             nd.get_id().split("-")[1]
#             for st in stages
#             for nd in b.get_nodes_by_stage_id(st)
#         ],
#     }
# )
# perfdf = pd.read_csv(f"{out_dir}/performances.tsv", sep="\t")
# perfdf = perfdf.merge(ndstdf, on = "module")

# perfdf_grouped = perfdf.groupby(['module', 'stage'])['s'].agg([
#     'mean',
#     'std',
#     'min',
#     'max',
#     'median',
#     ('q25', lambda x: x.quantile(0.25)),
#     ('q75', lambda x: x.quantile(0.75))
# ]).reset_index()


def get_module_timestamps(b: Benchmark) -> dict:
    """
    Get the timestamps of the module repositories used in the benchmark.
    Returns a dictionary with stage ids as keys and another dictionary as value,
    which has module ids as keys and timestamps as values.
    """

    repositories_dir = Path(".snakemake") / "repos"
    modules_repo_timestamps = {}
    for st in b.get_stage_ids():
        modules_repo_timestamps[st] = {}
        for node in b.get_nodes_by_stage_id(st):
            if node.module_id not in modules_repo_timestamps[st].keys():
                repo = node.get_repository()
                module_dir = repositories_dir / generate_unique_repo_folder_name(
                    repo["url"], repo["commit"]
                )
                if module_dir.is_dir():
                    timestamps = [
                        f.stat().st_mtime
                        for f in module_dir.iterdir()
                        if f.is_file() and f.stat().st_size > 0
                    ]
                    if len(timestamps) > 0:
                        timestamp = max(timestamps)
                    else:
                        timestamp = None
                else:
                    timestamp = None
                modules_repo_timestamps[st][node.module_id] = timestamp
    return modules_repo_timestamps


def get_exec_path_dict(b: Benchmark, modules_repo_timestamps: dict) -> dict:
    """
    Get the execution path dictionary for the benchmark.
    Returns a dictionary with execution path ids as keys and another dictionary as value,
    which has stage ids as keys and information about the nodes and their output files as values.
    """

    config = {
        "input": "",
        "output": "",
        "dataset": "[a-zA-Z0-9_.]*",
    }
    exec_paths = b.get_execution_paths()
    exec_path_dict = {}
    for exec_path_id in range(len(exec_paths)):
        exec_path_dict[exec_path_id] = {}
        exec_path = exec_paths[exec_path_id]

        # prepare init node
        node = exec_path[0]
        dataset = node.get_id().split("-")[1]

        tmp_config = config.copy()
        # tmp_config['input'] = os.path.commonpath(init_dirnames)
        tmp_config["dataset"] = dataset
        init_files = node.get_output_paths(tmp_config)
        init_files_exists = [
            Path(f).is_file() and Path(f).stat().st_size > 0 for f in init_files
        ]

        init_timestamps = [
            Path(f).stat().st_mtime if init_files_exists[i] else None
            for i, f in enumerate(init_files)
        ]
        max_init_timestamp = (
            datetime.now().timestamp() * 2
            if any([tst is None for tst in init_timestamps])
            else max([tst for tst in init_timestamps])
        )

        init_is_newer_file = init_files_exists

        init_timestamp_repo = modules_repo_timestamps[node.stage_id][node.module_id]
        if init_timestamp_repo is not None:
            init_is_newer_repo = [
                True if tst is None else tst > init_timestamp_repo
                for tst in init_timestamps
            ]
        else:
            init_is_newer_repo = [False for tst in init_timestamps]

        # dirnames of observed files for node
        init_dirnames = list(set([os.path.dirname(f) for f in init_files]))

        exec_path_dict[exec_path_id][node.stage_id] = {
            "node": node,
            "files": init_files,
            "exists": init_files_exists,
            "timestamps": init_timestamps,
            "dirnames": init_dirnames,
            "is_newer_file": [
                True for f in init_files
            ],  # shouldn't have any input files, thus should be true
            "is_newer_repo": init_is_newer_repo,
            "is_newer": init_is_newer_repo,
        }

        # iterate over rest of nodes in exec path
        for iter in range(len(exec_path) - 1):
            node2 = exec_path[iter + 1]
            tmp_config = config.copy()
            tmp_config["input"] = os.path.commonpath(init_dirnames)
            tmp_config["dataset"] = dataset
            matched_files = node2.get_output_paths(tmp_config)
            matched_files_exists = [
                Path(f).is_file() and Path(f).stat().st_size > 0 for f in matched_files
            ]

            matched_dirnames = list(set([os.path.dirname(f) for f in matched_files]))

            # check which dirnames match previous ones
            matched_timestamps = [
                Path(f).stat().st_mtime if matched_files_exists[i] else None
                for i, f in enumerate(matched_files)
            ]

            # compare timestamps with previous files from previous nodes
            # if any previous timestamp is None, consider all current timestamps as older
            max_init_timestamp = max(
                [
                    (
                        datetime.now().timestamp() * 2
                        if any([tst is None for tst in init_timestamps])
                        else max([tst for tst in init_timestamps])
                    ),
                    max_init_timestamp,
                ]
            )
            if all([iin for iin in init_is_newer_file if iin is not None]):
                matched_is_newer_file = [
                    None if tst is None else tst > max_init_timestamp
                    for tst in matched_timestamps
                ]
            else:
                matched_is_newer_file = [None for tst in matched_timestamps]

            matched_timestamp_repo = modules_repo_timestamps[node2.stage_id][
                node2.module_id
            ]
            if matched_timestamp_repo is not None:
                matched_is_newer_repo = [
                    True if tst is None else tst > matched_timestamp_repo and iin
                    for tst, iin in zip(matched_timestamps, init_is_newer_repo)
                ]
            else:
                matched_is_newer_repo = [False for tst in matched_timestamps]

            exec_path_dict[exec_path_id][node2.stage_id] = {
                "node": node2,
                "files": matched_files,
                "exists": matched_files_exists,
                "timestamps": matched_timestamps,
                "dirnames": matched_dirnames,
                "is_newer_file": matched_is_newer_file,
                "is_newer_repo": matched_is_newer_repo,
                "is_newer": [
                    i and j
                    for i, j in zip(matched_is_newer_file, matched_is_newer_repo)
                ],
            }
            init_dirnames = matched_dirnames
            init_timestamps = matched_timestamps
            init_is_newer_file = matched_is_newer_file
            init_timestamp_repo = matched_timestamp_repo
    return exec_path_dict


def get_filedict(b: Benchmark, exec_path_dict: dict) -> dict:
    """
    Get the file dictionary for the benchmark.
    Returns a nested dictionary with stage ids as keys, node ids as sub-keys,
    and information about observed, missing, and invalid files as values.
    """
    filedict = {
        **{
            st: {
                nd: {
                    "observed_files": set(),
                    "missing_files": set(),
                    "invalid_files_file": set(),
                    "invalid_files_repo": set(),
                    "invalid_files": set(),
                }
                for nd in b.get_nodes_by_stage_id(st)
            }
            for st in b.get_stage_ids()
        }
    }
    for st in b.get_stage_ids():
        for exec_path_id in exec_path_dict.keys():
            if st in exec_path_dict[exec_path_id].keys():
                de = exec_path_dict[exec_path_id][st]
                filedict[st][de["node"]]["observed_files"] = filedict[st][de["node"]][
                    "observed_files"
                ].union(
                    set(
                        [
                            f
                            for i, f in enumerate(de["files"])
                            if de["exists"][i] and de["is_newer"][i]
                        ]
                    )
                )
                filedict[st][de["node"]]["missing_files"] = filedict[st][de["node"]][
                    "missing_files"
                ].union(
                    set([f for i, f in enumerate(de["files"]) if not de["exists"][i]])
                )
                filedict[st][de["node"]]["invalid_files"] = filedict[st][de["node"]][
                    "invalid_files"
                ].union(
                    set(
                        [
                            f
                            for i, f in enumerate(de["files"])
                            if de["exists"][i] and (not de["is_newer"][i])
                        ]
                    )
                )
                filedict[st][de["node"]]["invalid_files_file"] = filedict[st][
                    de["node"]
                ]["invalid_files_file"].union(
                    set(
                        [
                            f
                            for i, f in enumerate(de["files"])
                            if not de["is_newer_file"][i]
                        ]
                    )
                )
                filedict[st][de["node"]]["invalid_files_repo"] = filedict[st][
                    de["node"]
                ]["invalid_files_repo"].union(
                    set(
                        [
                            f
                            for i, f in enumerate(de["files"])
                            if not de["is_newer_repo"][i]
                        ]
                    )
                )

    for st in b.get_stage_ids():
        for nd in b.get_nodes_by_stage_id(st):
            filedict[st][nd]["n_observed"] = len(filedict[st][nd]["observed_files"])
            filedict[st][nd]["n_missing"] = len(filedict[st][nd]["missing_files"])
            filedict[st][nd]["n_invalid"] = len(filedict[st][nd]["invalid_files"])
            filedict[st][nd]["n_invalid_file"] = len(
                filedict[st][nd]["invalid_files_file"]
            )
            filedict[st][nd]["n_invalid_repo"] = len(
                filedict[st][nd]["invalid_files_repo"]
            )
            filedict[st][nd]["n"] = (
                filedict[st][nd]["n_observed"]
                + filedict[st][nd]["n_missing"]
                + filedict[st][nd]["n_invalid"]
            )
    return filedict


# benchmark = "tests/data/Benchmark_001.yaml"
# out_dir = "out"
# b = validate_benchmark(benchmark, out_dir, echo=False)
# modules_repo_timestamps = get_module_timestamps(b)
# exec_path_dict = get_exec_path_dict(b, modules_repo_timestamps)


def print_exec_path_dict(
    exec_path_dict: dict,
    stages: list,
    threshold_n_missing: Union[int, None] = None,
    full: bool = False,
) -> str:
    def dict_encode(d, full: bool = False) -> str:
        if full:
            if any([not v for v in d["is_newer_file"]]):
                return "F"
            if any([not v for v in d["is_newer_repo"]]):
                return "R"
        if any([not v for v in d["exists"]]):
            return "M"
        return "."

    if threshold_n_missing is None:
        threshold_n_missing = len(stages) + 1

    outls = []
    nmissls = []
    for i in range(len(exec_path_dict)):
        tmps1 = ""
        nmiss = 0
        for st in stages:
            if st in exec_path_dict[i].keys():
                tmpcode = dict_encode(exec_path_dict[i][st], full=full)
                if tmpcode != ".":
                    nmiss += 1
                tmps1 += tmpcode
            else:
                tmps1 += " "
        nmissls.append(nmiss)
        filename = exec_path_dict[i][st]["files"][0]
        filename = filename.replace("default/", "")
        for st in stages:
            filename = filename.replace(f"{st}/", "")

        tmps1 += "    " + filename
        outls.append(tmps1)
    outls_final = []
    for out, nmiss in zip(outls, nmissls):
        if nmiss >= threshold_n_missing:
            outls_final.append(out)
    return "\n".join(outls_final)


def prepare_status(benchmark: str, out_dir: str, return_all: bool = False):
    """
    Prepare the status dictionary for the benchmark.
    If return_all is True, also return the filedict.
    """
    # benchmark = "tests/data/Benchmark_001.yaml"
    # out_dir = "out"
    b = validate_benchmark(benchmark, out_dir, echo=False)

    stages = list(b.get_stage_ids())
    modules_repo_timestamps = get_module_timestamps(b)
    exec_path_dict = get_exec_path_dict(b, modules_repo_timestamps)
    filedict2 = get_filedict(b, exec_path_dict)

    status_dict = {}
    status_dict["name"] = b.get_benchmark_name()
    status_dict["version"] = b.get_benchmark_version()
    aggr_results_files = list(b.get_metric_collector_output_paths())
    status_dict["results"] = {
        "observed_files": [
            f
            for f in aggr_results_files
            if Path(f).is_file() and Path(f).stat().st_size > 0
        ],
        "missing_files": [
            f
            for f in aggr_results_files
            if not Path(f).is_file() or Path(f).stat().st_size == 0
        ],
    }
    status_dict["stages"] = {
        st: {
            "n": sum([filedict2[st][nd]["n"] for nd in filedict2[st].keys()]),
            "n_observed": sum(
                [filedict2[st][nd]["n_observed"] for nd in filedict2[st].keys()]
            ),
            "n_missing": sum(
                [filedict2[st][nd]["n_missing"] for nd in filedict2[st].keys()]
            ),
            "n_invalid": sum(
                [filedict2[st][nd]["n_invalid"] for nd in filedict2[st].keys()]
            ),
            "n_invalid_file": sum(
                [filedict2[st][nd]["n_invalid_file"] for nd in filedict2[st].keys()]
            ),
            "n_invalid_repo": sum(
                [filedict2[st][nd]["n_invalid_repo"] for nd in filedict2[st].keys()]
            ),
            "n_nodes": len(filedict2[st].keys()),
            "n_modules": len(set([nd.module_id for nd in filedict2[st].keys()])),
        }
        for st in stages
    }

    status_dict["total"] = {
        "n": sum([status_dict["stages"][st]["n"] for st in stages]),
        "n_observed": sum([status_dict["stages"][st]["n_observed"] for st in stages]),
        "n_missing": sum([status_dict["stages"][st]["n_missing"] for st in stages]),
        "n_nodes": sum([status_dict["stages"][st]["n_nodes"] for st in stages]),
        "n_modules": sum([status_dict["stages"][st]["n_modules"] for st in stages]),
    }
    if return_all:
        return status_dict, filedict2, exec_path_dict
    else:
        return status_dict
