from datetime import datetime
import os
from pathlib import Path
from typing import Union
from omnibenchmark.benchmark import Benchmark, BenchmarkNode
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


def n2f(bl: Union[bool, None]) -> bool:
    return bl is True


def get_repo_timestamp(repo: dict) -> Union[float, None]:
    repositories_dir = Path(".snakemake") / "repos"
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
    return timestamp


def get_module_timestamps(b: Benchmark) -> dict:
    """
    Get the timestamps of the module repositories used in the benchmark.
    Returns a dictionary with stage ids as keys and another dictionary as value,
    which has module ids as keys and timestamps as values.
    """

    modules_repo_timestamps = {}
    for st in b.get_stage_ids():
        modules_repo_timestamps[st] = {}
        for node in b.get_nodes_by_stage_id(st):
            if node.module_id not in modules_repo_timestamps[st].keys():
                timestamp = get_repo_timestamp(node.get_repository())
                modules_repo_timestamps[st][node.module_id] = timestamp
    return modules_repo_timestamps


class ExecutionPathStageFile:
    def __init__(
        self,
        input_files: list[str],
        repo_timestamp: Union[float, None],
        output_file: str,
    ):
        self.output_file = Path(output_file)
        self.output_file_exist = (
            self.output_file.is_file() and self.output_file.stat().st_size > 0
        )
        self.timestamp = (
            self.output_file.stat().st_mtime if self.output_file_exist else None
        )

        self.input_files = [Path(f) for f in input_files]
        self.input_files_exist = [
            f.is_file() and f.stat().st_size > 0 for f in self.input_files
        ]
        self.all_input_files_exist = all(self.input_files_exist)
        self.input_files_timestamps = [
            f.stat().st_mtime if self.input_files_exist[i] else None
            for i, f in enumerate(self.input_files)
        ]
        self.max_input_file_timestamp = (
            datetime.now().timestamp() * 2 + 1e6
            if any([tst is None for tst in self.input_files_timestamps])
            and len(self.input_files_timestamps) > 0
            else max([tst for tst in self.input_files_timestamps if tst is not None])
            if len(self.input_files_timestamps) > 0
            else 0
        )

        self.input_file_is_newer = (
            self.max_input_file_timestamp > self.timestamp
            if self.timestamp is not None
            else None
        )

        self.repo_timestamp = repo_timestamp
        self.repo_is_newer = (
            self.repo_timestamp > self.timestamp
            if self.timestamp is not None and self.repo_timestamp is not None
            else None
        )

        self.any_is_newer = self.input_file_is_newer or self.repo_is_newer

        self.any_is_newer_cumulative: Union[bool, None] = None
        self.repo_is_newer_cumulative: Union[bool, None] = None
        self.input_file_is_newer_cumulative: Union[bool, None] = None
        self.input_file_exists_cumulative: Union[bool, None] = None

    def output_is_valid(
        self, type: str = "all", cumulative: bool = False
    ) -> Union[bool, None]:
        if type == "all":
            is_invalid = (
                self.any_is_newer_cumulative if cumulative else self.any_is_newer
            )
        elif type == "input_files_newer":
            is_invalid = (
                self.input_file_is_newer_cumulative
                if cumulative
                else self.input_file_is_newer
            )
        elif type == "input_files_exist":
            is_invalid = (
                self.all_input_files_exist_cumulative
                if cumulative
                else self.all_input_files_exist
            )
        elif type == "repo":
            is_invalid = (
                self.repo_is_newer_cumulative if cumulative else self.repo_is_newer
            )
        else:
            raise ValueError("type must be either 'all', 'input_file' or 'repo'")
        return None if is_invalid is None else not is_invalid

    def output_exists(self) -> bool:
        return self.output_file_exist

    def get_output_file(
        self, filter: str = "none", cumulative: bool = False
    ) -> Union[str, None]:
        if filter == "none":
            do_return = True
        elif filter == "observed":
            do_return = n2f(self.output_exists()) and n2f(
                self.output_is_valid(type="all", cumulative=cumulative)
            )
        elif filter == "missing":
            do_return = not n2f(self.output_exists())
        elif filter == "invalid":
            do_return = n2f(self.output_exists()) and not n2f(
                self.output_is_valid(type="all", cumulative=cumulative)
            )
        elif filter == "input_files_newer":
            do_return = n2f(self.output_exists()) and not n2f(
                self.output_is_valid(type="input_files_newer", cumulative=cumulative)
            )
        elif filter == "repo":
            do_return = n2f(self.output_exists()) and not n2f(
                self.output_is_valid(type="repo", cumulative=cumulative)
            )
        else:
            raise ValueError(
                "filter must be either 'none', 'observed', 'missing', 'invalid', 'input_files_newer' or 'repo'"
            )
        if do_return:
            return str(self.output_file)
        else:
            return None

    def set_cumulative(self, type: str, val: bool):
        if type == "all":
            self.any_is_newer_cumulative = val
        elif type == "input_files_newer":
            self.input_file_is_newer_cumulative = val
        elif type == "input_files_exist":
            self.all_input_files_exist_cumulative = val
        elif type == "repo":
            self.repo_is_newer_cumulative = val
        else:
            raise ValueError(
                "type must be either 'all', 'input_files_newer', 'input_files_exist' or 'repo'"
            )


class ExecutionPathStage:
    def __init__(
        self,
        node: BenchmarkNode,
        input_dir: str = "",
        repo_timestamp: Union[float, None] = None,
    ):
        self.node = node
        self.input_dir = input_dir
        self.stage_id = node.stage_id
        self.module_id = node.module_id
        self.param_id = node.param_id
        self.parameters = node.parameters

        if repo_timestamp is None:
            self.repo_timestamp = get_repo_timestamp(node.get_repository())
        else:
            self.repo_timestamp = repo_timestamp

        config = {
            "input": input_dir,
            "output": "",
            "dataset": self.module_id,
        }
        output_files = node.get_output_paths(config)

        input_files = []
        self.output_stage_files = [
            ExecutionPathStageFile(
                input_files=input_files,
                repo_timestamp=self.repo_timestamp,
                output_file=output_file,
            )
            for output_file in output_files
        ]

        self.output_dir = os.path.commonpath(
            list(set([os.path.dirname(f.output_file) for f in self.output_stage_files]))
        )

        stdout_log = Path(self.output_dir) / "stdout.log"
        self.stdout_log = stdout_log if stdout_log.is_file() else None
        stderr_log = Path(self.output_dir) / "stderr.log"
        self.stderr_log = stderr_log if stderr_log.is_file() else None

    def output_is_valid(
        self, type: str = "all", cumulative: bool = False, aggregate: str = "none"
    ) -> list[Union[bool, None]]:
        is_valid_list = [
            f.output_is_valid(type=type, cumulative=cumulative)
            for f in self.output_stage_files
        ]
        if aggregate == "none":
            return is_valid_list
        else:
            if any([v is None for v in is_valid_list]):
                return [None]
            if aggregate == "any":
                return [any(is_valid_list)]
            elif aggregate == "all":
                return [all(is_valid_list)]
            else:
                raise ValueError("aggregate must be either 'none', 'any' or 'all'")

    def output_exists(self, aggregate: str = "none") -> list[bool]:
        if aggregate == "none":
            return [f.output_exists() for f in self.output_stage_files]
        elif aggregate == "any":
            return [any(f.output_exists() for f in self.output_stage_files)]
        elif aggregate == "all":
            return [all(f.output_exists() for f in self.output_stage_files)]
        else:
            raise ValueError("aggregate must be either 'none', 'any' or 'all'")

    def set_cumulative(self, type: str, val: bool):
        for f in self.output_stage_files:
            f.set_cumulative(type=type, val=val)

    def get_output_files(
        self, filter: str = "none", cumulative: bool = False, remove_none: bool = True
    ) -> list[str] | list[Union[str, None]]:
        files = [
            f.get_output_file(filter=filter, cumulative=cumulative)
            for f in self.output_stage_files
        ]
        if remove_none:
            return [f for f in files if f is not None]
        else:
            return files

    def dict_encode(self, full: bool = False, cumulative: bool = False) -> str:
        output_exists = n2f(self.output_exists(aggregate="all")[0])
        if output_exists and self.output_is_valid(type="all", aggregate="all"):
            return "."
        if (
            not full
            and output_exists
            and not n2f(
                self.output_is_valid(
                    type="all", cumulative=cumulative, aggregate="all"
                )[0]
            )
        ):
            return "I"
        if (
            full
            and output_exists
            and not n2f(
                self.output_is_valid(
                    type="all", cumulative=cumulative, aggregate="all"
                )[0]
            )
        ):
            if not n2f(
                self.output_is_valid(type="input_files_newer", aggregate="all")[0]
            ):
                return "F"
            if not n2f(
                self.output_is_valid(type="input_files_newer", aggregate="all")[0]
            ):
                return "R"
        if not output_exists:
            return "M"
        else:
            return "?"


class ExecutionPath:
    def __init__(
        self, nodes: list[BenchmarkNode], repo_timestamps: Union[dict, None] = None
    ):
        self.nodes = nodes
        self.stages = [node.stage_id for node in nodes]

        # nodes = b.get_execution_paths()[0]
        self.exec_path = {}
        input_dir = ""
        for i, st in enumerate(self.stages):
            self.exec_path[st] = ExecutionPathStage(
                nodes[i],
                input_dir=input_dir,
                repo_timestamp=repo_timestamps[st] if repo_timestamps else None,
            )
            input_dir = self.exec_path[st].output_dir
        self.set_cumulative()

    def set_cumulative(self):
        for i in range(len(self.stages)):
            if i == 0:
                for t in ["all", "input_files_newer", "input_files_exist", "repo"]:
                    is_valid = self.exec_path[self.stages[i]].output_is_valid(
                        type=t, cumulative=False
                    )
                    is_valid_all = all(is_valid)
                    self.exec_path[self.stages[i]].set_cumulative(t, not is_valid_all)
            else:
                for t in ["all", "input_files_newer", "input_files_exist", "repo"]:
                    is_valid = self.exec_path[self.stages[i]].output_is_valid(
                        type=t, cumulative=False
                    )
                    is_valid_all = all(is_valid)
                    is_valid_cumulative = self.exec_path[
                        self.stages[i - 1]
                    ].output_is_valid(type=t, cumulative=True)
                    is_valid_cumulative_all = all(is_valid_cumulative)
                    self.exec_path[self.stages[i]].set_cumulative(
                        t, not (is_valid_all and is_valid_cumulative_all)
                    )

    def dict_encode(
        self, full: bool = False, stages: Union[list[str], None] = None
    ) -> tuple[str, int, Union[str, None]]:
        if stages is None:
            stages = self.stages
        tmps1 = ""
        nmiss = 0
        is_failed_stage = None
        for st in stages:
            if st in self.stages:
                tmpcode = self.exec_path[st].dict_encode(full=full)
                if tmpcode != ".":
                    if is_failed_stage is None:
                        is_failed_stage = st
                    nmiss += 1
                tmps1 += tmpcode
            else:
                tmps1 += " "
        return tmps1, nmiss, is_failed_stage


# repo_timestamps_all = get_module_timestamps(b)
# nodes = b.get_execution_paths()[0]
# repo_timestamps = {node.stage_id: repo_timestamps_all[node.stage_id][node.module_id] for node in nodes}
# ep = ExecutionPath(nodes=nodes, repo_timestamps = repo_timestamps)
# ep.exec_path['data'].output_is_valid(type="all", cumulative=False)
# ep.exec_path['process'].output_is_valid(type="all", cumulative=False)
# ep.exec_path['methods'].output_is_valid(type="all", cumulative=False)
# ep.exec_path['metrics'].output_is_valid(type="all", cumulative=False)


def get_exec_path_dict(b: Benchmark, modules_repo_timestamps: dict) -> dict:
    exec_paths = b.get_execution_paths()
    modules_repo_timestamps = get_module_timestamps(b)
    exec_path_dict = {}
    for exec_path_id in range(len(exec_paths)):
        nodes = b.get_execution_paths()[exec_path_id]
        repo_timestamps = {
            node.stage_id: modules_repo_timestamps[node.stage_id][node.module_id]
            for node in nodes
        }
        exec_path_dict[exec_path_id] = ExecutionPath(
            nodes=nodes, repo_timestamps=repo_timestamps
        )
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
                    "observed_output_files": set(),
                    "missing_output_files": set(),
                    "invalid_output_files_input_file_is_newer": set(),
                    "invalid_output_files_repo_is_newer": set(),
                    "invalid_output_files": set(),
                }
                for nd in b.get_nodes_by_stage_id(st)
            }
            for st in b.get_stage_ids()
        }
    }
    for st in b.get_stage_ids():
        for exec_path_id in exec_path_dict.keys():
            if st in exec_path_dict[exec_path_id].stages:
                de = exec_path_dict[exec_path_id].exec_path[st]
                filedict[st][de.node]["observed_output_files"] = filedict[st][de.node][
                    "observed_output_files"
                ].union(set(de.get_output_files("observed")))

                filedict[st][de.node]["missing_output_files"] = filedict[st][de.node][
                    "missing_output_files"
                ].union(set(de.get_output_files("missing")))
                filedict[st][de.node]["invalid_output_files"] = filedict[st][de.node][
                    "invalid_output_files"
                ].union(set(de.get_output_files("invalid")))
                filedict[st][de.node]["invalid_output_files_input_file_is_newer"] = (
                    filedict[
                        st
                    ][
                        de.node
                    ][
                        "invalid_output_files_input_file_is_newer"
                    ].union(set(de.get_output_files("input_files_newer")))
                )
                filedict[st][de.node]["invalid_output_files_repo_is_newer"] = filedict[
                    st
                ][de.node]["invalid_output_files_repo_is_newer"].union(
                    set(de.get_output_files("repo"))
                )

    for st in b.get_stage_ids():
        for nd in b.get_nodes_by_stage_id(st):
            filedict[st][nd]["n_observed"] = len(
                filedict[st][nd]["observed_output_files"]
            )
            filedict[st][nd]["n_missing"] = len(
                filedict[st][nd]["missing_output_files"]
            )
            filedict[st][nd]["n_invalid"] = len(
                filedict[st][nd]["invalid_output_files"]
            )
            filedict[st][nd]["n_invalid_input_file_is_newer"] = len(
                filedict[st][nd]["invalid_output_files_input_file_is_newer"]
            )
            filedict[st][nd]["n_invalid_repo_is_newer"] = len(
                filedict[st][nd]["invalid_output_files_repo_is_newer"]
            )
            filedict[st][nd]["n"] = (
                filedict[st][nd]["n_observed"]
                + filedict[st][nd]["n_missing"]
                + filedict[st][nd]["n_invalid"]
            )
    return filedict


benchmark = "tests/data/Benchmark_001.yaml"
out_dir = "out"
b = validate_benchmark(benchmark, out_dir, echo=False)
modules_repo_timestamps = get_module_timestamps(b)
exec_path_dict = get_exec_path_dict(b, modules_repo_timestamps)


def print_exec_path_dict(
    exec_path_dict: dict,
    stages: list,
    threshold_n_missing: Union[int, None] = None,
    full: bool = False,
    logs: bool = False,
) -> str:
    if threshold_n_missing is None:
        threshold_n_missing = len(stages) + 1

    outls = []
    logls = []
    nmissls = []
    for i in range(len(exec_path_dict)):
        outls2 = []
        tmps1, nmiss, is_failed_stage = exec_path_dict[i].dict_encode(
            full=full, stages=stages
        )
        outls2.append(tmps1)
        nmissls.append(nmiss)
        for st in stages:
            if st in exec_path_dict[i].stages:
                node = exec_path_dict[i].exec_path[st].node
                params_name = str(node.parameters).replace("'", "").replace('"', "")
                module_id = node.module_id
                outls2.append(
                    f"{module_id}{f"-{params_name}" if params_name != "None" else ""}"
                )

            else:
                outls2.append("")
        outls2.append(Path(exec_path_dict[i].exec_path[st].get_output_files()[0]).name)

        if logs and is_failed_stage is not None:
            tmps2 = ""
            if exec_path_dict[i].exec_path[st].stdout_log is not None:
                tmps2 += f"\n        STDOUT: {exec_path_dict[i].exec_path[st].stdout_log.as_posix()}"
            if exec_path_dict[i].exec_path[st].stderr_log is not None:
                tmps2 += f"\n        STDERR: {exec_path_dict[i].exec_path[st].stderr_log.as_posix()}"
            logls.append(tmps2)
        outls.append(outls2)

    tmplens = [
        max(
            [len(outls[i][j]) for i in range(len(outls))]
            + [0 if j < 1 or j > len(stages) else len(stages[j - 1])]
        )
        for j in range(len(stages) + 2)
    ]
    outls3 = []
    for out in outls:
        outls3.append(
            "   ".join([out[j].ljust(tmplens[j]) for j in range(len(stages) + 2)])
        )
    outls_final = []
    ascii_stages_header = "\n".join(
        [f"{'|' * (i + 1)}{st}" for i, st in enumerate(stages)]
    )
    ascii_stages_sep1 = "|" * len(stages)
    ascii_stages_sep2 = (
        "|" * len(stages)
        + "   "
        + "   ".join(
            [stages[j - 1].ljust(tmplens[j]) for j in range(1, len(stages) + 1)]
        )
    )
    outls_final.append(ascii_stages_header)
    outls_final.append(ascii_stages_sep1)
    outls_final.append(ascii_stages_sep2)
    for i, (out, nmiss) in enumerate(zip(outls3, nmissls)):
        if nmiss >= threshold_n_missing:
            tmpoutstr = out
            if logs and len(logls[i]) > 0:
                tmpoutstr += logls[i]
            outls_final.append(tmpoutstr)
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
        "observed_output_files": [
            f
            for f in aggr_results_files
            if Path(f).is_file() and Path(f).stat().st_size > 0
        ],
        "missing_output_files": [
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
            "n_invalid_input_file_is_newer": sum(
                [
                    filedict2[st][nd]["n_invalid_input_file_is_newer"]
                    for nd in filedict2[st].keys()
                ]
            ),
            "n_invalid_repo_is_newer": sum(
                [
                    filedict2[st][nd]["n_invalid_repo_is_newer"]
                    for nd in filedict2[st].keys()
                ]
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
