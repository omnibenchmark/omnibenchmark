from pathlib import Path
from typing import Union
from omnibenchmark.benchmark.execution_path import (
    ExecutionPathSet,
)
from omnibenchmark.benchmark import BenchmarkExecution


def print_exec_path_dict(
    exec_path_dict: dict,
    stages: list,
    threshold_n_missing: Union[int, None] = None,
    full: bool = False,
    logs: bool = False,
) -> str:
    """
    Print the execution path dictionary in a formatted way.
    Args:
        exec_path_dict (dict): The execution path dictionary
        stages (list): List of stage ids
        threshold_n_missing (int or None): Threshold for number of missing stages to print
        full (bool): Whether to include full status information
        logs (bool): Whether to include log file paths
    """
    if threshold_n_missing is None:
        threshold_n_missing = len(stages) + 1

    outls = []
    logls = []
    nmissls = []
    for i in range(len(exec_path_dict)):
        outls2 = []
        tmps1, nmiss, is_failed_stage = exec_path_dict[i].dict_encode(
            full=full, cumulative=True, stages=stages
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

        if logs:
            tmps2 = ""
            if is_failed_stage is not None:
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
            if logs and len(logls) > 0 and len(logls[i]) > 0:
                tmpoutstr += logls[i]
            outls_final.append(tmpoutstr)
    return "\n".join(outls_final)


def prepare_status(
    benchmark: str,
    out_dir: str,
    return_all: bool = False,
    cache_dir: Path = Path(".snakemake") / "repos",
):
    """
    Prepare the status dictionary for the benchmark.
    If return_all is True, also return the filedict.
    """
    # benchmark = "tests/data/Benchmark_001.yaml"
    # out_dir = "out"
    b = BenchmarkExecution(Path(benchmark), Path(out_dir))
    eps = ExecutionPathSet(b, cache_dir=cache_dir)

    stages = eps.stages
    exec_path_dict = eps.exec_path_dict
    filedict2 = eps.get_filedict(cumulative=True)

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
            "n_empty": sum(
                [filedict2[st][nd]["n_empty"] for nd in filedict2[st].keys()]
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
        "n_empty": sum([status_dict["stages"][st]["n_empty"] for st in stages]),
        "n_missing": sum([status_dict["stages"][st]["n_missing"] for st in stages]),
        "n_nodes": sum([status_dict["stages"][st]["n_nodes"] for st in stages]),
        "n_modules": sum([status_dict["stages"][st]["n_modules"] for st in stages]),
    }
    if return_all:
        return status_dict, filedict2, exec_path_dict
    else:
        return status_dict
