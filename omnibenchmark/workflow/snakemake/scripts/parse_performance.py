"""
Gathers performance benchmark.txt files and documents the parameters used for each stage/method, if applicable
"""

import csv
import glob
import os.path
from pathlib import Path
from typing import List

import pandas
import os.path as op


def write_combined_performance_file(out_dir: Path, performance_files: List[str]):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    fd = combine_performances(out_dir, performance_files)
    fd.to_csv(out_dir / "performances.tsv", sep="\t", index=False)


def combine_performances(
    out_dir: Path, performance_files: List[str]
) -> pandas.DataFrame:
    fd = pandas.DataFrame()

    for perf in performance_files:
        if os.path.exists(perf):
            curr = read_performance(perf)
            temp_df = pandas.DataFrame(curr, index=[1])
            temp_df["module"] = op.dirname(perf).split("/")[-2]
            temp_df["path"] = perf
            temp_df["params"] = read_params(out_dir, perf)
            fd = pandas.concat([fd, temp_df], ignore_index=True, axis=0)

    return fd


def read_performance(file_path: str):
    with open(file_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for record in reader:
            record.pop("h:m:s", None)
            for k, v in record.items():
                # Check for NA or invalid values and handle them
                if v in ["NA", "none", "", "null"]:
                    record[k] = 0
                else:
                    try:
                        record[k] = float(v)
                    except ValueError:
                        record[k] = 0
            yield record


def tokenize(output_path: Path, file_path: str):
    ## we get only after the 'out' directory
    # TODO(ben): be more careful here
    try:
        fp = file_path.split(f"{output_path}/")[1].split("/")
    except IndexError:
        return []
    ## and slice in stage/method/params triples
    return [x for x in zip(*(iter(fp),) * 3)]


def read_params(output_path: Path, file_path: str):
    triples = tokenize(output_path, file_path)
    res = ""
    parent = output_path
    for triple in triples:
        parent = parent / triple[0] / triple[1] / triple[2]
        if "default" not in triple[2]:
            param_file_path = parent / "parameters.json"
            with open(param_file_path) as fh:
                params = fh.read()
                res = "%s %s %s %s %s;" % (
                    res,
                    triple[0],
                    triple[1],
                    triple[2],
                    params.strip(),
                )

    return res


if __name__ == "__main__":
    files = glob.glob("**/*_performance.txt")
    write_combined_performance_file(Path("out"), files)
