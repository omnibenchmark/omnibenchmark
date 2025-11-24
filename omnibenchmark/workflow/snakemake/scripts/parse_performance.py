"""
Gathers performance benchmark.txt files and documents the parameters used for each stage/method, if applicable
"""

import csv
import glob
import os.path
from pathlib import Path
from typing import List, Dict, Any
import os.path as op


def write_combined_performance_file(out_dir: Path, performance_files: List[str]):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    rows = combine_performances(out_dir, performance_files)

    # Write to TSV file
    if rows:
        output_file = out_dir / "performances.tsv"
        with open(output_file, "w", newline="") as f:
            # Get all unique column names from all rows
            fieldnames = set()
            for row in rows:
                fieldnames.update(row.keys())
            fieldnames = sorted(fieldnames)  # Sort for consistent output

            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)


def combine_performances(
    out_dir: Path, performance_files: List[str]
) -> List[Dict[str, Any]]:
    rows = []

    for perf in performance_files:
        if os.path.exists(perf):
            for record in read_performance(perf):
                record["module"] = op.dirname(perf).split("/")[-2]
                record["path"] = perf
                record["params"] = read_params(out_dir, perf)
                rows.append(record)

    return rows


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
