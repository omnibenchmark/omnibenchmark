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
    """Write combined performance file from all performance files.

    Args:
        out_dir: Output directory to write the combined file
        performance_files: List of performance file paths to combine
    """
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    rows = combine_performances(out_dir, performance_files)

    # Write to TSV file
    output_file = out_dir / "performances.tsv"
    if rows:
        fieldnames = list(rows[0].keys())
        with open(output_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    else:
        # Create empty file to satisfy Snakemake output requirements
        with open(output_file, "w", newline="") as f:
            pass


def combine_performances(
    out_dir: Path, performance_files: List[str]
) -> List[Dict[str, Any]]:
    """Combine performance files into a list of dictionaries.

    Args:
        out_dir: Output directory path
        performance_files: List of performance file paths

    Returns:
        List of dictionaries containing combined performance data
    """
    combined_rows = []

    for perf in performance_files:
        if os.path.exists(perf):
            # read_performance is a generator, get first (and only) record
            for record in read_performance(perf):
                split_by_slash = op.dirname(perf).split("/")
                record["dataset"] = split_by_slash[2]
                record["module"] = split_by_slash[-2]
                record["lineage"] = "/".join(split_by_slash[3:-2])
                record["path"] = perf
                record["params"] = read_params(out_dir, perf)
                combined_rows.append(record)
                break  # Only process first record from each file

    return combined_rows


def read_performance(file_path: str):
    """Read performance data from a TSV file.

    Args:
        file_path: Path to performance file

    Yields:
        Dictionary of performance data
    """
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
    """Tokenize file path into stage/method/params triples.

    Args:
        output_path: Base output directory
        file_path: Full file path

    Returns:
        List of (stage, method, params) tuples
    """
    ## we get only after the 'out' directory
    # TODO(ben): catch the IndexError and log a more informative error message
    try:
        fp = file_path.split(f"{output_path}/")[1].split("/")
    except IndexError:
        return []
    ## and slice in stage/method/params triples
    return [x for x in zip(*(iter(fp),) * 3)]


def read_params(output_path: Path, file_path: str):
    """Read parameters from parameter files in the path.

    Args:
        output_path: Base output directory
        file_path: Full file path

    Returns:
        String containing parameter information
    """
    triples = tokenize(output_path, file_path)
    res = ""
    parent = output_path
    for triple in triples:
        parent = parent / triple[0] / triple[1] / triple[2]
        if "default" not in triple[2]:
            param_file_path = parent / "parameters.json"
            if param_file_path.exists():
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
    files = glob.glob("**/*_performance.txt", recursive=True)
    write_combined_performance_file(Path("out"), files)
