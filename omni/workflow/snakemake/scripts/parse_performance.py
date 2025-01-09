"""Found and parse snakemake benchmark file"""

# std import
import csv
import os
import re
import typing
import pandas
import glob


def write_combined_performance_file():
    fd = combine_performances()
    fd.to_csv("performances.tsv", sep="\t", index=False)


def combine_performances() -> pandas.DataFrame:
    fd = pandas.DataFrame()

    perfs = glob.glob("**/**performance.txt", recursive=True)

    for perf in perfs:
        curr = read_performance(perf)
        temp_df = pandas.DataFrame(curr, index=[1])
        temp_df["path"] = perf
        fd = pandas.concat([fd, temp_df], ignore_index=True, axis=0)

    return fd


def read_performance(file_path: str):
    with open(file_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for record in reader:
            record.pop("h:m:s", None)
            record = {k: float(v) for k, v in record.items()}
            yield record


if __name__ == "__main__":
    write_combined_performance_file()
