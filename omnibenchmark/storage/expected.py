"""Pure computation of the output files a benchmark is expected to produce."""

from pathlib import Path
from typing import List

from omnibenchmark.core import Benchmark
from omnibenchmark.storage.base import StorageOptions


def get_expected_benchmark_output_files(
    benchmark: Benchmark,
    storage_options: StorageOptions,
) -> List:
    object_names_to_keep = benchmark.get_output_paths()
    if storage_options.extra_files_to_version_not_in_benchmark_yaml:
        for (
            glob_expression
        ) in storage_options.extra_files_to_version_not_in_benchmark_yaml:
            # pathlib, not glob: `glob` skips names starting with a dot, and
            # every node's parameter directory is dotted (`.default`, `.<hash>`).
            for found_file in Path().glob(glob_expression):
                object_names_to_keep.add(str(found_file))
    return list(object_names_to_keep)
