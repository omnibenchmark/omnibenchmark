"""Pure computation of the output files a benchmark is expected to produce."""

from pathlib import Path
from typing import List

from omnibenchmark.core import Benchmark
from omnibenchmark.core._human_view import HUMAN_DIR
from omnibenchmark.storage.base import StorageOptions


def get_expected_benchmark_output_files(
    benchmark: Benchmark,
    storage_options: StorageOptions,
) -> List:
    object_names_to_keep = benchmark.get_output_paths()
    human_dirs = [Path(d) / HUMAN_DIR for d in storage_options.results_directories]
    if storage_options.extra_files_to_version_not_in_benchmark_yaml:
        for (
            glob_expression
        ) in storage_options.extra_files_to_version_not_in_benchmark_yaml:
            # pathlib, not glob: `glob` skips names starting with a dot, and
            # every node's parameter directory is dotted (`.default`, `.<hash>`).
            for found_file in Path().glob(glob_expression):
                # `human/` links to the results; version them once (007 §3.1.2).
                if any(found_file.is_relative_to(h) for h in human_dirs):
                    continue
                object_names_to_keep.add(str(found_file))
    return list(object_names_to_keep)
