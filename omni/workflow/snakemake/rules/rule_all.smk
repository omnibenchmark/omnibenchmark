import os
from pathlib import Path
from typing import List

from omni.workflow.snakemake.scripts.parse_performance import write_combined_performance_file


def create_all_rule(paths: List[str], aggregate_performance: bool = False):
    if not aggregate_performance:
        rule all:
            input: paths,
            # "out/data/D2/default/D2.data.ext",
            # "out/data/D2/default/process/P1/default/D2.txt.gz",
            # "out/data/D1/default/process/P2/default/methods/M2/default/D1.model.out.gz"
            # "out/data/D1/default/process/P2/default/methods/M2/default/m1/default/D1.results.txt"
    else:
        rule all:
            input: paths
            output:
                f"{benchmark.out_dir}/performances.tsv"
            run:
                result = glob_wildcards(str(benchmark.out_dir) + "/{path}/{dataset}_performance.txt")
                performances = expand(str(benchmark.out_dir) + "/{path}/{dataset}_performance.txt", path=result.path, dataset=result.dataset)

                output_dir = Path(str(os.path.commonpath(output)))
                if len(output) == 1:
                    output_dir = Path(os.path.dirname(output_dir))

                write_combined_performance_file(output_dir, performances)