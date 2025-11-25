import logging
import os
import shutil
import stat
from pathlib import Path

from omnibenchmark.benchmark import Benchmark
from omnibenchmark.workflow.snakemake import SnakemakeEngine

# TODO: no need to cleanup anything, just use tmp_dir
TO_CLEANUP = [".snakemake", "out", "Snakefile", "snakemake.log"]


def raise_if_file_not_found(file_path: str):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} does not exist.")


class SnakemakeSetup:
    def __init__(
        self,
        benchmark_file: Path,
        keep_files: bool = False,
        cwd: str = "",
        out_dir: Path = Path("out"),
    ) -> None:
        raise_if_file_not_found(benchmark_file.as_posix())
        self.benchmark_file = benchmark_file
        self.benchmark = Benchmark(benchmark_file, out_dir=out_dir)
        self.workflow = SnakemakeEngine()
        self.keep_files = keep_files
        self.cwd = cwd
        self._print_benchmark()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._cleanup_snakemake()

    def _cleanup_snakemake(self):
        if self.keep_files:
            logging.getLogger().warning("File cleanup skipped")
            return
        for file in TO_CLEANUP:
            file_path = os.path.join(self.cwd, file)
            if os.path.exists(file_path):
                if os.path.isfile(file_path):
                    os.remove(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path, onerror=self._remove_readonly)

    def _remove_readonly(self, func, path, _):
        "Clear the readonly bit and reattempt the removal"
        os.chmod(path, stat.S_IWRITE)
        func(path)

    def _print_benchmark(self):
        benchmark = self.benchmark.model
        print(self.benchmark.get_definition())

        stages = benchmark.get_stages()
        for stage_id in stages:
            stage = stages[stage_id]
            stage_name = stage.name
            print("Stage", stage_name if stage_name else stage_id)

            modules_in_stage = benchmark.get_modules_by_stage(stage)
            print("  ", stage_name, "with modules", modules_in_stage.keys(), "\n")
            print("  Implicit inputs:\n", benchmark.get_stage_implicit_inputs(stage))
            print(
                "  Explicit inputs:\n",
                [
                    benchmark.get_explicit_inputs(i)
                    for i in benchmark.get_stage_implicit_inputs(stage)
                ],
            )
            print("  Outputs:", benchmark.get_stage_outputs(stage))
            print("  All modules:", benchmark.get_modules().keys())
            print("------")

            for module_id in modules_in_stage:
                module = modules_in_stage[module_id]
                module_name = module.name
                print("  Module", module_name if module_name else module_id)
                print("    Repo:", benchmark.get_module_repository(module))
                print("    Excludes:", benchmark.get_module_excludes(module))
                print("    Params:", benchmark.get_module_parameters(module))
            print("------")

        print("------")

        nodes = self.benchmark.get_nodes()
        print("All nodes:", nodes)

        execution_paths = self.benchmark.get_execution_paths()
        print("All execution paths:", execution_paths)

        outputs_paths = sorted(self.benchmark.get_output_paths())
        print("All output paths:", outputs_paths)
