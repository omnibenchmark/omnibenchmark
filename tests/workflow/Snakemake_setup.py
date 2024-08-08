import os
import shutil
import stat
from pathlib import Path

from omni.benchmark import Benchmark
from omni.workflow.snakemake import SnakemakeEngine

TO_CLEANUP = [".snakemake", "out", "Snakefile", "snakemake.log"]


class SnakemakeSetup:
    def __init__(self, benchmark_file: Path) -> None:
        self.benchmark_file = benchmark_file

        if os.path.exists(benchmark_file):
            self.benchmark = Benchmark(benchmark_file)
            self.workflow = SnakemakeEngine()
            self._print_benchmark()
        else:
            raise FileNotFoundError(f"Benchmark file {benchmark_file} does not exist.")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._cleanup_snakemake()

    def _cleanup_snakemake(self):
        current_dir = os.getcwd()
        for file in TO_CLEANUP:
            file_path = os.path.join(current_dir, file)
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
        converter = self.benchmark.get_converter()
        print(self.benchmark.get_definition())

        stages = converter.get_stages()
        for stage_id in stages:
            stage = stages[stage_id]
            stage_name = stage["name"]
            print("Stage", stage_name if stage_name else stage_id)

            modules_in_stage = converter.get_modules_by_stage(stage)
            print("  ", stage_name, "with modules", modules_in_stage.keys(), "\n")
            print("  Implicit inputs:\n", converter.get_stage_implicit_inputs(stage))
            print(
                "  Explicit inputs:\n",
                [
                    converter.get_explicit_inputs(i)
                    for i in converter.get_stage_implicit_inputs(stage)
                ],
            )
            print("  Outputs\n", converter.get_stage_outputs(stage))
            print("------")

            for module_id in modules_in_stage:
                module = modules_in_stage[module_id]
                module_name = module["name"]
                print("  Module", module_name if module_name else module_id)
                print("    Repo:", converter.get_module_repository(module))
                print("    Excludes:", converter.get_module_excludes(module))
                print("    Params:", converter.get_module_parameters(module))
            print("------")

        print("------")

        nodes = self.benchmark.get_nodes()
        print("All nodes:", nodes)

        execution_paths = self.benchmark.get_execution_paths()
        print("All execution paths:", execution_paths)

        outputs_paths = sorted(self.benchmark.get_output_paths())
        print("All output paths:", outputs_paths)
