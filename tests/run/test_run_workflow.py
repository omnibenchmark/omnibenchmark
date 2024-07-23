import argparse
import os

from omni.benchmark.validation import ValidationError
from omni.workflow.snakemake import SnakemakeEngine
from omni.benchmark import Benchmark


def run_benchmark(benchmark_file):
    benchmark = Benchmark(benchmark_file)
    converter = benchmark.get_converter()
    print(benchmark.get_definition())

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

    nodes = benchmark.get_nodes()
    print("All nodes:", nodes)

    execution_paths = benchmark.get_execution_paths()
    print("All execution paths:", execution_paths)

    outputs_paths = sorted(benchmark.get_output_paths())
    print("All output paths:", outputs_paths)

    benchmark.plot_benchmark_graph()

    # Serialize workflow to Snakefile
    workflow = SnakemakeEngine()
    workflow.serialize_workflow(benchmark)
    # workflow.run_node_workflow(nodes[2], input_dir="out/data/D1/default", dataset="D1")


def test_run_workflow():
    benchmark_file = "../data/Benchmark_001.yaml"
    benchmark_file_path = os.path.join(os.path.dirname(__file__), benchmark_file)

    if os.path.exists(benchmark_file_path):
        try:
            run_benchmark(benchmark_file_path)
        except ValidationError as e:
            print("Validation failed: \n")
            for error in e.errors:
                print(error)
    else:
        raise RuntimeError(f"Benchmark file {benchmark_file} does not exist.")
