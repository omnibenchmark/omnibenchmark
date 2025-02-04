from omni.workflow.snakemake.scripts.parse_performance import write_combined_performance_file

onsuccess:
    print(benchmark.out_dir)
