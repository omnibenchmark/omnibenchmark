from omni.workflow.snakemake.scripts.parse_performance import write_combined_performance_file

onsuccess:
    write_combined_performance_file(benchmark.out_dir)
