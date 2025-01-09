from omni.workflow.snakemake.scripts import parse_performance

onsuccess:
    # shell('find out -name "*performance.txt" | sort | xargs head > performances.txt')
    write_combined_performance_file()
