#!/usr/bin/env python

"""
Conda-based software management, per benchmark yaml (and not per omniblock)
"""

def pin_conda_envs(benchmark_yaml):
    """
    If not benchmarking run mode: error
    """
    return('not implemented')

    
def check_node_conda(benchmark_yaml, stage_id, module_id):
    """
    Check if the env.yaml and the pin.txt for that conda env file exist locally.
    If so: valid, continue (Snakemake will install, if applicable).
    If not:
       - env.yaml not present: error
       - env.yaml present but pin.txt not present: error
    """
    return('not implemented')
