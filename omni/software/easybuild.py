#!/usr/bin/env python

"""
Easybuild-powered software management, mostly by omniblock. Also includes singularity image generation.
"""

def install_lmod():
    return('not implemented')

def check_easybuild_status():
    return('not implemented')

def easybuild_easyconfig(eb,
                   modulepath,
                   buildpath,
                   containerpath,
                   installpath,
                   repositorypath,
                   sourcepath,
                   robotpath,
                   lmodpath,
                   lmodinit,
                   lmodcmd,
                   containerize = FALSE,
                   container_build_image = FALSE):
    """
    easybuilds an easyconfig, potentially generating a (built) container image too
    """
    
    return('not implemented')

def easybuild_all_easyconfigs(benchmark_yaml):
    """
    Iterates over all easyconfigs from a benchmark yaml and easybuilds them all
    """
    return('not implemented')

def get_easyconfig_from_envmodule_name(envmodule):
    return('not implemented')

def get_envmodule_name_from_easyconfig(easyconfig):
    return('not implemented')

def check_envmodule_status(benchmark_yaml, stage_id, module_id):
    """
    Checks whether the envmodule is present for a given benchmarking node
    """
    return('not implemented')
