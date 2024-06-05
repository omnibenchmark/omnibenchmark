#!/usr/bin/env python

"""
Easybuild-powered software management, mostly by omniblock. Also includes singularity image generation.
"""

import subprocess
import os.path as op
# from easybuild.tools.module_naming_scheme.categorized_mns import det_full_module_name

def install_lmod():
    try:
        subprocess.call(["wget", "--help"])
    except FileNotFoundError:
        print("ERROR wget not found; please install it system-wise")
    try:
        subprocess.call(["make", "--help"])
    except FileNotFoundError:
        print("ERROR make (build tools) not found; please install them system-wise")    
    try:
        cmd = 'bash ./modules.sh'
        output = subprocess.check_output(
            cmd, stderr = subprocess.STDOUT, shell = True, timeout = 3,
            universal_newlines = True)
    except subprocess.CalledProcessError as exc:
        return("ERROR lmod install failed:", exc.returncode, exc.output)
    else:
        return("LOG lmod install output: \n{}\n".format(output))

def check_easybuild_status():
    try:
        subprocess.call(["eb", "--help"])
    except FileNotFoundError:
        return("ERROR easybuild not found")

def generate_default_easybuild_config_arguments(workdir, LMOD_VERS = "8.7"):
    HOME = op.expanduser("~")
    modulepath = op.join(workdir, 'easybuild', 'modules', 'all')
    buildpath = op.join(workdir, 'easybuild', 'build')
    containerpath = op.join(workdir, 'easybuild', 'containers')
    installpath = op.join(workdir, 'easybuild')
    repositorypath = op.join(workdir, 'easybuild', 'ebfiles_repo')
    # robotpath = op.join(workdir, 'easybuild', 'easyconfigs') ## let's use default's
    sourcepath = op.join(workdir, 'easybuild', 'sources')
    lmodpath=op.join(HOME, 'soft', 'lmod', LMOD_VERS, 'libexec'),
    lmodinit=op.join(HOME, 'soft', 'lmod', LMOD_VERS, 'init', 'bash'),
    lmodcmd=op.join(HOME, 'soft', 'lmod', LMOD_VERS, 'libexec', 'lmod')

    args = """--buildpath=%(buildpath)s  --installpath-modules %(modulepath)s \
              --containerpath %(containerpath)s --installpath %(installpath)s \
              --repositorypath %(repositorypath)s --sourcepath %(sourcepath)s""" %{
                  buildpath : buildpath,
                  modulepath : modulepath,
                  containerpath : containerpath,
                  installpath : installpath,
                  repositorypath : repositorypath,
                  sourcepath : sourcepath}
    return(args)
        
def easybuild_easyconfig(easyconfig,
                         workdir,
                         threads,
                         containerize = False,
                         container_build_image = False):
    """
    easybuilds an easyconfig, potentially generating a (built) container image too
    """

    args = generate_default_easybuild_config_arguments(workdir = workdir)
    
    cmd = """eb %(easyconfig)s --robot %(args)s  --parallel= %(threads)s \
              --detect-loaded-modules=unload --check-ebroot-env-vars=unset""" %{
                  easyconfig : easyconfig,
                  args : args,
                  threads : threads}

    if containerize or container_build_image:
        cmd += " --container-config bootstrap=localimage,from=example.sif --experimental"
        if container_build_image:
            cmd += " --container-build-image"
    try:
        output = subprocess.check_output(
            cmd, stderr = subprocess.STDOUT, shell = True, timeout = 3,
            universal_newlines = True)
    except subprocess.CalledProcessError as exc:
        return("ERROR easybuild failed:", exc.returncode, exc.output)
    else:
        return("LOG easybuild: \n{}\n".format(output))

def list_all_easyconfigs(benchmark_yaml):
    return('not_implemented')

def easybuild_all_easyconfigs(benchmark_yaml, workdir, threads, containerize = False,
                              container_build_image = False):
    """
    Iterates over all easyconfigs from a benchmark yaml and easybuilds them all
    """

    easyconfigs = list_all_easyconfigs()
    for eb in easyconfigs:
        easybuild_easyconfig(easyconfig = eb, workdir = workdir, threads = threads,
                             containerize = containerize, container_build_image = container_build_image)

def get_easyconfig_from_envmodule_name(envmodule):
    return('not implemented')

def get_envmodule_name_from_easyconfig(easyconfig):
    return('not implemented')

def check_envmodule_status(benchmark_yaml, stage_id, module_id):
    """
    Checks whether the envmodule is present for a given benchmarking node
    """
    return('not implemented')
