#!/usr/bin/env python

"""
Easybuild-powered software management, mostly by omniblock. Also includes singularity image generation.

Izaskun Mallona
Started 5th June 2024
"""

import subprocess
import os.path as op
import sys
from easybuild.tools.module_naming_scheme.mns import ModuleNamingScheme
from easybuild.tools.module_naming_scheme.utilities import det_full_ec_version, is_valid_module_name
from easybuild.framework.easyconfig.tools import det_easyconfig_paths, parse_easyconfigs
from easybuild.tools.options import set_up_configuration
from easybuild.tools.modules import get_software_root_env_var_name, modules_tool
# from warnings import deprecated

## shell-based stuff, partly to be replaced by direct eb API calls -------------------------------------

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

def export_lmod_env_vars(LMOD_VERS="8.7"):
    os.environ['PATH']= ':'.join([op.join(op.expanduser("~"), 'soft', 'lmod', LMOD_VERS, 'libexec'),
                                  os.environ['PATH']])
    os.system('/bin/bash -c "source %s"' % op.join(op.expanduser('~'), 'soft', 'lmod', LMOD_VERS, 'init', 'bash'))
    os.environ['LMOD_CMD'] =  op.join(op.expanduser('~'), 'soft', 'lmod', LMOD_VERS, 'libexec', 'lmod')
    
def check_easybuild_status():
    try:
        subprocess.call(["eb", "--help"])
    except FileNotFoundError:
        return("ERROR easybuild not found")

def generate_default_easybuild_config_arguments(workdir):
    HOME = op.expanduser("~")
    modulepath = op.join(workdir, 'easybuild', 'modules', 'all')
    buildpath = op.join(workdir, 'easybuild', 'build')
    containerpath = op.join(workdir, 'easybuild', 'containers')
    installpath = op.join(workdir, 'easybuild')
    repositorypath = op.join(workdir, 'easybuild', 'ebfiles_repo')
    # robotpath = op.join(workdir, 'easybuild', 'easyconfigs') ## let's use default's
    sourcepath = op.join(workdir, 'easybuild', 'sources')
    
    args = """--buildpath=%(buildpath)s  --installpath-modules=%(modulepath)s \
              --containerpath=%(containerpath)s --installpath=%(installpath)s \
              --repositorypath=%(repositorypath)s --sourcepath=%(sourcepath)s""" %{
                  'buildpath' : buildpath,
                  'modulepath' : modulepath,
                  'containerpath' : containerpath,
                  'installpath' : installpath,
                  'repositorypath' : repositorypath,
                  'sourcepath' : sourcepath}
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

## shell stuff end ------------------------------------------------------------------------------------

def list_all_easyconfigs(benchmark_yaml):
    return('not_implemented, perhaps import from workflow/schema?')

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

def parse_easyconfig(ec_fn, workdir):
    """
    Find and parse an easyconfig with specified filename,
    and return parsed easyconfig file (an EasyConfig instance).
    """
    opts, _ = set_up_configuration(args = generate_default_easybuild_config_arguments(workdir).split(),
                                   silent = True)    
    
    ec_path = det_easyconfig_paths([ec_fn])[0]
    
    # parse easyconfig file;
    # the 'parse_easyconfigs' function expects a list of tuples,
    # where the second item indicates whether or not the easyconfig file was automatically generated or not
    ec_dicts, _ = parse_easyconfigs([(ec_path, False)])
    
    # only retain first parsed easyconfig
    return ec_path, ec_dicts[0]['ec']

def get_easyconfig_full_path(easyconfig, workdir):
    ec_path, ec = parse_easyconfig(easyconfig, workdir)
    return(ec_path)

def get_envmodule_name_from_easyconfig(easyconfig, workdir):
    ec_path, ec = parse_easyconfig(easyconfig, workdir)
    return(os.path.join(ec['name'], det_full_ec_version(ec)))

def check_modules_tool():
    mod_tool = modules_tool()
    return("Modules tool: %s version %s" % (mod_tool.NAME, mod_tool.version))

def check_available_modules():
    return(modules_tool().available())

def check_envmodule_status(envmodule):
    mod_tool = modules_tool()    
    return(mod_tool.available(envmodule))


