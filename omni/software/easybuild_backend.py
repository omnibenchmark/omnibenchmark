#!/usr/bin/env python

"""
Easybuild-powered software management, mostly by omniblock. Also includes singularity image generation.

Izaskun Mallona
Started 5th June 2024
"""

import subprocess, os, sys
import os.path as op
from easybuild.tools.module_naming_scheme.mns import ModuleNamingScheme
from easybuild.tools.module_naming_scheme.utilities import det_full_ec_version, is_valid_module_name
from easybuild.framework.easyconfig.tools import det_easyconfig_paths, parse_easyconfigs
from easybuild.tools.options import set_up_configuration
from easybuild.tools.modules import get_software_root_env_var_name, modules_tool

from importlib import resources as impresources
from . import templates

HOME=op.expanduser("~")

## shell-based stuff, partly to be replaced by direct eb API calls -------------------------------------

def check_lmod_status():
    try:
        ret = subprocess.check_call(["module", "--version"])
    except FileNotFoundError:
        return("ERROR module not found")
    return(ret)

def check_singularity_status():
    try:
        ret =subprocess.check_call(["singularity", "--version"])
    except FileNotFoundError:
        return("ERROR singularity not found")
    return(ret)

def check_easybuild_status():
    try:
        ret = subprocess.check_call(["eb", "--version"])
    except FileNotFoundError:
        return("ERROR easybuild not found")
    return(ret)

# def install_lmod():
#     print('not implemented')
#     try:
#         subprocess.call(["wget", "--help"])
#     except FileNotFoundError:
#         print("ERROR wget not found; please install it system-wise")
#     try:
#         subprocess.call(["make", "--help"])
#     except FileNotFoundError:
#         print("ERROR make (build tools) not found; please install them system-wise")    
#     try:
#         cmd = 'bash ./modules.sh'
#         output = subprocess.check_output(
#             cmd, stderr = subprocess.STDOUT, shell = True,
#             universal_newlines = True)
#     except subprocess.CalledProcessError as exc:
#         return("ERROR lmod install failed:", exc.returncode, exc.output)
#     else:
#         return("LOG lmod install output: \n{}\n".format(output))

# def export_lmod_env_vars(LMOD_VERS="8.7"):
#     os.environ['PATH']= ':'.join([op.join(op.expanduser("~"), 'soft', 'lmod', LMOD_VERS, 'libexec'),
#                                   os.environ['PATH']])
#     os.system('/bin/bash -c "source %s"' % op.join(op.expanduser('~'), 'soft', 'lmod', LMOD_VERS, 'init', 'bash'))
#     os.environ['LMOD_CMD'] =  op.join(op.expanduser('~'), 'soft', 'lmod', LMOD_VERS, 'libexec', 'lmod')

def export_lmod_env_vars(export_modules_script = op.join(op.dirname(__file__), 'utils', 'export_modules.sh')):
    cmd = 'bash ' + export_modules_script
    output = subprocess.run(
        cmd.split(' '), stderr = subprocess.STDOUT, shell = False,
        universal_newlines = True)

# def export_lmod_env_vars():
#     from omni.software import utils
#     export_modules_script = op.join(str(utils.__path__[0]), 'export_modules.sh')
#     cmd = 'bash ' + str(export_modules_script)
#     output = subprocess.check_output(
#         cmd, stderr = subprocess.STDOUT, shell = True,
#         universal_newlines = True)
    
def generate_default_easybuild_config_arguments(workdir):
    modulepath = op.join(workdir, 'easybuild', 'modules', 'all')
    buildpath = op.join(workdir, 'easybuild', 'build')
    containerpath = op.join(workdir, 'easybuild', 'containers')
    installpath = op.join(workdir, 'easybuild')
    repositorypath = op.join(workdir, 'easybuild', 'ebfiles_repo')
    robotpath = op.join(workdir, 'easybuild', 'easyconfigs') ## let's use default's
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

# ## do not use without handling the lmod / module envs explicitly
# def easybuild_easyconfig(easyconfig,
#                          workdir,
#                          threads,
#                          containerize = False,
#                          container_build_image = False):
#     """
#     easybuilds an easyconfig, potentially generating a (built) container image too
#     """
#     cmd = build_easybuild_easyconfig_command(easyconfig = easyconfing,
#                                              workdir = workdir,
#                                              threads = threads,
#                                              containerize = containerize,
#                                              container_build_image = container_build_image)
    
#     try:
#         output = subprocess.check_output(
#             cmd, stderr = subprocess.STDOUT, shell = True,
#             universal_newlines = True)
#     except subprocess.CalledProcessError as exc:
#         return("ERROR easybuild failed:", exc.returncode, exc.output)
#     else:
#         return("LOG easybuild: \n{}\n".format(output))

## shell stuff end ------------------------------------------------------------------------------------

def list_all_easyconfigs(benchmark_yaml):
    return('not_implemented, perhaps import from workflow/schema?')

# def easybuild_all_easyconfigs(benchmark_yaml, workdir, threads, containerize = False,
#                               container_build_image = False):
#     """
#     Iterates over all easyconfigs from a benchmark yaml and easybuilds them all
#     """

#     easyconfigs = list_all_easyconfigs()
#     for eb in easyconfigs:
#         easybuild_easyconfig(easyconfig = eb, workdir = workdir, threads = threads,
#                              containerize = containerize, container_build_image = container_build_image)

# not needed, the yaml will have both
# def get_easyconfig_from_envmodule_name(envmodule):
#     return('not implemented')

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
    try:
        ec_path, ec = parse_easyconfig(easyconfig, workdir)
        return(ec_path)
    except:
        raise FileNotFoundError('ERROR: easyconfig not found.\n')
        # sys.exit('Premature exit')

    
def get_envmodule_name_from_easyconfig(easyconfig, workdir):
    export_lmod_env_vars()
    ec_path, ec = parse_easyconfig(easyconfig, workdir)
    return(os.path.join(ec['name'], det_full_ec_version(ec)))

# # 0 success, there is a module tool
# # 1 error, there is none
# def check_module_tool():
#     mod_tool = modules_tool()
#     if mod_tool.NAME is None:
#         return(1) #'Module tool: not installed')
#     else:
#         return(0) # "Module tool: %s version %s" % (mod_tool.NAME, mod_tool.version))

def check_available_modules():
    return(modules_tool().available())

def check_envmodule_status(envmodule):
    mod_tool = modules_tool()
    if len(mod_tool.available(envmodule)) == 0:
        return('not installed')
    else:
        return('installed')

# # beware containerpath hardcoded to default
# def get_toolchain_container_path(toochain,
#                                  containerpath = op.join(HOME, '.local', 'easybuild', 'containers')):
#     return(op.join(containerpath, toolchain + '.sif'))
                                      

# def pull_base_container():
#     cmd = "singularity build --sandbox debian_trixie_slim.sif docker://debian:trixie-slim"

#     try:
#         output = subprocess.check_output(
#             cmd, stderr = subprocess.STDOUT, shell = True,
#             universal_newlines = True)
#     except subprocess.CalledProcessError as exc:
#         return("ERROR pulling failed %s %s:" %(exc.returncode, exc.output))
#     else:
#         return("LOG pulling: \n{}\n".format(output))
    

# def build_easybuild_easyconfig_command(easyconfig,
#                                        workdir,
#                                        threads,
#                                        containerize = False,
#                                        container_build_image = False):

    
#     args = generate_default_easybuild_config_arguments(workdir = workdir)
    
#     cmd = """eb %(easyconfig)s --robot --parallel=%(threads)s \
#               --detect-loaded-modules=unload --check-ebroot-env-vars=unset""" %{
#                   'easyconfig' : easyconfig,
#                   # 'args' : args,
#                   'threads' : threads}
#     if containerize:
#         cmd += " --container-config bootstrap=localimage,from=example.sif --experimental"
#         if container_build_image:
#             cmd += " --container-build-image"        
#     return(cmd)

def create_definition_file(easyconfig, singularity_recipe, envmodule, nthreads):
    template = impresources.files(templates) / 'ubuntu_jammy.txt'
    with open(template, 'rt') as ubuntu, open(singularity_recipe, 'w') as sing:
        for line in ubuntu.read().split('\n'):
            if 'EASYCONFIG' in line:
                line = line.replace('EASYCONFIG', easyconfig)
            if 'ENVMODULENAME' in line:
                line = line.replace('ENVMODULENAME', envmodule)
            if 'EASYBUILDNTHREADSINT' in line:
                line = line.replace('EASYBUILDNTHREADSINT', nthreads)
            sing.write(line + '\n')

def singularity_build(easyconfig, singularity_recipe):
    image_name = op.basename(easyconfig) + '.sif'
    try:
        cmd = """ export PATH=/usr/sbin:$PATH; 
        singularity build --fakeroot """ + image_name + ' ' + singularity_recipe
        output = subprocess.check_output(
            cmd, stderr = subprocess.STDOUT, shell = True,
            universal_newlines = True)
    except subprocess.CalledProcessError as exc:
        return("ERROR singularity build failed:", exc.returncode, exc.output)
    else:
        return("LOG singularity build output: \n{}\n".format(output))
