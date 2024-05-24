#!/usr/bin/env python3

"""Easybuild-related stuff"""

import os.path as op

# import omni_schema.datamodel.omni_schema as model
# import yaml

# from ..utils import parse_instance

## basically from https://tutorial.easybuild.io/2021-lust/easybuild_library/

import sys

from easybuild.tools.filetools import remove_dir, which
from easybuild.tools.run import run_cmd
from easybuild.tools.options import set_up_configuration

    
opts, _ = set_up_configuration(args=[])

cmd = 'make'
cmd_path = which(cmd)
if cmd_path:
    print(">>> '%s' command found at %s" % (cmd, cmd_path))
else:
    sys.stderr.write("ERROR: '%s' command not found!\n" % cmd)
    sys.exit(1)

cmd = ' '.join(["make"] + sys.argv[1:])
out, ec = run_cmd(cmd)

print("\n>>> Output of '%s' (exit code %s):\n\n%s" % (cmd, ec, out))
