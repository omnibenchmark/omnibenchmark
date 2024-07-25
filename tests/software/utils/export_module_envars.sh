#!/bin/bash

# initialize environment for modules tool
if [ -f $HOME/moduleshome ]; then export MODULESHOME=$(cat $HOME/moduleshome); fi
source $(cat $HOME/mod_init)

# make sure 'eb' is available via $PATH, and that $PYTHONPATH is set (some tests expect that);
# also pick up changes to $PATH set by sourcing $MOD_INIT
WORKDIR=$GITHUB_WORKSPACE/easybuild-easyconfigs
export PATH=$(cat $HOME/path)
export PYTHONPATH=$WORKDIR
export EASYBUILD_MODULES_TOOL=Lmod

