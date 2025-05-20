#!/bin/bash
# This script creates a micromamba environment called "omnibenchmark"
# and installs the pinned system dependencies and an editable version of the
# omnibenchmark python package.
CMD=micromamba
ENV_FILE=conda-test-environment.yml
ENV_NAME=omnibenchmark

eval "$(micromamba shell hook --shell bash)"
$CMD create -n $ENV_NAME --yes
$CMD activate $ENV_NAME
$CMD install -f $ENV_FILE --yes
