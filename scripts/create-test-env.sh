#!/bin/bash
# This script creates a conda environment called "omnibenchmark"
# and installs the pinned system dependencies and an editable version of the
# omnibenchmark python package.
CMD=conda
ENV_FILE=omni-environment.yml
ENV_NAME=omnibenchmark

# Create environment
$CMD create -n $ENV_NAME --yes

# Activate environment
$CMD activate $ENV_NAME

# Install dependencies from environment file
$CMD env update -f $ENV_FILE
