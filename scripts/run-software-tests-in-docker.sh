#!/bin/bash

TEMP_WORKSPACE="/tmp/repo"
echo "Copying workspace to temporary directory with full user permissions..."
mkdir -p "$TEMP_WORKSPACE"

# Copy the entire directory, preserving hidden files (like .git)
cp -a /home/user/workspace/. "$TEMP_WORKSPACE/"
cd "$TEMP_WORKSPACE"

echo "Running tests from temporary workspace..."

set -ex

# Find the most recent wheel file in the dist directory
if [ -d "dist" ]; then
    WHEEL_FILE=$(ls -t dist/omnibenchmark-*.whl 2>/dev/null | head -1)
    if [ -n "$WHEEL_FILE" ]; then
        echo "Found wheel file: $WHEEL_FILE"
        # Replace the last line with the exact wheel file path
        sed -i "\$c\\      - $WHEEL_FILE" test-environment.yml
        echo "Updated test-environment.yml to use wheel file"
        cat test-environment.yml
    else
        echo "No wheel file found in dist directory"
    fi
else
    echo "No dist directory found, keeping source installation in test-environment.yml"
fi

# Setup micromamba environment
micromamba env create -f test-environment.yml -n omb --yes
eval "$(micromamba shell hook --shell bash)"

# Install extra dependencies
pip install -e ".[s3]"
pip install -e ".[test]"

# These seem a bit heavy for testing, we could use lighter alternatives in tests
pip install numpy pandas scikit-learn scipy

# Setup modules
source "$LMOD_PKG"/init/profile
export MODULEPATH="$HOME/.local/easybuild/modules/all:${REPO_ROOT}/tests/data/envs"
module use $MODULEPATH

# Validate modules
module spider 3.6.3-foss-2017b

# Set PYTHONPATH
export PYTHONPATH=${PYTHONPATH}:$LMOD_DIR/../init

# Run tests
cd tests/software
pytest -v -x --show-capture=stderr .

echo "Tests completed successfully!"
