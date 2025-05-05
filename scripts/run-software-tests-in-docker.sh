#!/bin/bash
# Usage: ./run-software-tests.sh [conda|envmodules|all]

# Initialize workspace and prevent permission issues
setup_workspace() {
    TEMP_WORKSPACE="/tmp/repo"
    echo "Copying workspace to temporary directory with full user permissions..."
    mkdir -p "$TEMP_WORKSPACE"
    cp -a /home/user/workspace/. "$TEMP_WORKSPACE/"
    cd "$TEMP_WORKSPACE"
    echo "Running tests from temporary workspace..."
    set -ex
}

# Update environment.yml with wheel if available
setup_wheel() {
    if [ -d "dist" ]; then
        WHEEL_FILE=$(ls -t dist/omnibenchmark-*.whl 2>/dev/null | head -1)
        if [ -n "$WHEEL_FILE" ]; then
            echo "Found wheel file: $WHEEL_FILE"
            sed -i "\$c\\      - $WHEEL_FILE" test-environment.yml
            echo "Updated test-environment.yml to use wheel file"
            cat test-environment.yml
        else
            echo "No wheel file found in dist directory"
        fi
    else
        echo "No dist directory found, keeping source installation in test-environment.yml"
    fi
}

# Setup conda/micromamba environment
setup_conda_env() {
    echo "Setting up conda environment..."
    micromamba env create -f test-environment.yml -n omb --yes
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate omb
}

# Install extra dependencies for both conda and envmodules
install_extra_dependencies() {
    echo "Installing extra dependencies..."
    pip install -e ".[s3]"
    pip install -e ".[test]"
    # These could be optional with a flag if they're too heavy
    pip install numpy pandas scikit-learn scipy
}

# Setup environment modules
setup_envmodules() {
    echo "Setting up environment modules..."
    source "$LMOD_PKG"/init/profile
    export MODULEPATH="/home/user/workspace/tests/data/envs"
    module use $MODULEPATH

    # Validate modules (using avail instead of spider to avoid paging)
    TERM=dumb module avail 3.6.3-foss-2017b

    # Set PYTHONPATH for module initialization
    export PYTHONPATH=${PYTHONPATH}:$LMOD_DIR/../init
}

# Run conda tests
run_conda_tests() {
    echo "Running conda tests..."
    cd tests/software
    pytest -v -x -s -m "conda"
}

# Run environment module tests
run_envmodule_tests() {
    echo "Running environment module tests..."
    cd tests/software
    pytest -v -x -s -m "envmodules"
}

# Run all software tests
run_all_tests() {
    echo "Running all software tests..."
    cd tests/software
    pytest -v -x -s
}

# Main function to orchestrate the testing process
main() {
    local test_type=${1:-"all"}

    setup_workspace
    setup_wheel

    case "$test_type" in
        "conda")
            setup_conda_env
            install_extra_dependencies
            run_conda_tests
            ;;
        "envmodules")
            setup_conda_env  # Still need conda for Python environment
            install_extra_dependencies
            setup_envmodules
            run_envmodule_tests
            ;;
        "all")
            setup_conda_env
            install_extra_dependencies
            setup_envmodules
            run_all_tests
            ;;
        *)
            echo "Invalid test type: $test_type"
            echo "Usage: $0 [conda|envmodules|all]"
            exit 1
            ;;
    esac

    echo "Tests completed successfully!"
}

# Execute main function with the provided argument
main "$@"
