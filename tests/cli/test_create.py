"""Tests for the create command CLI functionality"""

import os
import subprocess

import pytest

from .asserts import assert_in_output, assert_startswith
from .cli_setup import OmniCLISetup


@pytest.mark.short
def test_create_benchmark_help():
    """Test that create benchmark command shows help correctly"""
    with OmniCLISetup() as omni:
        run = omni.call(["create", "benchmark", "--help"])
        assert run.returncode == 0
        assert_in_output(run.stdout, "Create a new benchmark from a template")
        assert_in_output(run.stdout, "[PATH]")
        assert_in_output(run.stdout, "--no-input")


@pytest.mark.short
def test_create_module_help():
    """Test that create module command shows help correctly"""
    with OmniCLISetup() as omni:
        run = omni.call(["create", "module", "--help"])
        assert run.returncode == 0
        assert_in_output(run.stdout, "Create a new module from a template")
        assert_in_output(run.stdout, "PATH")
        assert_in_output(run.stdout, "--no-input")


@pytest.mark.short
def test_create_benchmark_stub(tmp_path):
    """Test create benchmark command stub functionality"""
    with OmniCLISetup() as omni:
        run = omni.call(["create", "benchmark", str(tmp_path)], input="\n\n\n\n\n")
        assert run.returncode == 0
        assert_startswith(run.stdout, "Creating new benchmark at")
        assert_in_output(run.stdout, "Benchmark scaffolding created successfully")
        assert_in_output(run.stdout, "Next steps:")


@pytest.mark.short
def test_create_benchmark_without_destination(tmp_path):
    """Test create benchmark command without destination starts questionnaire"""
    with OmniCLISetup() as omni:
        # Use echo to provide default values to copier prompts
        run = omni.call(["create", "benchmark"], input="\n\n\n\n\n", cwd=str(tmp_path))
        assert run.returncode == 0
        assert_in_output(run.stdout, "Starting copier questionnaire for new benchmark")
        assert_in_output(run.stdout, "Benchmark scaffolding created successfully")


@pytest.mark.short
def test_create_benchmark_non_tty_fallback(tmp_path):
    """
    Make sure that that when the stdin is not a TTY (for example, when the CLI is
    invoked with stdin redirected from /dev/null), we command fall back to
    copier defaults instead of failing with an InteractiveSessionError.
    """
    # Exec CLI directly, providing stdin from /dev/null to simulate a non-TTY execution
    result = subprocess.run(
        [
            "python",
            "-m",
            "omnibenchmark.cli.main",
            "create",
            "benchmark",
            str(tmp_path),
        ],
        stdin=open(os.devnull, "r"),
        capture_output=True,
        text=True,
    )

    combined = (result.stdout or "") + "\n" + (result.stderr or "")
    assert result.returncode == 0
    # Confirm scaffolding was created and that the fallback to defaults was used
    assert "Benchmark scaffolding created successfully" in combined
    assert (
        "Falling back to defaults" in combined or "Input is not a terminal" in combined
    )


@pytest.mark.short
def test_create_benchmark_no_input_requires_destination():
    """Test that --no-input requires destination to be specified"""
    with OmniCLISetup() as omni:
        run = omni.call(["create", "benchmark", "--no-input"])
        assert run.returncode == 1
        assert_in_output(run.stdout, "Path is required when using --no-input")


@pytest.mark.short
def test_create_module_stub(tmp_path):
    """Test create module command stub functionality"""
    with OmniCLISetup() as omni:
        run = omni.call(
            ["create", "module", str(tmp_path)], input="\n\n\n\n\n\n\n\n\n\n"
        )
        assert run.returncode == 0
        assert_startswith(run.stdout, "Creating new module at")
        assert_in_output(run.stdout, "Module scaffolding created successfully")
        assert_in_output(run.stdout, "Next steps:")


@pytest.mark.short
def test_create_module_with_no_input(tmp_path):
    """Test create module command with --no-input flag"""
    with OmniCLISetup() as omni:
        run = omni.call(["create", "module", str(tmp_path), "--no-input"])
        assert run.returncode == 0
        assert_in_output(run.stdout, "Creating new module at")
        assert_in_output(run.stdout, "Module scaffolding created successfully")


@pytest.mark.short
def test_create_module_without_destination(tmp_path):
    """Test create module command with destination starts questionnaire"""
    with OmniCLISetup() as omni:
        module_path = tmp_path / "my_module"
        run = omni.call(
            ["create", "module", str(module_path)],
            input="\n\n\n\n\n\n\n\n\n\n",
            cwd=str(tmp_path),
        )
        assert run.returncode == 0
        assert_in_output(run.stdout, "Creating new module at")


@pytest.mark.short
def test_create_with_no_input_flag(tmp_path):
    """Test create commands with --no-input flag"""
    with OmniCLISetup() as omni:
        # Test benchmark
        benchmark_path = tmp_path / "benchmark"
        run = omni.call(["create", "benchmark", str(benchmark_path), "--no-input"])
        assert run.returncode == 0

        # Test module
        module_path = tmp_path / "module"
        run = omni.call(["create", "module", str(module_path), "--no-input"])
        assert run.returncode == 0


@pytest.mark.short
def test_create_with_debug_flag(tmp_path):
    """Test create commands with debug flag"""
    with OmniCLISetup() as omni:
        run = omni.call(
            ["--debug", "create", "benchmark", str(tmp_path)],
            input="\n\n\n\n\n",
        )
        assert run.returncode == 0
        # Debug output should be visible in the logs
        assert_in_output(run.stdout, "Creating new benchmark at")


@pytest.mark.short
def test_create_benchmark_existing_directory_no_input(tmp_path):
    """Test that --no-input fails when target directory exists and is not empty"""
    # Create a non-empty directory
    existing_file = tmp_path / "existing.txt"
    existing_file.write_text("existing content")

    with OmniCLISetup() as omni:
        run = omni.call(["create", "benchmark", str(tmp_path), "--no-input"])
        assert run.returncode == 1
        assert_in_output(run.stdout, "already exists and is not empty")
        assert_in_output(run.stdout, "Cannot proceed with --no-input")


@pytest.mark.short
def test_create_benchmark_git_initialization(tmp_path):
    """Test that git repository is initialized and initial commit is created"""
    with OmniCLISetup() as omni:
        run = omni.call(["create", "benchmark", str(tmp_path), "--no-input"])
        assert run.returncode == 0
        assert_in_output(run.stdout, "Initialized Git repository")
        assert_in_output(run.stdout, "Created initial commit")

        # Verify git repository exists
        git_dir = tmp_path / ".git"
        assert git_dir.exists()

        # Verify files were committed
        import subprocess

        result = subprocess.run(
            ["git", "log", "--oneline"],
            cwd=str(tmp_path),
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "Initial commit" in result.stdout


@pytest.mark.short
def test_create_benchmark_existing_git_repo_no_input(tmp_path):
    """Test that command fails with existing git repo when using --no-input"""
    # Initialize git repo and add a file
    import subprocess

    subprocess.run(["git", "init"], cwd=str(tmp_path), check=True)
    existing_file = tmp_path / "existing.txt"
    existing_file.write_text("existing content")

    with OmniCLISetup() as omni:
        run = omni.call(["create", "benchmark", str(tmp_path), "--no-input"])
        assert run.returncode == 1
        assert_in_output(run.stdout, "Git repository already exists")
        assert_in_output(run.stdout, "Cannot proceed with --no-input")


@pytest.mark.short
def test_create_benchmark_empty_directory_with_git(tmp_path):
    """Test that command works with empty directory that has git repo"""
    # Initialize empty git repo
    import subprocess

    subprocess.run(["git", "init"], cwd=str(tmp_path), check=True)

    with OmniCLISetup() as omni:
        run = omni.call(["create", "benchmark", str(tmp_path), "--no-input"])
        assert run.returncode == 0
        assert_in_output(run.stdout, "Benchmark scaffolding created successfully")

        # Should not initialize git again but should create commit
        # (though in this case it might skip git init)
        benchmark_file = tmp_path / "benchmark.yaml"
        assert benchmark_file.exists()


@pytest.mark.short
def test_create_benchmark_version_information(tmp_path):
    """Test that version information is included in output and generated files"""
    with OmniCLISetup() as omni:
        run = omni.call(["create", "benchmark", str(tmp_path), "--no-input"])
        assert run.returncode == 0
        assert_in_output(run.stdout, "Benchmark scaffolding created successfully")

        # Check that version is in CLI output
        assert_in_output(run.stdout, "using OmniBenchmark v")
        assert_in_output(run.stdout, "successfully with OmniBenchmark v")

        # Check that version is in README
        readme_file = tmp_path / "README.md"
        assert readme_file.exists()
        readme_content = readme_file.read_text()
        assert "OmniBenchmark](https://omnibenchmark.org) v" in readme_content

        # Check that version is in git commit message
        import subprocess

        result = subprocess.run(
            ["git", "log", "--format=%B", "-n", "1"],
            cwd=str(tmp_path),
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "Generated with OmniBenchmark v" in result.stdout


@pytest.mark.short
def test_create_benchmark_non_interactive_missing_params():
    """Test that --non-interactive requires all mandatory parameters"""
    with OmniCLISetup() as omni:
        run = omni.call(["create", "benchmark", "/tmp/test", "--non-interactive"])
        assert run.returncode == 1
        assert_in_output(run.stdout, "all mandatory parameters are required")
        assert_in_output(run.stdout, "--name")
        assert_in_output(run.stdout, "--author-name")
        assert_in_output(run.stdout, "--author-email")


@pytest.mark.short
def test_create_benchmark_partial_params_requires_all():
    """Test that providing some parameters requires all mandatory ones"""
    with OmniCLISetup() as omni:
        run = omni.call(["create", "benchmark", "/tmp/test", "--name", "My Benchmark"])
        assert run.returncode == 1
        assert_in_output(run.stdout, "all mandatory parameters are required")
        assert_in_output(run.stdout, "--author-name, --author-email")


@pytest.mark.short
def test_create_benchmark_non_interactive_success(tmp_path):
    """Test successful benchmark creation with --non-interactive and all parameters"""
    with OmniCLISetup() as omni:
        run = omni.call(
            [
                "create",
                "benchmark",
                str(tmp_path),
                "--non-interactive",
                "--name",
                "Test Benchmark",
                "--author-name",
                "Test Author",
                "--author-email",
                "test@example.com",
                "--license",
                "Apache-2.0",
                "--description",
                "A test benchmark",
            ]
        )
        assert run.returncode == 0
        assert_in_output(run.stdout, "Creating new benchmark at")
        assert_in_output(run.stdout, "Benchmark scaffolding created successfully")

        # Check that files were created with correct content
        benchmark_yaml = tmp_path / "benchmark.yaml"
        assert benchmark_yaml.exists()


@pytest.mark.short
def test_create_benchmark_parameter_flags_auto_enable_non_interactive(tmp_path):
    """Test that providing parameter flags automatically enables non-interactive mode"""
    with OmniCLISetup() as omni:
        run = omni.call(
            [
                "create",
                "benchmark",
                str(tmp_path),
                "--name",
                "Test Benchmark",
                "--author-name",
                "Test Author",
                "--author-email",
                "test@example.com",
            ]
        )
        assert run.returncode == 0
        assert_in_output(run.stdout, "Benchmark scaffolding created successfully")


@pytest.mark.short
def test_create_module_non_interactive_success(tmp_path):
    """Test successful module creation with --non-interactive and all parameters"""
    with OmniCLISetup() as omni:
        run = omni.call(
            [
                "create",
                "module",
                str(tmp_path),
                "--non-interactive",
                "--name",
                "test-module",
                "--author-name",
                "Test Author",
                "--author-email",
                "test@example.com",
            ]
        )
        assert run.returncode == 0
        assert_in_output(run.stdout, "Creating new module at")
        assert_in_output(run.stdout, "Module scaffolding created successfully")

        # Check that files were created
        omnibenchmark_yaml = tmp_path / "omnibenchmark.yaml"
        assert omnibenchmark_yaml.exists()


@pytest.mark.short
def test_create_module_partial_params_requires_all():
    """Test that providing some module parameters requires all mandatory ones"""
    with OmniCLISetup() as omni:
        run = omni.call(["create", "module", "/tmp/test", "--name", "my-module"])
        assert run.returncode == 1
        assert_in_output(run.stdout, "all mandatory parameters are required")
        assert_in_output(run.stdout, "--author-name, --author-email")


@pytest.mark.short
def test_create_module_python_entrypoint_template(tmp_path):
    """Test that Python entrypoint gets the correct template content"""
    with OmniCLISetup() as omni:
        run = omni.call(
            [
                "create",
                "module",
                str(tmp_path),
                "--name",
                "python-module",
                "--author-name",
                "Test Author",
                "--author-email",
                "test@example.com",
                "--entrypoint",
                "run.py",
            ]
        )
        assert run.returncode == 0

        # Check that the Python script was created with correct template
        entrypoint_file = tmp_path / "run.py"
        assert entrypoint_file.exists()

        content = entrypoint_file.read_text()
        assert content.startswith("#!/usr/bin/env python3")
        assert "import argparse" in content
        assert "def parse_args():" in content
        assert "def main():" in content
        assert 'if __name__ == "__main__":' in content
        assert "from main import process_data" in content

        # Check that it's executable
        import stat

        assert bool(entrypoint_file.stat().st_mode & stat.S_IXUSR)

        # Check that src/main.py was created
        src_main = tmp_path / "src" / "main.py"
        assert src_main.exists()
        main_content = src_main.read_text()
        assert "def process_data(args):" in main_content
        assert "output_dir = Path(args.output_dir)" in main_content


@pytest.mark.short
def test_create_module_r_entrypoint_template(tmp_path):
    """Test that R entrypoint gets the correct template content"""
    with OmniCLISetup() as omni:
        run = omni.call(
            [
                "create",
                "module",
                str(tmp_path),
                "--name",
                "r-module",
                "--author-name",
                "Test Author",
                "--author-email",
                "test@example.com",
                "--entrypoint",
                "run.R",
            ]
        )
        assert run.returncode == 0

        # Check that the R script was created with correct template
        entrypoint_file = tmp_path / "run.R"
        assert entrypoint_file.exists()

        content = entrypoint_file.read_text()
        assert content.startswith("#!/usr/bin/env Rscript")
        assert "library(argparse)" in content
        assert "ArgumentParser(description=" in content
        assert "args <- parser$parse_args()" in content
        assert 'source("src/main.R")' in content

        # Check that it's executable
        import stat

        assert bool(entrypoint_file.stat().st_mode & stat.S_IXUSR)

        # Check that src/main.R was created
        src_main = tmp_path / "src" / "main.R"
        assert src_main.exists()
        main_content = src_main.read_text()
        assert "process_data <- function(args)" in main_content
        assert "dir.create(args$output_dir" in main_content


@pytest.mark.short
def test_create_module_bash_entrypoint_template(tmp_path):
    """Test that bash entrypoint gets the correct template content"""
    with OmniCLISetup() as omni:
        run = omni.call(
            [
                "create",
                "module",
                str(tmp_path),
                "--name",
                "bash-module",
                "--author-name",
                "Test Author",
                "--author-email",
                "test@example.com",
                "--entrypoint",
                "run.sh",
            ]
        )
        assert run.returncode == 0

        # Check that the bash script was created with correct template
        entrypoint_file = tmp_path / "run.sh"
        assert entrypoint_file.exists()

        content = entrypoint_file.read_text()
        assert content.startswith("#!/bin/bash")
        assert "set -euo pipefail" in content
        assert "while [[ $# -gt 0 ]]" in content
        assert "--output_dir)" in content
        assert "--name)" in content
        assert 'source "$(dirname "$0")/src/main.sh"' in content

        # Check that it's executable
        import stat

        assert bool(entrypoint_file.stat().st_mode & stat.S_IXUSR)

        # Check that src/main.sh was created
        src_main = tmp_path / "src" / "main.sh"
        assert src_main.exists()
        main_content = src_main.read_text()
        assert "process_data() {" in main_content
        assert 'mkdir -p "$OUTPUT_DIR"' in main_content


@pytest.mark.short
def test_create_module_default_entrypoint_with_no_input(tmp_path):
    """Test that --no-input creates default bash script with template"""
    with OmniCLISetup() as omni:
        run = omni.call(["create", "module", str(tmp_path), "--no-input"])
        assert run.returncode == 0

        # Check that the default run.sh script was created with template
        entrypoint_file = tmp_path / "run.sh"
        assert entrypoint_file.exists()

        content = entrypoint_file.read_text()
        assert content.startswith("#!/bin/bash")
        assert "set -euo pipefail" in content
        assert "process_data" in content

        # Check that it's executable
        import stat

        assert bool(entrypoint_file.stat().st_mode & stat.S_IXUSR)

        # Check that src/main.sh was created
        src_main = tmp_path / "src" / "main.sh"
        assert src_main.exists()
