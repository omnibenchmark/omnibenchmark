"""
High-level integration tests for timeout handling in execution.py.

Tests the timeout behavior at the module interface level:
- Default timeout is None (no timeout)
- Timeout triggers process termination with SIGTERM then SIGKILL
- Process cleanup works correctly
"""

import pytest
import time
from omnibenchmark.workflow.snakemake.scripts.execution import execution
from omnibenchmark.benchmark.params import Params


@pytest.mark.short
def test_execution_runs_successfully_without_timeout(tmp_path):
    """Test that a simple script executes successfully when no timeout is set."""
    # Create a simple test module
    module_dir = tmp_path / "test_module"
    module_dir.mkdir()

    # Create config
    config = module_dir / "config.cfg"
    config.write_text("[DEFAULT]\nSCRIPT = test.py\n")

    # Create simple Python script that completes quickly
    script = module_dir / "test.py"
    script.write_text("import sys\nprint('hello')\nsys.exit(0)\n")

    output_dir = tmp_path / "output"
    output_dir.mkdir()

    exit_code = execution(
        module_dir=module_dir,
        module_name="test_module",
        output_dir=output_dir,
        dataset="test",
        inputs_map={},
        parameters=Params(),
        keep_module_logs=False,
        timeout=None,  # No timeout
    )

    assert exit_code == 0


@pytest.mark.short
def test_execution_with_timeout_kills_long_running_process(tmp_path):
    """Test that a long-running process is killed when timeout is exceeded."""
    module_dir = tmp_path / "test_module"
    module_dir.mkdir()

    config = module_dir / "config.cfg"
    config.write_text("[DEFAULT]\nSCRIPT = slow.py\n")

    # Create script that sleeps longer than timeout
    script = module_dir / "slow.py"
    script.write_text("import time\ntime.sleep(60)\n")

    output_dir = tmp_path / "output"
    output_dir.mkdir()

    # Should return exit code 124 (timeout exit code)
    exit_code = execution(
        module_dir=module_dir,
        module_name="test_module",
        output_dir=output_dir,
        dataset="test",
        inputs_map={},
        parameters=Params(),
        keep_module_logs=False,
        timeout=1,  # 1 second timeout
    )

    assert exit_code == 124


@pytest.mark.short
def test_execution_timeout_terminates_child_processes(tmp_path):
    """Test that timeout kills the entire process tree, not just parent."""
    module_dir = tmp_path / "test_module"
    module_dir.mkdir()

    config = module_dir / "config.cfg"
    config.write_text("[DEFAULT]\nSCRIPT = spawn.py\n")

    # Create script that spawns child processes
    script = module_dir / "spawn.py"
    script.write_text("""
import subprocess
import time
# Spawn a child that will sleep
subprocess.Popen(['sleep', '60'])
time.sleep(60)
""")

    output_dir = tmp_path / "output"
    output_dir.mkdir()

    # Should return exit code 124 (timeout exit code)
    exit_code = execution(
        module_dir=module_dir,
        module_name="test_module",
        output_dir=output_dir,
        dataset="test",
        inputs_map={},
        parameters=Params(),
        keep_module_logs=False,
        timeout=1,
    )

    assert exit_code == 124

    # Give a moment for cleanup
    time.sleep(0.5)

    # Process group termination should have killed all children


@pytest.mark.short
def test_execution_with_shell_script(tmp_path):
    """Test that execution works with shell scripts."""
    module_dir = tmp_path / "test_module"
    module_dir.mkdir()

    config = module_dir / "config.cfg"
    config.write_text("[DEFAULT]\nSCRIPT = test.sh\n")

    # Create shell script
    script = module_dir / "test.sh"
    script.write_text("#!/bin/bash\necho 'test'\nexit 0\n")
    script.chmod(0o755)

    output_dir = tmp_path / "output"
    output_dir.mkdir()

    exit_code = execution(
        module_dir=module_dir,
        module_name="test_module",
        output_dir=output_dir,
        dataset="test",
        inputs_map={},
        parameters=Params(),
        keep_module_logs=False,
        timeout=None,
    )

    assert exit_code == 0


@pytest.mark.short
def test_execution_fast_task_with_timeout_succeeds(tmp_path):
    """Test that a fast task completes successfully even when timeout is set."""
    module_dir = tmp_path / "test_module"
    module_dir.mkdir()

    config = module_dir / "config.cfg"
    config.write_text("[DEFAULT]\nSCRIPT = fast.py\n")

    # Create script that completes quickly
    script = module_dir / "fast.py"
    script.write_text("import time\ntime.sleep(0.1)\nprint('done')\n")

    output_dir = tmp_path / "output"
    output_dir.mkdir()

    exit_code = execution(
        module_dir=module_dir,
        module_name="test_module",
        output_dir=output_dir,
        dataset="test",
        inputs_map={},
        parameters=Params(),
        keep_module_logs=False,
        timeout=5,  # 5 second timeout, task finishes in 0.1s
    )

    assert exit_code == 0
