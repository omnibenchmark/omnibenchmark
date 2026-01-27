# Test run benchmark with different timeouts
import shutil


from tests.cli.cli_setup import OmniCLISetup

from .path import get_benchmark_data_path

data = get_benchmark_data_path()
envs = data / "envs"
timeout_benchmark = "timeout_benchmark.yaml"


def test_benchmark_run_no_timeout_does_not_fail(tmp_path):
    """Test that benchmark run with no timeout (default None) does not fail"""

    shutil.copy(data / timeout_benchmark, tmp_path / timeout_benchmark)
    shutil.copytree(envs, tmp_path / "envs")

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                (tmp_path / timeout_benchmark).as_posix(),
            ],
            cwd=tmp_path,
        )
    # should not fail, just wait for completion (uses default timeout of None)
    assert result.returncode == 0


def test_benchmark_run_with_timeout_less_than_task_fails(tmp_path):
    """Test that benchmark run with a timeout lower than the task duration fails"""

    shutil.copy(data / timeout_benchmark, tmp_path / timeout_benchmark)
    shutil.copytree(envs, tmp_path / "envs")

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                (tmp_path / timeout_benchmark).as_posix(),
                "--task-timeout",
                "1s",
            ],
            cwd=tmp_path,
        )
    # should fail, timeout of 1s is less than the sleep passed to the task (2 seconds)
    assert result.returncode != 0


def test_benchmark_run_high_timeout(tmp_path):
    shutil.copy(data / timeout_benchmark, tmp_path / timeout_benchmark)
    shutil.copytree(envs, tmp_path / "envs")

    with OmniCLISetup() as omni:
        result = omni.call(
            [
                "run",
                (tmp_path / timeout_benchmark).as_posix(),
                "--task-timeout",
                "5s",
            ],
            cwd=tmp_path,
        )
    # should not fail (task sleep is 2 seconds, we're waiting for 5 seconds)
    assert result.returncode == 0
