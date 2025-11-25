from pathlib import Path
import os
import shutil
import pytest


from omnibenchmark.benchmark._node import BenchmarkNode
from tests.workflow.Snakemake_setup import SnakemakeSetup

from tests.workflow.path import data
from omnibenchmark.benchmark.execution_path import (
    ExecutionPathSet,
    ExecutionPathStageFile,
)


@pytest.fixture(
    scope="module",
    params=[
        pytest.param(
            {"current_dir": False, "keep_files": True, "test_case": "complete"},
            marks=pytest.mark.short,
            id="complete",
        ),
        pytest.param(
            {"current_dir": False, "keep_files": True, "test_case": "missing_files"},
            id="missing_files",
        ),
        pytest.param(
            {"current_dir": False, "keep_files": True, "test_case": "invalid_files"},
            id="invalid_files",
        ),
    ],
)
def snakemake_env(request):
    """Get environment from parametrization"""
    return request.param


@pytest.fixture(scope="module")
def benchmark_setup(tmp_path_factory, snakemake_env):
    """Fixture to prepare benchmark for the entire test class"""
    # Use tmp_path_factory for class-scoped fixtures
    if snakemake_env["current_dir"]:
        tmp_path = Path(os.getcwd()) / "tmp"
        tmp_path.mkdir(exist_ok=True)
    else:
        tmp_path = tmp_path_factory.mktemp("benchmark_test")

    benchmark_file = data / "Benchmark_001_trimmed.yaml"
    keep_files = snakemake_env["keep_files"]

    # Copy benchmark YAML to tmp_path so workflow can find it
    benchmark_file_in_tmp = tmp_path / "Benchmark_001_trimmed.yaml"
    with open(benchmark_file, "rb") as src, open(benchmark_file_in_tmp, "wb") as dst:
        dst.write(src.read())

    # Copy envs directory to tmp_path so environment files are available
    envs_src = data / "envs"
    envs_dst = tmp_path / "envs"
    if envs_src.exists():
        shutil.copytree(envs_src, envs_dst, dirs_exist_ok=True)

    with SnakemakeSetup(
        benchmark_file_in_tmp,
        keep_files=keep_files,
        cwd=tmp_path.as_posix(),
        out_dir=(tmp_path / "out").absolute(),
    ) as setup:
        benchmark = setup.benchmark
        assert benchmark.get_benchmark_name() == "Benchmark_001"

        # Run the whole workflow
        success = setup.workflow.run_workflow(
            benchmark,
            work_dir=tmp_path,
        )

        assert success is True

        if snakemake_env["test_case"] == "complete":
            pass  # All files should be present
        elif snakemake_env["test_case"] == "missing_files":
            # Remove some output files to simulate missing files
            for path in benchmark.get_execution_paths():
                for node in path:
                    if node.is_initial():
                        config = {
                            "input": "",
                            "output": "",
                            "dataset": node.module_id,
                        }
                        for output_file in node.get_output_paths(config):
                            if Path(output_file).exists():
                                Path(output_file).unlink()
        elif snakemake_env["test_case"] == "invalid_files":
            # Modify some output files to simulate invalid files
            for path in benchmark.get_execution_paths():
                for node in path:
                    if node.is_initial():
                        config = {
                            "input": "",
                            "output": "",
                            "dataset": node.module_id,
                        }
                        for output_file in node.get_output_paths(config):
                            if Path(output_file).exists():
                                with open(output_file, "a") as f:
                                    f.write("corrupted data")
        # return benchmark
        yield benchmark


# data = Path("tests/data")
# snakemake_env={"current_dir": True, "keep_files": True}
# benchmark = benchmark_setup(Path("."), snakemake_env = snakemake_env)
# benchmark_setup = benchmark


@pytest.mark.usefixtures("benchmark_setup")
class TestExecutionPathClasses:
    """Bundle all execution path tests together"""

    class TestExecutionPathSet:
        def test_ExecutionPathSet_initializes(self, benchmark_setup):
            eps = ExecutionPathSet(benchmark_setup)
            assert eps is not None
            assert eps.n_paths > 0
            assert eps.n_paths == 2
            assert len(eps.stages) > 0
            assert len(eps.stages) == 3

        def test_ExecutionPathSet_has_correct_stages(self, benchmark_setup):
            eps = ExecutionPathSet(benchmark_setup)
            expected_stages = list(benchmark_setup.get_stage_ids())
            assert eps.stages == expected_stages

        def test_ExecutionPathSet_creates_filedict(self, benchmark_setup):
            eps = ExecutionPathSet(benchmark_setup)
            filedict = eps.get_filedict(cumulative=False)
            assert filedict is not None
            assert len(filedict) > 0
            assert list(filedict.keys()) == eps.stages

        def test_ExecutionPathSet_filedict_structure(self, benchmark_setup):
            eps = ExecutionPathSet(benchmark_setup)
            filedict = eps.get_filedict(cumulative=False)
            for stage_id, stage_dict in filedict.items():
                assert isinstance(stage_dict, dict)
                assert isinstance(list(filedict[stage_id].keys())[0], BenchmarkNode)

                for node, stats in filedict[stage_id].items():
                    assert isinstance(stats, dict)
                    assert (
                        len(
                            set(filedict[stage_id][node].keys()).difference(
                                set(
                                    [
                                        "output_files",
                                        "valid_output_files",
                                        "observed_output_files",
                                        "empty_output_files",
                                        "missing_output_files",
                                        "invalid_output_files_input_file_is_newer",
                                        "invalid_output_files_repo_is_newer",
                                        "invalid_output_files",
                                        "n_valid",
                                        "n_observed",
                                        "n_empty",
                                        "n_missing",
                                        "n_invalid",
                                        "n_invalid_input_file_is_newer",
                                        "n_invalid_repo_is_newer",
                                        "n",
                                    ]
                                )
                            )
                        )
                        == 0
                    )

                    for n_test in [
                        "n_valid",
                        "n_observed",
                        "n_empty",
                        "n_missing",
                        "n_invalid",
                        "n_invalid_input_file_is_newer",
                        "n_invalid_repo_is_newer",
                        "n",
                    ]:
                        assert isinstance(filedict[stage_id][node][n_test], int)
                        assert filedict[stage_id][node][n_test] >= 0
                    for file_list in [
                        "output_files",
                        "valid_output_files",
                        "observed_output_files",
                        "empty_output_files",
                        "missing_output_files",
                        "invalid_output_files_input_file_is_newer",
                        "invalid_output_files_repo_is_newer",
                        "invalid_output_files",
                    ]:
                        assert isinstance(filedict[stage_id][node][file_list], set)
                        for file in filedict[stage_id][node][file_list]:
                            assert isinstance(file, str)
                            if file_list not in [
                                "output_files",
                                "missing_output_files",
                            ]:
                                assert (
                                    Path(file).is_file()
                                ), f"File {file} listed in {file_list} does not exist."
                            if file_list not in [
                                "output_files",
                                "missing_output_files",
                                "empty_output_files",
                            ]:
                                assert os.path.getsize(file) > 0

        def test_ExecutionPath_dict_encode(self, benchmark_setup):
            eps = ExecutionPathSet(benchmark_setup)
            for exec_path in eps.exec_path_dict.values():
                status, nmiss, failed_stage = exec_path.dict_encode(
                    full=True, cumulative=False
                )
                assert isinstance(status, str)
                assert isinstance(nmiss, int)
                assert failed_stage is None or isinstance(failed_stage, str)

        def test_ExecutionPath_stages_exist(self, benchmark_setup):
            eps = ExecutionPathSet(benchmark_setup)
            for exec_path in eps.exec_path_dict.values():
                assert len(exec_path.stages) > 0
                for stage in exec_path.stages:
                    assert stage in exec_path.exec_path

        def test_ExecutionPathStage_output_validity(self, benchmark_setup):
            eps = ExecutionPathSet(benchmark_setup)
            for exec_path in eps.exec_path_dict.values():
                for stage_id, stage in exec_path.exec_path.items():
                    validity = stage.output_is_valid(
                        type="all", cumulative=False, aggregate="none"
                    )
                    assert isinstance(validity, list)

        def test_get_module_timestamps(self, benchmark_setup):
            eps = ExecutionPathSet(benchmark_setup)
            assert eps.modules_repo_timestamps is not None
            assert isinstance(eps.modules_repo_timestamps, dict)

        def test_filedict_cumulative_flag(self, benchmark_setup):
            eps = ExecutionPathSet(benchmark_setup)
            eps.create_filedict(cumulative=False)
            assert not eps.filedict_is_cumulative
            eps.create_filedict(cumulative=True)
            assert eps.filedict_is_cumulative

    # @pytest.mark.usefixtures("benchmark_setup")
    class TestExecutionPath:
        """Test the ExecutionPath class"""

        def test_ExecutionPath_initializes(self, benchmark_setup):
            """Test that ExecutionPath can be initialized"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]
            assert exec_path is not None
            assert len(exec_path.stages) > 0
            assert len(exec_path.nodes) > 0
            assert len(exec_path.exec_path) > 0

        def test_ExecutionPath_has_correct_stages(self, benchmark_setup):
            """Test that ExecutionPath contains the expected stages"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]
            assert all(
                stage in benchmark_setup.get_stage_ids() for stage in exec_path.stages
            )

        def test_ExecutionPath_dict_encode_structure(self, benchmark_setup):
            """Test that dict_encode returns expected structure"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]
            status, nmiss, failed_stage = exec_path.dict_encode(
                full=False, cumulative=False
            )

            assert isinstance(status, str)
            assert isinstance(nmiss, int)
            assert nmiss >= 0
            assert failed_stage is None or isinstance(failed_stage, str)
            assert len(status) == len(exec_path.stages)

        def test_ExecutionPath_dict_encode_full(self, benchmark_setup):
            """Test dict_encode with full=True"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]
            status_full, nmiss_full, failed_full = exec_path.dict_encode(
                full=True, cumulative=False
            )
            status_simple, nmiss_simple, failed_simple = exec_path.dict_encode(
                full=False, cumulative=False
            )

            assert len(status_full) == len(status_simple)
            # Full status may provide more detail (F, R instead of I)
            assert set(status_full).issubset(set(".EFRIM "))
            assert set(status_simple).issubset(set(".EIOM "))

        def test_ExecutionPath_dict_encode_cumulative(self, benchmark_setup):
            """Test dict_encode with cumulative validity"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]
            status_cum, nmiss_cum, failed_cum = exec_path.dict_encode(
                full=True, cumulative=True
            )
            status_reg, nmiss_reg, failed_reg = exec_path.dict_encode(
                full=True, cumulative=False
            )

            # Cumulative should propagate failures forward
            assert isinstance(status_cum, str)
            assert isinstance(status_reg, str)

        def test_ExecutionPath_update(self, benchmark_setup):
            """Test that update refreshes the execution path state"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]

            # Should not raise an error
            exec_path.update()
            assert exec_path is not None

        def test_ExecutionPath_set_cumulative(self, benchmark_setup):
            """Test that set_cumulative propagates validity flags"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]

            # Should not raise an error
            exec_path.set_cumulative()

            # Check that cumulative flags are set for all stages
            for stage_id in exec_path.stages:
                stage = exec_path.exec_path[stage_id]
                for output_file in stage.output_stage_files:
                    # These should be set (not None)
                    assert (
                        output_file.any_is_newer_cumulative is not None
                        or output_file.any_is_newer_cumulative is None
                    )

    # @pytest.mark.usefixtures("benchmark_setup")
    class TestExecutionPathStage:
        """Test the ExecutionPathStage class"""

        def test_ExecutionPathStage_initializes(self, benchmark_setup):
            """Test that ExecutionPathStage can be initialized"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]
            stage = list(exec_path.exec_path.values())[0]

            assert stage is not None
            assert stage.stage_id is not None
            assert stage.module_id is not None
            assert stage.node is not None
            assert hasattr(stage, "output_stage_files")

        def test_ExecutionPathStage_has_output_files(self, benchmark_setup):
            """Test that stage has output files"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]
            stage = list(exec_path.exec_path.values())[0]

            assert len(stage.output_stage_files) > 0
            assert all(
                isinstance(f, ExecutionPathStageFile) for f in stage.output_stage_files
            )

        def test_ExecutionPathStage_output_exists(self, benchmark_setup):
            """Test output_exists method"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]
            stage = list(exec_path.exec_path.values())[0]

            # Test different aggregation modes
            exists_none = stage.output_exists(aggregate="none")
            exists_any = stage.output_exists(aggregate="any")
            exists_all = stage.output_exists(aggregate="all")

            assert isinstance(exists_none, list)
            assert isinstance(exists_any, list)
            assert isinstance(exists_all, list)
            assert len(exists_any) == 1
            assert len(exists_all) == 1

        def test_ExecutionPathStage_output_is_empty(self, benchmark_setup):
            """Test output_is_empty method"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]
            stage = list(exec_path.exec_path.values())[0]

            empty_none = stage.output_is_empty(aggregate="none")
            empty_any = stage.output_is_empty(aggregate="any")
            empty_all = stage.output_is_empty(aggregate="all")

            assert isinstance(empty_none, list)
            assert isinstance(empty_any, list)
            assert isinstance(empty_all, list)

        def test_ExecutionPathStage_output_is_valid_types(self, benchmark_setup):
            """Test output_is_valid with different types"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]
            stage = list(exec_path.exec_path.values())[0]

            for validity_type in ["all", "input_files_newer", "repo"]:
                validity = stage.output_is_valid(
                    type=validity_type, cumulative=False, aggregate="none"
                )
                assert isinstance(validity, list)
                assert all(v is None or isinstance(v, bool) for v in validity)

        def test_ExecutionPathStage_output_is_valid_aggregate(self, benchmark_setup):
            """Test output_is_valid with different aggregation modes"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]
            stage = list(exec_path.exec_path.values())[0]

            for agg_mode in ["none", "any", "all"]:
                validity = stage.output_is_valid(
                    type="all", cumulative=False, aggregate=agg_mode
                )
                assert isinstance(validity, list)

        def test_ExecutionPathStage_get_output_files(self, benchmark_setup):
            """Test get_output_files with different filters"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]
            stage = list(exec_path.exec_path.values())[0]

            filters = ["none", "valid", "observed", "empty", "missing", "invalid"]
            for filter_type in filters:
                files = stage.get_output_files(
                    filter=filter_type, cumulative=False, remove_none=True
                )
                assert isinstance(files, list)
                assert all(isinstance(f, str) for f in files)

        def test_ExecutionPathStage_dict_encode(self, benchmark_setup):
            """Test dict_encode returns valid status codes"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]
            stage = list(exec_path.exec_path.values())[0]

            code_simple = stage.dict_encode(full=False, cumulative=False)
            code_full = stage.dict_encode(full=True, cumulative=False)

            # Valid codes: '.', 'E', 'I', 'F', 'R', 'M', '?'
            assert code_simple in ".EIFRM?"
            assert code_full in ".EIFRM?"

        def test_ExecutionPathStage_has_repo_timestamp(self, benchmark_setup):
            """Test that stage has repo timestamp"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]
            stage = list(exec_path.exec_path.values())[0]

            # repo_timestamp can be None or a float
            assert stage.repo_timestamp is None or isinstance(
                stage.repo_timestamp, float
            )

        def test_ExecutionPathStage_update(self, benchmark_setup):
            """Test update method"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]
            stage = list(exec_path.exec_path.values())[0]

            # Should not raise an error
            stage.update()
            assert stage is not None

    # @pytest.mark.usefixtures("benchmark_setup")
    class TestExecutionPathStageFile:
        """Test the ExecutionPathStageFile class"""

        def test_ExecutionPathStageFile_initializes(self, benchmark_setup):
            """Test that ExecutionPathStageFile can be initialized"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]
            stage = list(exec_path.exec_path.values())[0]
            stage_file = stage.output_stage_files[0]

            assert stage_file is not None
            assert hasattr(stage_file, "output_file")
            assert hasattr(stage_file, "output_file_path")
            assert hasattr(stage_file, "input_files")
            assert hasattr(stage_file, "repo_timestamp")

        def test_ExecutionPathStageFile_has_attributes(self, benchmark_setup):
            """Test that stage file has expected attributes"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]
            stage = list(exec_path.exec_path.values())[0]
            stage_file = stage.output_stage_files[0]

            assert hasattr(stage_file, "output_file_exist")
            assert hasattr(stage_file, "output_file_empty")
            assert hasattr(stage_file, "timestamp")
            assert hasattr(stage_file, "input_files_exist")
            assert hasattr(stage_file, "all_input_files_exist")
            assert hasattr(stage_file, "input_file_is_newer")
            assert hasattr(stage_file, "repo_is_newer")

        def test_ExecutionPathStageFile_output_exists(self, benchmark_setup):
            """Test output_exists method"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]
            stage = list(exec_path.exec_path.values())[0]
            stage_file = stage.output_stage_files[0]

            exists = stage_file.output_exists()
            assert isinstance(exists, bool)

        def test_ExecutionPathStageFile_output_is_empty(self, benchmark_setup):
            """Test output_is_empty method"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]
            stage = list(exec_path.exec_path.values())[0]
            stage_file = stage.output_stage_files[0]

            is_empty = stage_file.output_is_empty()
            assert isinstance(is_empty, bool)

        def test_ExecutionPathStageFile_output_is_valid(self, benchmark_setup):
            """Test output_is_valid with different types"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]
            stage = list(exec_path.exec_path.values())[0]
            stage_file = stage.output_stage_files[0]

            for validity_type in ["all", "input_files_newer", "repo"]:
                validity = stage_file.output_is_valid(
                    type=validity_type, cumulative=False
                )
                assert validity is None or isinstance(validity, bool)

        def test_ExecutionPathStageFile_get_output_file(self, benchmark_setup):
            """Test get_output_file with different filters"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]
            stage = list(exec_path.exec_path.values())[0]
            stage_file = stage.output_stage_files[0]

            filters = [
                "none",
                "valid",
                "observed",
                "empty",
                "missing",
                "invalid",
                "input_files_newer",
                "repo",
            ]
            for filter_type in filters:
                result = stage_file.get_output_file(
                    filter=filter_type, cumulative=False
                )
                assert result is None or isinstance(result, str)

        def test_ExecutionPathStageFile_set_cumulative(self, benchmark_setup):
            """Test set_cumulative method"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]
            stage = list(exec_path.exec_path.values())[0]
            stage_file = stage.output_stage_files[0]

            # Test setting cumulative flags
            for validity_type in ["all", "input_files_newer", "repo"]:
                stage_file.set_cumulative(validity_type, True)
                stage_file.set_cumulative(validity_type, False)

            # Should not raise an error
            assert stage_file is not None

        def test_ExecutionPathStageFile_update(self, benchmark_setup):
            """Test update method"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]
            stage = list(exec_path.exec_path.values())[0]
            stage_file = stage.output_stage_files[0]

            # Should not raise an error
            stage_file.update()
            assert stage_file is not None

        def test_ExecutionPathStageFile_input_files_attributes(self, benchmark_setup):
            """Test input files related attributes"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]
            stage = list(exec_path.exec_path.values())[0]
            stage_file = stage.output_stage_files[0]

            assert isinstance(stage_file.input_files, list)
            assert isinstance(stage_file.input_file_paths, list)
            assert isinstance(stage_file.input_files_exist, list)
            assert isinstance(stage_file.all_input_files_exist, bool)

            if len(stage_file.input_files) > 0:
                assert len(stage_file.input_files) == len(stage_file.input_file_paths)
                assert len(stage_file.input_files_exist) == len(stage_file.input_files)

        def test_ExecutionPathStageFile_timestamps(self, benchmark_setup):
            """Test timestamp attributes"""
            eps = ExecutionPathSet(benchmark_setup)
            exec_path = list(eps.exec_path_dict.values())[0]
            stage = list(exec_path.exec_path.values())[0]
            stage_file = stage.output_stage_files[0]

            # Timestamp can be None or float
            assert stage_file.timestamp is None or isinstance(
                stage_file.timestamp, float
            )
            assert stage_file.repo_timestamp is None or isinstance(
                stage_file.repo_timestamp, float
            )

            if len(stage_file.input_files) > 0:
                assert isinstance(stage_file.input_files_timestamps, list)
                assert isinstance(stage_file.max_input_file_timestamp, float)
