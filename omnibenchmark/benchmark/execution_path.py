from datetime import datetime
import os
from pathlib import Path
from typing import Union
from omnibenchmark.benchmark import Benchmark, BenchmarkNode
from omnibenchmark.model.benchmark import Repository
from omnibenchmark.model.repo import get_repo_hash
from abc import ABCMeta, abstractmethod


def n2f(bl: Union[bool, None]) -> bool:
    """Convert None to False, keep True as True."""
    return bl is True


class ArtifactFile(metaclass=ABCMeta):
    """
    Class to represent an artifact file in the execution path.
    """

    def __init__(self, file_path: Union[str, Path]):
        """
        Initialize the ArtifactFile.
        Args:
            file_path (str): The path to the artifact file
        """
        self.file_path = file_path

    @abstractmethod
    def exists(self) -> bool:
        """Check if the artifact file exists."""
        pass

    @abstractmethod
    def is_empty(self) -> bool:
        """Check if the artifact file is empty."""
        pass

    def get_timestamp(self) -> Union[float, None]:
        """Get the modification timestamp of the artifact file."""
        pass


class LocalArtifactFile(ArtifactFile):
    def __init__(self, file_path: Union[str, Path]):
        if isinstance(file_path, str):
            file_path = Path(file_path)
        self.file_path = file_path

    def exists(self):
        return self.file_path.is_file()

    def is_empty(self):
        return self.exists() and self.file_path.stat().st_size == 0

    def get_timestamp(self) -> Union[float, None]:
        if self.exists():
            return self.file_path.stat().st_mtime
        else:
            return None


class ArtifactFileConfig:
    """Global configuration for ArtifactFile implementation."""

    _implementation = None

    @classmethod
    def set_implementation(cls, implementation_class):
        """Set the global ArtifactFile implementation."""
        cls._implementation = implementation_class

    @classmethod
    def get_implementation(cls):
        """Get the current ArtifactFile implementation."""
        if cls._implementation is None:
            return LocalArtifactFile  # Default implementation
        return cls._implementation


def create_artifact_file(file_path: Union[str, Path]) -> ArtifactFile:
    """Factory function to create ArtifactFile instances."""
    implementation = ArtifactFileConfig.get_implementation()
    return implementation(file_path)


def get_repo_timestamp(
    repo: Repository, cache_dir: Path = Path(".snakemake") / "repos"
) -> Union[float, None]:
    """
    Get the timestamp of the repository used in the module.
    Returns the latest modification time of files in the repository directory.
    If the repository directory does not exist, returns None.
    """
    module_dir = cache_dir / get_repo_hash(repo.url, repo.commit)
    if module_dir.is_dir():
        timestamps = [
            create_artifact_file(f).get_timestamp()
            for f in module_dir.iterdir()
            if create_artifact_file(f).exists()
            and not create_artifact_file(f).is_empty()
        ]
        timestamps = [t for t in timestamps if t is not None]
        if len(timestamps) > 0:
            timestamp = max(timestamps)
        else:
            timestamp = None
    else:
        timestamp = None
    return timestamp


def get_module_timestamps(
    b: Benchmark, cache_dir: Path = Path(".snakemake") / "repos"
) -> dict:
    """
    Get the timestamps of the module repositories used in the benchmark.
    Returns a dictionary with stage ids as keys and another dictionary as value,
    which has module ids as keys and timestamps as values.
    """

    modules_repo_timestamps = {}
    for st in b.get_stage_ids():
        modules_repo_timestamps[st] = {}
        for node in b.get_nodes_by_stage_id(st):
            if node.module_id not in modules_repo_timestamps[st].keys():
                timestamp = get_repo_timestamp(
                    node.get_repository(), cache_dir=cache_dir
                )
                modules_repo_timestamps[st][node.module_id] = timestamp
    return modules_repo_timestamps


class ExecutionPathStageFile:
    """
    Class to represent an output file in an execution path stage,
    and check its validity based on input files and repository timestamps.
    """

    def __init__(
        self,
        input_files: list[ArtifactFile],
        repo_timestamp: Union[float, None],
        output_file: ArtifactFile,
    ):
        """
        Initialize the ExecutionPathStageFile.
        Args:
            input_files (list): List of input file paths
            repo_timestamp (float or None): The repository timestamp
            output_file (str): The output file path
        """
        self.output_file = output_file
        self.output_file_path = output_file.file_path
        self.repo_timestamp = repo_timestamp
        self.input_file_paths = [f.file_path for f in input_files]
        self.input_files = input_files
        self.update()

    def update(self):
        self.output_file_exist = self.output_file.exists()
        self.output_file_empty = self.output_file.is_empty()
        self.timestamp = self.output_file.get_timestamp()
        self.input_files_exist = [
            f.exists() and not f.is_empty() for f in self.input_files
        ]
        self.all_input_files_exist = all(self.input_files_exist)
        self.input_files_timestamps = [
            f.get_timestamp() if self.input_files_exist[i] else None
            for i, f in enumerate(self.input_files)
        ]
        self.max_input_file_timestamp = (
            datetime.now().timestamp() * 2 + 1e6
            if any([tst is None for tst in self.input_files_timestamps])
            and len(self.input_files_timestamps) > 0
            else max([tst for tst in self.input_files_timestamps if tst is not None])
            if len(self.input_files_timestamps) > 0
            else 0
        )

        self.input_file_is_newer = (
            self.max_input_file_timestamp > self.timestamp
            if self.timestamp is not None
            else None
        )

        self.repo_is_newer = (
            self.repo_timestamp > self.timestamp
            if self.timestamp is not None and self.repo_timestamp is not None
            else None
        )

        self.any_is_newer = self.input_file_is_newer or self.repo_is_newer

        self.any_is_newer_cumulative: Union[bool, None] = None
        self.repo_is_newer_cumulative: Union[bool, None] = None
        self.input_file_is_newer_cumulative: Union[bool, None] = None
        self.all_input_files_exist_cumulative: Union[bool, None] = None

    def output_is_valid(
        self, type: str = "all", cumulative: bool = False
    ) -> Union[bool, None]:
        """
        Check if the output file is valid based on the specified type.
        Args:
            type (str): One of 'all', 'input_files_newer', 'input_files_exist', 'repo'
            cumulative (bool): Whether to use cumulative validity
        Returns:
            bool or None: True if valid, False if invalid, None if unknown
        """
        if type == "all":
            is_invalid = (
                self.any_is_newer_cumulative if cumulative else self.any_is_newer
            )
        elif type == "input_files_newer":
            is_invalid = (
                self.input_file_is_newer_cumulative
                if cumulative
                else self.input_file_is_newer
            )
        elif type == "input_files_exist":
            is_invalid = (
                self.all_input_files_exist_cumulative
                if cumulative
                else self.all_input_files_exist
            )
        elif type == "repo":
            is_invalid = (
                self.repo_is_newer_cumulative if cumulative else self.repo_is_newer
            )
        else:
            raise ValueError(
                "type must be either 'all', 'input_file_newer', 'input_files_exist' or 'repo'"
            )
        return not n2f(is_invalid)

    def output_exists(self) -> bool:
        return self.output_file_exist

    def output_is_empty(self) -> bool:
        return self.output_file_empty

    def get_output_file(
        self, filter: str = "none", cumulative: bool = False
    ) -> Union[str, None]:
        """
        Get the output file path based on the filter.
        Args:
            filter (str): One of 'none', 'valid', 'observed', 'empty', 'missing', 'invalid',
                'input_files_newer', 'repo'
            cumulative (bool): Whether to use cumulative validity
        """
        if filter == "none":
            do_return = True
        elif filter == "valid":
            do_return = (
                n2f(self.output_exists())
                and n2f(self.output_is_valid(type="all", cumulative=cumulative))
                and not n2f(self.output_is_empty())
            )
        elif filter == "observed":
            do_return = n2f(self.output_exists()) and n2f(
                self.output_is_valid(type="all", cumulative=cumulative)
            )
        elif filter == "empty":
            do_return = n2f(self.output_is_empty()) and n2f(
                self.output_is_valid(type="all", cumulative=cumulative)
            )
        elif filter == "missing":
            do_return = not n2f(self.output_exists())
        elif filter == "invalid":
            do_return = n2f(self.output_exists()) and not n2f(
                self.output_is_valid(type="all", cumulative=cumulative)
            )
        elif filter == "input_files_newer":
            do_return = n2f(self.output_exists()) and not n2f(
                self.output_is_valid(type="input_files_newer", cumulative=cumulative)
            )
        elif filter == "repo":
            do_return = n2f(self.output_exists()) and not n2f(
                self.output_is_valid(type="repo", cumulative=cumulative)
            )
        else:
            raise ValueError(
                "filter must be either 'none', 'observed', 'missing', 'invalid', 'input_files_newer' or 'repo'"
            )
        if do_return:
            return str(self.output_file_path)
        else:
            return None

    def set_cumulative(self, type: str, val: bool):
        """
        Set the cumulative validity flags.
        Args:
            type (str): One of 'all', 'input_files_newer', 'input_files_exist', 'repo'
            val (bool): Value to set
        """
        if type == "all":
            self.any_is_newer_cumulative = val
        elif type == "input_files_newer":
            self.input_file_is_newer_cumulative = val
        elif type == "input_files_exist":
            self.all_input_files_exist_cumulative = val
        elif type == "repo":
            self.repo_is_newer_cumulative = val
        else:
            raise ValueError(
                "type must be either 'all', 'input_files_newer', 'input_files_exist' or 'repo'"
            )


class ExecutionPathStage:
    """
    Class to represent an execution path stage,
    containing multiple output files and their validity checks.
    """

    def __init__(
        self,
        node: BenchmarkNode,
        dataset: str,
        input_dir: str = "",
        input_files: list[str] = [],
        repo_timestamp: Union[float, None] = None,
    ):
        """
        Initialize the ExecutionPathStage.
        Args:
            node (BenchmarkNode): The benchmark node for this stage
            input_dir (str): The input directory for the stage
            repo_timestamp (float or None): The repository timestamp for the stage
        """
        self.node = node
        self.dataset = dataset
        self.input_dir = input_dir
        self.stage_id = node.stage_id
        self.module_id = node.module_id
        self.param_id = node.param_id
        self.parameters = node.parameters

        if repo_timestamp is None:
            self.repo_timestamp = get_repo_timestamp(node.get_repository())
        else:
            self.repo_timestamp = repo_timestamp

        config = {
            "input": self.input_dir,
            "output": "",
            "dataset": self.dataset,
        }
        output_files = node.get_output_paths(config)

        # input_files = []
        self.output_stage_files = [
            ExecutionPathStageFile(
                input_files=[create_artifact_file(f) for f in input_files],
                repo_timestamp=self.repo_timestamp,
                output_file=create_artifact_file(output_file),
            )
            for output_file in output_files
        ]

        self.output_dir = os.path.commonpath(
            list(
                set(
                    [
                        os.path.dirname(f.output_file_path)
                        for f in self.output_stage_files
                    ]
                )
            )
        )

        # TODO: check if log file is newer than input files, can be both newer or older than output files
        stdout_log = Path(self.output_dir) / "stdout.log"
        self.stdout_log = stdout_log if stdout_log.is_file() else None
        stderr_log = Path(self.output_dir) / "stderr.log"
        self.stderr_log = stderr_log if stderr_log.is_file() else None

    def update(self):
        for f in self.output_stage_files:
            f.update()

    def output_is_valid(
        self, type: str = "all", cumulative: bool = False, aggregate: str = "none"
    ) -> list[Union[bool, None]]:
        """
        Check validity of output files based on the specified type and aggregation.
        Args:
            type (str): One of 'all', 'input_files_newer', 'input_files_exist', 'repo'
            cumulative (bool): Whether to check cumulative validity
            aggregate (str): One of 'none', 'any', 'all'
        """
        is_valid_list = [
            f.output_is_valid(type=type, cumulative=cumulative)
            for f in self.output_stage_files
        ]
        if aggregate == "none":
            return is_valid_list
        else:
            if any([v is None for v in is_valid_list]):
                return [None]
            if aggregate == "any":
                return [any(is_valid_list)]
            elif aggregate == "all":
                return [all(is_valid_list)]
            else:
                raise ValueError("aggregate must be either 'none', 'any' or 'all'")

    def output_exists(self, aggregate: str = "none") -> list[bool]:
        """
        Check existence of output files with specified aggregation.
        Args:
            aggregate (str): One of 'none', 'any', 'all'
        """
        if aggregate == "none":
            return [f.output_exists() for f in self.output_stage_files]
        elif aggregate == "any":
            return [any(f.output_exists() for f in self.output_stage_files)]
        elif aggregate == "all":
            return [all(f.output_exists() for f in self.output_stage_files)]
        else:
            raise ValueError("aggregate must be either 'none', 'any' or 'all'")

    def output_is_empty(self, aggregate: str = "none") -> list[bool]:
        """
        Check emptiness of output files with specified aggregation.
        Args:
            aggregate (str): One of 'none', 'any', 'all'
        """
        if aggregate == "none":
            return [f.output_is_empty() for f in self.output_stage_files]
        elif aggregate == "any":
            return [any(f.output_is_empty() for f in self.output_stage_files)]
        elif aggregate == "all":
            return [all(f.output_is_empty() for f in self.output_stage_files)]
        else:
            raise ValueError("aggregate must be either 'none', 'any' or 'all'")

    def set_cumulative(self, type: str, val: bool):
        """
        Set the cumulative validity flags for all output files.
        Args:
            type (str): One of 'all', 'input_files_newer', 'input_files_exist', 'repo'
            val (bool): Value to set
        """
        for f in self.output_stage_files:
            f.set_cumulative(type=type, val=val)

    def get_output_files(
        self, filter: str = "none", cumulative: bool = False, remove_none: bool = True
    ) -> list[str] | list[Union[str, None]]:
        """
        Get output files based on the filter and cumulative validity.
        Args:
            filter (str): One of 'none', 'observed','empty', 'missing', 'invalid',
                'input_files_newer', 'repo'
            cumulative (bool): Whether to use cumulative validity
            remove_none (bool): Whether to remove None values from the result
        """
        files = [
            f.get_output_file(filter=filter, cumulative=cumulative)
            for f in self.output_stage_files
        ]
        if remove_none:
            return [f for f in files if f is not None]
        else:
            return files

    def dict_encode(self, full: bool = False, cumulative: bool = False) -> str:
        """
        Encode the stage status as a single character code.
        Args:
            full (bool): Whether to include full status information
            cumulative (bool): Whether to use cumulative validity
        Returns:
            str: Single character code representing the stage status
                '.': all outputs valid
                'E': empty outputs
                'I': incomplete (some outputs invalid, 'F' or 'R')
                'F': input files newer than outputs
                'R': repo newer than outputs
                'M': missing outputs
        """
        output_exists = n2f(self.output_exists(aggregate="all")[0])
        output_is_empty = n2f(self.output_is_empty(aggregate="all")[0])
        if (
            output_is_empty
            and self.output_is_valid(
                type="all", cumulative=cumulative, aggregate="all"
            )[0]
        ):
            return "E"
        if (
            output_exists
            and self.output_is_valid(
                type="all", cumulative=cumulative, aggregate="all"
            )[0]
        ):
            return "."
        if (
            not full
            and output_exists
            and not n2f(
                self.output_is_valid(
                    type="all", cumulative=cumulative, aggregate="all"
                )[0]
            )
        ):
            return "I"
        if (
            full
            and output_exists
            and not n2f(
                self.output_is_valid(
                    type="all", cumulative=cumulative, aggregate="all"
                )[0]
            )
        ):
            if not n2f(
                self.output_is_valid(
                    type="input_files_newer", cumulative=cumulative, aggregate="all"
                )[0]
            ):
                return "F"
            if not n2f(
                self.output_is_valid(
                    type="input_files_newer", cumulative=cumulative, aggregate="all"
                )[0]
            ):
                return "R"
        if not output_exists:
            return "M"
        else:
            return "?"


class ExecutionPath:
    """
    Class to represent an execution path, consisting of multiple stages.
    """

    def __init__(
        self, nodes: list[BenchmarkNode], repo_timestamps: Union[dict, None] = None
    ):
        """
        Initialize the ExecutionPath.
        Args:
            nodes (list): List of BenchmarkNode objects representing the execution path
            repo_timestamps (dict or None): Dictionary with stage ids as keys and
                repository timestamps as values
        """
        self.nodes = nodes
        self.stages = [node.stage_id for node in nodes]

        # nodes = b.get_execution_paths()[0]
        self.exec_path = {}
        input_dir = ""
        input_files = []
        for i, st in enumerate(self.stages):
            self.exec_path[st] = ExecutionPathStage(
                nodes[i],
                dataset=nodes[0].module_id,  # use first node's module_id as dataset
                input_files=input_files,
                input_dir=input_dir,
                repo_timestamp=repo_timestamps[st] if repo_timestamps else None,
            )
            input_dir = self.exec_path[st].output_dir
            input_files = self.exec_path[st].get_output_files(filter="none")
        self.set_cumulative()

    def set_cumulative(self):
        """
        Set cumulative validity for all stages in the execution path.
        """
        for i in range(len(self.stages)):
            if i == 0:
                for t in ["all", "input_files_newer", "input_files_exist", "repo"]:
                    is_valid = self.exec_path[self.stages[i]].output_is_valid(
                        type=t, cumulative=False
                    )
                    is_valid_all = all(is_valid)
                    self.exec_path[self.stages[i]].set_cumulative(t, not is_valid_all)
            else:
                for t in ["all", "input_files_newer", "input_files_exist", "repo"]:
                    is_valid = self.exec_path[self.stages[i]].output_is_valid(
                        type=t, cumulative=False
                    )
                    is_valid_all = all(is_valid)
                    is_valid_cumulative = self.exec_path[
                        self.stages[i - 1]
                    ].output_is_valid(type=t, cumulative=True)
                    is_valid_cumulative_all = all(is_valid_cumulative)
                    self.exec_path[self.stages[i]].set_cumulative(
                        t, not (is_valid_all and is_valid_cumulative_all)
                    )

    def update(self):
        """
        Update the execution path by refreshing the status of all stages.
        """
        for st in self.stages:
            self.exec_path[st].update()
        self.set_cumulative()

    def dict_encode(
        self,
        full: bool = False,
        cumulative: bool = False,
        stages: Union[list[str], None] = None,
    ) -> tuple[str, int, Union[str, None]]:
        """
        Encode the execution path status as a string of single character codes.
        Args:
            full (bool): Whether to include full status information
            stages (list or None): List of stage ids to include, or None for all stages
        Returns:
            tuple: (status string, number of missing/invalid stages, first failed stage id or None)
        """
        if stages is None:
            stages = self.stages
        tmps1 = ""
        nmiss = 0
        is_failed_stage = None
        for st in stages:
            if st in self.stages:
                tmpcode = self.exec_path[st].dict_encode(
                    full=full, cumulative=cumulative
                )
                if tmpcode != ".":
                    if is_failed_stage is None:
                        is_failed_stage = st
                    nmiss += 1
                tmps1 += tmpcode
            else:
                tmps1 += " "
        return tmps1, nmiss, is_failed_stage

    def touch_outputs(self, type: str = "all", cumulative: bool = False):
        """
        Touch output files that are invalid based on the specified type.
        Args:
            type (str): One of 'all', 'input_files_newer', 'input_files_exist', 'repo'
            cumulative (bool): Whether to use cumulative validity
        """
        for st in self.stages:
            for f in self.exec_path[st].output_stage_files:
                if not f.output_is_valid(type=type, cumulative=cumulative):
                    f.output_file.touch()
                    self.exec_path[st].update()

    def remove_outputs(self, type: str = "all", cumulative: bool = False):
        """
        Remove output files that are invalid based on the specified type.
        Args:
            type (str): One of 'all', 'input_files_newer', 'input_files_exist', 'repo'
            cumulative (bool): Whether to use cumulative validity
        """
        for st in self.stages:
            for f in self.exec_path[st].output_stage_files:
                if not f.output_is_valid(type=type, cumulative=cumulative):
                    if f.output_file.is_file():
                        f.output_file.unlink()
                    self.exec_path[st].update()

    def remove_logs(self):
        """
        Remove stdout and stderr log files for all stages.
        """
        for st in self.stages:
            if (
                self.exec_path[st].stdout_log
                and self.exec_path[st].stdout_log.is_file()
            ):
                self.exec_path[st].stdout_log.unlink()
            if (
                self.exec_path[st].stderr_log
                and self.exec_path[st].stderr_log.is_file()
            ):
                self.exec_path[st].stderr_log.unlink()


# repo_timestamps_all = get_module_timestamps(b)
# nodes = b.get_execution_paths()[0]
# repo_timestamps = {node.stage_id: repo_timestamps_all[node.stage_id][node.module_id] for node in nodes}
# ep = ExecutionPath(nodes=nodes, repo_timestamps = repo_timestamps)
# ep.exec_path['data'].output_is_valid(type="all", cumulative=False)
# ep.exec_path['process'].output_is_valid(type="all", cumulative=False)
# ep.exec_path['methods'].output_is_valid(type="all", cumulative=False)
# ep.exec_path['metrics'].output_is_valid(type="all", cumulative=False)


class ExecutionPathSet:
    """
    A set of execution paths for a benchmark.
    """

    def __init__(self, b: Benchmark, cache_dir: Path = Path(".snakemake") / "repos"):
        """
        Initialize the ExecutionPathSet.
        Args:
            b (Benchmark): The benchmark object
        """
        self.b = b
        self.cache_dir = cache_dir
        self.modules_repo_timestamps = get_module_timestamps(
            b, cache_dir=self.cache_dir
        )
        self.exec_path_dict = self.get_exec_path_dict(b, self.modules_repo_timestamps)
        self.n_paths = len(self.exec_path_dict.keys())
        self.stages = list(b.get_stage_ids())

    def get_exec_path_dict(
        self, b: Benchmark, modules_repo_timestamps: Union[dict, None]
    ) -> dict:
        exec_paths = b.get_execution_paths()
        if modules_repo_timestamps is None:
            modules_repo_timestamps = get_module_timestamps(b, cache_dir=self.cache_dir)
        exec_path_dict = {}
        for exec_path_id in range(len(exec_paths)):
            nodes = b.get_execution_paths()[exec_path_id]
            repo_timestamps = {
                node.stage_id: modules_repo_timestamps[node.stage_id][node.module_id]
                for node in nodes
            }
            exec_path_dict[exec_path_id] = ExecutionPath(
                nodes=nodes, repo_timestamps=repo_timestamps
            )
        return exec_path_dict

    def create_filedict(self, cumulative: bool = False):
        """
        Get the file dictionary for the benchmark.
        Create a nested dictionary with stage ids as keys, node ids as sub-keys,
        and information about observed, missing, and invalid files as values.
        """
        self.filedict_is_cumulative = cumulative
        filedict = {
            st: {
                nd: {
                    "output_files": set(),
                    "valid_output_files": set(),
                    "observed_output_files": set(),
                    "empty_output_files": set(),
                    "missing_output_files": set(),
                    "invalid_output_files_input_file_is_newer": set(),
                    "invalid_output_files_repo_is_newer": set(),
                    "invalid_output_files": set(),
                    "n_valid": 0,
                    "n_observed": 0,
                    "n_empty": 0,
                    "n_missing": 0,
                    "n_invalid": 0,
                    "n_invalid_input_file_is_newer": 0,
                    "n_invalid_repo_is_newer": 0,
                    "n": 0,
                }
                for nd in self.b.get_nodes_by_stage_id(st)
            }
            for st in self.b.get_stage_ids()
        }
        for st in self.b.get_stage_ids():
            for exec_path_id in self.exec_path_dict.keys():
                if st in self.exec_path_dict[exec_path_id].stages:
                    de = self.exec_path_dict[exec_path_id].exec_path[st]
                    filedict[st][de.node]["output_files"] = filedict[st][de.node][
                        "output_files"
                    ].union(set(de.get_output_files("none", cumulative=cumulative)))
                    filedict[st][de.node]["valid_output_files"] = filedict[st][de.node][
                        "valid_output_files"
                    ].union(set(de.get_output_files("valid", cumulative=cumulative)))
                    filedict[st][de.node]["observed_output_files"] = filedict[st][
                        de.node
                    ]["observed_output_files"].union(
                        set(de.get_output_files("observed", cumulative=cumulative))
                    )
                    filedict[st][de.node]["empty_output_files"] = filedict[st][de.node][
                        "empty_output_files"
                    ].union(set(de.get_output_files("empty", cumulative=cumulative)))
                    filedict[st][de.node]["missing_output_files"] = filedict[st][
                        de.node
                    ]["missing_output_files"].union(
                        set(de.get_output_files("missing", cumulative=cumulative))
                    )
                    filedict[st][de.node]["invalid_output_files"] = filedict[st][
                        de.node
                    ]["invalid_output_files"].union(set(de.get_output_files("invalid")))
                    filedict[st][de.node][
                        "invalid_output_files_input_file_is_newer"
                    ] = filedict[st][de.node][
                        "invalid_output_files_input_file_is_newer"
                    ].union(
                        set(
                            de.get_output_files(
                                "input_files_newer", cumulative=cumulative
                            )
                        )
                    )
                    filedict[st][de.node]["invalid_output_files_repo_is_newer"] = (
                        filedict[
                            st
                        ][
                            de.node
                        ][
                            "invalid_output_files_repo_is_newer"
                        ].union(set(de.get_output_files("repo", cumulative=cumulative)))
                    )

        for st in self.b.get_stage_ids():
            for nd in self.b.get_nodes_by_stage_id(st):
                filedict[st][nd]["n"] = len(filedict[st][nd]["output_files"])
                filedict[st][nd]["n_valid"] = len(
                    filedict[st][nd]["valid_output_files"]
                )
                filedict[st][nd]["n_observed"] = len(
                    filedict[st][nd]["observed_output_files"]
                )
                filedict[st][nd]["n_empty"] = len(
                    filedict[st][nd]["empty_output_files"]
                )
                filedict[st][nd]["n_missing"] = len(
                    filedict[st][nd]["missing_output_files"]
                )
                filedict[st][nd]["n_invalid"] = len(
                    filedict[st][nd]["invalid_output_files"]
                )
                filedict[st][nd]["n_invalid_input_file_is_newer"] = len(
                    filedict[st][nd]["invalid_output_files_input_file_is_newer"]
                )
                filedict[st][nd]["n_invalid_repo_is_newer"] = len(
                    filedict[st][nd]["invalid_output_files_repo_is_newer"]
                )

        self.filedict = filedict

    def get_filedict(self, cumulative: bool = False) -> dict:
        """
        Get the file dictionary.
        """
        if not hasattr(self, "filedict"):
            self.create_filedict(cumulative=cumulative)
        if self.filedict_is_cumulative != cumulative:
            self.create_filedict(cumulative=cumulative)
        return self.filedict
