import os.path
from pathlib import Path
from abc import ABCMeta, abstractmethod

from omni_schema.datamodel.omni_schema import SoftwareBackendEnum

from omnibenchmark.benchmark import Benchmark, BenchmarkNode


class WorkflowEngine(metaclass=ABCMeta):
    """Interface for the workflow engine."""

    @abstractmethod
    def run_workflow(
        self,
        benchmark: Benchmark,
        cores: int = 1,
        update: bool = True,
        dryrun: bool = False,
        continue_on_error: bool = False,
        keep_module_logs: bool = False,
        backend: SoftwareBackendEnum = SoftwareBackendEnum.host,
        module_path: str = os.environ.get("MODULEPATH", None),
        work_dir: Path = Path(os.getcwd()),
        **kwargs,
    ) -> bool:
        """
        Serializes & runs benchmark workflow.

        Args:
            benchmark (Benchmark): benchmark to run
            cores (int): number of cores to run. Defaults to 1 core.
            update (bool): run workflow for non-existing outputs / changed nodes only. False means force running workflow from scratch. Default: True
            dryrun (bool): validate the workflow with the benchmark without actual execution. Default: False
            continue_on_error (bool): continue with independent jobs if a job fails. Default: False
            keep_module_logs (bool): keep module-specific log files after execution. Default: False
            backend (SoftwareBackendEnum): which software backend to use when running the workflow. Available: `host`, `docker`, `apptainer`, `conda`, `envmodules`. Default: `host`
            module_path (str): The path where the `envmodules` are located. This path will be searched during the workflow run using `envmodules` backend.
            work_dir (str): working directory. Default: current work directory
            **kwargs: keyword arguments to pass to the workflow engine

        Returns:
        - Status code (bool) of the workflow run.
        """
        raise NotImplementedError("Method not implemented yet")

    @abstractmethod
    def serialize_workflow(
        self, benchmark: Benchmark, output_dir: Path = Path(os.getcwd())
    ) -> Path:
        """
        Serializes a workflow file for the benchmark.

        Args:
            benchmark (Benchmark): benchmark to serialize
            output_dir (str): output directory for the workflow file

        Returns:
        - Workflow file path.
        """
        raise NotImplementedError("Method not implemented yet")

    @abstractmethod
    def run_node_workflow(
        self,
        node: BenchmarkNode,
        input_dir: Path,
        dataset: str,
        cores: int = 1,
        update: bool = True,
        dryrun: bool = False,
        continue_on_error: bool = False,
        keep_module_logs: bool = False,
        backend: SoftwareBackendEnum = SoftwareBackendEnum.host,
        module_path: str = os.environ.get("MODULEPATH", None),
        work_dir: Path = Path(os.getcwd()),
        **kwargs,
    ) -> bool:
        """
        Serializes & runs benchmark node workflow.

        Args:
            node (Benchmark): benchmark node to run
            input_dir (str): directory containing the inputs for the benchmark node
            dataset (str): file names corresponding to the dataset
            cores (int): number of cores to run. Defaults to 1 core.
            update (bool): run workflow for non-existing outputs / changed nodes only. False means force running workflow from scratch. Default: True
            dryrun (bool): validate the workflow with the benchmark without actual execution. Default: False
            continue_on_error (bool): continue with independent jobs if a job fails. Default: False
            keep_module_logs (bool): keep module-specific log files after execution. Default: False
            backend (SoftwareBackendEnum): which software backend to use when running the workflow. Available: `host`, `docker`, `apptainer`, `conda`, `envmodules`. Default: `host`
            module_path (str): The path where the `envmodules` are located. This path will be searched during the workflow run using `envmodules` backend.
            work_dir (str): working directory. Default: current work directory
            **kwargs: keyword arguments to pass to the workflow engine

        Returns:
        - Status code (bool) of the workflow run.
        """
        raise NotImplementedError("Method not implemented yet")

    @abstractmethod
    def serialize_node_workflow(
        self, node: BenchmarkNode, output_dir: Path = Path(os.getcwd())
    ) -> Path:
        """
        Serializes a workflow file for a benchmark node.

        Args:
            node (BenchmarkNode): benchmark node to serialize
            output_dir (str): output directory for the workflow file

        Returns:
        - Workflow file path.
        """
        raise NotImplementedError("Method not implemented yet")
