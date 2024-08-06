import os.path
from pathlib import Path

from omni.benchmark import Benchmark, BenchmarkNode


class WorkflowEngine:
    """Interface for the workflow engine."""

    def run_workflow(
        self,
        benchmark: Benchmark,
        cores: int = 1,
        update: bool = True,
        dryrun: bool = False,
        work_dir: Path = Path(os.getcwd()),
        **kwargs
    ) -> bool:
        """
        Serializes & runs benchmark workflow.

        Args:
            benchmark (Benchmark): benchmark to run
            cores (int): number of cores to run. Defaults to 1 core.
            update (bool): run workflow for non-existing outputs / changed nodes only. False means force running workflow from scratch. Default: True
            dryrun (bool): validate the workflow with the benchmark without actual execution. Default: False
            work_dir (str): working directory. Default: current work directory
            **kwargs: keyword arguments to pass to the workflow engine

        Returns:
        - Status code (bool) of the workflow run.
        """
        raise NotImplementedError("Method not implemented yet")

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

    def run_node_workflow(
        self,
        node: BenchmarkNode,
        input_dir: Path,
        dataset: str,
        cores: int = 1,
        update: bool = True,
        dryrun: bool = False,
        work_dir: Path = Path(os.getcwd()),
        **kwargs
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
            work_dir (str): working directory. Default: current work directory
            **kwargs: keyword arguments to pass to the workflow engine

        Returns:
        - Status code (bool) of the workflow run.
        """
        raise NotImplementedError("Method not implemented yet")

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
