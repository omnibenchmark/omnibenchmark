import os
import shutil
from typing import List, Optional

from typer.testing import CliRunner, Result

from omni.cli.main import cli

FILES_TO_CLEANUP = [".snakemake", "out", "Snakefile", "snakemake.log"]


class OmniCLISetup:
    def __init__(self):
        self.runner = CliRunner()

    def call(self, args: List[str], input: Optional = None) -> Result:
        return self.runner.invoke(cli, args, input=input)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._cleanup_run()

    def _cleanup_run(self):
        current_dir = os.getcwd()
        for file in FILES_TO_CLEANUP:
            file_path = os.path.join(current_dir, file)
            if os.path.exists(file_path):
                if os.path.isfile(file_path):
                    os.remove(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
