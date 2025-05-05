import subprocess
from typing import List, Optional


class OmniCLISetup:
    def call(
        self, args: List[str], input: Optional[str] = None, cwd: Optional[str] = None
    ) -> subprocess.CompletedProcess:
        """
        Call the CLI using subprocess.run.
        Args:
            args: List of arguments to pass to the CLI.
            input: Optional input to pass to the CLI's stdin.
        Returns:
            A subprocess.CompletedProcess object containing the result of the command.
        """
        # Construct the command to run
        command = ["python", "-m", "omnibenchmark.cli.main"] + args

        # Run the command using subprocess.run
        result = subprocess.run(
            command,
            input=input,
            text=True,  # Ensures input/output are treated as strings
            capture_output=True,  # Captures stdout and stderr
            cwd=cwd,  # Change working directory
        )

        return result

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass
