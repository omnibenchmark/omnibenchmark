## Izaskun Mallona
## 22th July 2024

import difflib
import subprocess
import os


def check_cmd_zero_exit(cmd_string):
    completed_process = subprocess.run(cmd_string.split(" "))
    assert completed_process.returncode == 0


def run_snmk_conda(Snakefile, produced, expected):
    completed_process = subprocess.run(
        [
            "snakemake",
            "-s",
            Snakefile,
            "--software-deployment-method",
            "conda",
            "--cores",
            "1",
        ]
    )
    assert completed_process.returncode == 0
    with open(expected, "r") as exp_fh, open(produced, "r") as prod_fh:
        expected_content = exp_fh.read()
        produced_content = prod_fh.read()

        if expected_content != produced_content:
            # Generate a unified diff format
            diff = difflib.unified_diff(
                expected_content.splitlines(),
                produced_content.splitlines(),
                fromfile="expected",
                tofile="produced",
                lineterm="",
            )
            diff_text = "\n".join(diff)
            assert False, f"Files differ:\n{diff_text}"
    # with open(expected, "r") as exp_fh, open(produced, "r") as prod_fh:
    #    assert cleaned(prod_fh.read()) == cleaned(exp_fh.read())


def run_snmk_apptainer(Snakefile, produced, expected):
    completed_process = subprocess.run(
        [
            "snakemake",
            "-s",
            Snakefile,
            "--software-deployment-method",
            "apptainer",
            "--cores",
            "1",
        ]
    )
    assert completed_process.returncode == 0
    with open(expected, "r") as exp_fh, open(produced, "r") as prod_fh:
        assert prod_fh.read() == exp_fh.read()


def run_snmk_envmodules(Snakefile, produced, expected):
    completed_process = subprocess.run(
        ["snakemake", "-s", Snakefile, "--use-envmodules", "--cores", "1"]
    )
    assert completed_process.returncode == 0
    with open(expected, "r") as exp_fh, open(produced, "r") as prod_fh:
        assert prod_fh.read() == exp_fh.read()


def run(Snakefile, produced, expected, method):
    if method == "conda":
        run_snmk_conda(Snakefile, produced, expected)
    elif method == "apptainer":
        run_snmk_apptainer(Snakefile, produced, expected)
    elif method == "envmodules":
        run_snmk_envmodules(Snakefile, produced, expected)
    else:
        print("not implemented")
    os.remove(produced)
