#!/usr/bin/env python
##
## Izaskun Mallona
## 22th July 2024

import pytest
import subprocess
import os

def run_snmk_conda(Snakefile, produced, expected):
    completed_process = subprocess.run(['snakemake', '-s', Snakefile, '--use-conda', '--cores', '1'])
    assert completed_process.returncode == 0
    with open(expected, 'r') as exp_fh, open(produced, 'r') as prod_fh:
        assert prod_fh.read()==exp_fh.read()


def run_snmk_apptainer(Snakefile, produced, expected):
    completed_process = subprocess.run(['snakemake', '-s', Snakefile, '--use-singularity', '--cores', '1'])
    assert completed_process.returncode == 0
    with open(expected, 'r') as exp_fh, open(produced, 'r') as prod_fh:
        assert prod_fh.read()==exp_fh.read()

def run(Snakefile, produced, expected, method):
    if method == 'conda':
        run_snmk_conda(Snakefile, produced, expected)
    elif method == 'apptainer':
        run_snmk_apptainer(Snakefile, produced, expected)
    else:
        print('not implemented')
    os.remove(produced)

