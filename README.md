## omnibenchmark

<p align="center">
<a href="https://github.com/omnibenchmark/omnibenchmark/tree/refs/heads/main"><img alt="tests result" src="https://github.com/omnibenchmark/omnibenchmark/actions/workflows/tests.yml/badge.svg?branch=main"></a>
<a href="https://github.com/omnibenchmark/omnibenchmark"><img alt="Coverage Status" src="./reports/coverage.svg"></a>
<a href="https://github.com/omnibenchmark/omnibenchmark/blob/main/LICENSE"><img alt="License: Apache 2.0" src="https://img.shields.io/badge/License-Apache_2.0-blue.svg"></a>
</p>

[![Linter: Ruff](https://img.shields.io/badge/Linter-Ruff-brightgreen?style=flat-square)](https://github.com/astral-sh/ruff)


[Omnibenchmark](https://omnibenchmark.org), a continuous benchmarking tool.

## Install

Recent versions of omnibenchmark should be installable [from pypi](https://pypi.org/project/omnibenchmark/):

```
pip install omnibenchmark
```

Do note that if you plan to use `conda` as a software execution backend, you will want to use omnibenchmark from within a conda environment manager. At the time of this writing, we recommend micromamba. See [the tutorial](https://omnibenchmark.org/tutorial/) for more details.

### Releases

See [our tags](https://github.com/omnibenchmark/omnibenchmark/tags).


### User documentation

[Live documentation](https://docs.omnibenchmark.org/latest) is avilable for every published release and the main branch.

Documentation includes a [CLI reference](https://docs.omnibenchmark.org/latest/reference/).

## Developer notes

Check [CONTRIBUTING.md](https://github.com/omnibenchmark/omnibenchmark/blob/main/CONTRIBUTING.md).

## Acknowledgements

Omnibenchmark incorporates great FOSS components, including but not limited to: [snakemake](https://snakemake.readthedocs.io/en/stable/), [easybuild](https://easybuild.io/), [apptainer](https://apptainer.org/), [lmod](https://lmod.readthedocs.io/en/latest/) and [git](https://git-scm.com/). Thank you!

## Pre-prints

- [Omnibenchmark (alpha) for continuous and open benchmarking in bioinformatics](https://arxiv.org/abs/2409.17038) (2024)
- [Building a continuous benchmarking ecosystem in bioinformatics](https://arxiv.org/abs/2409.15472) (2024)
