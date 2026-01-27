## omnibenchmark

<a href="https://github.com/omnibenchmark/omnibenchmark/actions/workflows/pixi.yml"><img alt="CI Pipeline" src="https://github.com/omnibenchmark/omnibenchmark/actions/workflows/pixi.yml/badge.svg?branch=main"></a>
<a href="https://pypi.org/project/omnibenchmark/"><img alt="PyPI" src="https://img.shields.io/pypi/v/omnibenchmark"></a>
<a href="https://codecov.io/gh/omnibenchmark/omnibenchmark"><img alt="codecov" src="https://codecov.io/gh/omnibenchmark/omnibenchmark/branch/main/graph/badge.svg"></a>
<a href="https://github.com/omnibenchmark/omnibenchmark/blob/main/LICENSE"><img alt="License: Apache 2.0" src="https://img.shields.io/badge/License-Apache_2.0-blue.svg"></a>
<a href="https://github.com/astral-sh/ruff"><img alt="Linter: Ruff" src="https://img.shields.io/badge/Linter-Ruff-brightgreen?style=flat-square"></a>


[Omnibenchmark](https://omnibenchmark.org), a continuous benchmarking tool.

## Install

Recent versions of omnibenchmark should be installable [from pypi](https://pypi.org/project/omnibenchmark/):

```
pip install omnibenchmark
```

Do note that if you plan to use `conda` as a software execution backend, you will want to use omnibenchmark from within a conda environment manager. At the time of this writing, we recommend miniforge. See [the tutorial](https://omnibenchmark.org/tutorial/) for more details.

### Releases

See [our tags](https://github.com/omnibenchmark/omnibenchmark/tags).


### User documentation

[Live documentation](https://docs.omnibenchmark.org/latest) is available for every published release and the main branch.

Documentation includes a [CLI reference](https://docs.omnibenchmark.org/latest/reference/).

### Configuration

Omnibenchmark uses a configuration file to store paths for Easybuild modules, dataset storage, and other settings. The configuration file is automatically created and managed by the system. See the [Configuration documentation](https://docs.omnibenchmark.org/latest/config/) for details.

## Developer notes

Check [CONTRIBUTING.md](https://github.com/omnibenchmark/omnibenchmark/blob/main/CONTRIBUTING.md).

## Acknowledgements

Omnibenchmark incorporates great FOSS components, including but not limited to: [snakemake](https://snakemake.readthedocs.io/en/stable/), [easybuild](https://easybuild.io/), [apptainer](https://apptainer.org/), [lmod](https://lmod.readthedocs.io/en/latest/) and [git](https://git-scm.com/). Thank you!

## Publications

- [Building a continuous benchmarking ecosystem in bioinformatics](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013658) (2025)
- [Omnibenchmark (alpha) for continuous and open benchmarking in bioinformatics](https://arxiv.org/abs/2409.17038) (2024)
