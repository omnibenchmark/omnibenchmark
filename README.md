## omnibenchmark

<p align="center">
<a href="ttps://github.com/omnibenchmark/omni-py"><img alt="Tests" src="./reports/tests.svg"></a>
<a href="https://github.com/omnibenchmark/omni-py/actions"><img alt="Actions Status" src="https://github.com/omnibenchmark/omni-py/workflows/Tests/badge.svg"></a>
<a href="ttps://github.com/omnibenchmark/omni-py"><img alt="Coverage Status" src="./reports/coverage.svg"></a>
<a href="https://github.com/omnibenchmark/omni-py/blob/main/LICENSE"><img alt="License: Apache 2.0" src="https://img.shields.io/badge/License-Apache_2.0-blue.svg"></a>
<a href="https://github.com/psf/black"><img alt="Code style: black" src="https://img.shields.io/badge/code%20style-black-000000.svg"></a>
</p>

[Omnibenchmark](https://omnibenchmark.org), a continuous benchmarking tool.

## Install

Recent versions of omnibenchmark should be installable via pip:

```
pip install omnibenchmark
```

Do note that if you plan to use `conda` as a software execution backend, you will want to use omnibenchmark from within a conda environment manager. At the time of this writing, we recommend micromamba. See [the tutorial](https://omnibenchmark.org/tutorial/) for more details.

### Releases

See [our tags](https://github.com/omnibenchmark/omnibenchmark/tags).


### User documentation

[Live documentation](https://omnibenchmark.github.io/omnibenchmark/) is built for every change to the dev branch.

Documentation includes a CLI reference. This reference is automated [via mkdocs-click](https://gitlab.renkulab.io/omnibenchmark/omni_site/-/blob/master/docs/reference.md?ref_type=heads) and extracts the current CLI commands from the built branch (dev, for now).

## Developer notes

Check the `CONTRIBUTING.md` doc.

## Acknowledgements

Omnibenchmark incorporates great FOSS components, including but not limited to: [snakemake](https://snakemake.readthedocs.io/en/stable/), [easybuild](https://easybuild.io/), [apptainer](https://apptainer.org/), [lmod](https://lmod.readthedocs.io/en/latest/) and [git](https://git-scm.com/). Thank you!

## Pre-prints

- [Omnibenchmark (alpha) for continuous and open benchmarking in bioinformatics](https://arxiv.org/abs/2409.17038) (2024)
- [Building a continuous benchmarking ecosystem in bioinformatics](https://arxiv.org/abs/2409.15472) (2024)
