# Change Log

This document records all notable changes to `omnibenchmark`.
This project adheres to [Semantic Versioning](https://semver.org/).

## [0.1.0-rc.2](https://github.com/omnibenchmark/omnibenchmark/releases/tag/v0.1.0-rc.2)
- example: added a realistic `tests/data/Clustering.yaml` defining a clustering benchmark 
- tests: replaced `tests/data/Benchmark_001.yaml` by `tests/data/Clustering.yaml` during benchmark execution testing (e.g. envmodules)
- tests: added dummy envmodule (lua) files to `./tests/data/envs` so envmodules are not only searched at and loaded from the easybuild modulepath
- tests: fixed several aspects of envmodule validations so envmodules are tested by searching their availability and load-ability from the modulepath
- remote storage: added filtering objects by version
- workflow: updated the materialized `Snakefile` storage policy: switched to temporary file
- workflow: added benchmark layout exporting to a graph or plot via `mermaid`
- workflow: decoupled workflow rule setup software backends (e.g. adding `conda`, `envmodules` etc to rules) from execution (e.g. `--use conda` during the snmk syscall)

## [0.1.0-rc.1](https://github.com/omnibenchmark/omnibenchmark/releases/tag/v0.1.0-rc.1) 
- Snakemake-based workflow dynamic generation
- Snakemake execution: benchmarks or modules
- Local and remote (s3) storage handling
- Easybuild, conda, apptainer software capabilities
- click CLI
- Not back-compatible with older omnibenchmark versions
- Renku-free

## [0.0.48](https://pypi.org/project/omnibenchmark/0.0.48/)
- Renku-based omnibenchmark (last release)

## [0.0.1](https://pypi.org/project/omnibenchmark/0.0.1/)
- Renku-based omnibenchmark (first release)

