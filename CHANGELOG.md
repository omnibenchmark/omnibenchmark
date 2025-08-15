# Changelog

This document records all notable changes to `omnibenchmark`.
This project uses [Semantic Versioning](https://semver.org/) and [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/).

## [0.3.0](https://github.com/omnibenchmark/omnibenchmark/releases/tag/v0.3.0) (Aug 15th 2025)

- feat: SLURM support (#147)
- feat(cli)!: --local argument has been renamed to --local-storage
- feat: Support passing of extra arguments from CLI -> Snakemake for `run` commands
- feat: add extra profiler to the snakemake execution (#151)
- fix: Restrict yaml output paths and hide implementation details (#138)
- refactor: remove the software command
- tests: disable snakemake-only tests that only exercise environemtes, under rework
- chore(pkg): use setuptools for builds instead of poetry

## [0.2.2](https://github.com/omnibenchmark/omnibenchmark/releases/tag/v0.2.2) (Aug 1st 2025)

- feat(cli): re-add temporarily disabled --task-timeout (#132)
- feat(cli): added `--out-dir` to `ob run benchmark` to configure a different output directory for the workflow
- feat: add general config file (#148)
- docs: clarify documentation on install methods (#144)
- docs: add design docs template

## [0.2.1](https://github.com/omnibenchmark/omnibenchmark/releases/tag/v0.2.1) (May 20th 2025)

- bugfix for `--continue-on-error` flag.

## [0.2.0](https://github.com/omnibenchmark/omnibenchmark/releases/tag/v0.2.0) (May 20th 2025)
- workflow parameters: Fix parameter instability, use human-friendly symlinks to reference stable, hashed parameter folders.
- workflow: added a metric collector capability linked to new YAML specs as defined in omni-schema version `94c312957f96d218369e5e1bcf7abaf976ed0fbc`
- cli: added `-k` or `--continue-on-error` to `ob run benchmark` to keep going on (rule) error
- storage: added archival capabilities
- fix: fixed race conditions at the module repository cloning step
- refactor: rename --threads parameter to --cores
- chore: make s3 dependencies optional
- tests: refactor CLI capture and execution tests
- tests: add bundled repositories for self-contained tests

### Known Bugs

- [#126](https://github.com/omnibenchmark/omnibenchmark/issues/126) `--continue-on-error` not working as intended. Patch release will follow.
- `--task-timeout` is disabled because it had effect on Snakemake's parellelization abilities; further investigation is needed before re-enabling it.

## [0.1.0](https://github.com/omnibenchmark/omnibenchmark/releases/tag/v0.1.0)
- cli: added dynamic versioning managed by poetry-dynamic-versioning , using a semvar style.
- cli: added `ob --version` which outputs the current version of the cli
- logging: added uniform logging + configuration and management across the whole codebase, no longer using click.echo, rather logger.log
- ci/cd: added a black pre-commit hook
- logging and workflow: for each rule execution, subprocess.stdout and subprocess.stderr, if not empty, are kept as logs next to the rule output(s)
- software: fixes path checks on `ob software conda prepare -b [yaml]`
- workflow: modifies rule generation to populate software (conda, envmodules, container) directives defensively to align with snakemake>8.25.2


## [0.1.0-rc.2](https://github.com/omnibenchmark/omnibenchmark/releases/tag/v0.1.0-rc.2)
- example: added a realistic `tests/data/Clustering.yaml` defining a clustering benchmark
- tests: replaced `tests/data/Benchmark_001.yaml` by `tests/data/Clustering.yaml` during benchmark execution testing (e.g. envmodules)
- tests: added dummy envmodule (lua) files to `./tests/data/envs` so envmodules are not only searched at and loaded from the easybuild modulepath
- tests: fixed several aspects of envmodule validations so envmodules are tested by searching their availability and load-ability from the modulepath
- remote storage: added filtering objects by version
- workflow: updated the materialized `Snakefile` storage policy: switched to temporary file
- workflow: added benchmark layout exporting to a graph or plot via `mermaid`
- workflow: decoupled workflow rule setup software backends (e.g. adding `conda`, `envmodules` etc to rules) from execution (e.g. `--use conda` during the snmk syscall)

## [0.1.0rc.1](https://github.com/omnibenchmark/omnibenchmark/releases/tag/v0.1.0-rc.1)
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
