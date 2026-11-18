# Glossary

The terms the documentation uses, as they apply from api 0.8.0.

## Benchmark

A study comparing methods on datasets with known outcomes, using shared
metrics. In omnibenchmark it is a [plan](#plan), the [modules](#module) the
plan points to, and the results of running them.

## Plan

The YAML file describing a benchmark: its identity, its
[software environments](#software-environment) and its [stages](#stage).
Also called the benchmark YAML. See the
[spec reference](spec-reference.md).

## API version

The version of the plan format, declared as `api_version`. It fixes the
fields a plan may use, the command line [modules](#module) receive and the
layout of the [output tree](#output-tree). See [API versions](api-versions.md).

## Stage

A step of the benchmark: a group of [modules](#module) that take the same
[inputs](#input) and write the same [outputs](#output), so any of them can
feed the next stage. Typical stages are data, methods and metrics.

## Module

One component of a stage (a dataset, a method, a metric): a git repository
pinned to a commit, with an [entrypoint](#entrypoint) omnibenchmark calls.
The same repository can back several modules.

## Entrypoint

The script a module runs, named in its `omnibenchmark.yaml`. It receives the
paths to read and write, and the module's [parameters](#parameter), as
command-line flags.

## Output

A file every module of a stage writes. A stage declares each one with an
`id`, which later stages refer to, and a `path`, the file name. A module
receives it as `--output id=PATH`.

## Input

An [output](#output) id of an earlier stage that a stage consumes. A module
receives it as `--<output id> PATH`.

## Parameter

A value passed to a module as `--key value`. A module with several parameter
sets runs once per set.

## Node

One unit of execution: a module, one parameter set, and the upstream nodes
it reads from. A plan expands into nodes; each node has one directory in the
[output tree](#output-tree).

## Lineage

The chain of nodes a result descends from, from the first stage down to the
node that wrote it.

- For an ordinary node, the directory path spells it out:
  `out/data.D1/methods.M2.9da5174a/` was computed by `M2` from `D1`.
- A [gather](#gather) collects results from many lineages, so its path cannot
  hold them. They are written to `lineage.json` in the gather's directory,
  one entry per collected node; the gather's module receives that file as
  `--from-files`.

## Label

A named value a stage attaches to its nodes with `provides`, carried down the
[lineage](#lineage). A later module with `requires` runs only where the label
has the value it asks for. See
[Run a module only on matching data](howto.md#run-a-module-only-on-matching-data).

## Output tree

The directory `ob run` writes to (`out/` by default), holding one directory
per [node](#node). How those directories are arranged is the
[layout](#layout).

## Layout

How node directories are arranged in the [output tree](#output-tree). The
[api version](#api-version) of the plan selects it.

- **Flat** (api 0.8.0): each node adds one directory,
  `stage.module[.params]`, inside the directory of the node it read from:
  `out/data.D1/methods.M2.9da5174a/`. No hidden directories or symlinks;
  `out/by-params/` is a derived view with readable parameter names.
- **Nested** (below api 0.8.0): each node adds three,
  `stage/module/.params`: `out/data/D1/.default/methods/M2/.9da5174a/`.
  **Deprecated**: `ob run` warns on such plans, and the nested layout is
  removed at api 0.10.0.

Moving a plan to the flat layout changes every output path, so the next run
recomputes everything. See
[Upgrade a benchmark to api 0.8.0](howto.md#output-layout).

## Run

One execution of `ob run` on an output tree. Each run gets a `run_id` (32
hex characters) and is recorded as one line of `out/.metadata/runs.jsonl`,
with the benchmark's identity and what was selected. With telemetry on, the
`run_id` is the trace id of the run's spans. See
[Keep an output directory tied to its benchmark](howto.md#keep-an-output-directory-tied-to-its-benchmark).

## Gather

A stage that runs one module over many upstream results at once, instead of
one at a time. With `group_by` it runs once per node of an upstream stage;
without, once for the whole run. Its module reads what was collected from
`lineage.json` (see [lineage](#lineage)). Replaces the metric collector of
older api versions. See [Gather results](tutorial.md#gather-results).

## Validator

A script, declared in the plan on a stage [output](#output), that checks one
produced file and exits 0 if it passes. Run with `ob validate outputs`. See
[Check outputs with validators](howto.md#check-outputs-with-validators).

## Software environment

A named, reproducible set of software a module runs in: a conda environment
file, an Apptainer image or an environment module. Declared in the plan under
`software_environments` and chosen per module.

## Software backend

How a benchmark's software environments are provided: `conda`, `apptainer`,
`docker`, `envmodules` or `host`. One per plan. See
[Software backends](howto.md#software-backends).

## obcommons

The community registry of omnibenchmark modules, at
[omnibenchmark.github.io/obcommons](https://omnibenchmark.github.io/obcommons/).
Each entry gives a module's repository and pinned commit, the stage it fits,
its inputs and outputs, and the api version it needs. Look there before
writing a module; add yours with a pull request to
[omnibenchmark/obcommons](https://github.com/omnibenchmark/obcommons).

## Benchmarker

The person who designs and maintains a benchmark: defines the stages, the
file formats between them and the validators.

## Contributor

Someone who adds a module (a dataset, method or metric) to a benchmark,
following the stage contract the benchmarker set.
