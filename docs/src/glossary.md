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

One run of one module: a module, one parameter set, and the upstream nodes
it reads from. A plan expands into nodes; each node has one directory in the
[output tree](#output-tree).

## Lineage

The chain of nodes a result descends from, from the first stage down to the
node that wrote it. A node's directory path spells out its lineage.

## Output tree

The directory `ob run` writes to (`out/` by default). Each node adds one
directory, `stage.module[.params]`, inside the directory of the node it read
from.

## Gather

A stage that runs one module over many upstream results at once, instead of
one at a time. With `group_by` it runs once per node of an upstream stage;
without, once for the whole run. Replaces the metric collector of older api
versions. See [Gather results](tutorial.md#gather-results).

## Member

A node whose output a [gather](#gather) collects.

## Manifest

The `lineage.json` file a gather module receives through `--from-files`. It
lists every [member](#member) with its files and its [lineage](#lineage).

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
