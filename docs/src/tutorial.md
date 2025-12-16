## Design a benchmark YAML

Benchmark specification files are written in YAML. They specify the formal
dependencies of benchmark components, as well as some metadata (e.g. the
repository containing their implementation, parameters to be run with, etc.).

Let's construct a simple example benchmark, shaped as follows:

- `D1` is a single (starting) dataset. (In real life, datasets will have meaningful
  names, e.g. `semisimulation_smith_2019`).
- `M1` and `M2` are methods. They process the dataset `D1` directly.
  (Similarly, methods would also have proper names, e.g. `limma` or `linreg`).
- `m1` and `m2` are metrics. They process the output of the methods `M1` and `M2`
  directly. (Again, naming is flexible, we're keeping it short for clarity.)

```mermaid
flowchart LR
  subgraph data
    D1
  end
  subgraph methods
    D1 --> M1
    D1 --> M2
  end
  subgraph metrics
    M1 --> m1
    M1 --> m2
    M2 --> m1
    M2 --> m2
  end
```

Benchmark specification files have a header and a body.

### Benchmark YAML header

Let's start with the header.

```yaml
---
## benchmark shortname
id: bench1

## benchmark description
description: a simple benchmark

## Benchmark version. `1.0`. This is our first attempt, so let's call it major version 1, minor version 0: `1.0`.
version: "1.0"

## Benchmark builder/contact person
benchmarker: "Mary the Benchmarker, mary@uzh.ch"

## Storage flavour for sharing results: currently only S3
storage_api: S3

## S3 endpoint to share our benchmark results.
##   `https://s3_object_storage_url.ch` does not exist, and we don't mind -
##   not sharing results from our benchmark yet.
storage: https://s3_object_storage_url.ch

## Benchmark YAML schema/specification version.
benchmark_yaml_spec: "0.0.1"

## The software backend used to run the benchmark
software_backend: apptainer

## Software environment recipes associated to this benchmark.
##  Suffice to say they are apptainer images in some ORAS-compatible registry.
software_environments:
  R:
    description: "R 4.3.3 with gfbf-2023 toolchain"
    apptainer: http://registry.ch/R_4.3.3-gfbf-2023b.sif
  python:
    description: "Python3.12.0 with gfbf-2023 toolchain"
    apptainer: http://registry.ch/python_vX-gfbf-2023b.sif
```

Hence, the header acts as a preamble, defining general attributes of the benchmark. The body contains the individual benchmark components (methods, metrics, etc) and their linking to each other.

### Benchmark YAML body

The benchmark body is structured in stages grouping benchmarking components that produce similarly shaped outputs and ingest similarly shaped inputs. That is:

```mermaid
flowchart LR
    classDef thing fill:#f96
    D1 -- produces --> image
    image:::thing -- "is ingested by" --> M1
    image:::thing -- "is ingested by" --> M2
    M1 -- produces --> matrix1
    M2 -- produces --> matrix2
    matrix1:::thing  -- "is ingested by" --> m1
    matrix1:::thing  -- "is ingested by" --> m2
    matrix2:::thing  -- "is ingested by" --> m1
    matrix2:::thing  -- "is ingested by" --> m2
```

In this example, `matrix1` and `matrix2` are similarly shaped, e.g. might be tab-separated files with some constraints, such as having a header and a rownames; and different from `image`, which might be a raster image in PNG format. We require `D1` to be part of a stage where modules *produce images, and ingest no inputs*; `M1` and `M2` to belong to a stage of *image-ingesting, matrix-producing* modules; and `m1` and `m2` to be part of a last stage of *matrix-ingesting* modules.

Let's start with the first stage, containing `D1`. We will call it `data` (naming is flexible).

```yaml
stages:
    ## the stage name
  - id: data
    ##  here we have a single data stage with one single module,
    ##  that outputs a single shape for the `data` stage.
    modules:
        ## unique module id
      - id: D1
        ## module name in a longer form
        name: "Dataset 1"
        ## software environment to run this module; maps to the header `software_environments`
        software_environment: "python"
        ## the git-compatible remote, and a particular pinned commit
        repository:
          url: https://github.com/omnibenchmark-example/data.git
          commit: 41aaa0a
    ## output file paths for this stage members. In this simple case, the output from D1.
    outputs:
        ## output id
      - id: data.image
        ## output path.
        path: "{dataset}.png"
```

Let's add another stage, for the modules `M1` and `M2`.
This stage is not initial: its modules have both inputs and outputs.

```yaml
    ## the stage name
  - id: methods
    ## a list of modules and their repositories, as above
    modules:
      - id: M1
        software_environment: "python"
        repository:
          url: https://github.com/omnibenchmark-example/method.git
          commit: 1004cdd
      - id: M2
        ## notice this method runs in a container offering some R capabilities
        software_environment: "R"
        repository:
          url: https://github.com/omnibenchmark-example/method2.git
          commit: 10sg4cdd
    ## input identifiers, refering to the `data stage` outputs
    inputs:
      - entries: data.image
    ## stage-specific outputs
    outputs:
      - id: methods.matrix
        ## output path.
        path: "{dataset}.matrix.tsv.gz"
```

You might be wondering: what does the wildcard `{dataset}` mean? It represents the dataset identifier that will be resolved at runtime. When writing the YAML, you only need to specify the filename template using simple wildcards like `{dataset}`. Omnibenchmark will automatically organize files in a hierarchical directory structure for provenance tracking. As a consequence, running module `D1` will generate files at:

```
./data/D1/default/D1.png
```

Hence, running modules `M1` and `M2` will produce files templated as `{dataset}.matrix.tsv.gz`, which, given there is only one dataset `D1` available, will result in:

```
./data/D1/default/methods/M1/default/D1.matrix.tsv.gz
./data/D1/default/methods/M2/default/D1.matrix.tsv.gz
```

Finally, we add the metrics stage containing modules `m1` and `m2`.

```yaml
    ## the stage name
  - id: metrics
    ## a list of modules and their repositories, as above
    modules:
      - id: m1
        software_environment: "python"
        repository:
          url: https://github.com/omnibenchmark-example/metric.git
          commit: 4504cdd
      - id: m2
        software_environment: "R"
        repository:
          url: https://github.com/omnibenchmark-example/metric2.git
          commit: 7sg4cdd
    ## input identifiers, refering to the `data stage` outputs
    inputs:
      - entries: methods.matrix
    ## stage specific-outputs
    outputs:
      - id: metrics.json
        ##  wildcard dataset: `D1` (here datasets refer to the initial stage above, not to the module name)
        path: "{dataset}.json"
```

Hence, running modules `m1` and `m2` will produce files templated as `{dataset}.json`; given there is only one dataset `D1` and two methods `M1` and `M2` available, will result in the following outputs:

```
./data/D1/default/methods/M1/default/metrics/m1/default/D1.json
./data/D1/default/methods/M2/default/metrics/m1/default/D1.json
./data/D1/default/methods/M1/default/metrics/m2/default/D1.json
./data/D1/default/methods/M2/default/metrics/m2/default/D1.json
```

The full benchmark YAML looks like this:

```yaml
---
## benchmark shortname
id: bench1

## benchmark description
description: a simple benchmark

## Benchmark version. `1.0`. This is our first attempt, so let's call it major version 1, minor version 0: `1.0`.
version: "1.0"

## Benchmark builder/contact person
benchmarker: "Mary the Benchmarker, mary@uzh.ch"

## Storage flavour for sharing results: currently only S3
storage_api: S3

## S3 endpoint to share our benchmark results.
##   `https://s3_object_storage_url.ch` does not exist, and we don't mind -
##   not sharing results from our benchmark yet.
storage: https://s3_object_storage_url.ch

## Benchmark YAML schema/specification version.
benchmark_yaml_spec: "0.0.1"

## License
# license: MIT # not yet part of the schema

software_backend: apptainer

## Software environment recipes associated to this benchmark.
## Suffice to say they are apptainer images in some ORAS-compatible registry.
software_environments:
  R:
    description: "R 4.3.3 with gfbf-2023 toolchain"
    apptainer: http://registry.ch/R_4.3.3-gfbf-2023b.sif
  python:
    description: "Python3.12.0 with gfbf-2023 toolchain"
    apptainer: http://registry.ch/python_vX-gfbf-2023b.sif

stages:
  - id: data
    modules:
      - id: D1
        name: "Dataset 1"
        software_environment: "python"
        repository:
          url: https://github.com/omnibenchmark-example/data.git
          commit: 41aaa0a
    outputs:
      - id: data.image
        path: "{dataset}.png"

  - id: methods
    modules:
      - id: M1
        software_environment: "python"
        repository:
          url: https://github.com/omnibenchmark-example/method.git
          commit: 1004cdd
      - id: M2
        software_environment: "R"
        repository:
          url: https://github.com/omnibenchmark-example/method2.git
          commit: 10sg4cdd
    inputs:
      - entries: data.image
    outputs:
      - id: methods.matrix
        path: "{dataset}.matrix.tsv.gz"

  - id: metrics
    modules:
      - id: m1
        software_environment: "python"
        repository:
          url: https://github.com/omnibenchmark-example/metric.git
          commit: 4504cdd
      - id: m2
        software_environment: "R"
        repository:
          url: https://github.com/omnibenchmark-example/metric2.git
          commit: 7sg4cdd
    inputs:
      - entries: methods.matrix
    outputs:
      - id: metrics.json
        path: "{dataset}.json"
```

### Metric collectors

The YAML stanzas above aim to scaffold a workflow by nesting inputs and outputs. Omnibenchmark automatically organizes files in a hierarchical directory structure where each module's outputs are stored in dedicated directories for provenance tracking. These files can be further processed by other modules, creating linear lineages with implicit provenance traceable by browsing the parent folder(s) of any folder and file. This can pose a challenge if multiple files (lineages) are meant to be gathered by a processing step.

An independent syntax allows collecting _multiple inputs across multiple folders and lineages_ to process them jointly. This usecase is typically needed when collecting metrics, that is, gathering all output files from some stage(s) to build a final aggregated report. Graphically, collection means adding the rightmost step  (`is collected by`) to the benchmarking workflow to produce `c1` (again, naming is flexible):

```mermaid
flowchart LR
    classDef thing fill:#f96
    D1 -- produces --> image
    image:::thing -- "is ingested by" --> M1
    image:::thing -- "is ingested by" --> M2
    M1 -- produces --> matrix1
    M2 -- produces --> matrix2
    matrix1:::thing  -- "is ingested by" --> m1
    matrix1:::thing  -- "is ingested by" --> m2
    matrix2:::thing  -- "is ingested by" --> m1
    matrix2:::thing  -- "is ingested by" --> m2
    m1 -- produces --> M1_m1
    m1 -- produces --> M2_m1
    m2 -- produces --> M1_m2
    m2 -- produces --> M2_m2
    M1_m1:::thing -- "is collected by\n (metric collector)" --> c1
    M1_m2:::thing -- "is collected by\n (metric collector)" --> c1
    M2_m1:::thing -- "is collected by\n (metric collector)" --> c1
    M2_m2:::thing -- "is collected by\n (metric collector)" --> c1
    c1 -- "renders" --> report
    report:::thing
```

The `is collected by` capability is specified within the benchmarking header (that is, before the `stages` stanzas) as a member of `metric_collectors` enumeration. (So multiple metric collectors can exist, if more than one stanza are added.) To specify a single metric collector reading all `metrics.json` outputs (while tracking their lineages, that is, their original dataset, methods, metric, parameters, etc):

```yaml
metric_collectors:
  - id: multiple_to_one
    name: "Single-backend (multiple) method outputs collector."
    software_environment: "python"
    repository:
      url: https://github.com/imallona/clustering_report
      commit: f1a5876
    inputs:
      - metrics.json
    outputs:
      - id: plotting.html
        path: "{input}/{name}/plotting_report.html"
```

Similarly to any other module, the associated code to run the processing is stored as a git-compatible remote (with `url` and `commit id`). In the example above, `multiple_to_one` only generates one output, a report named `plotting.html` collecting all computed `metrics.json` files.

## Validate a benchmark YAML

Let's save the benchmark above as a file named `benchmark_test.yaml`. Then we validate it with:

=== "Shell"

    ```shell
    ob run validate -b benchmark_test.yaml
    ```

=== "Output"


    ```
    Validating a benchmark yaml.
    Benchmark YAML file integrity check passed.
    ```


## Create a module suitable to be used in omnibenchmark

Any accesible git repository can host an omnibenchmark module. If it's convenient, you might want to push them to a remote in GitHub, Bitbucket, GitLab, etc. In reality, `omnibenchmark` just needs to be able to access the remote (clone or fetch), and be able to checkout your specified commit (so anything that works for your global git config should work for `omnibenchmark`).

We provide an example set of modules for the benchmark example file at [`tests/data/Benchmark_001.yaml`](https://github.com/omnibenchmark/omnibenchmark/blob/main/tests/data/Benchmark_001.yaml).

As shown below, module D1 points to the GitHub repository [example data](https://github.com/omnibenchmark-example/data.git) at the commit `63b7b36`. (Incidentally, in theory you should also be able to specify any valid dynamic git reference, like `HEAD` or a `tag` or `branch` name).

```yaml
stages:
  - id: data
    modules:
      - id: D1
        name: "Dataset 1"
        software_environment: "python"
        repository:
          url: https://github.com/omnibenchmark-example/data.git
          commit: 63b7b36
    outputs:
        ## output id
      - id: data.image
        ##   wildcard dataset: `D1` (module ids in initial stages - that is, the ones not ingesting inputs and only
        ##     generating outputs, are reused as `dataset` wildcards)
        path: "{dataset}.png"
```

Hence, the git repository implementing module `D1` doesn't have any input, but generates one output. In this case, the repository implementing `D1` has [an omnibenchmark.yaml file](https://github.com/omnibenchmark-example/data/blob/main/omnibenchmark.yaml) indicating the entrypoint is a python script named `entrypoint_data.py`:

```yaml
entrypoints:
  default: entrypoint_data.py
```

`entrypoint_data.py` uses the python library `argparse` to receive two arguments when called from the command line:

```python
parser.add_argument('--output_dir', type=str, help='output directory where dataset files will be saved.'))
parser.add_argument('--name', type=str, help='name of the dataset')
```

That is, the output directory where the `data.image` output is generated, and the dataset (`D1`) name.

Argument parsing aside, the `entrypoint_data.py` script structure is free: in this case, it materializes files with a dummy content.

Let's inspect another module, this time running in R and also receiving inputs.

```yaml
stages:
  - id: some_intermediate_step
    modules:
      - id: process
        exclude: [select_counts]
        software_environment: "R"
        repository:
          url: https://github.com/omnibenchmark-example/process.git
          commit: aeec1db
    inputs:
      - entries:
          - data.meta
          - data.counts
    outputs:
      - id: select_lognorm.selected
        path: "{dataset}.txt.gz"
```

So, in this case, the module `process` is likely to be implemented in R, receive three inputs, and produce one output. A dummy implementation is available at [https://github.com/omnibenchmark-example/process.git](https://github.com/omnibenchmark-example/process.git). There, the [omnibenchmark.yaml file](https://github.com/omnibenchmark-example/process/blob/main/omnibenchmark.yaml) indicates:

```yaml
entrypoints:
  default: entrypoint_process.R
```

so the script to be executed is named `entrypoint_process.R`. In this case, [the script](https://github.com/omnibenchmark-example/process/blob/aeec1db790542d447899d6ac4cb8564a9172b6e0/entrypoint_process.R#L3C1-L13C1) uses the R library `argparse` to provide a commandline interface:

```R
# Define argument parser
parser <- ArgumentParser(description="Process dataset files")

# Add arguments
parser$add_argument("--output_dir", "-o", dest="output_dir", type="character", help="output directory where files will be saved")
parser$add_argument("--name", "-n", dest="name", type="character", help="name of the dataset")
parser$add_argument("--data.counts", dest="data_counts", type="character", help="input file #1")
parser$add_argument("--data.meta", dest="data_meta", type="character", help="input file #2")
```

Notice these argument names **must match** the YAML's input ids: `data.counts` and `data.meta` are specified as inputs in the benchmark YAML; as before, `name` refers to the dataset name and `output_dir` to the path where outputs will be generated. As before, the script is free in structure - it implements some functionality, and can import other scripts as well, as long as it reads inputs and write outputs in a way compatible to the benchmark YAML specification.

## Run a benchmark

The benchmark [`tests/data/Benchmark_001.yaml`](https://github.com/omnibenchmark/omnibenchmark/blob/main/tests/data/Benchmark_001.yaml) above is a complex benchmark - but it runs quick enough. Let's try a dry run and inspect the rules that will be run:

=== "Shell"

    ```shell
    ob run tests/data/Benchmark_001.yaml  --cores 1 --dry
    ```

=== "Output"

    ```
    [snip]

    INFO:snakemake.logging:
    Job stats:
    job                   count
    ------------------  -------
    all                       1
    data_D1_default           1
    data_D2_default           1
    methods_M1_default        5
    methods_M2_param_0        5
    methods_M2_param_1        5
    metrics_m1_default       15
    metrics_m2_default       15
    metrics_m3_default       15
    process_P1_param_0        2
    process_P1_param_1        2
    process_P2_param_0        2
    process_P2_param_1        2
    total                    71

    [snip]
    ```
So it plans to run 71 jobs in total. Its methods are fast, so we can run it (it will take less than two minutes in most machines):

=== "Shell"

    ```shell
    ob run tests/data/Benchmark_001.yaml  --cores 1
    ```

=== "Output"

    ```
    [snip]

    resources: tmpdir=/home/imallona/tmp/eb-ge9tbg43

    INFO:snakemake.logging:
    [Thu Aug 29 10:23:23 2024]
    INFO:snakemake.logging:[Thu Aug 29 10:23:23 2024]
    Finished job 0.
    INFO:snakemake.logging:Finished job 0.
    71 of 71 steps (100%) done
    INFO:snakemake.logging:71 of 71 steps (100%) done
    Complete log: .snakemake/log/2024-08-29T102204.875104.snakemake.log
    WARNING:snakemake.logging:Complete log: .snakemake/log/2024-08-29T102204.875104.snakemake.log
    Benchmark run has finished successfully.
    ```

## Run an initial module

The benchmark [`tests/data/Benchmark_001.yaml`](https://github.com/omnibenchmark/omnibenchmark/blob/main/tests/data/Benchmark_001.yaml) contains several initial steps which generate datasets and don't receive any input. To run these, we have to use the `ob run [benchmark.yaml] --module [MODULE_ID]` verb.

=== "Shell"

    ```shell
    ob run tests/data/Benchmark_001.yaml --module D1
    ```

=== "Output"

    ```
    Running module on a dataset provided in a custom directory.
    Benchmark YAML file integrity check passed.
    Found 1 workflow nodes for module D1.
    Running module benchmark...
    Assuming unrestricted shared filesystem usage.
    WARNING:snakemake.logging:Assuming unrestricted shared filesystem usage.
    Building DAG of jobs...
    WARNING:snakemake.logging:Building DAG of jobs...
    Using shell: /usr/bin/bash
    WARNING:snakemake.logging:Using shell: /usr/bin/bash
    Provided cores: 1 (use --cores to define parallelism)
    WARNING:snakemake.logging:Provided cores: 1 (use --cores to define parallelism)
    Rules claiming more threads will be scaled down.
    WARNING:snakemake.logging:Rules claiming more threads will be scaled down.
    Job stats:
    job                count
    ---------------  -------
    all                    1
    data_D1_default        1
    total                  2

    reason: Missing output files: out/data/D1/default/D1.txt.gz, out/data/D1/default/D1_params.txt, out/data/D1/default/D1.meta.json; Code has changed since last execution
    resources: tmpdir=/home/imallona/tmp/eb-w7lf3kqk
    INFO:snakemake.logging:localrule data_D1_default:

    [snip]

    INFO:snakemake.logging:
    [Fri Sep  6 12:26:23 2024]
    INFO:snakemake.logging:[Fri Sep  6 12:26:23 2024]
    Finished job 0.
    INFO:snakemake.logging:Finished job 0.
    2 of 2 steps (100%) done
    INFO:snakemake.logging:2 of 2 steps (100%) done
    Complete log: .snakemake/log/2024-09-06T122622.173281.snakemake.log
    WARNING:snakemake.logging:Complete log: .snakemake/log/2024-09-06T122622.173281.snakemake.log
    Module run has finished successfully.
    ```

## Run a module specifying the inputs

The benchmark [`tests/data/Benchmark_001.yaml`](https://github.com/omnibenchmark/omnibenchmark/blob/main/tests/data/Benchmark_001.yaml) contains some data processing steps (e.g. `P1`) which take data inputs and produce outputs. To run only module `P1` only on data inputs already available locally at `out/data/D1/default/`, so results will be generated at `out/data/D1/default/process/P1/params/`, first double check the inputs are already where expected:

```
$ ls out/data/D1/default/
D1.meta.json  D1_params.txt  D1.txt.gz
```

If not, run the whole benchmark first (with [`ob run`](https://omnibenchmark.org/tutorial/#run-a-benchmark)). Once the input files are at `out/data/D1/default/`,
run `ob run [benchmark.yaml] --module [MODULE_ID]` with:

=== "Shell"

    ```shell
    ob run tests/data/Benchmark_001.yaml --module P1 --input-dir out/data/D1/default
    ```

=== "Output"

    ```
    Running module on a dataset provided in a custom directory.
    Benchmark YAML file integrity check passed.
    Found 2 workflow nodes for module P1.
    Running module benchmark...
    Assuming unrestricted shared filesystem usage.

    input: out/data/D1/default/D1.txt.gz, out/data/D1/default/D1.meta.json
    [snip]

    localrule all:
    input: out/data/D1/default/D1.txt.gz, out/data/D1/default/D1.meta.json, out/data/D1/default/process/P1/param_0/D1.txt.gz
    jobid: 0
    reason: Input files updated by another job: out/data/D1/default/process/P1/param_0/D1.txt.gz
    resources: tmpdir=/home/imallona/tmp/eb-unlssiuj
    INFO:snakemake.logging:localrule all:
    input: out/data/D1/default/D1.txt.gz, out/data/D1/default/D1.meta.json, out/data/D1/default/process/P1/param_0/D1.txt.gz
    jobid: 0
    reason: Input files updated by another job: out/data/D1/default/process/P1/param_0/D1.txt.gz
    resources: tmpdir=/home/imallona/tmp/eb-unlssiuj

    INFO:snakemake.logging:
    [Fri Sep  6 12:35:15 2024]

    INFO:snakemake.logging:[Fri Sep  6 12:35:15 2024]
    Finished job 0.
    INFO:snakemake.logging:Finished job 0.
    2 of 2 steps (100%) done
    INFO:snakemake.logging:2 of 2 steps (100%) done
    Complete log: .snakemake/log/2024-09-06T123513.568197.snakemake.log
    WARNING:snakemake.logging:Complete log: .snakemake/log/2024-09-06T123513.568197.snakemake.log
    Module run has finished successfully.
    ```


## Remote storage - S3 (AWS or MinIO)

To restrict access to a dedicated bucket an access key with a specific policy have to generated.

### Create policy

Create new policy with
```
ob remote policy create --benchmark tests/data/Benchmark_001.yaml
```
The output of this command needs to be added to either MinIO or AWS as described below.

### Create new access key

#### MinIO

In the MinIO Console navigate to 'Access Keys' and click 'Create access key'. Set 'Restrict beyond user policy' to 'ON'. Replace the displayed policy with the output of the above command.
Optionally enter a name and a description. Click on `Create` and copy the access key and secret key.

#### AWS

Create a new user. Create a new policy with the output of the above command. Attach policy to user. Create access key for user.

### Save access key information locally (Optional)

Save the access key and secret key in a `<CONFIG>.json` file somewhere with the following format:

```
{"access_key": "<ACCESS_KEY>", "secret_key": "<SECRET_KEY>"}
```

### Usage

To use the credentials to write to the remote storage the access key and secret key are passed to omnibenchmark with environment variables. If the credentials have been stored as described above the environment variable `OB_STORAGE_S3_CONFIG` can be set which is the name of the config file. For example:

```
OB_STORAGE_S3_CONFIG=<CONFIG>.json ob run tests/data/Benchmark_003.yaml
```

alternatively `OB_STORAGE_S3_ACCESS_KEY` and `OB_STORAGE_S3_SECRET_KEY` can be set. For example:

```
OB_STORAGE_S3_ACCESS_KEY=<ACCESS_KEY> OB_STORAGE_S3_SECRET_KEY=<SECRET_KEY> ob run tests/data/Benchmark_003.yaml
```

### Versioning

To version currently stored data in the remote (i.e. make it read-only) run the following:
```
ob remote version create -b tests/data/Benchmark_003.yaml
```

A second execution will result in an error since this version now already exists. To create a new version, first update the version in the Benchmark.yaml file and then rerun the above command.
