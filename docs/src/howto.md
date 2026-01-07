## Installation

Omnibenchmark is a pip-installable python package ([PyPI](https://pypi.org/project/omnibenchmark/), [source code](https://github.com/omnibenchmark/omnibenchmark)).

### Supported platforms

Even if Omnibenchmark installs, there are limitations of running benchmarks on some operating systems.

| Backend      | Linux | MacOS | Windows |
|--------------|-------|-------|---------|
| Conda        | ✅    | ✅    | ⚠️       |
| Apptainer    | ✅    | ❌    | ❌       |
| Easybuild    | ✅    | ⚠️    | ❌       |
| Lmod         | ✅    | ⚠️    | ❌       |


Similarly, the system architecture matters. Conda packages built for amd64 do not run on arm64 machines.

### Quick start using conda

```shell
# Install Miniforge and git if not already installed
# See: https://conda-forge.org/miniforge/ and https://git-scm.com/

curl -sSL https://raw.githubusercontent.com/omnibenchmark/omnibenchmark/main/omni-environment.yml -o omni-environment.yml

conda create -n omnibenchmark python=3.12 -y
conda activate omnibenchmark
conda env update -f omni-environment.yml

ob --version
```

Start a new benchmark:

```bash
ob create benchmark ~/my_benchmark
```

Creat a new module:

```bash
ob create module ~/my_new_module
```


For detailed instructions, see below.

### Installation via conda

This is the **recommended** way to install Omnibenchmark because it also enables using conda-managed workflows. Similarly, we provide a conda environment YAML to help installing other dependencies, such as `lmod` or `easybuild`.

First, you need to install a Conda-based Python3 distribution. The recommended choice is [Miniforge](https://conda-forge.org/download/).

*Note: Omnibenchmark expects a `conda` command to be available in the PATH, root environment or in the same environment as omnibenchmark itself.*

First, install [miniforge](https://github.com/conda-forge/miniforge).

=== "Shell"

    ```shell
    conda --version
    ```

=== "Output"

    ```
    conda 24.9.2
    ```

Then, download the omnibenchmark environment file and install it in a new conda environment.

=== "Shell"

    ```shell
    curl -sSL https://raw.githubusercontent.com/omnibenchmark/omnibenchmark/main/omni-environment.yml -o omni-environment.yml

    conda init "$(basename "${SHELL}")"
    conda create -n omnibenchmark python=3.12 -y
    conda activate omnibenchmark
    conda env update -f omni-environment.yml
    ```

=== "Output"

    ```
    Empty environment created at prefix: /home/user/miniforge3/envs/omnibenchmark
    ...
    Successfully built omnibenchmark easybuild easybuild-easyblocks easybuild-easyconfigs easybuild-framework snakedeploy
    Installing collected packages: throttler, sortedcontainers, pytz, fastjsonschema, easybuild-framework, easybuild-easyconfigs, easybuild-easyblocks, connection_pool, async, appdirs, wrapt, urllib3, tzdata, typing-extensions, traitlets, tqdm, tabulate, spdx-license-list, smmap, six, rpds-py, reretry, PyYAML, pytrie, pyparsing, pyjwt, pygments, pycparser, pulp, psutil, propcache, pluggy, platformdirs, pillow, packaging, numpy, networkx, multidict, MarkupSafe, kiwisolver, iniconfig, immutables, idna, humanfriendly, hbreader, frozenlist, fonttools, filelock, easybuild, dpath, docutils, cycler, configargparse, click, charset_normalizer, certifi, attrs, argparse-dataclass, annotated-types, aiohappyeyeballs, yte, yarl, typing-inspection, snakemake-interface-common, smart-open, requests, referencing, rdflib, python-dateutil, pytest, pydot, pydantic-core, jupyter-core, jsonasobj2, json-flattener, jinja2, gitdb, deprecated, contourpy, conda-inject, cffi, aiosignal, snakemake-interface-storage-plugins, snakemake-interface-report-plugins, snakemake-interface-logger-plugins, snakemake-interface-executor-plugins, pytest-logging, pynacl, pydantic, pandas, matplotlib, jsonschema-specifications, GitPython, cryptography, aiohttp, prefixcommons, jsonschema, curies, pygithub, prefixmaps, nbformat, snakemake, snakedeploy, linkml-runtime, omni-schema, omnibenchmark
    ...
    Successfully installed GitPython-3.1.44 MarkupSafe-3.0.2 PyYAML-6.0.2 aiohappyeyeballs-2.6.1 aiohttp-3.12.14 aiosignal-1.4.0 annotated-types-0.7.0 appdirs-1.4.4 argparse-dataclass-2.0.0 async-0.6.2 attrs-25.3.0 certifi-2025.7.14 cffi-1.17.1 charset_normalizer-3.4.2 click-8.2.1 conda-inject-1.3.2 configargparse-1.7.1 connection_pool-0.0.3 contourpy-1.3.2 cryptography-45.0.5 curies-0.10.19 cycler-0.12.1 deprecated-1.2.18 docutils-0.21.2 dpath-2.2.0 easybuild-5.1.1 easybuild-easyblocks-5.1.1 easybuild-easyconfigs-5.1.1 easybuild-framework-5.1.1 fastjsonschema-2.21.1 filelock-3.18.0 fonttools-4.59.0 frozenlist-1.7.0 gitdb-4.0.12 hbreader-0.9.1 humanfriendly-10.0 idna-3.10 immutables-0.21 iniconfig-2.1.0 jinja2-3.1.6 json-flattener-0.1.9 jsonasobj2-1.0.4 jsonschema-4.25.0 jsonschema-specifications-2025.4.1 jupyter-core-5.8.1 kiwisolver-1.4.8 linkml-runtime-1.9.4 matplotlib-3.8.0 multidict-6.6.3 nbformat-5.10.4 networkx-3.5 numpy-1.26.4 omni-schema-0.0.6 omnibenchmark-0.4.0 packaging-25.0 pandas-2.3.1 pillow-11.3.0 platformdirs-4.3.8 pluggy-1.6.0 prefixcommons-0.1.12 prefixmaps-0.2.6 propcache-0.3.2 psutil-7.0.0 pulp-3.2.1 pycparser-2.22 pydantic-2.11.7 pydantic-core-2.33.2 pydot-4.0.1 pygithub-2.6.1 pygments-2.19.2 pyjwt-2.10.1 pynacl-1.5.0 pyparsing-3.2.3 pytest-8.4.1 pytest-logging-2015.11.4 python-dateutil-2.9.0.post0 pytrie-0.4.0 pytz-2025.2 rdflib-7.1.4 referencing-0.36.2 requests-2.32.4 reretry-0.11.8 rpds-py-0.26.0 six-1.17.0 smart-open-7.3.0.post1 smmap-5.0.2 snakedeploy-0.11.0 snakemake-9.8.1 snakemake-interface-common-1.20.2 snakemake-interface-executor-plugins-9.3.8 snakemake-interface-logger-plugins-1.2.3 snakemake-interface-report-plugins-1.1.1 snakemake-interface-storage-plugins-4.2.1 sortedcontainers-2.4.0 spdx-license-list-3.27.0 tabulate-0.9.0 throttler-1.2.2 tqdm-4.67.1 traitlets-5.14.3 typing-extensions-4.14.1 typing-inspection-0.4.1 tzdata-2025.2 urllib3-2.5.0 wrapt-1.17.2 yarl-1.20.1 yte-1.9.0
    ```

Check Omnibenchmark has successfully installed.

=== "Shell"

    ```shell
    ob --version
    ```

=== "Output"

    ```
    S3 storage not available. You might want to install the 'minio' and 'boto3' packages.
    OmniBenchmark CLI, version 0.4.0
    ```

### Installation via pip

Omnibenchmark requires `python==3.12`. You might want to configure a [virtual env](https://docs.python.org/3/library/venv.html#module-venv).

You can install omnibenchmark as a python package with `pip`.

=== "Shell"

    ```shell
    pip install omnibenchmark==0.4.0
    ```

### Installation from source

If you want to become a contributor, then you need to install omnibenchmark from source. For more details check out [CONTRIBUTING.md](https://github.com/omnibenchmark/omnibenchmark/blob/main/CONTRIBUTING.md).

## Install Additional Dependencies

Omnibenchmark aims to faciltate benchmarking using different software backends, including `conda`, `apptainer` (formerly `singularity`), or `easybuild`-built envmodules. Hence, extra steps to install some requirements (e.g., `apptainer`, `envmodules`, and so on) are required.

### Installation on Linux

To support custom backends like `conda`, `apptainer` (formerly `singularity`), or `easybuild`-built environment modules, you'll need a few system-wide dependencies. These are all readily available on most Linux distributions.

#### 1. Install Apptainer

Follow the [installation guide for Apptainer](https://apptainer.org/docs/admin/main/installation.html).

#### 2. Install Required System Dependencies

- **debootstrap**
  Required for building Debian-based containers — this is needed even on non-Debian Linux distributions.

- **fakeroot**
  Allows users to simulate root privileges without actual root access. This is especially useful for unprivileged container builds.

After installing `fakeroot`, configure it for apptainer with [`apptainer config fakeroot`](https://docs.sylabs.io/guides/3.5/admin-guide/user_namespace.html#config-fakeroot) to allow non-root users to simulate root privileges while managing containers.

```shell
sudo apt install lua5.2 liblua5.2-dev lua-filesystem lua-posix tcl tcl-dev wget debootstrap software-properties-common
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update
sudo apt install openmpi-bin libopenmpi-dev apptainer
```

Check everything works with:

=== "Shell"

    ```shell
    conda --version
    apptainer --version
    eb --version
    module --version
    ```

=== "Output"


    ```
    conda 24.9.2
    apptainer version 3.11.4
    This is EasyBuild 5.1.1 (framework: 5.1.1, easyblocks: 5.1.1)

    Modules based on Lua: Version 8.7.53 2024-10-12 19:57 -05:00
        by Robert McLay mclay@tacc.utexas.edu

    ```


### Installation on Mac

Mac devices only support `conda` and `easybuild`-built environment modules as custom backends.

#### 1. Install required prerequisites

=== "Shell"

    ```shell
    brew upgrade
    brew install coreutils
    brew install gcc
    brew install lmod
    if [ -f /usr/local/opt/lmod/init/profile ]; then
        source /usr/local/opt/lmod/init/profile
    fi
    if [ -f /opt/homebrew/opt/lmod/init/profile ]; then
        source /opt/homebrew/opt/lmod/init/profile
    fi

    brew install wget
    brew reinstall cmake
    ```

Check everything works with:

=== "Shell"

    ```shell
    conda --version
    eb --version
    module --version
    ```

=== "Output"

    ```
    conda 24.9.2
    This is EasyBuild 5.1.1 (framework: 5.1.1, easyblocks: 5.1.1)

    Modules based on Lua: Version 8.7.53 2024-10-12 19:57 -05:00
        by Robert McLay mclay@tacc.utexas.edu

    ```

#### 2. Persist Lmod setup

To ensure `lmod` is available in every terminal session, add the following to your shell profile:

- For `bash`, edit `~/.bash_profile` or `~/.bashrc`
- For `zsh` (default on modern macOS), edit `~/.zprofile` or `~/.zshrc`

```bash
# Intel-based Macs
if [ -f /usr/local/opt/lmod/init/profile ]; then
    source /usr/local/opt/lmod/init/profile
fi

# Apple Silicon Macs
if [ -f /opt/homebrew/opt/lmod/init/profile ]; then
    source /opt/homebrew/opt/lmod/init/profile
fi
```

## Create a new benchmark

Create a benchmark scaffold with all necessary files:

```bash
ob create benchmark ~/my_benchmark
```

This creates a benchmark directory with:
- `benchmark.yaml` - main configuration file
- `CITATION.cff` - citation metadata
- `envs/` - software environment definitions
- `.git/` - initialized git repository

## Create a new module

Create a standalone module:

```bash
ob create module ~/my_module
```

### Create a module for a specific stage

Generate a module template with pre-configured inputs for a benchmark stage using `--for-stage`:

```bash
ob create module ~/my_method \
  --benchmark ~/my_benchmark/benchmark.yaml \
  --for-stage methods
```

This generates:
- `omnibenchmark.yaml` - module configuration
- Entrypoint script with CLI parsing for stage-specific inputs
- `CITATION.cff` - citation metadata

The `--for-stage` option automatically:
- Reads input/output requirements from the specified stage
- Generates CLI argument parsing code
- Creates appropriate file I/O boilerplate

Example for a methods stage with inputs from data stage:

```bash
# Creates a module that expects data.counts and data.meta as inputs
ob create module ~/clustering_method \
  --benchmark ~/clustering/benchmark.yaml \
  --for-stage methods
```

## Add parameters to a benchmark YAML

Module `P1` is parametrized and will run twice, once with `-a 0 -b 0.1` and second with `-a 1 -b 0.1`.


```yaml
[snip]
stages:
    [snip]
  - id: process
    modules:
      - id: P1
        software_environment: "R"
        parameters:
          - values: ["-a 0", "-b 0.1"]
          - values: ["-a 1", "-b 0.1"]
        repository:
          url: https://github.com/omnibenchmark-example/process.git
          commit: 706edb9
[snip]

```

## Exclude certain module-module chains

Module `M1` won't use inputs from module `D2`.

```yaml
[snip]
stages:
- id: methods
    modules:
      - id: M1
        software_environment: "python"
        exclude:
          - D2
        repository:
          url: https://github.com/omnibenchmark-example/method.git
          commit: 1004cdd
   [snip]
```

## Use a custom apptainer container to run methods

We recommend building apptainer containers using apptainer. Still, it is possible to use any apptainer container from an ORAS-compatible registry (could be a GitLab registry), or available locally as a SIF file.

```yaml
---
id: bench1

[snip]

software_environments:
  remote_custom_container:
    description: "An apptainer container from a registry"
    ## update the path to an ORAS-compatible registry
    apptainer: oras://registry.renkulab.io/izaskun.mallona/sing
  local_custom_container:
    description: "A singularity container - locally available as a SIF"
    ## local path to a local SIF file
    apptainer: /home/user/singularity_image.sif
```

## Simplify benchmark YAMLs using branch names instead of commit names

The `commit` field within a module stanza can accept branch names or git tags. Using commits or tags is recommended.

```yaml
stages:
  - id: data
    modules:
      - id: D1
        name: "Dataset 1"
        software_environment: "python"
        repository:
          url: https://github.com/omnibenchmark-example/data.git
          commit: main # pointing to the latest commit in branch main
    outputs:
      - id: data.image
        path: "{dataset}.png"
```


## Use local module repositories

Local Git repositories can be referenced directly in the `url` field of benchmarking YAML manifestos, without the need to push to GitHub, GitLab, Bitbucket, or any other remote. To ensure changes are tracked, remember to stage them with `git add` and commit them with `git commit` in the local repository. It is recommended to specify full paths (beginning with `/`) to the local Git repository (i.e. `/home/user/repos/data` in the example below).



```yaml
stages:
  - id: data
    modules:
      - id: D1
        name: "Dataset 1"
        software_environment: "python"
        repository:
          url: /home/user/repos/data # full path to the local git repository
          commit: 41aaa0a            # note the commit is still needed
    outputs:
      - id: data.image
        path: "{dataset}.png"
```
