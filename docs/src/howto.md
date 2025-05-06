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
          commit: ac5365e
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

## Use a custom singularity container to run methods

We recommend building singularity containers using easybuild and `ob software singularity build --easyconfig [easyconfig]`. Still, it is possible to use any singularity container from an ORAS-compatible registry (could be a GitLab registry), or available locally as a SIF file.

```yaml
---
id: bench1

[snip]

software_environments:                                 
  remote_custom_container:
    description: "A singularity container from a registry"
    ## update the path to an ORAS-compatible registry
    apptainer: oras://registry.renkulab.io/izaskun.mallona/sing
  local_custom_container:
    description: "A singularity container - locally available as a SIF"
    ## local path to a local SIF file
    apptainer: /home/user/singularity_image.sif
```

## Choosing the right software backend (lmod, singularity, conda)

Omnibenchmark is written in python and depends only on pypi packages. Some of its dependencies, namely those related to (reproducible) software stack management, are OS-specific.


| OS                 |  `envmodules`           | `singularity`          |      `conda`          |      `cernvmfs`        |
| -----------        |  -------------------    | ----------------       |  ----------------     | ---------------        |
| `Linux`            |  :material-check-all:   | :material-check-all:   |  :material-check-all: | (:material-check-all:) |
| `MacOS`            |  :material-check:       | :material-close:       |  :material-check-all: | (:material-check-all:) |
| `Windows`          |  :material-close:       | :material-close:       |    :material-help:    |   :material-help:     |



On Linux, software can be managed:

- Using the host's binaries. If relevant interpreters/software are in your $PATH (perhaps using a virtual environment, or directly), they're accessible to omnibenchmark.
- Using conda. For that [mamba](https://github.com/mamba-org/mamba) is required. We provide an `environment.yaml` to help building the environment. 
- Using singularity. For that, apptainer is needed.
- Using environment modules (lmod). For that, lmod is needed.

Similarly, on MacOS:

- Using the host's binaries. If relevant interpreters/software are in your $PATH (perhaps using a virtual environment, or directly), they're accessible to omnibenchmark.
- Using conda. For that [mamba](https://github.com/mamba-org/mamba) is required. We provide an `environment.yaml` to help building the environment. 
- Using singularity. It won't work unless using a virtual machine to provide a Linux-friendly host to singularity.
- Using environment modules (lmod). For that, lmod needs to be installed.

We haven't fully tested omnibenchmark on Windows, but we would recommend using the  Windows Subsystem for Linux (WSL).



## Enabling networking in singularity containers

(Updating)
