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

Omnibenchmark itself is a python package. Some of its dependencies, namely those related to (reproducible) software stack management, are OS-specific.

Software can be managed:

- Using the host's binaries. If relevant interpreters/software are in your $PATH (perhaps using a virtual environment, or directly), they're accessible to omnibenchmark.
- Using conda. For that `conda` is required. We provide a conda environment YAML to help installing all dependencies. 
- Using apptainer. For that, apptainer is needed.
- Using environment modules (lmod). For that, lmod is needed.
