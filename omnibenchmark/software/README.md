# Design


# Use cases

- Omniblock = omnibenchmark module
- Module = lmod software module, as configured by easybuild
- Environment = conda env
- Singularity = benchmark-specific image(s) reproducibly built with easybuild (apptainer export)
- Apptainer = singularity

Depending on the CLI use:

1. Not contributing, just downloading / checking benchmarks
   - Do nothing (software-wise)
2. Running an omniblock with just its direct inputs
   - (Download the inputs)
```
if using system binaries
    run
else if using conda
    if conda not installed
        error - conda needs to be installed by the user
    else:
        ## pin files are:
        ## https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#freezing-environments-to-exactly-pinned-packages
        download the benchmark's yaml conda env YAMLs and pin files
        run with --use-conda flag (following the pins)
else if using singularity
    if singularity not installed
        error - the user - singularity needs to be installed
    else
        pull relevant image(s)
        run with --use-singularity flag
else
    break we don't offer easybuild for module contribs
```
3. Running an omniblock and some ancestor omniblocks
   - get the module list sorted from initial to terminal
     try `running an omniblock with just its direct inputs` for each omniblock
4. Pinning conda versions (benchmark- and version-specific)
```
if conda not installed
    error - conda needs to be installed by the user
else:
    for each omniblock
        build the conda environment
        conda list --explicit the pins
        store the pins on the remote storage, so people can conda-run stuff with conda-pinned versions
        
```
5. Building and pushing the apptainer images (benchmark- and version-specific)
```
if lmod not installed
    install lmod

for each omniblock:
    easybuild the easyconfig, build the singularity image
    push image to the registry
```

6. Running a benchmark
```
if using system binaries
    disable writing to the remote but allow running (locally)
else if using conda
    if conda not installed
        error - conda needs to be installed by the user
    else
        warn the user the benchmark is using conda!
        run with --use-conda using the "pinned" (explicit) versions (pin file)
        log this has been run with conda
else if using easyconfig
    if lmod not installed
        install lmod
    for each omniblock
        easybuild the easyconfig
else if using apptainer:
    if apptainer not installed:
        error - apptainer needs to be installed by the user
    else:
        pull the matching apptainer image(s)
        run with --use-singularity
```


