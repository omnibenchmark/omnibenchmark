# Logics

Omniblock = omnibenchmark module
Module = lmod software module

Depending on the CLI use:

1. Not contributing, just downloading / checking benchmarks
   - Do nothing.
2. Running something
   - Install lmod
3. Running an omniblock with just its direct inputs
   - (Download the inputs)
```
if using system binaries:
    run
else if using singularity and singularity installed
    run with singularity
else if using singularity and singularity not installed
    install singularity
    pull image
    run with singularity
else if using easybuild
parse the benchmark yaml to get the module name
    if module not installed
    easybuild the easyconfig
     ```
4. Running an omniblock and some ancestor modules
   - get the module list from earlier to late
     try `running a module with just its direct inputs` for each omniblock
5. Running a benchmark
   - list all modules and install all easyconfigs
