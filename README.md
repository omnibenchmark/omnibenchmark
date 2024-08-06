## omnibenchmark

<p align="center">
<a href="ttps://github.com/omnibenchmark/omni-py"><img alt="Tests" src="./reports/tests.svg"></a>
<a href="https://github.com/omnibenchmark/omni-py/actions"><img alt="Actions Status" src="https://github.com/omnibenchmark/omni-py/workflows/Tests/badge.svg"></a>
<a href="ttps://github.com/omnibenchmark/omni-py"><img alt="Coverage Status" src="./reports/coverage.svg"></a>
<a href="https://github.com/omnibenchmark/omni-py/blob/main/LICENSE"><img alt="License: Apache 2.0" src="https://img.shields.io/badge/License-Apache_2.0-blue.svg"></a>
<a href="https://github.com/psf/black"><img alt="Code style: black" src="https://img.shields.io/badge/code%20style-black-000000.svg"></a>
</p>


The entrypoint to omnibenchmark. It contains a cli and multiple utility functions to support module development via requirements checks, validation and container handling. 

# How to install

0. Install singularity and debootstrap.
1. Install micromamba, mamba or conda, activate with `micromamba activate`.
2. Create an environment. If named omnibenchmark, `micromamba create -n omnibenchmark`.
3. Activate it with `micromamba activate omnibenchmark`.
4. Clone this repository and get to its root.
5. Install with `micromamba install -f test-environment.yml`.
