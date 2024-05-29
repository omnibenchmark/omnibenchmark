## omnibenchmark

<p align="center">
<a href="ttps://github.com/omnibenchmark/omni-py"><img alt="Tests" src="./reports/tests.svg"></a>
<a href="https://github.com/omnibenchmark/omni-py/actions"><img alt="Actions Status" src="https://github.com/omnibenchmark/omni-py/workflows/Tests/badge.svg"></a>
<a href="ttps://github.com/omnibenchmark/omni-py"><img alt="Coverage Status" src="./reports/coverage.svg"></a>
<a href="https://github.com/omnibenchmark/omni-py/blob/main/LICENSE"><img alt="License: Apache 2.0" src="https://img.shields.io/badge/License-Apache_2.0-blue.svg"></a>
<a href="https://github.com/psf/black"><img alt="Code style: black" src="https://img.shields.io/badge/code%20style-black-000000.svg"></a>
</p>


The entrypoint to omnibenchmark. It contains a cli and multiple utility functions to support module development via requirements checks, validation and container handling. 

# How to run

## Installing poetry too

Create a virtualenv and install poetry, with a compilled (altinstall) python3.12

```
mkdir -p ~/soft/python
cd $_
wget https://www.python.org/ftp/python/3.12.3/Python-3.12.3.tgz
tar xzvf Python-3.12.3.tgz
cd Python-3.12.3
./configure --enable-optimizations
make -j 8
sudo make altinstall # to /usr/local/bin/python3.12

```

```
mkdir -p ~/virtenvs
/usr/local/bin/python3.12 -m venv ~/virtenvs/poetry
source ~/virtenvs/poetry/bin/activate
curl -sSL https://install.python-poetry.org | python3 -
poetry --version
```

Clone the repo

```
mkdir -p ~/src
cd ~/src
git clone git@github.com:omnibenchmark/omni-py.git
```

```
cd ~/src/omni-py
poetry shell
poetry install
ob --help
deactivate
```

## Directly

Clone the repo

```
mkdir -p ~/src
cd ~/src
git clone git@github.com:omnibenchmark/omni-py.git
```

Install

```
pip install omni-py
```

To develop

```
pip install poetry
poetry shell
ob --help
ob fetch --help
```

# Once installed

```
source ~/virtenvs/poetry/bin/activate 
poetry shell
```

# How tu contribute

Lorem ipsum
