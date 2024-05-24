## omnibenchmark

The entrypoint to omnibenchmark. It contains a cli and multiple utility functions to support module development via requirements checks, validation and container handling. 

# How to run

## Installing poetry too

Create a virtualenv and install poetry

```
mkdir -p ~/virtenvs
python3 -m venv ~/virtenvs/poetry
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

# How tu contribute

Lorem ipsum
