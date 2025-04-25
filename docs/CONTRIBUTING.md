# Contributing to omnibenchmark

The main project documentation lives at [https://omnibenchmark.org](https://omnibenchmark.org).

Here you wil find information for developers that want to contribute to the codebase.

## Setting up a Virtual Environment for development

Currently, the project is using poetry for builds, but the declaration of dependencies also supports other package managers.
This assumes you have
[uv](https://docs.astral.sh/uv/getting-started/installation/) available in your
system (uv is a quite fast manager compatible with pip), but you should be able
to adapt to any other python environment manager of your choice:

```bash
uv venv --python 3.12
uv pip install -r pyproject.toml --all-extras
uv pip install -e .
```


### Using contraints to lock ranges or pin dependencies

The base dependencies are not strictly reproducible. The policy is that they should
be specified as loosely as possible (so that we can keep evolving the codebse,
and above all, adapt to API changes in our dependencies as soon as possible).

If you have specific constrains for whatever reason (temporary pinnings due to
un-addressed issues in third-party libraries, for instance), you can specify
them in the `constraints.txt` file.

We're generating two forms of pinned dependencies, using the `pin-reqs` script:

* A `requirements.txt` file that takes the current `constraints.txt` into account.
* A `poetry.lock` file (to be phased out).

### Activate the Environment

As usual:

```bash
source .venv/bin/activate
```

## Using `pre-commit-config.yaml`

This repository uses [pre-commit](https://pre-commit.com/) to manage Git hooks for enforcing code quality and standards.

### Setup

You should have `pre-commit` installed in your virtualenv, since we installed all extras:

1. **Install Hooks**
   Run this in the root of the repository to set up the hooks:
   ```bash
   pre-commit install
   ```

2. **Run Hooks Manually (Optional)**
   To check all files manually:
   ```bash
   pre-commit run --all-files
   ```

### Updating Hooks

To update to the latest versions of the hooks:
```bash
pre-commit autoupdate
pre-commit install
```

### Example Workflow

1. Make changes and stage them (`git add`).
2. Commit your changes (`git commit`). Hooks will run automatically.
3. Fix any issues reported by the hooks and re-commit.

## Testing

### Overview: The testing pyramid

While contributing to the codebase, you're encouuraged to think about the concept of a [test pyramid](https://martinfowler.com/articles/practical-test-pyramid.html).

In practical terms, this means:

- **Unit Tests**: Fast and isolated tests for individual components (modules). As many as possible, covering all possible scenarios. They should run quickly and in isolation.
- **Integration Tests**: Tests that verify interactions between components. In our case, this can be checking how software provisioning interacts with the workflows.
- **End-to-End Tests**: Tests that simulate real-world usage scenarios. Here we want to check almost from the user perspective: we input a benchmark file, some environment definitions, and we expect the full benchmark execution. Tests that depend on real services, like dummy containers for external storage, belong here. They should be few and strongly bounded in time.

### Running only a category of tests

As an example, we use the `short` tag for tests that we know _should_ run fast. Unit tests, smoke tests etc.

```bash
pytest -m short
```

### Running tests in parallel

Use pytest `-n` flag to increase the level of paralellism. Here we use 6 workers:

```bash
pytest -s -v -n 6 tests/cli/test_run_benchmark.py
```

### Increasing verbosity during tests execution

If you need to debug tests, you might want to increase verbosity while capturing library logs:

```bash
pytest -s -v --log-cli-level=DEBUG tests/cli/test_run_benchmark.py::test_local
```
