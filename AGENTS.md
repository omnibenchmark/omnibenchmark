# AGENTS.md — OmniBenchmark

Guidance for AI coding agents working on this repository.

## Project Overview

**OmniBenchmark** is an open-source CLI tool (`ob`) for automated scientific benchmarking. It orchestrates complex benchmark studies that evaluate computational methods on datasets with known outcomes, using a declarative YAML-based format and Snakemake as the workflow engine.

- **Docs**: https://docs.omnibenchmark.org
- **Repo**: https://github.com/omnibenchmark/omnibenchmark
- **License**: Apache 2.0
- **Python**: `>=3.12, <3.14`

---

## Repository Layout

```
omnibenchmark/          # Main package
  cli/                  # Click-based CLI subcommands
  benchmark/            # Execution orchestration, DAG management
  model/                # Pydantic data models for benchmark YAML
  workflow/snakemake/   # Snakemake workflow engine
  dag/                  # Custom SimpleDAG (replaces NetworkX)
  remote/               # S3/MinIO remote storage integration
  versioning/           # Benchmark version management
  archive/              # Artifact archiving
  git/                  # Git clone utilities
  templates/            # Copier/Jinja2 templates for scaffolding
  config.py             # ConfigAccessor (~/.config/omnibenchmark)
  constants.py          # Global constants
tests/
  benchmark/            # Benchmark layer unit/integration tests
  cli/                  # CLI tests
  model/                # Model layer tests
  workflow/             # Workflow tests
  remote/               # Remote storage tests
  e2e/                  # End-to-end scenarios
  data/bundles/         # Git bundles for offline test repos
docs/
  src/                  # Markdown sources for MkDocs
  design/               # Architecture Decision Records (ADR-style)
```

---

## Architecture

### Layer Separation (Important)

The codebase enforces a strict separation:

| Layer | Location | Responsibility |
|-------|----------|----------------|
| **Model** | `omnibenchmark/model/` | Pure data models (Pydantic), path-agnostic, YAML ↔ Python |
| **Benchmark** | `omnibenchmark/benchmark/` | Execution context, file system awareness, DAG construction |
| **Workflow** | `omnibenchmark/workflow/` | Engine abstraction; `SnakemakeEngine` is the only implementation |
| **CLI** | `omnibenchmark/cli/` | User interaction layer only; thin wrappers over benchmark/model |

Do not mix concerns between layers. Model classes should not know about file paths; CLI commands should not contain business logic.

### Key Classes

- `BenchmarkExecution` (`benchmark/benchmark.py`): Central orchestrator; loads YAML, builds DAG, drives execution.
- `Benchmark` (`model/benchmark.py`): Top-level Pydantic model for the benchmark YAML spec.
- `SnakemakeEngine` (`workflow/snakemake/snakemake.py`): Translates benchmark definition to Snakemake rules and runs them.
- `SimpleDAG` (`dag/simple_dag.py`): Lightweight DAG implementation (no NetworkX dependency) with topological sort, cycle detection.
- `DAGBuilder` (`benchmark/_dag_builder.py`): Constructs the execution graph from the benchmark model.
- `BenchmarkVersionManager` (`versioning/manager.py`): File-locked version tracking.
- `MinIOStorage` (`remote/MinIOStorage.py`): S3/MinIO backend implementation.

---

## Development Commands

### Running Tests

```bash
pixi run test-short                             # short tests (excluding e2e)
pixi run test-e2e                               # end-to-end tests
pixi run pytest -m "not short and not e2e_s3"   # integration tests
pixi run pytest -k "test_run"                   # by name filter
pixi run pytest tests/cli/ -v                   # specific suite
```

### Linting & Formatting

```bash
pixi run lint                            # lint (ruff)
```

### Type Checking

```bash
pixi run typecheck                       # pyright
```

- **Strict** type checking: `omnibenchmark/dag/`, `omnibenchmark/model/`
- **Standard** type checking: all other modules
- `omnibenchmark/io/` has type checking relaxed (known debt — needs fixes)

### Documentation

```bash
cd docs/ && mkdocs serve                 # live preview at localhost:8000
```

---

## Coding Conventions

### Python Style

- **Formatter**: Ruff (Black-compatible, line length 88)
- **Import sorting**: isort with Black profile
- Pre-commit hooks run Ruff automatically; run `pre-commit install` when setting up.

### Pydantic Models

- All YAML-facing data structures use Pydantic v2.
- Use field validators (`@field_validator`) for semantic validation (file existence, SPDX license IDs, URLs).
- Line number tracking is implemented in `model/_parsing.py` for user-friendly error messages — preserve this when adding new parsed types.

### CLI Commands

- All subcommands use Click and live in `omnibenchmark/cli/{command}.py`.
- Register new commands by importing and adding them to the group in `omnibenchmark/cli/main.py`.
- Keep CLI thin — delegate logic to `benchmark/` or `model/` layers.

### Workflow Engine

- `WorkflowEngine` in `workflow/workflow.py` is the interface; all engines must implement it. **Deprecated — scheduled for removal** ([#201](https://github.com/omnibenchmark/omnibenchmark/issues/201)).
- Snakemake is pinned to `>=9.0,<9.10`. The `<9.10` upper bound is a workaround: Snakemake ≥9.10 detects duplicate rule names, which our current implementation can produce — a module with the same parameters processing different input collections (e.g., `methods-M2-.dc8f7f3b-after_data` and `methods-M2-.dc8f7f3b-after_process`) both generate the rule name `methods_M2_.dc8f7f3b`. The fix is to include the `after` field in rule name generation in `rule_node.smk`. The upper bound will be relaxed once this is resolved ([#269](https://github.com/omnibenchmark/omnibenchmark/issues/269)).

### Error Handling

- Use `omnibenchmark/cli/error_formatting.py` for pretty-printing user-facing errors.
- Custom exceptions live close to their domain (e.g., `remote/exception.py`, `versioning/exceptions.py`).

---

## Adding a New Feature — Checklist

1. **Model changes**: Edit `omnibenchmark/model/benchmark.py` (or `module.py`).
2. **Execution logic**: Edit `omnibenchmark/benchmark/benchmark.py` or add a helper in `benchmark/`.
3. **CLI exposure**: Add a file in `omnibenchmark/cli/` and register in `main.py`.
4. **Tests**: Add tests in `tests/{module}/test_{feature}.py`. Mirror the package structure.
5. **Type check**: `pyright omnibenchmark/` — zero new errors.
6. **Format/lint**: `ruff check --fix && ruff format`.

---

## Testing Notes

- Test repositories are stored as **Git bundles** in `tests/data/bundles/` and extracted by `GitBundleManager` (`tests/git_bundle.py`) to avoid network calls.
- Use the `git_bundle_manager` and `data_repo` pytest fixtures (defined in `tests/conftest.py`) for tests that need real Git repos.
- End-to-end tests in `tests/e2e/` are numbered (`test_00_`, `test_01_`, …) and run in order — be careful about ordering dependencies.
- Remote storage tests (`tests/remote/`) use `testcontainers` to spin up a real MinIO instance; they require Docker.

---

## Key File Paths (Quick Reference)

| Purpose | Path |
|---------|------|
| CLI entry point | `omnibenchmark/cli/main.py` |
| Benchmark execution | `omnibenchmark/benchmark/benchmark.py` |
| Benchmark Pydantic model | `omnibenchmark/model/benchmark.py` |
| Snakemake engine | `omnibenchmark/workflow/snakemake/snakemake.py` |
| DAG implementation | `omnibenchmark/dag/simple_dag.py` |
| DAG builder | `omnibenchmark/benchmark/_dag_builder.py` |
| Remote storage factory | `omnibenchmark/remote/storage.py` |
| Config accessor | `omnibenchmark/config.py` |
| Global constants | `omnibenchmark/constants.py` |
| Test fixtures | `tests/conftest.py`, `tests/fixtures.py` |
| Package metadata | `pyproject.toml` |
| Pixi environments | `pixi.toml` |