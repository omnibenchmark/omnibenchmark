# Architecture

> **Generated** by `scripts/gen_architecture.py` from package `__init__` docstrings and actual imports. Do not edit by hand; run `python scripts/gen_architecture.py` to refresh.

Packages are layered top (user-facing) to bottom (foundation). Imports only ever point downward — the graph is acyclic.

```mermaid
flowchart TB
  subgraph cli[cli]
    cli
  end
  subgraph backend[backend]
    remote
  end
  subgraph service[service]
    archive
    backend
    versioning
  end
  subgraph interface[interface]
    storage
  end
  subgraph engine[engine]
    core
  end
  subgraph infrastructure[infrastructure]
    git
  end
  subgraph foundation[foundation]
    dag
    model
  end

  archive --> core
  archive --> git
  archive --> model
  archive --> storage
  backend --> core
  backend --> git
  backend --> model
  cli --> archive
  cli --> backend
  cli --> core
  cli --> dag
  cli --> model
  cli --> remote
  core --> dag
  core --> git
  core --> model
  git --> model
  remote --> archive
  remote --> core
  remote --> storage
  remote --> versioning
  storage --> core
  versioning --> model
```

**Import cycles:** none — every package is its own strongly-connected component (clean DAG).

## Packages

| Package | Layer | Represents | Depends on |
|---|---|---|---|
| `cli` | cli (top — nothing imports it) | the user-facing ``ob`` commands; each module maps to a subcommand and wires user input to the layers below | `archive`, `backend`, `core`, `dag`, `model`, `remote` |
| `remote` | backend (implementation) | actually moving bytes to/from object storage so results can be shared and versioned outside the local working tree | `archive`, `core`, `storage`, `versioning` |
| `archive` | service | assembling the parts of a benchmark (configuration, code, software environments, results) into a single archive, for local-only use or upload | `core`, `git`, `model`, `storage` |
| `backend` | service | the compilation pipeline that turns resolved nodes into an executable Snakemake workflow | `core`, `git`, `model` |
| `versioning` | service | the single source of truth for parsing, validating, comparing, and incrementing benchmark versions, plus version-manager lifecycle and locking | `model` |
| `storage` | interface (port) | the *interface* (port) every storage backend implements, the value objects that describe what gets tracked, and the pure rules for which files a benchmark is expected to produce | `core` |
| `core` | engine (central domain) | the engine that turns a validated benchmark definition into a run — builds the DAG of stages/modules and drives it to completion | `dag`, `git`, `model` |
| `git` | infrastructure | fetching and caching the module git repositories a benchmark references, via a two-tier scheme | `model` |
| `dag` | foundation (generic) | backend-agnostic graph operations (topological sort, simple paths, node attributes) — a small standalone toolkit, not benchmark-specific | — |
| `model` | foundation (data model) | what a benchmark *is*, as passive Pydantic structures (stages, modules, parameters, storage config, ...) plus the rules that validate them | — |

