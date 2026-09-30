# 012: Typed Ports and a Node Editor for Plans

[![Status: Draft](https://img.shields.io/badge/Status-Draft-yellow.svg)](https://github.com/omnibenchmark/docs/design)
[![Version: 1](https://img.shields.io/badge/Version-1-blue.svg)](https://github.com/omnibenchmark/docs/design)

| | |
|---|---|
| **Authors** | ben |
| **Date** | 2026-09-25 |
| **Status** | Draft |
| **Version** | 1 |
| **Supersedes** | N/A |
| **Reviewed-by** | TBD |
| **Related Issues** | TBD |
| **Related designs** | [004](004-yaml-specification.md) — plan syntax; 009 — named outputs (`kind`); [010](010-gather.md) — gather |

## Changes

| Version | Date | Description | Author |
|---------|------|-------------|--------|
| 1 | 2026-09-25 | Initial draft | ben |

## 1. Problem Statement

We want to assemble a plan visually, Galaxy-style: pick datasets, methods,
metrics and collectors from a catalog, wire them on a canvas, and hand the
result to an agent that fills in whatever is missing.

The canvas is not the hard part. obeditor already loads the real `Benchmark`
model in the browser, and obcommons is already a module catalog. What blocks
it is that **ports are untyped**:

- `IOFile` is `id` + `path`. An output id means something only inside the
  plan that declares it.
- obcommons records `inputs: [methods.mapping, data.meta]`: stage ids of one
  specific benchmark.

So no tool can answer "which catalog modules can consume this output?" across
benchmarks. Galaxy can, because every port carries a datatype.

## 2. Design Goals

- One optional field that makes an output's content recognisable outside its
  plan.
- A catalog query of the form *role + input formats → modules*.
- A canvas that is a view of the plan YAML, not a second source of truth.
- An agent handoff that needs no new schema: `ob validate plan` going green is
  the completion condition.

### Non-Goals

- Format subtyping or hierarchy (EDAM has one; we compare strings).
- Format converters (Galaxy's implicit conversion).
- Checking file contents against a format. That belongs to output validators
  ([002](002-module-artifact-validation.md)).
- Any change to execution. `format` is metadata; the run loop ignores it.

## 3. Proposed Solution

### 3.1 `format` on outputs

```yaml
outputs:
  - id: labels
    path: "{module.id}.labels.tsv"
    format: ob:cluster-labels-tsv     # or edam:format_3475
```

- Optional. Value is a CURIE: `edam:format_NNNN` or `ob:<kebab-name>`.
- Outputs only. An input references an output id, so it inherits that
  output's format.
- ob checks the CURIE pattern and nothing else. It holds no vocabulary.
- Orthogonal to 009's `kind`: `kind` is the container (file, zip), `format`
  is the content.
- Gated on the api version that ships 009 (0.8.0).

### 3.2 The `ob:` vocabulary

One file, `formats.yaml`, in obcommons. Each entry has `id`, `description`
and an example path. Use an EDAM term when one fits; mint an `ob:` term only
for formats EDAM lacks (cluster labels, score tables). Adding a term is a
pull request, like adding a module.

### 3.3 obcommons entry changes

```yaml
role: metric                 # data | method | metric | collector
inputs:
  - {id: labels, format: ob:cluster-labels-tsv}
  - {id: truth,  format: ob:cluster-labels-tsv}
outputs:
  - {id: scores, format: ob:scores-json}
presets:                     # mostly for data fetchers
  - {name: iris, params: {dataset: iris}}
```

- `role` sits next to the free-string `stage`. `stage` still names the stage
  in a given benchmark; `role` is benchmark-independent.
- The current string form of `inputs`/`outputs` stays valid and is untyped.
- A dataset is a fetcher module plus a preset (e.g. `omni-huggingface` with
  `{repo: …}`), not one module per dataset.
- A **metric catalog** is the query `role: metric` filtered by input formats.
  It needs no separate repository. Seed it for clustering from
  `clustbench_metrics`.

### 3.4 Canvas mapping (obeditor)

| Canvas | Plan |
|---|---|
| Box | stage; its ports are the stage's output contracts |
| Chip inside a box | module, with its parameter grid |
| Edge | an entry in the downstream stage's input collection |
| Box with a many-in port | gather stage (010) |

Boxes are stages, not modules. An ob edge means every module of stage B
consumes every output of stage A. One box per module would make the user draw
n×m edges that say nothing more than one edge does.

An edge is allowed when both ends have the same format, or when either end is
untyped; the untyped case shows a warning. Box positions go in a top-level
`x-layout:` key. The model ignores unknown keys, so the plan still validates
and runs unchanged.

### 3.5 Agent handoff

The editor exports a plain plan YAML. A catalog module the user did not pick
is a **stub**: a module entry with a `description` stating its intent and no
`repository`. The stage it sits in supplies its input and output contract.

The plan does not validate while stubs remain, and that is intended. For each
stub the agent:

1. queries obcommons by role and formats, and uses a match if one exists;
2. otherwise scaffolds with `ob create module` against the stage (the
   `omnibenchmark-module` skill) or wraps a named tool
   (`omnibenchmark-wrapping`);
3. writes the `repository` back.

The agent is done when `ob validate plan` passes.

## 4. Alternatives Considered

- **Boxes are modules.** Rejected: n×m redundant edges (§3.4).
- **EDAM only.** Rejected: it has no terms for most of our intermediate
  formats (cluster labels, score tables). It stays the preferred namespace
  where a term exists.
- **Structured type system (schemas, generics).** Rejected: nobody has asked
  for anything beyond equality checks.
- **Separate metric catalog repo.** Rejected: it would duplicate obcommons.
- **Export to Galaxy or CWL workflows.** Rejected: neither models stage-level
  fan-out or gather.
- **Stubs in a sidecar file.** Rejected: the plan would validate while the
  work is unfinished. A missing `repository` makes validation fail until the
  agent is done.

## 5. Implementation Plan

1. `IOFile.format` (optional, CURIE pattern), 004 section, api 0.8.0 gate.
2. obcommons: `role`, typed `inputs`/`outputs`, `presets`, `formats.yaml`;
   seed clustering metrics.
3. obeditor: render an existing plan as a read-only graph (§3.4).
4. obeditor: editing, catalog palette, stub export (§3.5).

### Testing Strategy

- Model: `format` accepted/rejected by pattern; rejected below 0.8.0.
- obcommons: the Zod schema covers the new fields; `npm test` covers the
  role/format filter.
- obeditor: round-trip a plan through the graph and compare YAML, ignoring
  `x-layout`.

## 6. References

1. [EDAM ontology, format branch](https://edamontology.org/format_1915)
2. [Galaxy datatypes](https://docs.galaxyproject.org/en/latest/dev/data_types.html)
