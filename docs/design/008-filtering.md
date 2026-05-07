# 008: Filtering and gating mechanisms

[![Status: Draft](https://img.shields.io/badge/Status-Draft-yellow.svg)](https://github.com/omnibenchmark/docs/design)
[![Version: 0.1](https://img.shields.io/badge/Version-0.1-blue.svg)](https://github.com/omnibenchmark/docs/design)

**Authors**: btraven00
**Date**: 2026-05-07
**Status**: Draft
**Version**: 0.1
**Supersedes**: N/A
**Reviewed-by**: TBD
**Related Issues**: [#331](https://github.com/omnibenchmark/omnibenchmark/issues/331), [#330](https://github.com/omnibenchmark/omnibenchmark/pull/330) (provenance metadata)

## Changes

| Version | Date | Description | Author |
|---------|------|-------------|--------|
| 0.1 | 2026-05-07 | Initial draft | btraven00 |

## 1. Problem Statement

A benchmark plan describes the full universe of stages, modules, and parameter
combinations. In practice users want to run a *subset* of that universe for
several distinct reasons:

- **Capability gating** — a module needs a GPU, large memory, or a network
  resource that is not present on the current host (#331).
- **Partial pipelines** — during development a user wants to run only up to a
  given stage and inspect intermediates before paying for the rest.
- **Slicing a large benchmark** — a "child" run that is a stable, addressable
  subset of a "parent" plan (e.g. one dataset out of ten, two methods out of
  twenty), produced by tooling such as the `obeditor` web wizard.
- **Ad-hoc skipping** — a known-broken module that should not block the rest
  of a run.

Today, the only mechanism that approaches this is `excludes:` (and `requires:`)
on modules. `excludes:` is *defined by absence*: to exclude a combination, you
must enumerate it. This is brittle for any non-trivial benchmark and does not
help with the four cases above.

## 2. Design Goals

- **One selector vocabulary** shared by every gating feature, so that users
  learn one mental model.
- **Transparent specs**: filter definitions must be human-readable and
  diffable. Opaque blobs are acceptable as a *transport* for tooling output but
  not as the source of truth.
- **Provenance integration**: every applied filter is recorded in the run's
  provenance metadata (#330) so that a child run can be reproduced from
  *parent + filter*.
- **Additive, not subtractive**: prefer mechanisms that say "include if X"
  over "exclude when Y".
- **Cheap to start**: this design is rolled out feature-by-feature; no
  big-bang migration.

### Non-Goals

- A general filter language. We are not building XPath/JMESPath; the selector
  vocabulary stays small and typed.
- Free-form `tags:` on modules. Tags ossify into a parallel DSL that the model
  layer cannot type-check. We prefer typed fields (`requires_capabilities`,
  parameter predicates) for every use case we can foresee.
- The `.obfilter` opaque blob as a v1 surface. It remains a possible *output*
  of `obeditor`, but the underlying spec must be transparent.
- Replacing the cartesian-expansion engine. `provides`/`requires` already
  exists and is the seed of an additive model; we grow it, we do not rewrite.

## 3. Proposed Solution

### 3.1 Selector vocabulary

A *selector* identifies a set of nodes in the resolved DAG. The vocabulary in
v1 is restricted to four predicates, all typed:

| Selector | Matches |
|---|---|
| `stage_id: <id>` | every node whose `stage_id` equals `<id>` |
| `module_id: <id>` | every node whose `module_id` equals `<id>` |
| `capability: <name>` | every node whose module declares `<name>` in `requires_capabilities` |
| `param: {key: <k>, value: <v>}` | every node whose `parameters[k] == v` |

A *selector set* is a list of selectors interpreted as a union. There is no
boolean algebra in v1; if combinations are needed they are expressed as
multiple selector sets (see §3.5).

### 3.2 `--until <stage>`

Stops the resolved DAG at and including `<stage>`. Exactly one stage may be
named.

- Validate that `<stage>` exists; error otherwise.
- Compute `cutoff = stage_order.index(<stage>)` and prune `resolved_nodes` to
  those whose `stage_id` is in `stage_order[:cutoff+1]`.
- Skip metric-collector resolution if any of its declared inputs reference
  pruned outputs.

### 3.3 `--capability <name>` (repeatable)

Declares a capability available on the current host. Modules may declare
required capabilities:

```yaml
modules:
  - id: gpu_pca
    requires_capabilities: [gpu]
```

Resolution rule: a module whose `requires_capabilities` is not a subset of the
provided set is *silently pruned* — its outputs are never registered in
`output_to_nodes`, so downstream stages naturally lose those branches without
explicit exclusion logic. One info-level log line per pruned module.

Rationale for silent (not strict) pruning: typos in capability names surface
loudly via "no nodes resolved for stage X." A `--strict` mode would be
redundant.

### 3.4 Capabilities are *not* tags

`requires_capabilities` is deliberately scoped to host capabilities. Any field
that begs for free-form values ("experimental", "old", "phase2") is a code
smell — it should be modeled either as a parameter, a stage split, or a
separate benchmark.

### 3.5 Filter spec (deferred to v2)

Sketch only; do not implement in v1.

```yaml
filter:
  include:
    - stage_id: data
      module_id: D1
    - stage_id: methods
      module_id: M1
  exclude:
    - capability: gpu
```

Rules: `include` is a positive selector set; if absent, all nodes are
included. `exclude` is removed after `include`. The wizard may emit this spec
verbatim or transport it as base64-JSON (`.obfilter`); both forms must
round-trip.

### 3.6 Provenance hooks

Every applied filter (including `--until` and `--capability`) is recorded
under a `provenance.filters` block in the run metadata (#330). A child run
declares its parent benchmark by canonical URL plus the applied filter spec.
This is what makes a sliced run reproducible.

Concretely:

```yaml
provenance:
  parent: https://example.org/benchmarks/foo@v1.2.3
  filters:
    until: methods
    capabilities: [gpu]
```

## 4. Alternatives Considered

### Alternative 1: free-form `tags:` on modules
- **Description**: arbitrary string tags + `--skip <tag>`.
- **Pros**: maximally flexible.
- **Cons**: opaque to the model, no type checking, attracts every gating
  question into one over-loaded field.
- **Reason for rejection**: see §2 non-goals. Capabilities cover the concrete
  use case in #331 with stronger semantics.

### Alternative 2: replace cartesian + excludes with provides/wants slots
- **Description**: full slot-and-signal mechanism, declarative wiring.
- **Pros**: cleanest model, additive by construction.
- **Cons**: large engine change; pushback on perceived complexity.
- **Reason for rejection (for now)**: `requires:`/`provides:` already exist
  on stages and modules and resolve through `_satisfies_requires`. We extend
  rather than replace.

### Alternative 3: `.obfilter` opaque blob as primary surface
- **Description**: wizard emits a signed/encoded blob; tooling consumes it.
- **Pros**: tamper-evident, easy to embed in URLs.
- **Cons**: unreviewable in PRs; drifts silently when the parent plan
  evolves; teaches users a thing they can't read.
- **Reason for rejection**: kept as a *transport* for the underlying
  transparent spec, not as the source of truth.

## 5. Implementation Plan

### Phase 1 — `--until <stage>`
- CLI flag, single stage, validation against `stage_order`.
- Prune resolved nodes; skip metric collectors that reference pruned outputs.
- Tests: until=initial / middle / terminal / unknown / stage with multiple
  providers.

### Phase 2 — `--capability` + `requires_capabilities`
- Model: `requires_capabilities: list[str]` on `Module` (and optionally
  `Stage`).
- CLI: repeatable `--capability` flag.
- Resolution: prune modules whose required set ⊄ available set; do not
  register their outputs.
- Tests: missing capability prunes module; pruning cascades; provided
  capabilities log line; interaction with `--until`.

### Phase 3 — Provenance plumbing
- After #330 lands, emit applied filters into the provenance block.
- Round-trip test: serialize → load → re-run produces the same node set.

### Phase 4 — Filter spec (deferred)
- Implement §3.5 once the v1 features have shipped and we understand the
  ergonomics.

### Testing Strategy
- Unit: selector predicates, `--until` boundary cases, capability pruning
  cascades.
- Integration: an e2e fixture that exercises capability-gated modules across
  a multi-stage DAG.
- Provenance round-trip in a separate test once #330 lands.

## 6. References

1. [Issue #331 — conditional execution of modules](https://github.com/omnibenchmark/omnibenchmark/issues/331)
2. [PR #330 — provenance metadata tracking](https://github.com/omnibenchmark/omnibenchmark/pull/330)
3. [Snakemake `--until` semantics](https://snakemake.readthedocs.io/en/stable/executing/cli.html)
