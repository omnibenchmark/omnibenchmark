# 008: Filtering and gating mechanisms

[![Status: Draft](https://img.shields.io/badge/Status-Draft-yellow.svg)](https://github.com/omnibenchmark/docs/design)
[![Version: 0.1](https://img.shields.io/badge/Version-0.1-blue.svg)](https://github.com/omnibenchmark/docs/design)

**Authors**: btraven00, atchox
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
multiple selector sets (see §3.6).

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

> **TODO (follow-up):** support `requires_capabilities` at the *stage* level
> as well, with **union** semantics — a module's effective set =
> `stage.requires_capabilities ∪ module.requires_capabilities`. Override-down
> (a module silently dropping a stage-level requirement) must not be allowed,
> because it would let a module fake-pass a host gate the stage author
> declared. Defer until a real benchmark asks: the DRY case is narrower than
> it looks (if every module in a stage needs the same capability, that often
> reads cleaner as a separate stage or benchmark). Tracked separately.

> **TODO (follow-up):** emit a diagnostic when a stage ends up with zero
> nodes after capability pruning (or after any other filter). Today the
> resolved DAG silently cascades-empty downstream, which makes a typo in
> a capability name look like "nothing happened." A single
> `logger.warning("stage X pruned to empty by --capability …")` line at
> graph-build time would surface this. Out of scope for the v1 cut.

> **TODO (follow-up):** make capabilities settable once per host so that
> operators don't have to retype `--capability` every run. Two layered
> sources, both strictly additive (host facts only ever accumulate, never
> revoke):
>
> - **`~/.omnibenchmarkrc`** (or equivalent user config) — the right home
>   for *persistent* host facts: capabilities don't change between shell
>   sessions, they're a property of the box. A workstation operator
>   declares `capabilities: [gpu, large_mem]` once.
> - **`OMNIBENCHMARK_CAPABILITIES=gpu,large_mem`** env var — useful for
>   *ephemeral / scoped* contexts: CI jobs that want per-job control
>   without writing files, `direnv`-managed project shells, container
>   entrypoints.
>
> Precedence is union, not override: effective set = rc-file ∪ env ∪ CLI.
> Each layer can only add; nothing can subtract. (If you reach for
> subtraction, see §3.6 / the rejected-`!gpu` note.) Cheap to add when
> someone asks.

#### Considered and rejected: negative capabilities (`--capability !gpu`)

A negation syntax was considered and deliberately not adopted in v1.
`--capability` is scoped to *host facts*: "what does this machine have?"
A `!gpu` form would conflate that with *node selection*: "which modules do
I want to exclude?" Two reasons this is the wrong place:

- **Two semantics on one flag is a UX tax.** Host-fact declaration and
  module-set filtering are orthogonal axes; mixing them invites the user to
  reason about which interpretation applies in each case.
- **Most "I want to skip gpu modules" cases collapse to "don't pass
  `--capability gpu`."** The remaining real case — "skip gpu modules on a
  host that does have gpu" — belongs in the v2 filter spec (§3.6) as
  `exclude: [{capability: gpu}]`, where it sits alongside other exclusions
  in one shared selector grammar.
- **Slippery slope.** Once `!gpu` is accepted, the next requests are
  `gpu|tpu`, `gpu&large_mem`, regex, glob. We do not want to negotiate a
  mini-language on this flag.

If you find yourself reaching for negation, that is the signal to use the
filter spec instead.

#### What capabilities are good for

Capabilities are *host facts*: things that are true (or not) about the
machine running the benchmark. They are static for the duration of a run.
Useful categories:

| Category | Example labels | Meaning |
|---|---|---|
| Hardware accelerator | `gpu`, `tpu`, `cuda_12` | a specific class of accelerator is present |
| Memory tier | `large_mem` | the host has enough RAM for the heavy method |
| Compute tier (laptop → cluster) | `compute_xs`, `compute_sm`, `compute_md`, `compute_lg`, `compute_xl` | coarse profile of available compute (cores × RAM × time budget) |
| Network / data access | `internet`, `s3_egress` | the host can reach remote resources at run time |

The compute-tier vocabulary is a *suggested convention*, not enforced by the
schema. A module declaring `requires_capabilities: [compute_md]` is asking
the operator to assert "this is a medium-tier machine"; it does not measure
RAM. The benchmark author chooses the granularity. Standardizing the names
in the suggested set lets the wizard surface them as fixed checkboxes
instead of free-form strings.

### 3.4 Capabilities are *not* tags, and *not* lineage labels

`requires_capabilities` is deliberately scoped to **host facts**, not to:

- **Free-form tags** (`experimental`, `old`, `phase2`). Any value that begs
  for free-form text is a code smell — it should be modeled as a parameter,
  a stage split, or a separate benchmark.
- **Lineage labels** (`dataset_size: lg`, `treatment: ctrl`). Those describe
  the *data* flowing through an execution path, not the *host*. They are
  expressed via the existing `provides:` (on stages) / `requires:` (on
  modules) mechanism — see §3.5 below.

This distinction matters: a host with `--capability lg` should not be
allowed to fake-run a "large dataset only" method on small data. Host gates
and lineage gates compose independently and a module may need both.

### 3.5 Lineage-derived gates (`provides` / `requires`)

Independent of capabilities, a module can be gated on properties of its
upstream lineage. Stages emit labels via `provides:`; modules opt in via
`requires:`. The values flow with the data, not with the host.

```yaml
stages:
  - id: data
    provides: [dataset_size]
    modules:
      - id: small_set
        parameters: [{dataset_size: sm}]
      - id: huge_set
        parameters: [{dataset_size: lg}]
        requires_capabilities: [large_mem]   # host gate at the source

  - id: methods
    modules:
      - id: cheap_method
        requires: {dataset_size: sm}         # lineage gate
      - id: scalable_method
        requires_capabilities: [gpu]         # host gate
        # no `requires` → runs on every dataset_size
```

Two axes, composed by AND:

- **Host axis** (`requires_capabilities` ↔ `--capability`): is this box
  capable?
- **Lineage axis** (`requires` ↔ `provides`): is this execution path
  carrying the right kind of data?

A module is included only when both gates pass for the candidate node.

#### Suggested canonical labels

| Label scope | Example labels | Where declared |
|---|---|---|
| Lineage / dataset size | `xs`, `sm`, `md`, `lg`, `xl` (as values of a `dataset_size` provides-label) | parameter on the dataset module + `provides` on the stage |
| Lineage / domain | `treatment: ctrl|drug`, `species: human|mouse` | stage `provides:` |
| Host / compute tier | `compute_xs`, `compute_sm`, `compute_md`, `compute_lg`, `compute_xl` | module `requires_capabilities:` + run-time `--capability` |
| Host / hardware | `gpu`, `tpu`, `large_mem`, `internet` | module `requires_capabilities:` |

These are conventions, not enforced names. Encouraging consistency lets the
filter wizard offer fixed checkboxes and a benchmarker can pick a tier
without reinventing one per benchmark.

### 3.6 Filter spec (deferred to v2)

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

### 3.7 Provenance hooks

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
- Implement §3.6 once the v1 features have shipped and we understand the
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
