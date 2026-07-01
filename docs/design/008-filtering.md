# 008: Filtering and gating mechanisms

[![Status: Draft](https://img.shields.io/badge/Status-Draft-yellow.svg)](https://github.com/omnibenchmark/docs/design)
[![Version: 4](https://img.shields.io/badge/Version-4-blue.svg)](https://github.com/omnibenchmark/docs/design)

| Field           | Value                                                                                                               |
|-----------------|---------------------------------------------------------------------------------------------------------------------|
| Authors         | btraven00, atchox                                                                                                   |
| Date            | 2026-05-07                                                                                                          |
| Status          | Draft                                                                                                               |
| Version         | 4                                                                                                                   |
| Supersedes      | N/A                                                                                                                 |
| Reviewed-by     | TBD                                                                                                                 |
| Related Issues  | [#360](https://github.com/omnibenchmark/omnibenchmark/issues/360) (umbrella: filtering & slicing of the execution graph), [#331](https://github.com/omnibenchmark/omnibenchmark/issues/331) (conditional execution), [#330](https://github.com/omnibenchmark/omnibenchmark/pull/330) (provenance metadata) |

## Changes

| Version | Date       | Description   | Author    |
|---------|------------|---------------|-----------|
| 1       | 2026-05-07 | Initial draft | btraven00 |
| 2       | 2026-06-12 | Two-level label resolution (param-name fallback rejected); reserved builtin labels; diagnostics promoted to phase scope; phases reordered — lineage gating lands before capability gating | btraven00 |
| 3       | 2026-06-29 | Add #360 as the umbrella issue; develop reusable `.obfilter` filters and the drifting-parent re-application case (§3.6); rebut "Snakemake already filters" (§4, Alt 4) | btraven00 |
| 4       | 2026-06-30 | §3.6 rephrased for simplicity; changed drift behaviour to match the prototype in #363; §3.7 marked planned-not-built and grounded on the actual `out/.metadata/manifest.json` | btraven00 |

## 1. Problem Statement

A benchmark plan describes the full universe of stages, modules, and parameter
combinations. In practice users might want to run a *subset* of that universe for
several distinct reasons:

- **Capability gating** — a module needs a GPU, large memory, or a network
  resource that is not present on the current host (#331).
- **Partial pipelines** — during development a user wants to run only up to a
  given stage and inspect intermediates before paying for the rest.
- **Slicing a large benchmark** — a "child" run that is a stable, addressable
  subset of a "parent" plan (e.g. one dataset out of ten, two methods out of
  twenty), produced by tooling such as the `obeditor` web wizard.
- **Ad-hoc skipping** — a known-broken module that should not block the rest
  of a run (just commenting out is a way of getting this, but then we have a different hash
  of the plan and we need to diff what was left out vs. previous versions).

Today, the only mechanism that approaches this is `excludes:` (and `requires:`)
on modules. To exclude a combination, you must enumerate it. This is noisy and brittle
for any non-trivial benchmark and does not help with the four cases above.

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
  over "exclude when Y" (needs changes in less places as the benchmark evolves).
- **No silent absence**: every pruned node must be observable: per-prune log
  lines plus an end-of-expansion summary. In a benchmark, a silently missing
  result cell is worse than a crash: it corrupts the conclusion without
  failing the run.
- **Easy to start**: this design is rolled out feature-by-feature; no
  big migration needed.

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

Stops the resolved DAG at `<stage>` and its lineage. Exactly one stage may be
named.

- Validate that `<stage>` exists; error otherwise.
- Keep `<stage>` plus its **transitive ancestors** — computed from the declared
  cross-stage topology (`stage_adjacency`, the same source the renderers use),
  *not* from YAML declaration order. A benchmark may declare an upstream stage
  after the one that consumes it; declaration-order truncation would silently
  drop a required ancestor (or keep an unrelated downstream branch). Pruning
  happens *before* module resolution, so modules in pruned stages are never
  checked out or environment-resolved.
- Skip metric-collector resolution if any of its declared inputs reference
  pruned outputs.
- Rejected up front when combined with `-m/--module`, which already truncates
  at the target module's stage.

> **TODO (follow-up):** narrow the git prefetch as well.
> `populate_git_cache` still fetches every repository referenced by the
> benchmark; only resolution (checkout + environment setup) is currently
> scoped to the selected stages. Cheap to add once the selected-stages set
> is threaded into the prefetch step.

### 3.3 `--capability <name>` (repeatable)

Declares a capability available on the current host. Modules may declare
required capabilities:

```yaml
modules:
  - id: gpu_pca
    requires_capabilities: [gpu]
```

Resolution rule: a module whose `requires_capabilities` is not a subset of the
provided set is pruned — its outputs are never registered in
`output_to_nodes`, so downstream stages naturally lose those branches without
explicit exclusion logic. One info-level log line per pruned module, and the
pruned count appears in the end-of-expansion summary (§3.5, Diagnostics).

Interaction with `-m/--module` (dev mode): explicit module selection wins.
When `-m` names a module, host gates are bypassed for the whole sub-graph,
with a log line noting the bypass. Silently pruning the very module the user
asked for is never acceptable; if the host truly cannot run it, the module
fails loudly at execution time.

> **TODO (follow-up):** support `requires_capabilities` at the *stage* level
> as well, with **union** semantics — a module's effective set =
> `stage.requires_capabilities ∪ module.requires_capabilities`. Override-down
> (a module silently dropping a stage-level requirement) must not be allowed,
> because it would let a module fake-pass a host gate the stage author
> declared. Defer until a real benchmark asks: the DRY case is narrower than
> it looks (if every module in a stage needs the same capability, that often
> reads cleaner as a separate stage or benchmark). Tracked separately.

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
        provides: {dataset_size: sm}
      - id: huge_set
        provides: {dataset_size: lg}
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

#### Resolution chain (api ≥ 0.6.0)

For each label declared in `Stage.provides`, the per-node value is
resolved in order, most specific first:

1. **`Module.provides[label]`** — explicit benchmark-local binding. This is
   the recommended form: it keeps label routing out of the module's CLI
   parameter contract, so the same module can be reused across benchmarks
   that use different label vocabularies.
2. **`module_id`** — silent default. The zero-config pattern for "label =
   module identity" (e.g. a data stage where each module *is* a dataset
   and you gate downstream by dataset name).

Downstream modules with `requires: {label: value}` run only on upstream
nodes whose resolved label matches by exact string equality.

#### Considered and rejected: parameter-name fallback

An intermediate resolution rule was prototyped and cut before shipping:
when a module had a *parameter literally named* like the label
(`parameters: [{dataset_size: lg}]`), its value would become the label
value. Rejected because it silently couples two contracts with different
owners:

- A module's parameter list is its **CLI surface** — `--dataset_size lg`
  is passed to the entrypoint; it belongs to the module repository.
- A label is **routing metadata** — it belongs to the benchmark.

With the fallback in place, a module author renaming a CLI flag (an
internal refactor, from their perspective) would silently reroute the
benchmark's DAG: label resolution falls through to the module-id default,
downstream `requires:` stops matching, and the benchmark loses cells with
no error at either end. That is action-at-a-distance *across
repositories* — invisible at the point of cause, silent at the point of
effect.

The one case the fallback handled that `Module.provides` cannot: a label
varying *per parameter expansion* of a single module (one module declaring
`parameters: [{dataset_size: sm}, {dataset_size: lg}]`). The answer for
that case is to split the module entry in two — which reads better anyway,
since the label claims to describe the data. If a real benchmark hits this
and splitting is genuinely painful, the fallback (or a templated
`Module.provides: {dataset_size: "{params.dataset_size}"}` form) can be
added back **backward-compatibly**. Removing it after specs depend on it
cannot. That asymmetry decides it.

#### Reserved labels: `name` and `dataset`

The runtime auto-populates two labels on every node: `dataset` (the root
dataset identity, inherited down the lineage) and `name` (the current
module's own id, never inherited — see PR #323). Declaring either in
`Stage.provides` is a **parse-time error**: the runtime would silently
clobber the user's value on every node, and a downstream
`requires: {name: …}` would gate on an accident of propagation rather
than a declared contract.

#### Diagnostics

Pruning must never be invisible (§2, "No silent absence"). Ships with the
lineage phase:

- **Empty-stage warning** — when a stage expands to zero nodes after any
  filter (lineage gate, capability gate, exclusion), emit
  `logger.warning("stage X pruned to empty by …")` at graph-build time.
  Without it, a typo in a label value cascades-empty downstream and looks
  like "nothing happened."
- **Pruned-combinations summary** — one end-of-expansion line:
  `N combinations pruned: X by requires, Y by exclude, Z by capability`.
  This is also the seed of the provenance `filters` block (§3.7): the same
  record, persisted.

#### Working example

```yaml
api_version: "0.6.0"
stages:
  - id: data
    provides: [dataset_size]            # stage advertises this label
    modules:
      - id: small
        provides: {dataset_size: sm}    # explicit binding
        parameters: [{evaluate: "1+1"}]
      - id: huge
        provides: {dataset_size: lg}
        parameters: [{evaluate: "9+9"}]
    outputs:
      - {id: data.raw, path: "{dataset}_data.json"}

  - id: methods
    inputs: [data.raw]
    modules:
      - id: M_lg_only
        requires: {dataset_size: lg}    # lineage gate
        parameters: [{evaluate: "input+10"}]
    outputs:
      - {id: methods.result, path: "{dataset}_method.json"}
```

Result: `M_lg_only` runs only on the `huge` execution path; the `small`
branch is pruned at DAG-construction time. No exclusion list, no
cartesian-product subtraction.

The "label = module identity" case is even shorter — no `Module.provides`
needed:

```yaml
stages:
  - id: data
    provides: [source]
    modules:
      - id: iris
      - id: penguins
      - id: atlas_v2

  - id: methods
    modules:
      - id: cheap
        requires: {source: iris}        # match by module id
```

(The builtin `dataset` label already behaves this way for data stages; a
custom label is needed only when the identity pattern applies to a
non-root stage, or when you want a name that survives the eventual
deprecation of the `dataset` builtin.)

#### Scope and naming

This is the *list-of-labels* form of `Stage.provides:` plus the
*label → value* form of `Module.provides:`. A separate *map* form
(`Stage.provides: {label: output_id}`) was proposed for gather contracts
(see PR #291) but is intentionally out of scope here. If gather lands
later, its binding will use a different keyword (e.g. `gather_outputs:`,
`exports:`) to keep the namespaces clean.

> **TODO (deferred — add only if needed):** *stage-level partition syntax*
> for the case where many modules share a label value. Today, declaring
> `Module.provides: {dataset_size: lg}` on each of N modules is verbose.
> A partition form would let the stage author group modules by label
> value once:
>
> ```yaml
> # Proposed; NOT implemented.
> provides:
>   dataset_size:
>     sm: [iris, penguins, mnist_subset]
>     lg: [atlas_v2, big_brain]
> ```
>
> Resolution would slot in after `Module.provides` and before the
> module-id fallback. Defer until a real benchmark hits this volume — for
> the small-N case, per-module `provides:` is fine and reads more locally.

> **TODO (follow-up):** a `provides` consistency check at parse time. When
> a stage declares `provides: [X]`, warn if a module in the stage has no
> `Module.provides` entry for X — i.e. it would silently fall through to
> module_id. This catches the foot-gun where the author forgot to declare
> the value and a downstream `requires` silently prunes the module. Cannot
> warn unconditionally at the producer side (the "label = module_id"
> pattern is legitimate); needs to look at sibling modules in the stage
> and the consumer side to be precise.

#### Suggested canonical labels

| Label scope | Example labels | Where declared |
|---|---|---|
| Lineage / dataset size | `xs`, `sm`, `md`, `lg`, `xl` (as values of a `dataset_size` provides-label) | `Module.provides` binding + `provides` label on the stage |
| Lineage / domain | `treatment: ctrl|drug`, `species: human|mouse` | stage `provides:` |
| Host / compute tier | `compute_xs`, `compute_sm`, `compute_md`, `compute_lg`, `compute_xl` | module `requires_capabilities:` + run-time `--capability` |
| Host / hardware | `gpu`, `tpu`, `large_mem`, `internet` | module `requires_capabilities:` |

These are conventions, not enforced names. Encouraging consistency lets the
filter wizard offer fixed checkboxes and a benchmarker can pick a tier
without reinventing one per benchmark.

### 3.6 Filter spec

A filter says which parts of a benchmark to run. There are two ways to write
one, and they suit different needs.

**A pick list — name the exact parts to keep.** For each stage, you list which
modules to keep, and for each module which parameter combinations to keep (`all`
of them, just the `first`, or specific ones). Writing `"*"` in place of a module
name means "every module in this stage". This is what the obeditor web wizard
exports today, and what a native `ob run --filter` reads. It is exact and
reproducible. Its trade-off: it is tied to the exact benchmark it was made from
(see "Reusing a saved filter"). This is the form being prototyped first.

**A rule-based filter — describe the parts to keep (sketch only, not built).**
Instead of naming exact parts, you match them by rule:

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

`include` lists what to keep (if omitted, everything is kept); `exclude` then
removes from that result. This style keeps working as the benchmark grows —
"every GPU module" still means the right thing after new modules are added — but
it is less exact than a pick list. It is a possible future addition, not part of
the first version.

#### Reusing a saved filter

A filter is a value you can save to a file, commit to git, and apply again
later — even after the benchmark it came from has changed. This is the third use
case in #360: keeping a stable slice ("just my two methods on the small
dataset") as the parent benchmark keeps evolving, instead of hand-editing a
forked copy every time the parent gains a stage, module, or parameter.

obeditor does this today: the wizard exports a saved filter (a small compressed
text token), and a later run applies it to the current benchmark. See the
[obeditor implementation][obfilter-impl]. When native support lands,
omnibenchmark must read the same format so the two stay interchangeable —
obeditor produces it, omnibenchmark consumes it.

[obfilter-impl]: https://github.com/btraven00/obeditor

**What happens when the benchmark has changed since the filter was made.** A
saved pick list names exact parts, and also records a short fingerprint of the
benchmark it was made from. When you apply it, three things can happen:

- *Benchmark unchanged* (same fingerprint): it runs exactly as picked.
- *Benchmark changed, but everything the filter names still exists*: it runs the
  slice and prints one note saying the benchmark has moved on.
- *Something the filter names is gone* (a stage or module was renamed or
  removed, or a chosen parameter combination no longer exists): it stops with an
  error listing exactly what is missing — never a silently smaller run. You can
  pass an override to run the parts that still match instead.

Two things keep this honest: the fingerprint, which flags any change so nothing
slips by unnoticed, and the `"*"` form, which keeps matching modules that are
added to a stage later.

### 3.7 Provenance hooks (planned — not yet built)

This section describes intended behaviour, not current behaviour. It depends on
the provenance-metadata work in #330.

Today, each run writes a `.metadata/` folder next to its outputs
(`out/.metadata/`): `manifest.json` (a run ID, a timestamp, and host details),
a copy of the benchmark YAML, and the list of resolved modules. None of these
records which filter was applied.

The plan is to add, to that same `manifest.json`, a record of every filter that
was applied to the run (`--until`, `--capability`, or a saved pick list) plus a
fingerprint of the benchmark it ran against. Together these would let the same
slice be reproduced later: take the parent benchmark, apply the recorded filter,
get the same run.

The intended shape:

```yaml
provenance:
  parent: https://example.org/benchmarks/foo@v1.2.3
  parent_fingerprint: 3fd5c5ca…          # of the parent benchmark
  filters:
    until: methods
    capabilities: [gpu]
    pick_list: <saved filter token>      # the obeditor-style pick list, if used
```

This is a small, isolated addition to the manifest writer; the fingerprint is
already computed at run time, so it does not block on #330 if we want it sooner.

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

### Alternative 4: rely on Snakemake's native filtering

- **Description**: Snakemake already exposes `--until` / `--omit-from` (by
  rule name), explicit file targets, and wildcard/target-pattern selection.
  Why not let users filter at the Snakemake layer instead of building a
  selector vocabulary in omnibenchmark?
- **Reason for rejection**: it cannot express the cases in §1, for two
  structural reasons rooted in *how omnibenchmark uses Snakemake today*:

  1. **The stage / module / contract hierarchy is gone by the time it is a
     Snakefile.** Omnibenchmark resolves a plan (stages → modules → parameter
     expansions, plus the `provides` / `requires` lineage contracts) into a
     flat list of concrete rules. The semantic layer a filter needs to select
     *on* — "this rule belongs to the `methods` stage", "this node carries
     `dataset_size: lg`", "this is module `M1`" — is not encoded in the
     generated Snakefile; only rule names and output paths survive.
     Snakemake's `--until RULE` operates on that residue, so the natural
     selectors (`stage_id`, `capability`, lineage label) have nothing to bind
     to at the Snakemake layer. Filtering must happen *above* the lowering
     step, where the hierarchy still exists — which is precisely the layer
     this design targets.

  2. **The generated Snakefiles are explicit, not wildcard-driven.** Snakemake's
     filtering leverage (DAG pruning by wildcard, target patterns, inference
     back from a requested output) assumes rules parameterized by wildcards
     that Snakemake expands. Omnibenchmark currently emits one fully-expanded
     rule per resolved node with materialized paths; there are no wildcards
     for Snakemake's pattern machinery to act on. So even the filtering
     Snakemake *could* do has no surface to grip. Until/unless the generated
     Snakefiles move to wildcard expansion, native filtering is not an
     available lever — and even then, reason (1) keeps the semantic selectors
     above the Snakemake boundary.

  Net: Snakemake filtering and omnibenchmark filtering are at different
  altitudes. We resolve and prune the benchmark DAG *before* lowering to
  Snakemake; Snakemake then schedules whatever survived. The two compose; one
  does not replace the other.

## 5. Implementation Plan

Phases land as separate PRs, in this order. Lineage gating moved ahead of
capability gating (v2 of this doc): it is the pressing replacement for
enumerate-everything `excludes:`, while capabilities serve a narrower
host-heterogeneity case.

### Phase 1 — `--until <stage>` (this PR)
- CLI flag, single stage, existence validation.
- Keep `<stage>` + transitive ancestors via declared topology, not YAML order.
- Stage selection happens *before* module resolution (no checkout/env cost
  for pruned stages); git prefetch narrowing is a follow-up (§3.2 TODO).
- Skip metric collectors that reference pruned outputs.
- Reject combination with `-m/--module`.
- Tests: until=initial / middle / terminal / unknown / ancestor declared
  after target / unrelated branch pruned / diamond lineage.

### Phase 2 — lineage gates (`Stage.provides` / `Module.provides` / `requires`)
- Model: `Stage.provides: list[str]`, `Module.provides: dict[str, str]`;
  both gated on api ≥ 0.6.0. **Mints `APIVersion.V0_6_0`.**
- Two-level resolution chain (§3.5); no parameter-name fallback.
- Reserved-label parse error for `name` and `dataset`.
- Diagnostics: empty-stage warning + pruned-combinations summary (§3.5).
- Tests: unit tests assert on resolved node sets (not generated-Snakefile
  text); e2e fixture for the lineage gate across a multi-stage DAG.

### Phase 3 — `--capability` + `requires_capabilities`
- Model: `requires_capabilities: list[str]` on `Module` (extends the 0.6
  api gate).
- CLI: repeatable `--capability` flag.
- Resolution: prune modules whose required set ⊄ available set *before*
  module resolution; do not register their outputs.
- `-m` bypasses host gates with a log line (§3.3).
- Tests: missing capability prunes module; pruning cascades; provided
  capabilities log line; interaction with `--until` and `-m`.

### Phase 4 — Provenance plumbing
- After #330 lands, emit applied filters into the provenance block,
  reusing the pruned-picks record from Phase 2.
- Round-trip test: serialize → load → re-run produces the same node set.

### Phase 5 — Filter spec (deferred)
- Implement §3.6 once the earlier phases have shipped and we understand the
  ergonomics.

### Testing Strategy
- Unit: selector predicates, `--until` boundary cases, label resolution
  precedence, capability pruning cascades — all asserting on resolved node
  sets.
- Integration: e2e fixtures exercising lineage-gated and capability-gated
  modules across a multi-stage DAG.
- Provenance round-trip in a separate test once #330 lands.

## 6. References

1. [Issue #360 — filtering of the execution graph (umbrella ticket)](https://github.com/omnibenchmark/omnibenchmark/issues/360)
2. [Issue #331 — conditional execution of modules](https://github.com/omnibenchmark/omnibenchmark/issues/331)
3. [PR #330 — provenance metadata tracking](https://github.com/omnibenchmark/omnibenchmark/pull/330)
4. [PR #323 — `{name}` resolves to the current module id](https://github.com/omnibenchmark/omnibenchmark/pull/323)
5. [Snakemake `--until` semantics](https://snakemake.readthedocs.io/en/stable/executing/cli.html)
6. [obeditor — `.obfilter` reference implementation](https://github.com/btraven00/obeditor)
7. [obflow — Go chat-ops automation over benchmark plans](https://github.com/btraven00/obflow)
