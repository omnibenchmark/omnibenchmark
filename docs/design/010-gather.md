# 010: Generic Gather

[![Status: Draft](https://img.shields.io/badge/Status-Draft-yellow.svg)](https://github.com/omnibenchmark/docs/design)
[![Version: 0.1](https://img.shields.io/badge/Version-0.1-blue.svg)](https://github.com/omnibenchmark/docs/design)

**Authors**: btraven00
**Date**: 2026-07-07
**Status**: Draft
**Version**: 0.1
**Supersedes**: N/A
**Reviewed-by**: TBD
**Related Issues**: [#289](https://github.com/omnibenchmark/omnibenchmark/issues/289) (multi-stage outputs), [#291](https://github.com/omnibenchmark/omnibenchmark/pull/291) (earlier gather proposal, held back)

## Changes

| Version | Date | Description | Author |
|---------|------|-------------|--------|
| 0.1 | 2026-07-07 | Initial draft | btraven00 |

## 1. Problem Statement

A stage cannot consume outputs from *many* upstream nodes at once. The only
fan-in construct is `metric_collectors`, which is privileged in every way that
matters:

- **Global**: it gathers every matching output in the whole benchmark, with no
  way to group (e.g. "one summary per dataset").
- **Terminal**: its outputs are not registered in the output registry, so
  nothing can consume them downstream.
- **Outside the lineage system**: synthetic stage id, no parent, no template
  context, no participation in `requires`/`exclude`.
- **Hacky**: list-valued inputs are smuggled through enumerated
  `input_0..input_N` keys.

Users asking for "multi-stage input collection" (#289) want an *ordinary*
stage that can collect outputs produced by several other stages, be grouped
along a chosen axis, and feed further stages. Note that multi-stage inputs
along a single lineage already work (ancestor-chain resolution); what is
missing is fan-in.

The structural obstacle is the internal representation: a node's identity is
its lineage *chain* (`parent_id` chain, prefix-composed node ids), and a
fan-in node has no single chain. This document defines what a gather node is,
what its lineage means, and how the metadata that used to live on the chain is
represented once the chain is cut.

## 2. Design Goals

- **One new construct** (`gather:`), reusing vocabularies that already exist:
  output ids for selection, lineage labels (design 008) for filtering and
  grouping.
- **Gather nodes are ordinary nodes**: they run a module, they have outputs,
  downstream stages chain off them normally (scatter after gather).
- **No silent absence** (008 §2): every dropped member and every empty group
  is logged; impossible references are plan-time errors.
- **Provenance survives**: the full set of contributing upstream nodes is
  recorded, even though gating lineage is deliberately reduced.
- **Deprecation path** for `metric_collectors`, which becomes the degenerate
  case (global gather, terminal by choice).

### Non-Goals

- A type system for outputs. Selection is by shared output id; there is no
  `shape:`/`exports:` field (see Alternatives).
- Format compatibility enforcement between producers of a shared output id.
  That is the job of the planned output validators; until then it is a
  convention between benchmark authors.
- Dynamic fan-out (Snakemake checkpoints, unknown cardinality). Membership is
  fully known at plan time.
- Removing `metric_collectors` in this release (deprecation only).

## 3. Proposed Solution

### 3.1 Selection: shared output ids

An output id declared by more than one stage is a *contract*: the files are
interchangeable for consumers that reference that id.

```yaml
stages:
  - id: method_a
    outputs:
      - id: clustering              # shared id — not "method_a.result"
        path: "{name}_assignments.tsv"

  - id: method_b
    outputs:
      - id: clustering
        path: "{name}_labels.tsv"
```

Resolution rules for an output id with multiple producers:

1. **Plain (chain) input** — nearest ancestor on the lineage wins, and a
   warning is emitted when shadowing occurs. (The current implementation
   already resolves duplicates nearest-ancestor-wins, silently; the warning is
   the new part.)
2. **Gather input** — *all* producers, from any stage, post-pruning (§3.4).

### 3.2 DSL

```yaml
stages:
  - id: metrics
    gather:
      - from: clustering            # required: output id to collect
        where: {dataset_size: lg}   # optional: member filter, label vocabulary (008)
        group_by: [dataset]         # optional: labels held fixed; omitted → one global gather
    inputs:
      - some.reference              # optional: plain chain inputs may coexist*
    modules:
      - id: summarizer
        software_environment: env
        repository: {url: ..., commit: ...}
    outputs:
      - id: metrics.summary
        path: "{dataset}_summary.tsv"   # {dataset} valid because dataset ∈ group_by
```

- `from` — an output id. Every node producing it (any stage) is a candidate
  member. Referencing an id with zero producers is a plan-time error.
- `where` — label predicate applied per member, identical semantics to
  `Module.requires` (exact string equality against the member's resolved
  labels). Members that fail are dropped from the list, with a log line.
- `group_by` — list of label names. One gather node is created per distinct
  combination of label *values* among surviving members. Omitted or empty →
  a single gather node over all members (collector semantics).

\* Plain `inputs:` on a gather stage are deferred to a later phase (see §5);
v1 rejects mixing them so the initial semantics stay small.

Grouping is by label **value**, not by minting node: `group_by: [dataset]`
puts members from `data-D1.p1` and `data-D1.p2` (two parameter expansions of
the same dataset module) into the *same* group, because both carry
`dataset=D1`. If per-parameter grouping is needed, the label must encode it.

### 3.3 Lineage: what the chain becomes

A gather is a **quotient** over the lineage tree: everything not named in
`group_by` is deliberately forgotten; what remains is a tuple of label
bindings. This splits what the chain used to carry into two distinct records:

**Gating lineage** — what `requires`, `exclude`, and path templates reason
over. For a gather node this is exactly the group key:

- `TemplateContext.provides` = the `group_by` bindings (plus the builtin
  `name` = the gather module's own id). Downstream `requires:` matches against
  these labels and nothing else.
- Output path templates may reference `{label}` only for labels in
  `group_by`. Referencing anything else (e.g. `{dataset}` downstream of a
  global gather) is a **plan-time error**, never an empty substitution.
- `exclude` rules pairing a pre-gather module with a post-gather module have
  no chain to be transitive over. They degrade to **membership filtering**:
  the excluded lineage's file drops out of the gathered list. A rule that can
  never fire across a gather boundary is reported at plan time.

**Provenance lineage** — which concrete nodes fed this node. This is a new
explicit field:

```
ResolvedNode.gathered_from: List[str]     # member node ids, plan order
```

Full ancestry of any post-gather node is recovered by walking `parent_id`
chains as today and, at a gather node, fanning out through `gathered_from`
into each member's chain. The run metadata (`out/.metadata/`) records
`gathered_from` alongside the resolved-modules list, so a result file's
complete upstream closure — every module id, commit, and parameter hash that
contributed — remains reconstructable even though no single chain encodes it.

**Identity.** Gather node ids are composed from the group key instead of a
parent chain, and downstream ids prefix-compose off them as usual. The group
key is encoded as the ordered label *values* (not `key=value` pairs — `=` is
not path-safe and the label order is fixed by `group_by`), joined with the
same separators used for scatter node ids:

```
metrics-summarizer-D1.default                 # group_by: [dataset], value D1
metrics-summarizer-D1.default-report-R1.default
```

For a global gather (no `group_by`) the key segment is empty:
`metrics-summarizer.default`. This reuses the id-separator reservation from
the layout-v2 cleanup (no dots in stage/module ids; see §3.8): because ids
cannot contain the separators, a group value is unambiguous in the id even
though it is a bare value rather than a `key=value` pair.

`parent_id` is `None` for gather nodes: with grouping by label value there is
in general no single ancestor node to point at. Chains resume immediately
downstream — the plan becomes alternating scatter segments and gather cut
points (series-parallel), never a general DAG.

### 3.4 Ordering of pruning and membership

Membership is computed **after** all pruning:

```
excludes → requires/label gates (008) → capability gates (008) → gather membership → where → group_by partition
```

Consequences, all logged per 008's diagnostics rules:

- A pruned lineage simply never appears in any gather list.
- `where:` drops are counted per group.
- A group that ends up **empty** is a warning (and the gather node for that
  key is not created); *all* groups empty is an error.

### 3.5 Internal representation changes

| Where | Change |
|---|---|
| `model/benchmark.py` | `Stage.gather: Optional[List[GatherSpec]]`; `GatherSpec {from_: str, where: Optional[Dict[str,str]], group_by: Optional[List[str]]}`. Gated on api ≥ 0.7.0. Duplicate output ids across stages become legal (currently an implicit-uniqueness assumption, now a documented contract). |
| `model/resolved.py` | `ResolvedNode.gathered_from: List[str]`; `ResolvedNode.group_key: Dict[str, str]`; `inputs` values may be `List[str]` for gather entries (removes the collector's `input_N` enumeration). `is_gather` already exists — it finally gets set. |
| `cli/run.py` expansion loop | Gather stages do not take the cartesian product over `previous_stage_nodes`. Instead: collect producers of `from` from `output_to_nodes`, apply `where`, partition by `group_by` values (read from each member's `template_context.provides`), emit one node per group × parameter combination. Gather node outputs *are* registered in `output_to_nodes`, so downstream stages consume them normally. |
| `backend/snakemake.py` | `_write_gather_shell` already exists (currently only reachable via collectors); it gains named list inputs — Snakemake native — instead of enumerated keys. Shell contract: `--<from_id> path1 path2 …` (space-separated, as collectors do today). |

### 3.6 Module CLI contract for gather modules

A gather module receives, per gather entry, one flag named after the output
id with all member paths:

```
<entrypoint> --output_dir … --name summarizer \
    --clustering /abs/p1.tsv /abs/p2.tsv /abs/p3.tsv \
    [param flags]
```

Under api ≥ 0.6 named-output rules (design 009), `--output` flags are emitted
as for any other node. Member order is plan order (stage declaration order,
then module order, then parameter order) and is deterministic.

### 3.7 Metric collector deprecation

`metric_collectors` is re-expressible as a gather stage with no `group_by`.
Plan:

1. When gather lands, the collector resolver is reimplemented *on top of*
   gather nodes (`is_collector` kept for output-layout compat).
2. A deprecation warning on `metric_collectors:` points at the equivalent
   `gather:` stanza.
3. Removal in a later major api version, not scheduled here.

### 3.8 Output layout and id separators (rides along)

Gather forces a path-layout decision — a gather node has no parent chain to
extend, so its output directory cannot be built the current way (append
`<stage>/<module>/.<hash>/` to the parent's path). Rather than special-case
gather, this design adopts the **layout-v2** cleanup for api ≥ 0.7.0
benchmarks and lets gather fall out of it:

- **One directory per node segment**, collapsing today's three components:
  `stage.module.hash/` instead of `stage/module/.hash/`. A gather node's
  segment is just `stage.module.<groupvalue>.hash/` — no parent, no special
  case. Lineage depth becomes N+1 for N stages instead of 3N.
- **Reserve `.` as the segment separator** by validating that stage and
  module ids contain no dots (parse-time error, same api 0.7 gate). 004 §3.1
  already declares ids alphanumeric+underscore but nothing enforces it; this
  makes the id-in-path and group-value-in-id encodings (§3.3) collision-proof
  by construction. Output ids keep their dots (`data.raw`) — they are never
  path segments. `:` was considered and rejected (illegal on Windows paths,
  needs shell/URL quoting).
- Old-api benchmarks keep the existing nested layout — the change is gated,
  so there is zero back-compat risk to the on-disk trees of existing runs.
- Path *parsers* to audit when this lands: the v0.4 `--name` extraction in
  `snakemake.py _write_shell`, `_node.py get_benchmark_path` (commonpath),
  `ob collect` / dashboard tree walking, and e2e `expected.json` globs.

This is documented in full as a standalone cleanup
(`scratch/cleanups_layout_v2_and_run_refactor.md`); it is folded into 010
because gather is the feature that makes it load-bearing rather than merely
nice. 007 (output layout) gets a v2 section when Phase 1 lands.

## 4. Alternatives Considered

### Alternative 1: a shape/type field on outputs (`exports:` / `shape:`)
- **Description**: outputs declare a semantic type; `gather.from` selects by
  type name; stage-local output ids remain unique.
- **Pros**: a file can have both a local name and a contract name; precise
  chain references when two producers share a lineage.
- **Cons**: introduces a second "this offers X" keyword next to `provides:` —
  a near-synonym pair whose distinction (node-level label vs file-level type)
  is hard to teach; one more namespace to validate and document.
- **Reason for rejection**: the shared-output-id form covers the same cases
  with zero new vocabulary. If a benchmark ever needs both names for one
  file, a type field can be added **back-compatibly** later; removing one
  after specs depend on it cannot. Same asymmetry argument as 008's
  parameter-name-fallback cut.

### Alternative 2: map-form `Stage.provides: {label: output_id}` (PR #291)
- **Description**: bind lineage labels to output ids; gather selects by label.
- **Reason for rejection**: already rejected in 008 ("Scope and naming") — it
  overloads the label namespace, and a label attaches to a *node*, which
  underdetermines *which file* to collect.

### Alternative 3: topology-directed gather (collapse an axis under a common prefix)
- **Description**: gather = "all descendants of a prefix node at stage S";
  grouping key = shared lineage prefix.
- **Reason for rejection**: selection must be type-directed — producers of the
  same contract may live in different stages and branches, and once gathers
  chain, members need not share any meaningful prefix. The common-prefix case
  falls out of `group_by` when the grouping labels happen to be prefix-minted
  (e.g. `dataset`); it is not the model.

### Alternative 4: Snakemake checkpoints / `directory()` fan-in
- **Reason for rejection**: membership is fully known at plan time; the
  explicit-Snakefile backend has no wildcard machinery for Snakemake to
  leverage (same altitude argument as 008 Alt 4 and 009 Alt 1).

## 5. Implementation Plan

### Dependencies

This design builds on two in-flight pieces and assumes they merge first:

- **Named outputs (design 009, PR #329)** — required before Phase 1.
  Not strictly for membership (`output_to_nodes` is already id-keyed in
  main), but for the representation gather stands on: 009 turns
  `ResolvedNode.outputs` into an id-keyed dict with `output_name_mapping`,
  replacing the index-matching resolution gather would otherwise have to
  build against and immediately rewrite. Both changes also touch the same
  files (`resolved.py`, the expansion-loop output block, `snakemake.py`,
  the collector resolver), and 009 mints api `0.6.0`, which this design's
  `0.7.0` sequences after.
- **Lineage labels (design 008 Phase 2, PR #354)** — required before
  Phases 2–3. `group_by` reads label values from
  `TemplateContext.provides` and `where` reuses the `requires` matching;
  until labels are populated, only the builtin `dataset` label can group,
  and `where` has nothing to match.

Net landing order: 008 doc (#353) → labels (#354) → named outputs (#329) →
this design, Phase 1 onward.

Phased; each lands separately, smallest working piece first.

### Phase 0 — layout v2 + id validation (prerequisite, api 0.7.0)
- Collapse the per-node path segment to one directory (§3.8); reserve `.` by
  rejecting dots in stage/module ids at parse time. Gate on api 0.7.0.
- One production site (`run.py` output-path block); audit the path parsers
  listed in §3.8; add the 007 v2 section.
- Lands first because it defines the path/id scheme gather nodes reuse, and
  it is self-contained (no gather concepts yet). ~1 day.
- Tests: v0.7 layout is flat-segment; v0.4/0.5/0.6 fixtures unchanged; dotted
  id is a parse error under 0.7.

### Phase 1 — global gather (collector parity, api 0.7.0)
- Requires 009 (#329) merged and Phase 0; see Dependencies.
- Model: `GatherSpec` with `from` only; duplicate output ids legalized;
  nearest-ancestor shadowing warning for chain inputs.
- Expansion: implement as a **sibling** `_expand_gather_stage()`, not a new
  branch inside the existing 7-deep loop (see refactor note below). Gather
  nodes with `group_by` absent (single group), `gathered_from`, list-valued
  inputs, outputs registered downstream.
- Backend: named list inputs in `_write_gather_shell`.
- Tests: two stages sharing an output id; gather collects both; downstream
  stage consumes the gather output; zero-producer `from` errors.

### Phase 2 — `group_by`
- `group_by` on custom labels requires 008 Phase 2 (#354); with only #329
  merged, grouping is limited to the builtin `dataset` label.
- Partition by label values; group-key node ids; `TemplateContext.provides`
  = group key; plan-time error for out-of-key template variables;
  empty-group warning / all-empty error.
- Tests: per-dataset grouping; template error cases; group counts in the
  expansion summary.

### Phase 3 — `where` + cross-boundary `exclude` semantics
- Depends on 008 Phase 2 (labels) being in main.
- Member filtering; unsatisfiable cross-gather exclude rules reported.

### Phase 4 — provenance + collector deprecation
- `gathered_from` in `out/.metadata/`; collectors reimplemented over gather;
  deprecation warning.

### Expansion-loop refactor (accompanies Phase 1, no behavior change)

`run.py` `_generate_explicit_snakefile` expands stages in a ~7-deep loop with
helper functions defined inside the innermost body. Phase 1 must not deepen
it. Approach:

1. In Phase 1, add gather as a top-level `_expand_gather_stage()` reached by a
   single `if stage.gather:` dispatch — do **not** inline it.
2. Immediately after (mechanical), extract the existing scatter body into
   `_expand_scatter_stage()` and hoist the two loop-local helpers
   (`get_output_ids_for_node`, `get_ancestor_nodes`) to module level, passing
   their closed-over state as arguments. Target main loop:
   `nodes = expand_gather(...) if stage.gather else expand_scatter(...)`.
   Existing e2e suite is the regression net.

Deferred (TODOs, add when a benchmark asks):
- Plain `inputs:` mixed with `gather:` on one stage (chain input resolved
  against… which chain? needs the group key to pin it; defer).
- Multiple `gather:` entries with *different* `group_by` on one stage
  (cross-product of groups; unclear demand).
- `group_by` by minting node rather than label value.
- **Wildcard-based Snakefile emission** — an alternative lowering (one rule
  per stage×module, concrete targets only in `rule all`) that shrinks
  Snakefiles from O(nodes) to O(stages×modules) rules at scale. Feasible
  behind the same `ResolvedNode` list (resolver keeps all pruning), and
  flat-segment layout (§3.8) makes the wildcard patterns clean — but it is a
  backend-emitter rewrite that loses the materialized/diffable Snakefile
  artifact (007), and it removes no resolver complexity, so it buys nothing
  for gather. Explicitly **not** coupled here; separate proposal for when node
  counts make explicit Snakefiles hurt.

### Testing Strategy

- Unit: membership selection, `where` filtering, `group_by` partitioning,
  group-key id composition — asserting on resolved node sets, not Snakefile
  text (008's convention).
- Integration: e2e fixture with two method stages sharing an output id, a
  grouped gather, and a downstream consumer (scatter → gather → scatter).
- Provenance: walk `gathered_from` from a terminal node and assert the full
  upstream closure matches the plan.

## 6. References

1. [Issue #289 — support for multi-stage outputs](https://github.com/omnibenchmark/omnibenchmark/issues/289)
2. [Design 008 — filtering and gating mechanisms](008-filtering.md)
3. [Design 009 — named outputs](008-named-outputs.md) (renumbering to 009 pending)
4. [Köster et al., Snakemake paper, Fig 8a scatter-gather](https://f1000research.com/articles/10-33/v3#f8)
5. [PR #291 — earlier gather proposal](https://github.com/omnibenchmark/omnibenchmark/pull/291)
