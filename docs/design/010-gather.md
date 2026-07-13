# 010: Generic Gather

[![Status: Draft](https://img.shields.io/badge/Status-Draft-yellow.svg)](https://github.com/omnibenchmark/docs/design)
[![Version: 1](https://img.shields.io/badge/Version-1-blue.svg)](https://github.com/omnibenchmark/docs/design)

TODO make metadata a table, this renders badly

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
| 1       | 2026-07-07 | Initial draft | btraven00 |

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
    prefix: metrics                 # required: filesystem prefix for the cut chain (§3.3)
    gather:
      - from: clustering            # required: output id to collect
        group_by: data              # required: STAGE id to group members by
    modules:
      - id: summarizer
        software_environment: env
        repository: {url: ..., commit: ...}
    outputs:
      - id: metrics.summary
        path: "{data}_summary.tsv"      # {data} = the group value (which `data` node)
```

- `from` — an output id. Every node producing it (any stage) is a candidate
  member. Referencing an id with zero producers is a plan-time error.
- `group_by` — **required**, a single **stage id** (not a tuple, not a label).
  Members are partitioned by the ancestor module of that stage they descend
  from: one gather node per distinct ancestor module. This is *structural*
  grouping over the existing lineage chain — it needs no `provides` labels. The
  stage id doubles as the template variable (`{data}` above) bound to the group
  value.
- `prefix` (stage-level, **required**) — the filesystem prefix under which the
  gathered chain is written. A gather cuts the parent chain, so there is no
  parent directory to extend; the stage must name where the fresh chain
  starts. See §3.3.

**Deferred** (later phases, see §5): `where` (member filtering), label-value
`group_by` (grouping by a `provides` label rather than a stage — the general
form of the structural rule here), and plain `inputs:` mixed with `gather:`.

Grouping merges parameter expansions: with `group_by: data`, members descending
from `data-D1.p1` and `data-D1.p2` (two parameter expansions of the same `data`
module `D1`) land in the *same* group, because the group value is the ancestor
**module** id (`D1`), not its node id. If per-parameter grouping is needed,
that is a future label-based `group_by`.

### 3.3 Lineage: what the chain becomes

A gather is a **quotient** over the lineage tree: everything below the
`group_by` stage is deliberately forgotten; what remains is the single group
value (the ancestor module of that stage). This splits what the chain used to
carry into two distinct records:

**Gating lineage** — what `requires`, `exclude`, and path templates reason
over. For a gather node this is exactly the group key:

- `TemplateContext.provides` binds one label: the `group_by` stage id → the
  group value (plus the builtin `name` = the gather module's own id).
  Downstream `requires:` matches against this and nothing else.
- Output path templates may reference `{<group_by stage>}` (the group value)
  and `{name}`/`{params.*}`. Referencing any other label is a **plan-time
  error** (raised in `substitute()`), never an empty substitution.
- `exclude` rules pairing a pre-gather module with a post-gather module have
  no chain to be transitive over. (Cross-boundary exclude degradation is a
  later phase, with `where`; the MVP does not filter members.)

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

**New filesystem prefix (required).** Because the parent chain is cut, a gather
node has no parent directory to extend — the current layout builds every
node's output dir by appending `<stage>/<module>/.<hash>/` to its parent's
path, and a gather node has no parent. The gather stage must therefore declare
an explicit `prefix:` (§3.2) that names where the fresh chain begins on disk.
The group-key segment is written *under* that prefix:

```
<root>/metrics.summarizer.D1.default/…            # root=metrics, group_by=[dataset], value D1
<root>/metrics.summarizer.D1.default/report.R1.default/…
```

This is the minimal, gather-local version of the path decision; it does not
require migrating scatter nodes to the layout-v2 scheme (§3.8). Reserving `.`
as the segment separator (reject dots in stage/module ids) is still needed so
the group value is unambiguous in the segment.

**Provenance now that nesting is broken (open).** With the chain cut, the
directory tree no longer encodes the full upstream closure — the on-disk path
of a post-gather result no longer *is* its lineage. `gathered_from` records the
member node ids in the run metadata, but the layout itself is lossy across the
cut. Candidate mechanism (TBD): a **sidecar file at the gather prefix** — in the
spirit of the existing `parameters.json` — recording, per gather node, the
resolved member set (module id, commit, parameter hash) so the closure is
reconstructable from the filesystem alone, not only from the run's metadata db.
Deciding the sidecar's name, location (per-root vs per-node), and relationship
to `out/.metadata/` is deferred to Phase 4 (provenance); the requirement here
is only that the root exists and is declared.

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
| `model/benchmark.py` | `Stage.gather: Optional[List[GatherSpec]]`; `Stage.prefix: Optional[str]` (**required when `gather` is set** — the filesystem prefix for the cut chain, §3.3); `GatherSpec {from_: str (alias `from`), group_by: str (a stage id)}`. Both fields required; `where` deferred. Gated on api ≥ 0.6.0 (rides with `provides`/named outputs; no new api bump). `group_by` validated to name a real stage. Duplicate output ids across stages are already legal (the `get_outputs()` dict collapses them, so the uniqueness check never fires) — now documented as the intended contract. |
| `model/resolved.py` | `ResolvedNode.gathered_from: List[str]` (member node ids, plan order); `is_gather` already exists — it finally gets set. Gather inputs reuse the collector's enumerated `input_N` keys + `input_name_mapping`, which `_write_gather_shell` already groups back by flag name. |
| `cli/run.py` | Sibling `_expand_gather_stage()` reached by an `if stage.gather:` dispatch (not threaded through the scatter loop). Collect producers of `from` from `output_to_nodes`; partition by the ancestor module of the `group_by` stage (walk `parent_id` via `nodes_by_id`); emit one node per group × module × parameter. `parent_id=None`; outputs written under `prefix/<group>/…` and registered in `output_to_nodes` so downstream stages consume them. |
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
gather, this design adopts the **layout-v2** cleanup for api > 0.6
benchmarks and lets gather fall out of it:

- **One directory per node segment**, collapsing today's three components:
  `stage.module.hash/` instead of `stage/module/.hash/`. A gather node's
  segment is just `stage.module.<groupvalue>.hash/` — no parent, no special
  case. Lineage depth becomes N+1 for N stages instead of 3N.
- **Reserve `.` as the segment separator** by validating that stage and
  module ids contain no dots (parse-time error, same api > 0.6 gate). 004 §3.1
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

### 3.9 Node identity, rule labels, and multi-parent lineage

Today a node's `id` is overloaded: it is simultaneously (a) the **identity**,
(b) the **human-readable label** (the Snakemake rule name is
`_sanitize_rule_name(node.id)`, and the log file is `<rule>.log`), and (c) the
**lineage**, encoded as a prefix-composed chain (`A.p-C2.p-E.p`) that
`get_ancestor_nodes` recovers by id-prefix matching. Overloading (c) into the id
string is *why* a node can have only one parent — a single prefix path — which
forces both gather and the diamond join (#289) into `parent_id`/`gathered_from`
workarounds.

**Decided (used by the join + gather work):**

- **Lineage becomes explicit edges.** `ResolvedNode.parents: List[str]` holds the
  direct parent ids — `[]` root, `[p]` linear, the full producer set for a
  fan-in (join or gather). Ancestry recovery walks `parents` (unioned with the
  id-prefix walk for the linear spine), which removes the "downstream of a
  fan-in sees only the primary branch" limitation. The full record is also
  emitted to a **`lineage.json` sidecar** (a file for now; a DB later).
- **Identity vs label are separated.** The id stays the readable label. Linear
  nodes are unchanged (full lineage-readable id). Fan-in nodes, having no single
  prefix chain, get a `stage-module`-stemmed id plus a short disambiguating hash
  — readable stem, guaranteed unique. The rule name still derives from the id.

**Deferred (capture only — a later concern, not part of the #289/gather fix):**

- **Author-controlled rule labels.** Reuse `TemplateContext` (the same
  `{provides}` / `{module.*}` / `{params.*}` vocabulary as output paths) to let
  a benchmark name its rules, e.g. `label: "cluster_{dataset}_{method}"`.
  Snakemake does not care what rules are called, so the label is free — *except*
  two hard constraints it does enforce: **uniqueness** (resolve → sanitize →
  detect collisions → append a short hash only to the colliders; never trust the
  template to be unique) and **validity** (sanitize to a Python identifier).
  Enabling hook: a `label` field on `ResolvedNode` defaulting to `node.id`, with
  `snakemake.py` using `node.label` for the rule name — then this feature only
  has to populate `label`. Scope (open): per-stage template with an optional
  per-module override (mirrors how `outputs.path` templates are stage-scoped).
- **Content-addressed identity.** A fully stable `key = hash(sorted(parents) +
  module + params)`, order-independent and stable across runs, for caching / a
  results DB / provenance dedup. Separate from the readable label; not needed for
  execution correctness because Snakemake incrementality here keys on **output
  file paths**, not rule names. Fold the short fan-in hash (above) forward into
  this when caching/DB actually lands.

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

The **MVP (Phase 1, global gather)** is deliberately self-contained and does
**not** depend on 008 or full 009:

- Membership uses `output_to_nodes`, already id-keyed in main.
- The only 009-flavoured change it needs is **legalizing duplicate output ids
  across stages** (the shared-id contract, §3.1) — not the full
  `output_name_mapping` rewrite.
- It consumes **no lineage labels**: with `group_by`/`where` deferred, there is
  nothing to group or filter, so `provides` (008) is out of scope. The builtin
  `dataset` label is being deprecated in favour of explicit `provides` naming;
  the MVP relies on neither.

Later phases pick up the heavier pieces when they land:

- **Named outputs (design 009, PR #329)** — the id-keyed `ResolvedNode.outputs`
  dict with `output_name_mapping` replaces the index-matching resolution and
  shares files with gather (`resolved.py`, the expansion-loop output block,
  `snakemake.py`, the collector resolver). Folded in when the representation
  is rewritten, not required to prove fan-in.
- **Lineage labels (design 008, `provides`)** — required before `group_by`
  (Phase 2) and `where` (Phase 3): both read label values from
  `TemplateContext.provides`. Grouping keys are `provides` labels; the
  deprecated `dataset` builtin is **not** a grouping axis.

Landing order: this design's Phase 1 (global gather) can land on main directly;
009 and 008-`provides` sequence in ahead of Phases 2–3. Gather gates on api
`≥ 0.6.0` — it rides with `provides`/named outputs rather than minting a new
api version.

Phased; each lands separately, smallest working piece first.

### Phase 0 — layout v2 + id validation (DEFERRED — not an MVP prerequisite)
Originally a prerequisite so gather nodes could reuse a flat-segment path
scheme. The MVP does **not** need it: a gather node writes under its declared
`prefix:` as `<prefix>/<stage>/<module>/<param>/…`, reusing the existing nested
path builder with the prefix as the new root. The full layout-v2 cleanup
(collapse every node to one segment, reserve `.`, re-audit all path parsers,
§3.8) is a self-contained future change on its own later api bump; gather does
not block on it.

### Phase 1 — structural gather (api ≥ 0.6.0) — **the MVP**
- Lands on current main; needs neither 009 nor 008 (see Dependencies).
- Model: `GatherSpec {from, group_by}` (both required); `Stage.prefix`
  required with `gather`; `group_by` validated to name a real stage; duplicate
  output ids legalized.
- Expansion: sibling `_expand_gather_stage()` reached by an `if stage.gather:`
  dispatch (not inside the 7-deep scatter loop). Partition members by the
  ancestor module of the `group_by` stage (walk `parent_id`); one node per
  group × module × parameter; `parent_id=None`; `gathered_from`; outputs under
  `prefix/<group>/…`, registered downstream.
- Backend: reuses `_write_gather_shell` (enumerated `input_N` + name mapping).
- Tests: two datasets × two methods sharing an output id → one gather node per
  dataset; downstream consumes it; group value in path/template; zero-producer
  `from` and unknown `group_by` stage error.

### Phase 2 — label-based `group_by` (generalization)
- Group by a `provides` label (008) instead of a stage — the general form of
  Phase 1's structural rule. Reads label values from `TemplateContext.provides`.
- Multiple `group_by` labels (tuple); per-parameter grouping.

### Phase 3 — `where` + cross-boundary `exclude` semantics
- Depends on 008 labels. Member filtering; unsatisfiable cross-gather exclude
  rules reported; empty-group warning / all-empty error.

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
