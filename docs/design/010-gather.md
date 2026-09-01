# 010: Generic Gather and Joins on Stage Output Contracts

[![Status: Draft](https://img.shields.io/badge/Status-Draft-yellow.svg)](https://github.com/omnibenchmark/docs/design)
[![Version: 4](https://img.shields.io/badge/Version-4-blue.svg)](https://github.com/omnibenchmark/docs/design)

| | |
|---|---|
| **Authors** | btraven00 |
| **Date** | 2026-08-25 |
| **Status** | Draft |
| **Version** | 4 |
| **Supersedes** | N/A |
| **Reviewed-by** | daninci, atchox |
| **Related Issues** | [#289](https://github.com/omnibenchmark/omnibenchmark/issues/289) (multi-stage outputs), [#291](https://github.com/omnibenchmark/omnibenchmark/pull/291) (earlier gather proposal, held back) |
| **Related designs** | [004](004-yaml-specification.md) §3.11–3.12 — the YAML surface; [008](008-filtering.md) — filtering and gating (landed); 009 — named outputs (not landed) |

## Changes

| Version | Date | Description | Author |
|---------|------|-------------|--------|
| 1       | 2026-07-07 | Initial draft | btraven00 |
| 2       | 2026-08-13 | §3.1 chain resolution split into *shadowing* and *alternatives*; define *producer* / *maximal*; alternatives ungated on api_version | btraven00 |
| 3       | 2026-08-25 | Simplify: syntax moves to 004 §3.11–3.12, mechanism moves to §5, alternatives cut to one line each. | btraven00 |
| 4       | 2026-08-28 | Drop `prefix:`; a gather's output tree roots at the stage id (§3.3, §4.5) | btraven00 |

## 1. Problem Statement

At the current state (v0.6.0), benchmark provenance is only captured via the
nested tree model. Each module adds three levels of nesting (stage, module and
parameter set) to the tree. This allows us to express only trivial fan-outs,
with combinatorial pruning in the form of excludes
([004 §3.8](004-yaml-specification.md)) and metric-collectors as a special case
of terminal leaf gather.

Two related cases are needed to express real-world benchmarks, both fan-in
variants:

1. **The join** — a stage cannot consume outputs from stages on *divergent
   branches* (no producer is an ancestor of the other, even when they share the
   same named output). The resolver linearises each stage onto one lineage, so
   one branch silently falls out at run time ("Could not resolve input"). PR
   #367 added a regression fixture and made plan validation reject the diamond
   up front; this design makes it legal at api ≥ 0.7.0.
2. **The gather** — a stage cannot consume outputs from *many* upstream nodes at
   once: N files fanned into one list-valued input, e.g. "all clustering results
   for dataset D". A metric module wanting that shape has no way to express it,
   yet it is the efficient one — a single pass over the N files, read once,
   computed vectorised, written back per-row or aggregated.

#367 also made stage expansion topological, which allows out-of-order
declarations on a single chain. Still, an output id with several producers
resolves to a single one silently: nearest-ancestor for chain inputs,
first-declared-stage for topology edges. The second is why stage-level diagrams
draw fewer edges than the YAML declares. Both come from one single-producer
lookup, so both are fixed in one place — `get_stages_by_output`, which
`stage_adjacency` and input resolution now share. Every topology consumer
inherits the fix: the mermaid and dot exporters, and **obeditor**, which loads
omnibenchmark as a Pyodide wheel and calls `stage_adjacency` directly. No
obeditor-side change is needed, only a wheel bump.

The only fan-in construct at v0.5 is the original `metric_collector`, which is
privileged in every way that matters:

- **Global** — gathers every matching output in the benchmark, with no way to
  group (e.g. "one summary per dataset").
- **Terminal** — its outputs are not registered, so nothing can consume them.
- **Outside the lineage system** — synthetic stage id, no parent, no template
  context, no participation in `requires`/`exclude`.
- **Hacky** — list-valued inputs smuggled through enumerated `input_0..input_N`
  keys.

Users asking for "multi-stage input collection" (#289) want an *ordinary* stage
that collects outputs from several other stages, groups them along a chosen
axis, and feeds further stages. Multi-stage inputs along a single lineage
already work; what is missing is fan-in.

The obstacle is the lineage *chain*, which constrains two things at once.
Internally, a node's identity **is** its chain (`parent_id`, prefix-composed
ids), and a fan-in node has no single chain. On disk the folder tree is the
chain too — every node's directory extends its parent's — so the path doubles
as the provenance record, coupling storage and provenance.

Both fan-in shapes break that, differently:

- **A gather has no parent directory to extend**, because the chain is cut. Its
  fresh tree roots at the stage id, which is unique by construction. A join
  keeps its parents and still extends its deepest input's directory.
- **A join's path stops identifying it.** Two joins sharing a deepest input but
  differing in another parent land on the same path, which Snakemake rejects as
  an ambiguous rule. The directory segment therefore carries a digest of the
  parent set.
- **Neither path encodes its members.** A gather's records only the group key, a
  join's only its deepest input, so the rest of the closure survives only if
  written down — `gathered_from`/`parents` in the plan, `lineage.json` on disk.

This document defines what a gather node is, what its lineage means, and how the
metadata that used to live on the chain is represented once the chain is cut.
The YAML surface itself is specified in
[004 §3.11–3.12](004-yaml-specification.md).

## 2. Design Goals

- **Multi-stage joins behave as expected**: any stage pair connected by a named
  output is connected in the plan and in the rendered topology.
- **One new construct** (`gather:`), reusing existing vocabulary — output ids
  now, lineage labels (008) for later phases.
- **Gather nodes are ordinary nodes**: they run a module, have outputs, and
  downstream stages chain off them (scatter after gather).
- **No silent absence** (008 §2): every dropped member and empty group is
  logged; impossible references are plan-time errors.
- **Provenance survives** the cut, even though gating lineage is reduced.
- **Deprecation path** for `metric_collectors`, the degenerate case: a global
  gather, terminal by choice.

### Non-Goals

- A type system for outputs. Selection is by shared output id; no
  `shape:`/`exports:` field (§4).
- Format compatibility between producers of a shared id. That belongs to the
  planned output validators; until then it is a convention between authors.
- Dynamic fan-out (checkpoints, unknown cardinality). Membership is known at
  plan time.
- Removing `metric_collectors` in this release (deprecation notices only).

## 3. Model

### 3.1 Selection: shared output ids

An output id declared by more than one stage is a *contract*: the files are
interchangeable for consumers referencing that id. A stage declaring an output
id is a **producer** of it.

Output ids and `provides` labels are separate namespaces and never interact: an
id names a file contract, a label names a property of a node's lineage. Binding
one to the other was considered and rejected (§4).

When a consumer references an id with several producers, one question decides
what happens: **does one producer run downstream of another?**

```
chain:       data ──▶ norm ──▶ cluster      norm and cluster both declare `result`
                                            → cluster shadows norm

branches:    data ──▶ stage_a               stage_a and stage_b both declare `result`
                  ╰─▶ stage_b              → alternatives
```

- **Yes — one is downstream of the other.** The downstream producer wins; the
  upstream one is **shadowed**. Its file is still written, but consumers bind to
  the nearer producer.
- **No — they are on parallel branches.** They are **alternatives**, and all of
  them run: each becomes its own expansion base, so the consumer produces one
  node per producer. There is no "nearer" between parallel branches, and picking
  one silently strands the other — yet the files are interchangeable by the
  contract above, so both belong in the plan.

With more than two producers, apply it pairwise: drop any producer that has
another producer downstream of it; every survivor is an alternative.
"Downstream" is read off the resolved nodes — the `parent_id` chain plus fan-in
`parents` edges — not off the declared stage DAG.

> The survivors are the **maximal** producers: those with no other producer
> below them under the ancestry order. That is the term the implementation uses
> and the compact way to state the rule for any number of producers — *every
> maximal producer is an expansion base*. Shadowed producers are exactly the
> non-maximal ones.

**A gather ignores shadowing.** `gather.from` collects every producer of the id,
shadowed ones included. Shadowing decides which producer a chain input binds to;
it is not a rule about membership.

**Alternatives are same-id only.** Two producers of *different* ids that a
consumer needs together are a fan-in, not alternatives: the consumer needs both
at once, so it resolves to one node with several parents (§3.3). Treating them
as alternatives would emit that join once per branch.

**Not gated on api_version**, unlike `gather` and the join, because it is a bug
fix rather than a new construct. Shadowing is a deliberate rule and is
unchanged. What changes is the parallel case, where there is no "nearest" and
the resolver was picking a producer by topological index — an implementation
detail with no meaning in the spec. Nothing worth preserving sits behind a gate.
Not gating is also safe: with one producer per id there is nothing to shadow or
alternate, so existing plans are identical.

### 3.2 The gather stage

Syntax and validation rules: [004 §3.12](004-yaml-specification.md). In short, a
stage declares `gather: [{from, group_by}]` in place of `inputs:`.

Members are every node producing `from`. They are partitioned by the ancestor
module of the `group_by` stage, one gather node per group. This is *structural*
grouping over the existing chain: it needs no `provides` labels.

`group_by` is optional. Omitted, every producer lands in one node and the
output path carries no group segment — that is the global form, and it is what
a `metric_collector` is (§3.6).

Grouping merges parameter expansions: members descending from `data-D1.p1` and
`data-D1.p2` land in the same group, because the key is the ancestor **module**
id (`D1`), not its node id. Per-parameter grouping needs a label-based
`group_by` (§6, Phase 2).

**Group population**, for a stage with several `gather` entries:

- **Every entry populated** — the normal case, one node.
- **No entry populated** — the group does not apply. Warn, emit no node.
- **Some entries populated, others not** — plan-time error. The node's CLI
  would be missing a flag the module declares, so it fails inside argparse at
  run time, far from the cause. Both origins are rejected: no producer of that
  id under the group at all, and every producer excluded for this module.

**Deferred**: `where` (member filtering), label-value `group_by`, and plain
`inputs:` mixed with `gather:`.

### 3.3 The cut: what lineage becomes

A gather **partitions** its members by their `group_by` ancestor and keeps only
the partition key. Everything below that stage is deliberately forgotten.

```
members (producers of `clustering`)      group (ancestor `data` module)
D1-method_a ─┐
D1-method_b ─┴─▶ metrics-summarizer-D1     ({data} = D1)
D2-method_a ─┐
D2-method_b ─┴─▶ metrics-summarizer-D2     ({data} = D2)
```

That splits what the chain used to carry into two records.

**Gating lineage** — what `requires`, `exclude` and path templates see:

- A gather node carries exactly one label: the `group_by` stage id bound to the
  group value, plus the builtin `name`. Downstream `requires:` matches that and
  nothing else. Output templates may reference it and `{name}`/`{params.*}`;
  any other label is a plan-time error, never an empty substitution.
- `exclude` rules pairing the gather's own module with a member's lineage drop
  that member from that module's gather, so one excluded member never poisons
  the whole group. Beyond that the cut holds: an `exclude` pairing a pre-gather
  module with a post-gather one has no chain to be transitive over.

A **join** does not cut. It keeps every parent, and its labels are the union
over all of them — a gate downstream of a join may name a label from any
branch. Labels cannot conflict in that union because a label is owned by
exactly one stage (008 §3.5).

**Provenance lineage** — which concrete nodes fed this one. `gathered_from`
holds the member node ids in the resolved plan; `lineage.json`, written beside
the outputs, holds the same record on disk (§5.2). Full ancestry of any
post-gather node is recovered by walking `parent_id` as usual and, at a gather,
fanning out through the members.

**Identity and paths.** A gather node's id is composed from the group key
(`metrics-summarizer-D1.default`) rather than a parent chain, and `parent_id` is
`None`. Downstream ids prefix-compose off it as usual, so chains resume
immediately: the plan alternates scatter segments and gather cut points, never a
general DAG. Outputs land under `<stage>/<group>/<module>/<param>/…`, rooted at
the stage id — no author-supplied path, so nothing to validate and no way for
two gathers to share a tree — with the existing nested builder resuming below
it.

### 3.4 Membership and pruning order

Membership is computed **after** all pruning:

```
excludes → requires → capability gates → gather membership → where → group_by partition
```

- A pruned lineage never appears in any gather list. An `exclude` naming the
  gather module itself additionally drops matching members from that module's
  gather.
- A group that ends up empty is a warning, and no node is created for that key.
  All groups empty is also a warning; promoting it to an error should land with
  `where`, when silent emptiness becomes likelier.

### 3.5 Module CLI contract

A gather module receives one flag per gather entry, named after the output id,
with every member path:

```
<entrypoint> --output_dir … --name summarizer \
    --clustering /abs/p1.tsv /abs/p2.tsv /abs/p3.tsv \
    [param flags]
```

Output handling is unchanged from any other node and follows whatever the
declared api version specifies — the api ≥ 0.5 contract today, design 009's
named outputs when it lands. There is no gather-specific rule.

Member order is plan order: topological by stage, then declaration order within
a stage, then parameter expansion. Arbitrary but deterministic: identical YAML
must produce byte-identical Snakefiles, or Snakemake reruns unchanged work.
Modules must not read meaning from position; identity comes from the paths.

### 3.6 Metric collectors are a global gather

A `metric_collector` is a gather with no `group_by`, declared terminal by
convention rather than by the engine. The two forms line up field for field:

| `metric_collectors:` | `gather:` |
|---|---|
| `inputs: [methods.result]` | `gather: [{from: methods.result}]` |
| implicit global scope | no `group_by` |
| implicit output location | rooted at the stage id |
| terminal — outputs unregistered | outputs registered; downstream stages consume them |
| — | add `group_by: data` for one report per dataset |

The last row is the reason to migrate: a collector cannot produce one summary
per dataset, and that is the common request.

Deprecation is a **notice, not a rewrite**. Reimplementing the collector
resolver on top of gather nodes rewires a working path for no user-visible
gain; a warning pointing at the equivalent `gather:` stanza does the actual
work of moving people, and the code deletes itself once nobody declares
`metric_collectors:`. Removal waits for a later major api version.

### 3.7 Worked example: a report module

The global form is what a "render a report over everything" module wants:

```yaml
  - id: report
    gather:
      - from: metrics.summary        # every metric node in the benchmark
    modules:
      - id: rmd
        software_environment: r
        repository: {url: …, commit: …}   # the .Rmd lives here, pinned
    outputs:
      - id: report.html
        path: report.html
```

The template belongs in the **module repository**, pinned by commit like any
other module asset. A module is already "code and files at a commit", so
nothing new is needed to make the rendered report reproducible; a
benchmark-level asset mechanism would be a second way to say the same thing.

The part that looks like it needs a mechanism is labelling. A report over N
files must know which run produced each one, or every row is anonymous. That is
already covered: `lineage.json` is written into the node's output directory
*before* the entrypoint runs, so the module can read it and get each member's
stage, module, commit, parameter hash and directory:

```
render.R          # reads $OUTPUT_DIR/lineage.json and the --metrics.summary paths
template.Rmd
```

So the sidecar is not only a provenance artifact for humans — for a gather
module it is part of the input contract (§5.2).

## 4. Alternatives Considered

1. **A shape/type field on outputs (`exports:`/`shape:`)** — outputs declare a
   semantic type and `gather.from` selects by it. Rejected: shared output ids
   cover the same cases with no new vocabulary, and a second "this offers X"
   keyword next to `provides:` is hard to teach. A type field can be added
   back-compatibly later; removing one cannot.
2. **Map-form `Stage.provides: {label: output_id}`** (PR #291) — bind labels to
   output ids and gather by label. Rejected in 008: it overloads the label
   namespace, and a label attaches to a *node*, which underdetermines which file
   to collect.
3. **Topology-directed gather** — "all descendants of a prefix node at stage S".
   Rejected: selection must follow the contract, not the tree. Producers of one
   contract may live on different branches, and once gathers chain, members need
   not share a prefix at all.
4. **Snakemake checkpoints / `directory()` fan-in** — rejected: membership is
   known at plan time, and the explicit-Snakefile backend has no wildcard
   machinery for Snakemake to leverage.
5. **An author-supplied `prefix:` for the cut chain** (in versions 1–3 of this
   document) — rejected: a user-controlled path field owes three answers (is it
   relative to the output directory, is `..` allowed, may two stages share
   one?), and the stage id answers all three by construction — it is relative,
   inert, and unique. An explicit root can be added back if someone wants two
   gathers under one tree; the reverse is not possible.

## 5. Implementation Notes

### 5.1 Representation

| Where | Change |
|---|---|
| `model/benchmark.py` | `Stage.gather: Optional[List[GatherSpec]]`, `GatherSpec {from_ (alias `from`), group_by}`. Gated on api ≥ 0.7.0; `group_by` validated against real stage ids. Duplicate output ids across stages were already legal (`get_outputs()` collapses them) — now the intended contract. |
| `model/resolved.py` | `gathered_from: List[str]` (members, plan order) and `parents: List[str]` (direct parent ids of any fan-in). `is_gather` finally gets set. Gather inputs reuse the collector's enumerated `input_N` keys plus `input_name_mapping`. |
| `core/_expand.py` | `expand_gather_stage()`, a sibling of the scatter path reached by an `if stage.gather:` dispatch. Members come from `output_to_nodes`; partition by the `group_by` ancestor; one node per group × module × parameter; `parent_id=None`; outputs registered so downstream stages consume them normally. |
| `core/_lineage.py` | `select_input_bundles()` pairs producers on divergent branches whose ancestries agree at every shared stage; one node per bundle. `inherited_provides()` unions a bundle's labels. |
| `backend/snakemake.py` | `_write_gather_shell` reused unchanged (it predates gather, reachable only via collectors). Shell contract `--<from_id> path1 path2 …`. Migrating to Snakemake-native list inputs is a possible later cleanup. |

### 5.2 Node identity and multi-parent lineage

A node's `id` was doing three jobs: identity, human-readable label, and lineage
(a prefix-composed chain recovered by string matching). The third is why a node
could have only one parent, which is what forced both fan-in shapes into
workarounds.

- **Lineage is explicit edges.** `parents` holds the direct parent ids of a
  fan-in; root and linear nodes leave it empty, their single parent already
  being in `parent_id`. Ancestry recovery walks `parents` unioned with the
  linear spine, which removes the "downstream of a fan-in sees only the primary
  branch" limitation. This traversal answers *provenance* questions; gating
  deliberately stops at a gather cut (§3.3).
- **Identity is separate from label.** Linear nodes keep their
  lineage-readable id. Fan-in nodes, having no prefix chain, get a
  `stage-module` stem plus a short digest of the parent set — readable and
  unique. The rule name still derives from the id.
- **`lineage.json`** is written beside the outputs of any fan-in node,
  recording each member's node id, stage, module, commit, parameter hash and
  directory. This is what makes the closure recoverable from the filesystem
  alone once a path no longer encodes it.

## 6. Implementation Plan

### Phase 1 — structural gather + fan-in join (api ≥ 0.7.0) — **landed**

Model, expansion, and backend per §5.1, plus `lineage.json`, which was pulled
forward from Phase 4 because a fan-in path is unreadable without it. The scatter
body was extracted alongside so gather is a sibling rather than a deeper nesting
level.

Phase 1 is self-contained: grouping is structural, membership uses the
already-id-keyed `output_to_nodes`, and neither 009 nor a layout change is a
prerequisite.

### Phase 2 — label-based `group_by` — unblocked (008 landed)

Group by a `provides` label instead of a stage: the general form of Phase 1's
structural rule, reading values from the node's labels. Enables tuple keys and
per-parameter grouping. The deprecated `dataset` builtin is not a grouping axis.

### Phase 3 — `where` + cross-boundary `exclude` — unblocked (008 landed)

Member filtering; unsatisfiable cross-gather exclude rules reported;
empty-group warning promoted to an all-groups-empty error.

### Phase 4 — provenance + collector deprecation

Persist the member record into the run metadata (`out/.metadata/`) alongside the
on-disk `lineage.json`; reimplement collectors on top of gather; deprecation
warning on `metric_collectors:`.

### Landing order

010 landed first in code, then 008; both share api `0.7.0`. Only 009 is still
outstanding, and its named-output representation folds in when it lands.

### Deferred

- Plain `inputs:` mixed with `gather:` on one stage — which chain resolves the
  plain input? It needs the group key to pin it.
- Multiple `gather:` entries with different `group_by` axes on one stage
  (cross-product of groups; unclear demand, make it illegal for now).
- **Author-controlled rule labels** — a `label:` template over the same
  `{provides}`/`{module.*}`/`{params.*}` vocabulary as output paths. Snakemake
  does not care what rules are called, so the label is free except for
  uniqueness (resolve, sanitize, hash only the colliders) and validity (a Python
  identifier). Hook: a `label` field defaulting to `node.id`.
- **Content-addressed identity** — a stable `hash(sorted(parents) + module +
  params)` for caching or a results DB. Not needed for execution, since
  Snakemake keys on output paths rather than rule names. The fan-in digest folds
  forward into this if it ever lands.
- **A flat output layout** (one segment per node instead of three) would make a
  gather a non-special case, but nothing here depends on it. It is a future
  evolution belonging to a 007 revision, on its own api bump. The one piece
  gather would want from it is reserving the id separator, so group values stay
  unambiguous inside node ids.
- **Wildcard-based Snakefile emission** — a backend rewrite that loses the
  diffable Snakefile artifact (007). See `scratch/011-wildcard-lowering.md`.

### Testing

- Unit: membership selection, `group_by` partitioning, group-key id composition
  — asserting on resolved node sets, not Snakefile text.
- Integration: `tests/e2e/test_13_gather.py` — two method stages sharing an
  output id, a grouped gather, a downstream consumer, executed.
- Provenance: walk the member record from a terminal node and assert the closure
  matches the plan.

## 7. References

1. [Issue #289 — multi-stage outputs](https://github.com/omnibenchmark/omnibenchmark/issues/289)
2. [004 — YAML specification](004-yaml-specification.md) §3.11–3.12, the gather and fan-in surface
3. [008 — filtering and gating](008-filtering.md) — landed
4. 009 — named outputs — not landed; [PR #329](https://github.com/omnibenchmark/omnibenchmark/pull/329)
5. [Köster et al., Snakemake paper, Fig 8a scatter-gather](https://f1000research.com/articles/10-33/v3#f8)
6. [PR #291 — earlier gather proposal](https://github.com/omnibenchmark/omnibenchmark/pull/291)
