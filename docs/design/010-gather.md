# 010: Generic Gather

[![Status: Draft](https://img.shields.io/badge/Status-Draft-yellow.svg)](https://github.com/omnibenchmark/docs/design)
[![Version: 1](https://img.shields.io/badge/Version-1-blue.svg)](https://github.com/omnibenchmark/docs/design)

| | |
|---|---|
| **Authors** | btraven00 |
| **Date** | 2026-07-07 |
| **Status** | Draft |
| **Version** | 1 |
| **Supersedes** | N/A |
| **Reviewed-by** | daninci (TBD) |
| **Related Issues** | [#289](https://github.com/omnibenchmark/omnibenchmark/issues/289) (multi-stage outputs), [#291](https://github.com/omnibenchmark/omnibenchmark/pull/291) (earlier gather proposal, held back) |
| **Related designs** | 008 — filtering and gating ([PR #354](https://github.com/omnibenchmark/omnibenchmark/pull/354), in flight), 009 — named outputs ([PR #329](https://github.com/omnibenchmark/omnibenchmark/pull/329), in flight); referenced throughout — they land *after* 010 (see §6) |

## Changes

| Version | Date | Description | Author |
|---------|------|-------------|--------|
| 1       | 2026-07-07 | Initial draft | btraven00 |

## 1. Problem Statement

Two related things are missing, both flavours of fan-in:

1. **The join**: a stage cannot consume outputs from stages on *divergent
   branches* (no producer is an ancestor of the other). The resolver
   linearises each stage onto one lineage, so one branch silently falls out at
   run time ("Could not resolve input"). PR #367 (daninci) added the
   regression fixture and made plan validation reject the diamond up front;
   this design makes it legal at api ≥ 0.7.0 — the same fixture now drives
   both `test_validate_plan_accepts_diamond_input_collection` (0.7, accepted)
   and `test_validate_plan_rejects_diamond_below_api_0_7` (0.5, the old
   failure mode) in `tests/cli/test_validate.py`.
2. **The gather**: a stage cannot consume outputs from *many* upstream nodes
   at once — N files fanned into one list-valued input, e.g. "all clustering
   results for dataset D".

#367 also made stage expansion topological, which *allows* out-of-order
declarations on a single chain — but an output id with several producers is
still resolved to a single one silently (nearest-ancestor for chain inputs;
first-declared-stage for topology edges, fixed here by
`get_stages_by_output`, see `test_get_stages_by_output_returns_all_producers_in_order`).

The only fan-in construct at v0.5 is the original `metric_collector`
concept, which is privileged in every way that matters:

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

The structural obstacle is the lineage *chain*, which is load-bearing in two
places at once. Internally, a node's identity **is** its chain (`parent_id`
chain, prefix-composed node ids), and a fan-in node has no single chain. On
disk, the folder tree **is** the chain too — every node's output directory
extends its parent's directory, so the path doubles as the provenance record,
and a fan-in node has no parent directory to extend (hence `prefix:`, §3.3)
and no path that encodes its members (hence `gathered_from`/the sidecar,
§3.3). This document defines what a gather node is, what its lineage means,
and how the metadata that used to live on the chain is represented once the
chain is cut.

## 2. Design Goals

- **One new construct** (`gather:`), reusing vocabulary defined elsewhere:
  output ids for selection now, lineage labels (design 008, in flight) for
  filtering and grouping in later phases.
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

1. **Plain (chain) input** — nearest ancestor on the lineage wins. (The
   current implementation already resolves duplicates nearest-ancestor-wins,
   silently; emitting a shadowing warning is planned work, not in this PR.)
2. **Gather input** — *all* producers, from any stage, post-pruning (§3.4).

### 3.2 DSL

```yaml
stages:
  - id: metrics
    prefix: aggregated              # required: filesystem prefix for the cut chain (§3.3)
    gather:
      - from: clustering            # required: output id to collect
        group_by: data              # required: STAGE id to group members by
    modules:
      - id: summarizer
        software_environment: env
        repository: {url: ..., commit: ...}
    outputs:
      - id: metrics.summary
        path: "{data}_summary.tsv"      # {data} = the group value (which `data` module)
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

A gather **partitions** the members by their `group_by` ancestor and keeps
only the partition key: everything below the `group_by` stage is deliberately
forgotten; what remains is the single group value (the ancestor module of
that stage). This splits what the chain used to carry into two distinct
records:

With the §3.2 example — two `data` modules (`D1`, `D2`), two method stages
producing `clustering`, and `group_by: data` — the four members partition
into two groups, one gather node each:

```
members (producers of `clustering`, ids abbreviated)   group (ancestor `data` module)
D1-method_a.default ─┐
D1-method_b.default ─┴─▶ metrics-summarizer-D1.default   ({data} = D1)
D2-method_a.default ─┐
D2-method_b.default ─┴─▶ metrics-summarizer-D2.default   ({data} = D2)
```

Each gather node keeps only its group key: `{data} = D1` is usable in
templates and matched by downstream `requires:`, while *which methods* fed
the node is forgotten from gating — it survives only in `gathered_from`
(below).

**Gating lineage** — what `requires`, `exclude`, and path templates reason
over. Note the label *system* itself is 008 work, scoped for api 0.7 and not
landed with this PR. The MVP only binds the group key into the existing
internal `TemplateContext.provides` — the same container that carries the
builtin `dataset`/`name` today — which templates and the already-wired
`requires:` match against; the rest below is the intended semantics once 008
lands. For a gather node the gating record is exactly the group key:

- `TemplateContext.provides` binds one label: the `group_by` stage id → the
  group value (plus the builtin `name` = the gather module's own id).
  Downstream `requires:` matches against this and nothing else.
- Output path templates may reference `{<group_by stage>}` (the group value)
  and `{name}`/`{params.*}`. Referencing any other label is a **plan-time
  error** (raised in `substitute()`), never an empty substitution.
- `exclude` rules pairing the gather stage's own module with a member's
  lineage drop that member from that module's gather — transitive exclude at
  the member level, so one excluded member never poisons the whole group.
  Beyond that boundary the cut holds: an `exclude` pairing a pre-gather
  module with a *post*-gather module has no chain to be transitive over
  (richer cross-boundary semantics arrive with `where`, Phase 3).

**Provenance lineage** — which concrete nodes fed this node. This is a new
explicit field:

```
ResolvedNode.gathered_from: List[str]     # member node ids, plan order
```

Full ancestry of any post-gather node is recovered by walking `parent_id`
chains as today and, at a gather node, fanning out through `gathered_from`
into each member's chain. In the MVP this record exists only in the resolved
plan (and in the Snakefile's input lists); *persisting* it is Phase 4: the
run metadata (`out/.metadata/`) will record `gathered_from` alongside the
resolved-modules list, so a result file's complete upstream closure — every
module id, commit, and parameter hash that contributed — stays
reconstructable even though no single chain encodes it.

**Identity.** Gather node ids are composed from the group key instead of a
parent chain — `metrics-summarizer-D1.default`, the bare group *value* in the
id — and downstream ids prefix-compose off them as usual
(`metrics-summarizer-D1.default-report-R1.default`). `parent_id` is `None`:
a fan-in has no single ancestor node to point at. Chains resume immediately
downstream — the plan becomes alternating scatter segments and gather cut
points (series-parallel), never a general DAG. Bare group values stay
unambiguous in ids only while ids cannot contain the separators — the
id-separator reservation, deferred with layout v2 (§3.8); the full identity
story (labels, fan-in hashes) is §3.9.

**New filesystem prefix (required).** Because the chain is cut, a gather node
has no parent directory to extend, so the stage must declare `prefix:` (§3.2)
naming where the fresh chain begins; outputs land under
`<prefix>/<group>/<stage>/<module>/<param>/…`, the existing nested builder
resuming below the prefix. See §3.8 for why this local answer suffices and
layout v2 stays orthogonal.

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

Membership is computed **after** all pruning (stages marked 008/`where` are
later phases; in the MVP only `excludes` prune):

```
excludes → requires/label gates (008) → capability gates (008) → gather membership → where → group_by partition
```

Consequences, all logged per 008's diagnostics rules:

- A pruned lineage simply never appears in any gather list; an `exclude`
  naming the gather module itself additionally drops the matching members
  from that module's gather (member-level filtering, MVP).
- `where:` drops are counted per group (Phase 3).
- A group that ends up **empty** is a warning (and the gather node for that
  key is not created). In the MVP *all* groups empty is also a warning;
  promoting it to an error lands with `where` (Phase 3), when silent
  emptiness becomes likelier.

### 3.5 Internal representation changes

| Where | Change |
|---|---|
| `model/benchmark.py` | `Stage.gather: Optional[List[GatherSpec]]`; `Stage.prefix: Optional[str]` (**required when `gather` is set** — the filesystem prefix for the cut chain, §3.3); `GatherSpec {from_: str (alias `from`), group_by: str (a stage id)}`. Both fields required; `where` deferred. Gated on api ≥ 0.7.0 (mints `0.7.0`, shared with 008 and 009; 0.6 is left unminted). `group_by` validated to name a real stage. Duplicate output ids across stages are already legal (the `get_outputs()` dict collapses them, so the uniqueness check never fires) — now documented as the intended contract. |
| `model/resolved.py` | `ResolvedNode.gathered_from: List[str]` (member node ids, plan order); `is_gather` already exists — it finally gets set. Gather inputs reuse the collector's enumerated `input_N` keys + `input_name_mapping`, which `_write_gather_shell` already groups back by flag name. |
| `cli/run.py` | Sibling `_expand_gather_stage()` reached by an `if stage.gather:` dispatch (not threaded through the scatter loop). Collect producers of `from` from `output_to_nodes`; partition by the ancestor module of the `group_by` stage (walk `parent_id` via `nodes_by_id`); emit one node per group × module × parameter. `parent_id=None`; outputs written under `prefix/<group>/…` and registered in `output_to_nodes` so downstream stages consume them. |
| `backend/snakemake.py` | No change: `_write_gather_shell` already exists (previously only reachable via collectors) and gather nodes reuse it as-is, enumerated `input_N` keys included. Shell contract: `--<from_id> path1 path2 …` (space-separated, as collectors today). Migrating to Snakemake-native named list inputs is a possible later cleanup. |

### 3.6 Module CLI contract for gather modules

A gather module receives, per gather entry, one flag named after the output
id with all member paths:

```
<entrypoint> --output_dir … --name summarizer \
    --clustering /abs/p1.tsv /abs/p2.tsv /abs/p3.tsv \
    [param flags]
```

Outputs follow the existing api ≥ 0.5 contract: no `--output` flags are
emitted; the module writes its declared output filenames into `--output_dir`
(exactly what `_write_gather_shell` produces for collectors today). Design
009's named-output rules (explicit `--output` flags) are not yet merged; when
009 lands at api 0.7, its rules apply to gather nodes as to any other node —
same gate, no compatibility split.
Member order is plan order — the order the planner created the member nodes:
stage *topological* order (declaration order breaks ties between independent
stages, the sort being stable), then module declaration order, then
input-lineage × parameter expansion order. The order is arbitrary but
deterministic: it is baked into the emitted command line, and identical YAML
must produce byte-identical Snakefiles or Snakemake reruns unchanged work.
Modules must not read meaning from positions — identity comes from the paths.

### 3.7 Metric collector deprecation

`metric_collectors` is re-expressible as a *global* gather — a gather with
no `group_by`, one node collecting every producer. `group_by` is required in
the MVP, so this form arrives with the reimplementation itself. Plan:

1. When gather lands, the collector resolver is reimplemented *on top of*
   gather nodes (`is_collector` kept for output-layout compat).
2. A deprecation warning on `metric_collectors:` points at the equivalent
   `gather:` stanza.
3. Removal in a later major api version, not scheduled here.

### 3.8 Output layout: `prefix` is enough (layout v2 is orthogonal)

Gather forces a path-layout decision — the current builder appends
`<stage>/<module>/.<hash>/` to the parent's directory, and a gather node has
no parent directory to extend. The MVP answers it locally: the required
`prefix:` names a fresh root, the group value is the first segment under it,
and the existing nested builder resumes unchanged:
`<prefix>/<group>/<stage>/<module>/<param>/…`. Nothing else about the layout
changes; existing trees are untouched.

The **layout-v2** un-nesting — one flat `stage.module.hash/` segment per node
instead of today's three components, depth N+1 instead of 3N, stage-to-stage
nesting untouched (an unapproved cleanup note in `scratch/`, not a design; 007
remains the layout spec of record) — is orthogonal: it would make gather a
non-special case (a parentless segment is just another segment), but nothing in this
design depends on it, and the overly deep nesting is a pre-existing cost, not
a gather cost. It remains a standalone cleanup on its own api bump; 007
(output layout) gets a v2 section if and when it lands. The one piece gather
will eventually want from it is the id-separator reservation (no dots in
stage/module ids) so group values stay unambiguous inside node ids — noted in
§3.3 and deferred with the rest.

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

- **Lineage becomes explicit edges.** `ResolvedNode.parents: List[str]` holds
  the direct parent ids of a fan-in (join or gather) — the full producer set.
  Root and linear nodes leave it empty: their single parent already lives in
  `parent_id` and the id prefix. Ancestry recovery walks `parents` (unioned with the
  id-prefix walk for the linear spine), which removes the "downstream of a
  fan-in sees only the primary branch" limitation. Persisting the full record
  — a **`lineage.json` sidecar** (a file for now; a DB later) — is deferred to
  Phase 4 with the rest of provenance; the MVP lands only the in-memory field.

  Formally the lineage is a directed **hypergraph** (a B-graph): each
  derivation is one hyperedge `{parents} → child`, consumed *together*, not a
  bundle of independent pairwise edges. Storing it per-node loses nothing
  because of an invariant the planner maintains: **one node = one
  derivation**. A node never has alternative input bundles — the expansion
  materializes every combination as its own node (`_select_input_bundles`
  emits one node per bundle, gather one node per group), so
  `(node.parents, node.id)` *is* the hyperedge and upward traversal is
  unambiguous: everything reachable genuinely contributed. Two caveats keep
  this sound: which parent supplied which file lives in
  `input_name_mapping`/`gathered_from` (not in the edge set), and the linear
  spine is walked via `parent_id` (not id-prefix matching), so ids containing
  the separator cannot corrupt ancestry. Note this traversal answers *provenance*
  questions only; gating deliberately does not walk through a gather cut
  (§3.3 — the partition forgets, so e.g. one excluded member among hundreds
  gathered must not poison every downstream node).
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

Phased; each phase lands separately, smallest working piece first.

### Phase 1 — structural gather + fan-in join (api ≥ 0.7.0) — **this PR**

- **Model**: `GatherSpec {from, group_by}` (both required); `Stage.prefix`
  required with `gather`; `group_by` must name a real stage; duplicate output
  ids across stages legalized (`get_stages_by_output` replaces the
  single-producer lookup).
- **Expansion**: `_expand_gather_stage()`, a sibling of the scatter path
  reached by one `if stage.gather:` dispatch. Members = producers of `from`;
  partition by the ancestor module of the `group_by` stage; one node per
  group × module × parameter; `parent_id=None`, `parents`/`gathered_from`
  set; outputs under `prefix/<group>/…`, registered for downstream stages.
- **Join** (§3.9) rides along: `_select_input_bundles` pairs producers on
  divergent branches that share a lineage root; one node per bundle.
- **Backend**: reuses the collector's `_write_gather_shell` unchanged
  (enumerated `input_N` keys + name mapping).
- **Refactor, mechanical**: the scatter body is extracted to
  `_expand_scatter_stage()` and its loop-local helpers hoisted to module
  level, so gather is a sibling rather than a deeper nesting level. No
  behaviour change; the e2e suite is the regression net.

Phase 1 is deliberately self-contained — it needs neither 008 nor 009.
Grouping is structural (walks `parent_id`; no `provides` labels), membership
uses the already-id-keyed `output_to_nodes`, and the only 009-flavoured piece
is legalizing duplicate output ids, not the `output_name_mapping` rewrite.
Layout v2 is not a prerequisite either: `prefix:` answers the path question
locally (§3.8).

### Phase 2 — label-based `group_by` (needs 008)

Group by a `provides` label instead of a stage — the general form of Phase
1's structural rule, reading values from `TemplateContext.provides`. Enables
tuple keys and per-parameter grouping. The deprecated `dataset` builtin is
not a grouping axis.

### Phase 3 — `where` + cross-boundary `exclude` (needs 008)

Member filtering; unsatisfiable cross-gather exclude rules reported;
empty-group warning, all-groups-empty error.

### Phase 4 — provenance + collector deprecation

Persist `gathered_from` (`out/.metadata/`, `lineage.json` sidecar, §3.3/§3.9);
reimplement collectors on top of gather; deprecation warning on
`metric_collectors:`.

### Landing order

The priority inversion: 010 lands first in *code*; 008 and 009 sequence in
behind it, all sharing api `0.7.0` (this PR mints it). 009's named-output
representation is folded in when it lands; 008's labels unlock Phases 2–3.

### Deferred (add when a benchmark asks)

- Plain `inputs:` mixed with `gather:` on one stage (which chain resolves the
  plain input? needs the group key to pin it).
- Multiple `gather:` entries with *different* `group_by` axes on one stage
  (cross-product of groups; unclear demand).
- Wildcard-based Snakefile emission (O(stages×modules) rules instead of
  O(nodes)) — a backend-emitter rewrite that loses the diffable Snakefile
  artifact (007) and buys nothing for gather; separate proposal if node
  counts ever make explicit Snakefiles hurt.

### Testing Strategy

- Unit: membership selection, `group_by` partitioning, group-key id
  composition — asserting on resolved node sets, not Snakefile text (008's
  convention). `where` filtering joins in Phase 3.
- Integration: e2e fixture with two method stages sharing an output id, a
  grouped gather, and a downstream consumer (scatter → gather → scatter).
- Provenance: walk `gathered_from` from a terminal node and assert the full
  upstream closure matches the plan.

## 6. References

010 lands ahead of 008 and 009 (they sequence in with Phases 2–3, see §5), so
the relative links below dangle on main until their PRs merge; the PR links
are the live copies meanwhile.

1. [Issue #289 — support for multi-stage outputs](https://github.com/omnibenchmark/omnibenchmark/issues/289)
2. [Design 008 — filtering and gating mechanisms](008-filtering.md) — in flight on [PR #354](https://github.com/omnibenchmark/omnibenchmark/pull/354)
3. [Design 009 — named outputs](009-named-outputs.md) — in flight on [PR #329](https://github.com/omnibenchmark/omnibenchmark/pull/329)
4. [Köster et al., Snakemake paper, Fig 8a scatter-gather](https://f1000research.com/articles/10-33/v3#f8)
5. [PR #291 — earlier gather proposal](https://github.com/omnibenchmark/omnibenchmark/pull/291)
