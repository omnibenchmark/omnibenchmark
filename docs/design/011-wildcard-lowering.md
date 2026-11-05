# 011: Wildcard lowering and its prerequisites

[![Status: Draft](https://img.shields.io/badge/Status-Draft-yellow.svg)](https://github.com/omnibenchmark/docs/design)
[![Version: 1](https://img.shields.io/badge/Version-1-blue.svg)](https://github.com/omnibenchmark/docs/design)

| Field          | Value |
|----------------|-------|
| Authors        | btraven00 |
| Date           | 2026-08-13 |
| Status         | Draft |
| Version        | 1 |
| Supersedes     | N/A |
| Reviewed-by    | TBD |
| Related Issues | [#289](https://github.com/omnibenchmark/omnibenchmark/issues/289) (multi-stage input collection / diamonds), [#360](https://github.com/omnibenchmark/omnibenchmark/issues/360) (filtering & slicing), [#201](https://github.com/omnibenchmark/omnibenchmark/issues/201) (per-backend rule duplication) |

## Changes

| Version | Date       | Description   | Author    |
|---------|------------|---------------|-----------|
| 1       | 2026-08-13 | Initial draft | btraven00 |
| 2       | 2026-11-05 | Tracked in the repo; §8: what 013 §4.1 implements of §3.6 | ben |

## 1. Problem Statement

`ob run` lowers a benchmark plan to a **fully materialised Snakefile: one rule per
resolved node**, with every path concrete. `SnakemakeGenerator.generate_snakefile`
(`omnibenchmark/backend/snakemake.py:72-78`) writes one `_write_node_rule` per
`ResolvedNode`, and the node list itself is produced by the expansion loop at
`omnibenchmark/cli/run.py:1236-1458`.

That loop is a **reimplementation of Snakemake's DAG builder in Python**: it walks
lineages, materialises the cartesian product (`itertools.product` over
`(upstream nodes × parameter combinations)`, `run.py:1289`), selects a single
parent per node, and gathers ancestor outputs. Every structural limitation we
carry lives in that reimplementation, not in Snakemake:

- **Emission cost is O(nodes).** The 30-module scRNA plan produces ~26k rules
  (recorded in `scratch/wildcard_backend_and_labels.md`). Rule text, node ids,
  and log filenames all scale with the full expansion.
- **Ancestor resolution is ours to get right.** The `get_ancestor_nodes`
  2^depth walk was a real hang (`fix/ancestor-walk-exponential`); the O(N) node
  scans were a real regression. Neither is a problem Snakemake has.
- **No joins.** `select_input_nodes` (`core/_lineage.py:10-65`) collapses a
  stage's inputs to a single deepest parent stage, and
  `model/validation.py:181-236` exists solely to *reject* the diamonds the tree
  representation cannot express (#289).
- **No native partial execution.** `_run_snakemake` passes no targets at all
  (`run.py:440-465`); the run relies entirely on `rule all` +
  `default_target: True`. Incremental / staged submission (the parked SLURM
  direction) needs `snakemake <target>` to mean something.

**Filtering is explicitly *not* a motivation.** Design 008 §4 Alternative 4 gives
two reasons for keeping semantic filtering above the lowering step; this document
changes only the second (*"the generated Snakefiles are explicit, not
wildcard-driven … there are no wildcards for Snakemake's pattern machinery to act
on"*). The first reason — that `stage_id`, `capability` and lineage labels have
nothing to bind to below the lowering step — stands unchanged, and pruning stays
in the resolver under every option here.

There is also a documentation inconsistency to settle: **`004-yaml-specification.md`
§3.8 "Wildcard System" is still normative and still describes a wildcard
implementation that no longer exists** (wildcard naming from output templates,
module ids as wildcard values, propagation creating cartesian products). The spec
says wildcards; the code says explicit rules.

## 2. Prior art: the wildcard backend we deleted

This is not a greenfield proposal. A wildcard-based generator existed and was
removed in `1a88b0e1` (2026-03-30, *"feat: explicit Snakefile generation"*), which
deleted `omnibenchmark/workflow/` wholesale. The source is still recoverable:

```
git show 1a88b0e1^:omnibenchmark/workflow/snakemake/rules/rule_node.smk
git show 1a88b0e1^:omnibenchmark/workflow/snakemake/scripts/run_module.py
```

It used broadly the shape proposed in §3.6 below — wildcards `{input}` (the whole
upstream prefix), `{stage}`, `{module}`, `{params}`, `{dataset}`, pinned per rule
with `wildcard_constraints`, and inputs supplied by an input function
(`formatter.format_input_templates_to_be_expanded(benchmark, wildcards)`).

It was, in the words of the person who replaced it, *very buggy*. The causes are
worth recording precisely, because **none of them is "wildcards"**:

1. **`script:` instead of `shell:` — the dominant bug source.** Every rule ran
   `script: run_module.py`, a 119-line shim (plus a 224-line `execution.py`)
   executing *inside* the Snakemake job. Per job it performed:

   | What it did at job runtime | Why that is wrong |
   |---|---|
   | `clone_module(repository_url, commit_hash)` | a **git checkout per job**, under concurrency |
   | `extract_module_name(snakemake.rule)` | module identity by string surgery on the rule name |
   | `build_inputs_map(template, snakemake.input)` | input wiring rebuilt from a template at runtime |
   | `output_dir = os.path.commonpath(snakemake.output)` | output directory *inferred* rather than known |

   Live Python objects (`Params` instances, maps, lambdas) crossed the `params:`
   boundary. Every row above re-derives, N times and remotely, something the
   resolver already knew exactly once.

2. **Rule bodies quadruplicated per software backend.** Because `conda:` /
   `envmodules:` / `container:` are directives rather than values, each rule was
   written out four times — ~400 lines of near-identical `rule:` blocks, with a
   `TODO #201` sitting on top of every copy.

3. **A dynamic `.smk` executing Python inside Snakemake**, not an emitted file.
   There was nothing to read, diff, archive, or run standalone. 007 §3 later made
   *"Fully materialised Snakefile; all wildcards resolved to concrete paths. Can
   be run standalone"* an explicit product promise — a direct reaction to this.

4. **The path grammar re-implemented as a regex.** Collector inputs were selected
   by `_compile_regex_pattern_for_collectors_input`, which rebuilt the output-path
   layout as a regular expression. That duplication still exists today in a
   different form (§3.1).

### 2.1 Constraints inherited from that failure

These are stated as **hard constraints on any future wildcard lowering**, not
preferences:

- Resolution — git checkout, entrypoint detection, environment resolution —
  happens **once, before Snakemake starts**, as it does today via `ModuleResolver`.
- Rules use `shell:` only. Never `script:`.
- `params:` carries strings. Never live Python objects.
- Identity is never re-derived at job runtime from `snakemake.rule` or from
  `commonpath(snakemake.output)`.
- The generator **emits a text file**. Rules are not constructed by executing
  Python inside Snakemake.

### 2.2 The observation that makes this tractable

**The current `shell:` block is already wildcard-ready.** `_write_shell`
(`backend/snakemake.py:238-327`) is written entirely in terms of
`{params.module_dir}`, `{params.entrypoint}`, `{params.cli_args}`,
`{params.output_dir}`, `{input.*}`, `{log}` and `{rule}`. **Nothing in the shell
text names a node.** The only per-node literals in a rule are the *values*
assigned in `params:`, plus the `parameters.json` echo and the human-name symlink.

So the thesis of this document is narrow:

> Keep emitting a text file of `shell:` rules. Change only the **granularity**
> (node → stage×module) and turn the per-node literals into **wildcard-keyed
> lookups**.

Wildcards were never the problem. Runtime re-derivation, dynamic rule
construction, and a missing path abstraction were.

## 3. Proposed Solution

The proposal in this version is **the five prerequisites (§3.1–§3.5)**, each
landable on its own and each with a payoff that does not depend on wildcards ever
shipping. §3.6 sketches what the lowering looks like once they are in place; it is
a payoff description, not a commitment, and is deliberately not scheduled.

### 3.0 Why prerequisites at all

A wildcard rule *is* the output-path grammar written as a pattern. Today that
grammar is an f-string in one loop and a positional parser in five other places:

| Concern | Where it lives today | Why it blocks a pattern |
|---|---|---|
| Output path construction | `cli/run.py:1388-1406` — `f"{base_path}/{stage.id}/{module_id}/{param_id}/{tmpl}"` | grammar exists only as an inline f-string |
| Positional re-parsing | `backend/snakemake.py:144,306` (`parts[1]` = dataset); `core/_paths.py:109` `normalize_dataset_path`; `backend/_metric_collector.py:403`; `core/_node.py:107` | five parsers to keep in sync with any pattern change |
| Three dirs per lineage step | `<stage>/<module>/.<hash>/`, depth 3N | a prefix wildcard cannot be split back into segments without knowing N |
| Single parent | `select_input_nodes` + `base_path` = parent of the *deepest* input (`run.py:1372-1378`) | the pattern is chain-derived, not stage-local |
| Per-node text | `cli_args`, `parameters.json`, `--name`, env, resources baked into each rule | O(nodes) text; needs wildcard-keyed lookup |
| `--name` from path position | `backend/snakemake.py:292-314` (API ≤ 0.4) | hard-codes layout v1 into generated shell text |
| Rule name as job identity | `sanitize_rule_name(node.id)`, consumed by `telemetry/emitter.py:263`, `core/_paths.py:66-76`, and the progress parser at `run.py:507` | wildcard rules share one name across many jobs |

### 3.1 S1 — One canonical layout module

New `omnibenchmark/core/layout.py`, the single owner of the output-path grammar:

```python
def segment(stage_id: str, module_id: str, param_hash: str) -> str: ...
def node_path(lineage: Sequence[Segment], filename: str) -> str: ...
def parse(path: str) -> list[Segment]: ...        # inverse of node_path
```

One production site to convert (`run.py:1388-1406`). Every positional parser in
the table above becomes a `parse()` call.

*Payoff without wildcards:* the layout stops being tribal knowledge spread over
six files, and layout v2 (§3.2) becomes a one-file change instead of an audit.

**This is the prerequisite for everything else.** A wildcard pattern is then just
one more consumer — `pattern(stage_id)` alongside `node_path`.

### 3.2 S2 — Layout v2 (one directory per lineage step)

Already scoped in `scratch/cleanups_layout_v2_and_run_refactor.md`, gated on api
0.7.0: collapse `<stage>/<module>/.<hash>/` (three components per step) to a
single `stage.module.hash/` segment — depth N+1 instead of 3N.

*Why it matters here:* it makes the prefix wildcard **self-delimiting**.
`lineage.split("/")` recovers the ancestor record directly, so an input function
can resolve any ancestor's output from the path alone. Under the three-directory
layout the same split requires knowing the depth in advance.

*Payoff without wildcards:* shorter paths, which reduces reliance on
`truncate_path_filename` (`core/_paths.py:100`) and the collisions it papers over.

### 3.3 S3 — Parameters as an opaque wildcard plus a lookup table

`param_id` is already an opaque `.{hash8}` (`run.py:1321`); parameter *values*
reach the module only as `cli_args` and `parameters.json`. Move those into a table
emitted once at the top of the Snakefile:

```python
PARAMS = {
    ("methods", "M1"): {
        "317a5066": {"cli_args": "--a 0 --b 0.1", "json": '{"a": "0", "b": "0.1"}',
                     "human": "a-0_b-0.1"},
        ...
    },
}
```

read back as `params: cli_args=lambda w: PARAMS[("methods","M1")][w.phash]["cli_args"]`.

*Payoff without wildcards:* Snakefile text drops from O(nodes) to O(distinct
parameter sets), and the `parameters.json` escaping dance plus the human-name
symlink (`backend/snakemake.py:270-278`) leave the shell string, where they are
currently brace-escaped by hand.

### 3.4 S4 — Stop deriving identity from path position

Delete the API ≤ 0.4 `--name` workaround (`backend/snakemake.py:292-314`), which
recovers the dataset name from `first_output.split("/")[1]`, and let design 009's
named outputs make stage output filenames lineage-independent: **the directory
carries identity, the filename does not.**

*Payoff without wildcards:* removes the last consumer of positional path parsing
from generated shell text, and unblocks S2 (which changes exactly the position
that workaround indexes into).

### 3.5 S5 — Job identity is the output path, not the rule name

Rework the three consumers that key on `sanitize_rule_name(node.id)` —
`telemetry/emitter.py:263`, the `describe status` log lookup at
`core/_paths.py:66-76`, and the `rule (\w+):` progress parser at `run.py:507` — to
key on the first output path instead.

*Payoff without wildcards:* node ids grow with lineage depth, so rule names and
log filenames already hit `truncate_filename`'s 255-byte cap and its
hash-suffix collision handling; an output path is a stable identity that does not.

**This is the one prerequisite whose absence silently breaks a wildcard backend**,
since many jobs would then share a single rule name. It is listed last because it
is the only one motivated more by the destination than by present-day pain.

### 3.6 What the lowering would then look like (sketch — not scheduled)

One rule per (stage, module): its own segment literal, the upstream prefix and its
own parameter hash as wildcards.

```python
rule methods_M1:
    input:  lambda w: inputs_for("methods", "M1", w)   # splits w.lineage via layout.parse
    output: "{lineage}/methods.M1.{phash}/result.json"
    wildcard_constraints:
        phash  = "|".join(PARAM_HASHES["methods", "M1"]),
        lineage = LINEAGE_RE["methods"]
    params:
        module_dir  = ".modules/methods/3b17081",
        entrypoint  = "entrypoint.py",
        cli_args    = lambda w: PARAMS["methods", "M1"][w.phash]["cli_args"],
        output_dir  = lambda wildcards, output: os.path.dirname(output[0]) or ".",
    conda: "envs/methods.yaml"
    shell: """ ...unchanged from today... """
```

Properties worth stating explicitly:

- **Rule count is O(stages × modules)** — tens, not tens of thousands. The module's
  own segment stays literal, which is also what keeps rules unambiguous.
- **All pruning stays in the resolver.** `exclude`, `requires`, capability gating,
  `--until`, and saved filters continue to run above the lowering step; the result
  is expressed as the **target list in `rule all`**, and Snakemake demand-drives
  the intermediates. Nothing semantic moves below the boundary, so 008 §4 Alt 4
  reason (1) is untouched.
- **The target list is where `exclude` and `provides`/`requires` differ.** With
  lineage labels (008 Phase 2) the surviving set is a product of independent
  label axes and could be written as an `expand()`. With `exclude` — a negative,
  pairwise, *cross-field* constraint — it must stay enumerated, because
  `wildcard_constraints` are per-wildcard regexes and cannot express a joint
  condition across two wildcards. See `scratch/wildcard_backend_and_labels.md`.
- **Diamonds (#289) are not fixed for free.** The wildcard form makes joins
  natural to *express*, but `select_input_nodes` still picks one parent; the
  resolver change is needed either way.

## 4. Alternatives Considered

### Alternative 1: keep explicit lowering, do nothing
- **Description**: leave the one-rule-per-node generator in place.
- **Pros**: the Snakefile stays a fully materialised, diffable, standalone
  artifact (007 §3); obeditor and the dashboard keep a concrete job list; zero
  risk.
- **Cons**: we own the DAG builder, with its performance characteristics and its
  #289 limitation, and partial execution stays unavailable.
- **Reason for rejection**: not rejected as a *default*. This is the honest status
  quo, and §3.1–§3.5 improve it whether or not §3.6 ever happens.

### Alternative 2: wildcard lowering without the prerequisites
- **Description**: write a wildcard emitter directly against the current layout.
- **Cons**: the emitter would have to embed the path grammar a seventh time, work
  around the 3N depth, and collide with rule-name-keyed telemetry — which is
  approximately how the `1a88b0e1` version accumulated its bugs.
- **Reason for rejection**: it converts five small, independently valuable changes
  into one large risky one.

### Alternative 3: a different engine (Nextflow)
- **Description**: lower to Nextflow, whose channel model is closer to the
  label-bag representation. Sketched in `scratch/nextflow-backend.md`.
- **Reason for rejection**: orthogonal and far larger. The prerequisites here
  benefit any backend, since they are about our own representation.

### Alternative 4: wildcards as an additional lowering behind a flag
- **Description**: keep both emitters, select with `--emit explicit|wildcard`.
- **Pros**: preserves the 007 artifact guarantee for provenance while allowing the
  wildcard form at scale.
- **Status**: the likely shape if §3.6 is ever scheduled; recorded here, not
  decided.

## 5. Implementation Plan

Only §3.1–§3.5 are proposed for implementation. Each is a separate PR.

1. **S1 — `core/layout.py`** and conversion of the six sites. Pure refactor; e2e
   fixtures are the safety net. **Blocks the rest.**
2. **S2 — layout v2**, api 0.7.0 gate (rides along with the gather work per
   `scratch/cleanups_layout_v2_and_run_refactor.md`).
3. **S3 — parameter lookup table** in the emitted Snakefile.
4. **S4 — drop the API ≤ 0.4 `--name` path workaround** (coordinate with 009).
5. **S5 — output-path job identity** for telemetry, status, and progress.

S2, S3 and S4 are independent of one another once S1 has landed.

### Open questions to settle before scheduling §3.6

- **Provenance.** 007 §3 promises a materialised, standalone Snakefile; a
  wildcard Snakefile cannot be one. The natural mitigation is to move that role to
  a `plan.json` — `ResolvedNode.to_dict()` already exists and is unused for
  persistence, and `save_metadata` currently writes only `benchmark.yaml` and
  `modules.txt`. Not proposed here.
- **Gather and collectors.** Superseded by §7 — design 010 landed, `is_gather` is
  now set, and gather turns out to be the easy case rather than an open question.
  Collector lowering follows gather's once collectors are reimplemented on top of
  it (010 §3.6).
- **Per-backend duplication (#201).** The four-way env-directive split that
  bloated the deleted implementation is handled today by writing one directive per
  concrete rule; a wildcard rule must still choose one at emission time.
- **Test surface.** `tests/backend/test_snakemake_gen*.py` and
  `tests/cli/test_run_benchmark.py:254,340` assert on Snakefile text — the latter
  asserting that braces are *absent*, precisely because they would be read as
  wildcards.

### Regression guard

Any future wildcard emitter that reaches for `script:`, performs git operations
inside a job, or infers identity from `snakemake.rule` or
`commonpath(snakemake.output)` is reproducing the `1a88b0e1` failure and should be
rejected in review on that ground alone (§2.1).

### Testing Strategy

- **S1**: unit tests for `node_path` / `parse` round-tripping, including the
  `.default` parameter case and templated filenames; existing e2e fixtures must
  produce byte-identical output paths.
- **S2**: api-gated — 0.6 benchmarks keep layout v1, 0.7 benchmarks get v2; assert
  both in the same test run.
- **S3**: assert the emitted Snakefile contains one `PARAMS` entry per distinct
  parameter set and that `cli_args` resolves identically to today.
- **S5**: assert telemetry and `describe status` resolve a job from its output
  path when two nodes share a rule name.

## 6. References

1. `008-filtering.md` §4 Alternative 4 — why semantic filtering stays above the lowering step (not yet in-tree; branch `feat/filtering-phase1`)
2. [007-output-layout.md](./007-output-layout.md) §3 — the materialised-Snakefile guarantee
3. [004-yaml-specification.md](./004-yaml-specification.md) §3.8 — the wildcard system as still specified
4. `009-named-outputs.md` — lineage-independent output filenames (not yet in-tree; branch `feat/named-outputs`)
5. `scratch/wildcard_backend_and_labels.md` — why `exclude` is wildcard-hostile and `provides`/`requires` is native
6. `scratch/cleanups_layout_v2_and_run_refactor.md` — layout v2 and the parked backend note
7. Commit `1a88b0e1` — *"feat: explicit Snakefile generation"*, which deleted the wildcard backend
8. [Issue #289](https://github.com/omnibenchmark/omnibenchmark/issues/289) — multi-stage input collection

## 7. After design 010 landed: what changes, and the simplest way forward

Written after the gather + fan-in join work merged. Three things this note
assumed are now settled, and they move the recommendation.

### 7.1 Gather lowers cleanly. The join does not.

The intuition was backwards. Gather looks like the hard construct because it cuts
the chain — but a gather node's path is
`<prefix>/{group}/<stage>/<module>/{phash}/`, `{group}` is a real wildcard, and
its members are recoverable from it: every member descends from the `group_by`
ancestor, so that ancestor's segment appears in each member's path. Gather inputs
are `expand()` over (stage, module, phash) with `{group}` pinned. That is the
idiom, not a workaround. The deliberate forgetting in 010 §3.3 also stops being a
rule we enforce: a gather node's path carries only `{group}`, so downstream rules
can only unify on `{group}` anyway.

The **join** is the one that resists. Its path carries only the *deepest*
branch's prefix; the other parents appear nowhere in it, only as the `join_hash`
digest on the segment. So `inputs_for(w)` cannot derive the partner paths from
wildcards and must look them up by `w.jhash`. Snakemake would be calling back into
our resolver for exactly the construct that motivated the lowering — the
unification never happens.

This splits §3.6 in two. Under the `{lineage}`-as-opaque-prefix sketch above,
joins stay table-driven and only scatter and gather benefit. Under the label-bag
model (`scratch/wildcard_backend_and_labels.md`), both branches' labels are *in
the path*, the partner path is derivable, and a join becomes what that note
claims: one rule with two inputs matched on shared wildcards. **010's join is an
argument for the label-bag model over the conservative lowering**, and that choice
should be made before §3.6 is scheduled rather than during it.

`join_hash` is transitional in exactly that sense: it is standing in for the
labels that are not in the path yet. It exists because identity is still
path-position-derived, which is S4's premise. When labels move into the path, the
digest has nothing left to disambiguate. A related prerequisite already landed: a
join node now inherits the union of its branches' labels rather than only the
anchor's, so its label bag is complete enough to be addressed by a pattern
spanning both branches.

### 7.2 S1 got more valuable, and it is the recommended next step

The table in §3 counts one output-path construction site. There are now three,
all in `core/_expand.py`: the nested scatter path, the flat scatter path, and the
gather's `f"{stage.prefix}/{gval}/{stage.id}/{module_id}/{param_id}/{tmpl}"`. Each
new fan-in shape adds one, and `expansion_segment` has also made the parameter
segment stop being a single axis — a join emits `.{phash}-{jhash}`, which S3
assumes is purely the parameter hash.

So the recommendation is unchanged in content but sharper in priority: **do S1,
now, before anything else in this note.** It is the only item whose cost grows
with every construct added to the DSL, and 010 just added two.

S2's stated coupling is stale: it says layout v2 "rides along with the gather
work", but gather shipped without it — 010 §3.3 answers the path question locally
with `prefix:` and keeps layout v2 orthogonal. S2 is now free-standing and can be
scheduled on its own merits.

### 7.3 Deprecating `exclude`: helps, but is not on the critical path

Worth stating precisely, because the answer is smaller than it looks.

§3.6 already puts all pruning above the lowering: the surviving set is expressed
as the target list in `rule all`. So `exclude` does **not** block a wildcard
backend. What it costs is that the target list stays *enumerated* instead of
collapsing to an `expand()` over independent label axes — the Snakefile shrinks
from O(nodes) rules to O(stages×modules) rules either way, but with `exclude` in
play it keeps an O(nodes) target list, and the resolver still enumerates every
node internally. That last part is where the actual performance win lives, and it
is the only part `exclude` withholds.

Fan-in makes it slightly worse: exclusions now union over every branch of a join
(`lineage_module_ids` walks `parents`), so a cross-branch rule is a joint
constraint across more wildcards than before. `wildcard_constraints` are
per-wildcard regexes and still cannot express it.

**Recommendation: demote, do not deprecate.**

1. Document `provides`/`requires` as the primary conditional-execution axis and
   `exclude` as the escape hatch for genuine incompatibilities ("metric M is only
   defined for clustering C"), which are real and survive any representation.
2. Do not schedule a removal. Full deprecation is not reachable today anyway:
   008 has no negative gate, so "everything except X" has to be enumerated as
   N−1 positive `requires`, which is worse YAML than the `exclude` it replaces.
   A `requires_not` would change that, and is not designed.
3. **Measure before committing.** The pruned-combinations summary 008 already
   ships (`prune_counts`, one line per run: *N pruned, X by requires, Y by
   exclude*) is the instrument. If Y trends to zero across real benchmarks once
   labels are available, deprecation becomes a cheap follow-up decided on data.
   If it does not, the escape hatch was load-bearing and the question is settled
   the other way.

Nothing here needs to happen before S1–S5.

## 8. Implemented so far (2026-11-05)

Design 013 §4.1 implements the §3.6 sketch for the subtree below a `foreach`
producer: one rule per (chain of stages, stage, module, params), each chain
position a `{s_<stage>}` wildcard constrained to `<stage>\.[^/]+`
(`backend/_wildcard.py`). Status of the prerequisites it rests on:

- S2 (layout v2) landed in 007 §3.1 and is what makes the chain addressable.
- S3 is not done: parameters stay in the rule (`cli_args` per rule), so a
  module with several parameter sets keeps one rule per set.
- S5 holds for grouped rules: the log is `.logs/<rule>.<wildcards>.log`, so
  job identity is the instantiated path, not the rule name. Status still
  looks logs up by node id and does not find these (013 §5).
- The `{module.id}` file-name template is the one thing that forces a lookup
  table (`_OB_IN`): an input whose file name depends on an ancestor's module
  is not one pattern across the group. File names equal to output ids would
  remove it; that is a plan-format decision (api bump), not a backend one.
- Pruning never enters the Snakefile: `rule all` enumerates the surviving
  leaves, so `exclude` and `requires` lower identically. §7.3's "demote, do
  not deprecate" stands; dropping `exclude` is a DSL question.
- Joins, gathers, collectors and nodes outside a producer stay one rule per
  node, as §7.1 predicts for joins.
