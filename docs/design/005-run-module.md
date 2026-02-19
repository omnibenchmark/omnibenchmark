# 005: Module-Scoped Run (`ob run -m`)

[![Status: Draft](https://img.shields.io/badge/Status-Draft-yellow.svg)](https://github.com/omnibenchmark/docs/design)
[![Version: 0.1](https://img.shields.io/badge/Version-0.1-blue.svg)](https://github.com/omnibenchmark/docs/design)

**Authors**: ben
**Date**: 2026-02-19
**Status**: Draft
**Version**: 0.1
**Supersedes**: N/A
**Reviewed-by**: TBD
**Related Issues**: TBD

## Changes

| Version | Date | Description | Author |
|---------|------|-------------|--------|
| 0.1 | 2026-02-19 | Initial draft | ben |

## 1. Problem Statement

`ob run benchmark.yaml` always runs the full benchmark DAG: every stage, every
module, every parameter expansion, every input combination. This is appropriate
for production but wasteful during **module development**, where a developer
wants to verify that a single module works correctly before the full benchmark
is executed.

The common developer workflow:

1. Write or modify a module (e.g. a new clustering method `M3` in stage `methods`).
2. Run it against *some* representative upstream data to check it produces valid
   outputs.
3. Iterate until satisfied, then submit for a full benchmark run.

Running the full benchmark for this costs unnecessary compute and time, and may
not even be possible if other modules or environments are unavailable locally.

## 2. Design Goals

- **Fast feedback loop**: Generate and run the minimum DAG sub-graph needed to
  exercise a single module.
- **Reproducible path structure**: Paths generated in module-run mode must be
  identical to paths that would be generated in full-run mode for the same
  first-expansion combination, so outputs can be reused.
- **No new YAML keys**: This is purely a CLI/compiler concern; the benchmark
  YAML is unchanged.
- **Simple, predictable heuristic**: The user should be able to reason about
  exactly which nodes get included.

### Non-Goals

- Selecting *which* upstream data to use (see Future Work).
- Automatically download missing upstream data from a remote source.
- Running a module against *all* upstream combinations (use `ob run` for that).
- Cross-module comparison within a single invocation.
- Gathering/aggregation stages when the target module belongs to a gather stage
  (not supported in v1 — the gather node requires all provider outputs to
  exist).

## 3. Proposed Solution

### 3.1 CLI Change

We reuse the option `-m / --module` in `ob run`, as in the `0.4.0` release:

```
ob run benchmark.yaml -m <module_id>
```

`<module_id>` must exactly match the `id` field of a module declared in the
benchmark YAML. The option is mutually exclusive with nothing (it can be
combined with `--dry`, `--dirty`, `--cores`, etc.).

### 3.2 DAG Pruning Heuristic

When `-m <module_id>` is given, the compiler applies two pruning steps before
Snakefile generation.

#### Step 1 — Stage ceiling

Find the stage `S` that contains module `<module_id>`.

Include **only** stages whose index in the benchmark YAML is ≤ index(`S`).

Stages after `S` are dropped entirely. This means downstream consumers of `S`'s
outputs are not built.

**Rationale**: The developer wants to verify their module, not its consumers.

**Limitation**: If `S` is a gather stage (`stage.is_gather_stage() == True`)
module-run mode is not supported and the command exits with an error. Gather
nodes require *all* provider outputs to be present (they are the `rule all`
target), so running a single gather module in isolation is not meaningful
without a separate "bless a dataset" mechanism (see Future Work).

#### Step 2 — First expansion only

Within the included stages, each stage is expanded as normal *except* that only
the **first** input-node × parameter combination is kept for each module.

Concretely, for any module `m` in stage `S'` (where `S'` ≤ `S`):

- Take `input_node_combinations[0]` (the first upstream node, in document
  order).
- Take `params_list[0]` (the first expanded parameter combination, in document
  order).
- Keep exactly that one `(input_node, params)` pair; discard the rest.

For the **target module** (`<module_id>`) in stage `S`:
- Apply the same first-expansion rule.
- Additionally, keep only the node whose `input_node` is the first node
  produced by the immediately preceding included stage.

**Rationale**: The goal is a *smoke test*, not an exhaustive cross-product. The
first combination is deterministic (document order), requires no user
configuration, and produces a minimal but runnable sub-graph. Paths produced
are identical to the paths full-run would have generated for the same
combination, so if the developer inspects or reuses the output it is consistent.

**Assumption**: "First" means index 0 in document order for both input nodes
and parameter expansions. This is stable and predictable.

### 3.3 Implementation Sketch

```python
# In _generate_explicit_snakefile():

if module_filter:                          # set by -m CLI option
    # Identify target stage index
    target_stage_idx = None
    for idx, stage in enumerate(benchmark.model.stages):
        if any(m.id == module_filter for m in stage.modules):
            target_stage_idx = idx
            break
    if target_stage_idx is None:
        sys.exit(f"Module '{module_filter}' not found in benchmark")
    target_stage = benchmark.model.stages[target_stage_idx]
    if target_stage.is_gather_stage():
        sys.exit("ob run -m is not supported for gather stages (see docs/design/005)")

    # Prune: only stages[0..target_stage_idx]
    stages_to_expand = benchmark.model.stages[:target_stage_idx + 1]
else:
    stages_to_expand = benchmark.model.stages

# During node expansion, for each stage:
for stage in stages_to_expand:
    for module in stage.modules:
        ...
        # Prune to first expansion when in module-filter mode
        if module_filter:
            node_combinations = node_combinations[:1]
        ...
```

### 3.4 Path Consistency

Because the first-expansion heuristic picks the same `(input_node, params)` as
full-run would, and because path construction is deterministic, the output paths
for the selected node are **byte-for-byte identical** to what full-run would
produce. This means:

- Outputs written during module-run can be reused by a subsequent full-run
  without recomputation (Snakemake's up-to-date check will skip them).
- Log files, metadata, and benchmark output are stored in standard locations.

## 4. Alternatives Considered

### Alternative 1: Snakemake target-rule execution

Pass `-t <rule_name>` to Snakemake instead of pruning the DAG at the compiler
level.

- **Pros**: No compiler changes needed; reuses Snakemake's own dependency
  resolution.
- **Cons**: The compiler still generates the full DAG; Snakemake then needs to
  read and parse all rules. For large benchmarks this is slow. Also, rule names
  in the generated Snakefile are hashed node IDs — not user-friendly.
- **Reason for rejection**: Compiler-level pruning is faster, more transparent,
  and gives a smaller Snakefile for better debuggability.

### Alternative 2: User-selected upstream combination

Allow the user to specify which upstream module to use as input (e.g.
`-m M3 --input-dataset D2`).

- **Pros**: More control.
- **Cons**: Introduces new CLI surface area and requires the user to understand
  the wildcard structure. The first-expansion heuristic is simpler and covers
  the dominant use case (smoke testing).
- **Reason for rejection**: Out of scope for v1; see Future Work.

### Alternative 3: Stub/mock upstream outputs

Instead of running upstream stages, synthesise placeholder inputs.

- **Pros**: Fastest possible feedback; no upstream execution at all.
- **Cons**: Module may behave differently with synthetic vs real inputs; harder
  to implement correctly.
- **Reason for rejection**: Not representative enough to be useful.

## 5. Implementation Plan

1. **Phase 1** (this PR): Add `-m / --module` CLI option; implement stage-ceiling
   and first-expansion pruning in `_generate_explicit_snakefile()`; add error
   path for gather stages.
2. **Phase 2** (future): Tests — unit tests for pruning logic, integration test
   with a minimal benchmark.

### Testing Strategy

- Unit test: given a benchmark with 3 stages (data × 2, methods × 2,
  metrics × 2), `run -m M1` produces exactly 3 nodes (data_D1, methods_M1, and
  the single combined node).
- Integration test: generate Snakefile for module run, verify output paths match
  a full run's first combination.

## 6. Future Work

### 6.1 Blessed test dataset

The first-expansion heuristic is a pragmatic default. A richer future feature
would allow benchmark maintainers to designate one upstream combination as the
**blessed dataset** — a curated, small, representative input that is explicitly
intended for module development and CI smoke tests.

```yaml
# Hypothetical future syntax
stages:
  - id: data
    modules:
      - id: D1
        test: true   # this is the canonical dev/test dataset
      - id: D2
```

With a blessed dataset, `ob run -m` would automatically select the blessed
upstream node rather than document-order-first. This decouples the "what to
smoke-test against" decision from "which happens to appear first in the YAML",
making the selection more intentional and maintainable.

### 6.2 Multi-module run

`ob run -m M1 -m M3` — run multiple modules in a single pass, building their
shared upstream sub-graph once.

### 6.3 Gather module support

To support `ob run -m <gather_module>`, the compiler would need to run all
provider stages (with their first expansions) so the gather node has inputs.
This is straightforward but was deferred to keep the initial implementation
simple.

### 6.4 Cross-stage smoke test

`ob run --stage <stage_id>` — run one expansion of *every* module in a given
stage, against the first upstream combination. Useful for stage-level
integration testing.

## 7. References

1. [Design 004: YAML Specification](004-yaml-specification.md)
