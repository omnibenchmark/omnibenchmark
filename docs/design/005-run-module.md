# 005: Module-Scoped Run (`ob run -m`)

[![Status: Draft](https://img.shields.io/badge/Status-Draft-yellow.svg)](https://github.com/omnibenchmark/docs/design)
[![Version: 0.2](https://img.shields.io/badge/Version-0.2-blue.svg)](https://github.com/omnibenchmark/docs/design)

**Authors**: ben
**Date**: 2026-02-18
**Status**: Draft
**Version**: 2
**Supersedes**: N/A
**Reviewed-by**: daniel
**Related Issues**: TBD

## Changes

| Version | Date | Description | Author |
|---------|------|-------------|--------|
| 1       | 2026-02-19 | Initial draft | ben |
| 2       | 2026-03-31 | Add pluggable `--strategy` option; promote first-expansion to `first`; add `full-upstream` and `blessed` strategies | ben |

## 1. Problem Statement

`ob run benchmark.yaml` always runs the full benchmark DAG: every stage, every
module, parameter expansion, and input combination. This is appropriate
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
- **Extensible expansion strategy**: The choice of "which upstream combination
  to use" is a separate, swappable concern; adding new strategies should not
  require changes to the stage-ceiling or path-construction logic.

### Non-Goals

- Automatically download missing upstream data from a remote source.
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
combined with `--dry`, `--dirty`, `--dev`, `--cores`, etc.).

A `--strategy` option (see §3.2) is planned for Phase 2. In the current
implementation the expansion strategy is fixed to `first` and cannot be
changed at the command line.

### 3.2 DAG Pruning

When `-m <module_id>` is given, the compiler applies two steps before
Snakefile generation: a fixed **stage ceiling** (step 1, always applied) and a
**strategy-dependent expansion filter** (step 2, controlled by `--strategy`).

#### Step 1 — Stage ceiling (invariant)

Find the stage `S` that contains module `<module_id>`.

Include **only** stages whose index in the benchmark YAML is ≤ index(`S`).

Stages after `S` are dropped entirely. This means downstream consumers of `S`'s
outputs are not built.

**Rationale**: The developer wants to verify their module, not its consumers.

**Limitation**: If `S` is a gather stage (`stage.is_gather_stage() == True`)
module-run mode is not supported and the command should exit with an error. Gather
nodes require *all* provider outputs to be present (they are the `rule all`
target), so running a single gather module in isolation is not meaningful
(see §6.3).

#### Step 2 — Expansion strategy

The expansion strategy determines which input-node × parameter combinations are
kept for each module in the pruned sub-graph. Three strategies are defined:

---

##### `--strategy=first` (default)

Each module in every included stage keeps only the **first** input-node ×
parameter combination (index 0 in document order).

Concretely, for any module `m` in stage `S'` (where `S'` ≤ `S`):

- Take `input_node_combinations[0]` (the first upstream node, in document order).
- Take `params_list[0]` (the first expanded parameter combination, in document order).
- Keep exactly that one `(input_node, params)` pair; discard the rest.

**Rationale**: Deterministic, requires no configuration, produces a minimal but
runnable sub-graph. Paths are byte-for-byte identical to what a full run would
produce for the same combination, so outputs are reusable (see §3.4).

**Use case**: Fast local smoke test during module development.

---

##### `--strategy=full-upstream`

All included stages (upstream and target) are expanded fully — all input-node ×
parameter combinations for all modules are kept. The only difference from a
full `ob run` is the stage ceiling: stages after `S` are still dropped.

In practice this means: run all upstream work exactly as a full run would, then
run the target module against every upstream output × every parameter combination.

**Rationale**: Useful when the developer wants to verify their module works
correctly across the full range of upstream inputs, not just a representative
sample. This is the right choice before a PR, or when correctness (not just
smoke-testing) is the goal.

**Use case**: Pre-submission validation; CI pipelines for a single module.

---

##### `--strategy=blessed`

Upstream modules that carry a `test: true` annotation are preferred over
document-order-first. Each included stage selects the **blessed** module (if
one exists) and applies the first-expansion rule to it; if no module is blessed
in a stage the strategy falls back to `first` for that stage.

```yaml
# Example benchmark YAML with blessed annotation
stages:
  - id: data
    modules:
      - id: D1
        test: true   # canonical dev/test dataset
      - id: D2
  - id: methods
    modules:
      - id: M1
        test: true   # reference method for metric development
      - id: M2
```

This strategy requires YAML-side coordination by the benchmark maintainer.
The `test: true` flag is a new optional key on the module object; it carries no
semantics outside of `ob run -m --strategy=blessed`.

**Rationale**: Decouples "which input to smoke-test against" from "which
happens to appear first in the YAML", making the selection intentional and
maintainable across benchmark evolution.

**Use case**: Project-level convention for developer onboarding and CI; the
maintainer curates a small, representative input once and all module developers
benefit automatically.

**Status**: `blessed` requires the new `test` YAML key and is deferred to a
future release. The CLI should reject `--strategy=blessed` with a clear error
until it is implemented.

### 3.3 Implementation Sketch

The strategy is represented as a callable (or small protocol) that receives the
full list of `(input_node, params)` combinations for a given module and returns
the filtered subset.

```python
from typing import Protocol, Sequence

class ExpansionStrategy(Protocol):
    def filter_combinations(
        self,
        combinations: list[tuple],   # all (input_node, params) pairs
        module,                       # the module being expanded
        stage,                        # the stage it belongs to
    ) -> list[tuple]: ...


class FirstStrategy:
    """Keep only the first combination (document order)."""
    def filter_combinations(self, combinations, module, stage):
        return combinations[:1]


class FullUpstreamStrategy:
    """Keep all combinations (no pruning)."""
    def filter_combinations(self, combinations, module, stage):
        return combinations


# In _generate_explicit_snakefile():

if module_filter:
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

    stages_to_expand = benchmark.model.stages[:target_stage_idx + 1]

    strategy_map = {
        "first": FirstStrategy(),
        "full-upstream": FullUpstreamStrategy(),
        # "blessed": BlessedStrategy(),  # future
    }
    strategy = strategy_map[expansion_strategy]  # from --strategy CLI arg
else:
    stages_to_expand = benchmark.model.stages
    strategy = None

# During node expansion, for each stage:
for stage in stages_to_expand:
    for module in stage.modules:
        ...
        if strategy is not None:
            node_combinations = strategy.filter_combinations(
                node_combinations, module, stage
            )
        ...
```

The `BlessedStrategy` class is a future addition; it would scan `module.test`
flags and fall back to `FirstStrategy` when none are set.

### 3.4 Path Consistency

Because every strategy selects only combinations that full-run would also
generate, and because path construction is deterministic, the output paths for
any selected node are **byte-for-byte identical** to what full-run would
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
  and gives a smaller Snakefile for better debugging.

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
- **Reason for rejection**: Not representative enough to be useful. Additionally, this can be better supported when we land the validation feature (valid data ranges and types).

## 5. Implementation Plan

1. **Phase 1** (this PR): Add `-m / --module` CLI option; implement
   stage-ceiling logic and the `first` expansion (hardcoded) in
   `_generate_explicit_snakefile()`; add error path for gather stages.
2. **Phase 2** (future): Add `--strategy` CLI option; implement
   `full-upstream` strategy; reject `--strategy=blessed` with a clear error.
   Unit tests for both strategies, integration test with a minimal benchmark.
3. **Phase 3** (future): `blessed` strategy — add `test: true` YAML key,
   implement `BlessedStrategy`, update YAML schema.

### Testing Strategy

- Unit test (`first`): given a benchmark with 3 stages (data × 2, methods × 2,
  metrics × 2), `run -m M1 --strategy=first` produces exactly 2 rules
  (data_D1, methods_M1).
- Unit test (`full-upstream`): same benchmark, `run -m M1 --strategy=full-upstream`
  produces rules for all data modules plus methods_M1 × all data outputs × all parameter combinations.
- Integration test: generate Snakefile for module run, verify output paths match
  the corresponding paths from a full run.

## 6. Future Work

### 6.1 Multi-module run

`ob run -m M1 -m M3` — run multiple modules in a single pass, building their
shared upstream sub-graph once.

### 6.2 Gather module support

To support `ob run -m <gather_module>`, the compiler would need to run all
provider stages (with their first expansions) so the gather node has inputs.
This is straightforward but is deferred to keep the initial implementation
simple.

### 6.3 Cross-stage smoke test

`ob run --stage <stage_id>` — run one expansion of *every* module in a given
stage, against the first upstream combination. Useful for stage-level
integration testing.

## 7. References

1. [Design 004: YAML Specification](004-yaml-specification.md)
