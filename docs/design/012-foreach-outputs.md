# 012: `foreach` outputs (runtime cardinality)

[![Status: Draft](https://img.shields.io/badge/Status-Draft-yellow.svg)](https://github.com/omnibenchmark/docs/design)
[![Version: 1](https://img.shields.io/badge/Version-1-blue.svg)](https://github.com/omnibenchmark/docs/design)

**Authors**: ben
**Date**: 2026-10-23
**Status**: Draft
**Version**: 1
**Supersedes**: N/A
**Reviewed-by**: TBD
**Related Issues**: 009 §3.1 (`kind: zip`), 009 §4 Alt. 1 (`directory()`), 010 (gather), 011 (wildcard lowering, parked)

## Changes

| Version | Date | Description | Author |
|---------|------|-------------|--------|
| 1 | 2026-10-23 | Initial draft, implemented on api 0.8.0 | ben |

## 1. Problem Statement

A stage output is one file per node. A stage whose file count is only known
after it runs (download every sample of a 10x study, split a library by
barcode) cannot be expressed: the plan enumerates every node before Snakemake
starts.

When the count is known up front, parameter values cover it: one node per
value, and `{a..b}` ranges keep long lists short. This design covers the case
where it is not.

## 2. Design Goals

- One stage output yields N files, N chosen by the module at run time.
- Each file feeds its own downstream node, so Snakemake runs them in parallel.
- A gather collects them, unchanged in syntax.
- Everything else stays explicit: every node exists at plan time; only one
  path segment is left to Snakemake.

### Non-Goals

- Nested `foreach` (one below another).
- Gating (`exclude`, `requires`, obfilter) on captured values.
- A gather grouped by a captured value.
- The general wildcard backend (011).

## 3. Model

Spec: [004 §3.14](004-yaml-specification.md).

```yaml
- id: data
  outputs:
    - id: counts
      path: "{sample}/counts.h5"
      foreach: sample
- id: process
  inputs: [counts]
  outputs: [{id: processed, path: "{module.id}.h5ad"}]
- id: per_dataset
  gather: [{from: processed, group_by: data}]
```

Plan-time errors:

1. `foreach` is an identifier other than `module` and `params`.
2. `path` starts with `{<foreach>}/` and uses it nowhere else. The first
   segment sits directly under the node directory, and a downstream node's
   directory extends its input's, so every descendant path carries the value.
3. At most one `foreach` output per stage; none on metric collectors.
4. No `foreach` on a stage whose nodes already carry one.
5. Every wildcard a node's inputs carry must reach its outputs. This rejects
   `nesting_strategy: flat` and a join whose deepest input is outside the
   `foreach` branch.
6. api_version ≥ 0.8.0.

### 3.1 Module contract

The producer gets the pattern with the name left in:

```
--output counts=/abs/out/data/D1/.a1b2/{sample}/counts.h5
```

It writes one file per value, substituting `{sample}` itself. Values match
`[^/]+`. Zero files is valid: the downstream nodes then have no instances.
Downstream modules see an ordinary single input path.

### 3.2 What a captured value is

Not a lineage label. It names a directory and nothing else: it is not in
`TemplateContext.provides`, not visible to `requires`, and a gather does not
bind it.

## 4. Lowering

Invariant: after `TemplateContext.substitute`, a resolved path holds no `{`
except a kept `foreach` name. A `{` in a resolved path therefore means "below a
`foreach`", and the generator derives everything from that.

```python
wildcard_constraints:
    sample="[^/]+",

def _ob_foreach(ckpt, pattern, name, template):
    ckpt.get()
    found = glob_wildcards(pattern.replace("{%s}" % name, "{%s,[^/]+}" % name))
    return expand(template, **{name: sorted(getattr(found, name))})

checkpoint data_D1_a1b2:
    output:
        "data/D1/.a1b2/.counts.foreach",
    shell:
        """
        FOREACH_DIR=$(cd $(dirname {output[0]}) && pwd)
        rm -f $FOREACH_DIR/*/counts.h5
        OUTPUT_FILE_0=$FOREACH_DIR/{{sample}}/counts.h5
        ...
        --output counts=$OUTPUT_FILE_0
        touch $FOREACH_DIR/$(basename {output[0]})
        """

rule data_D1_a1b2_process_M_default:
    input:
        counts="data/D1/.a1b2/{sample}/counts.h5",
    output:
        "data/D1/.a1b2/{sample}/process/M/.default/M.h5ad",
    log:
        ".logs/data_D1_a1b2_process_M_default.{sample}.log"

rule per_dataset_G_D1_default:
    input:
        input_0=lambda wc: _ob_foreach(checkpoints.data_D1_a1b2,
            "data/D1/.a1b2/{sample}/counts.h5", "sample",
            "data/D1/.a1b2/{sample}/process/M/.default/M.h5ad"),
```

- **Producer**: written as `checkpoint`. Its rule output for the `foreach` id
  is a sentinel `.<id>.foreach` in the node directory, touched after the module
  succeeds; the checkpoint itself has no wildcard. Files matching the pattern
  are removed before the module runs, so a re-run that yields fewer values
  leaves nothing stale to glob. Downstream outputs of a dropped value stay on
  disk, unreferenced. No `directory()`: downstream outputs live under the same
  directory, and Snakemake rejects outputs inside another rule's directory
  output.
- **Wildcard node**: any node with `{` in a path. Unchanged except the log
  name, which gains `.{sample}` (output, log and benchmark must share
  wildcards).
- **Consumer**: any rule that references a path carrying a wildcard its own
  outputs lack — a gather, a metric collector, or `rule all` — gets the path
  through `_ob_foreach`. The producer is the one whose node directory plus
  `/{name}/` starts the path; node directories are unique, so the match is.
- `checkpoint` is a Snakefile keyword, hence `ckpt`. Glob order is unstable
  and becomes a gather's CLI order, hence `sorted`.

Because `rule all` is wrapped the same way, a wildcard leaf with no gather
still runs, and `--until` a stage below the `foreach` works.

## 5. Costs

| Area | Behaviour |
|------|-----------|
| `ob describe`, mermaid, dot | One node per (module, params); paths show `{sample}` |
| status, archive, telemetry node counts | See the pattern, not instances |
| `snakemake -n` | Lists jobs below the checkpoint only after it has run |
| `exclude` / `requires` / obfilter | Work on modules and params; cannot select values |
| gather `group_by` | Any stage; values always merge into the group, as parameter expansions do (010 §3.2) |

## 6. Alternatives Considered

1. **Plan-time values (`{a..b}` ranges)**: known cardinality only.
2. **`kind: zip` (009 §3.1)**: one tracked file, no DAG change; the next
   module loops over samples itself, so no per-sample parallelism.
3. **Two-phase run (manifest read at plan time)**: the plan depends on an
   earlier run's output; needs staged execution.
4. **`directory()` + checkpoint (009 §4 Alt. 1)**: rejected there for DAG
   complexity. This design is that alternative without `directory()`, confined
   to one segment.
5. **General wildcard backend (011)**: subsumes this, far larger.
6. **Manifest written by the module** instead of ob globbing: gives status and
   provenance a record, at the cost of a module contract. Deferred; the plan
   syntax does not change if it is added.

## 7. Implementation Notes

- `IOFile.foreach` and rules 1–3, 6: `model/benchmark.py`. The `/` deprecation
  warning on output paths skips the leading `{name}/`.
- `TemplateContext.substitute(keep=…)` and `ResolvedNode.foreach = (name,
  output id)` on producers: `model/resolved.py`. Descendants carry no flag.
- Rules 4–5: `core/_expand.py`, scatter expansion, via `path_wildcards`
  (`core/_paths.py`).
- Lowering: `backend/snakemake.py` (`_rule_outputs`, `_path_ref`,
  `_write_foreach_header`, producer shell).
- Tests: `tests/model/test_foreach.py`, `test_foreach_lowering` in
  `tests/backend/test_snakemake_gen.py`, `tests/e2e/test_15_foreach.py` (fan
  out + gather; a re-run dropping a value; nested rejection). The e2e module is
  a local git repo built in `tmp_path`, not a bundle.

Deferred: status / archive expansion of patterns, manifest, gather grouping by
captured value, nested `foreach`.
