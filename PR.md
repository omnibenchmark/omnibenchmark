# feat: fan-in — joins and gather (#289)

Makes fan-in expressible. A benchmark can now consume outputs from several
upstream nodes at once, whether they are alternatives, branches to join, or a
set to collect.

Closes #289.

## What changes

Three ways a stage can reference outputs produced by more than one upstream
stage. Each was broken differently before:

| producers declare | consumer declares | before | after |
|---|---|---|---|
| `method_a` → `outputs: [{id: clustering}]`<br>`method_b` → `outputs: [{id: clustering}]` | `inputs: [clustering]` | **accepted and silently wrong** — one consumer node; `method_b` ran and nothing consumed its output | one consumer node **per producer** |
| `method_a` → `outputs: [{id: a.out}]`<br>`method_b` → `outputs: [{id: b.out}]`<br>on divergent branches | `inputs: [a.out, b.out]` | **rejected at load**: *"…on divergent branches… not supported yet"* | one node with **two parents** (a join) |
| any stages declaring `clustering` | `gather: [{from: clustering}]`<br>plus `prefix:` | **not expressible** | one node **per group**, members passed as a list |

`inputs:` syntax is unchanged in the first two rows. What decides the outcome is
whether the producers declare the **same** output id: same id means the files
are interchangeable, so the consumer runs once per producer; different ids mean
the consumer needs both at once, so they join into one node.

The first row is a resolver fix and is **not api-gated** — with a single
producer per id the plan is byte-identical, so the only benchmarks whose plans
change are the ones that were already losing work. The join and `gather` are
gated on **api 0.7.0**.

## Gather

```yaml
- id: metrics
  prefix: aggregated              # required: root of the cut chain
  gather:
    - from: clustering            # every producer of this id
      group_by: data              # one node per ancestor `data` module
  outputs:
    - id: metrics.summary
      path: "{data}_summary.tsv"
```

`group_by` is optional — omitted, every producer collects into one node, which
is what a `metric_collector` is. Gather nodes are ordinary nodes: downstream
stages chain off them normally.

## Other behaviour changes

- **`provides` labels are owned by one stage.** A second stage declaring the
  same label is now a parse-time error. This is what makes a join's label union
  well-defined, and it removes an undesigned silent overwrite on plain chains.
- **A join inherits labels from every branch**, so `requires:` downstream of a
  join can name a label from any of them. `exclude` already worked this way; the
  two agree now.
- **Stage-level topology drew one edge per output id** regardless of how many
  stages produced it, so a diagram showed a consumer depending on one of its
  producers. Fixed in `get_stages_by_output`, which mermaid, dot and obeditor
  all reach through `stage_adjacency`.

## Provenance

A fan-in node's path no longer *is* its lineage — a gather's records only the
group key, a join's only its deepest input. Two additions cover that: the
directory segment carries a digest of the parent set (without it two joins
sharing a branch collide, and Snakemake raises `AmbiguousRuleException`), and
every fan-in node writes a `lineage.json` sidecar naming each contributing node.
The sidecar is written before the module runs, so a gather module can read it to
label its output.

## Also in here

- The planner moved out of `cli/run.py` into `core/_expand.py`,
  `core/_prune.py` and `core/_lineage.py`. `run.py` drops 2195 → 1313 lines and
  keeps the CLI command, the Snakemake subprocess and the resolve/report loop.
  No behaviour change.
- Design docs updated to match: 004 v6 (the YAML surface), 007 v2, 010 v3.

## Testing

Full suite green. New e2e fixture runs scatter → gather → scatter for real and
checks the aggregate against what the upstream modules computed, so collecting
the wrong members fails, not just the wrong number. Unit coverage for
membership, grouping, bundle selection and topology edges — all now marked
`short`, since most of the gather tests were unmarked and `tox -m short` was
skipping them.

36 files, +4798 / −1451, 20 commits.
