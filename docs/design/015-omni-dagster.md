# 015 — omni_dagster: the plan model with a Python face

Status: sketch, 2026-10-25. Depends on 014 (roles, vectorized modules, one node per module).

## 1. Idea

A small Python library that wraps every obcommons module as an asset factory.
The registry already carries what a factory needs: stage, inputs by output id,
outputs, environment, parameters.

```python
from omni_dagster import data, method, metric, report, nextflow, plan

ds  = data("omni-sc-datasets", ids=["pbmc3k", "lung_atlas"])            # partitions = datasets
m1  = method("scanpy-combat", on=ds.h5ad, params={"batch_col": ["batch"]})
m2  = method("op-bbknn", on=ds.h5ad)                                    # OP component, via 013
m3  = nextflow("nf-core/scdownstream@825cf27", block="INTEGRATE/scvi", on=ds.h5ad)  # via nf2ob
sc  = metric("scib-metrics", on=[m1, m2, m3], exclude=["op-bbknn"])      # exclude = partition subset
rep = report("omni-quarto-report", gather=sc, template="report.Rmd")

defs = plan("batch_integration", [ds, m1, m2, m3, sc, rep])             # Definitions + plan.yaml
```

Each factory returns a partitioned `@asset`: the op body is the module's
entrypoint with the ob CLI contract, the partitions definition is the lineage
axis the stage adds (upstream partitions × its own parameter grid), the
`group_name` is the role. `plan()` builds the Dagster `Definitions` and emits
`plan.yaml` through the reverse mapping; the factories *are* the plan model, so
the compiler is a projection and nothing has to be recovered.

## 2. Why roles make it clean

`role` is the only typing a factory carries, and it lets the builder enforce the
stage contract the DSL cannot yet: a method takes a data or method output, a
metric takes a method output, a report gathers metrics, data takes nothing.
That is the "matching types" the node-editor question (012) wanted, as Python
type hints over roles instead of UI code. It becomes a slot-level check the
moment output contracts land, because the factory reads the module's declared
outputs. Dagster's UI then renders the plan-level DAG with per-partition status,
so the dashboard comes along.

## 3. Mapping

| ob | Dagster |
|---|---|
| stage × module | one asset; `group_name` = role |
| resolved node (lineage) | one partition; key = lineage path |
| `parameters` grid | a partition dimension (multi-partitions: upstream × params) |
| `exclude` / `requires` | partition subset |
| gather / report | an unpartitioned asset depending on all partitions (`AllPartitionMapping`) |
| vectorized module | `BackfillPolicy.single_run()`: one op run for all partitions, `--nodes` records |
| per-node module | default backfill policy: one op run per partition |
| software environment | op resource (conda, apptainer, host); Nextflow blocks bring their own containers |
| output layout | `IOManager` writing the ob paths, so `plan.yaml` and Dagster runs are interchangeable |

Nextflow comes in through 014's module kind: `nextflow(pipeline, block)` lowers
the block with nf2ob (block source, channel shapes, vals) and returns the same
kind of asset as `method()`, backed by the generic runner. nf-core modules
proper are the degenerate block (one process call). So the three sources of
modules (obcommons, OpenProblems via Viash, nf-core via Nextflow) share one
factory shape and one plan.

## 4. What you get for free, and what you don't

Free: the DAG and per-partition status UI, backfills and retries per partition,
asset lineage and run history, sensors and schedules, materialising a single
method for a single dataset from the UI.

Not free: the plan-level semantics stay ob's. Lineage paths, `exclude`,
`provides`/`requires` and the output layout are computed by ob's resolver, not
re-implemented in Dagster; the factories call it. Reproducibility stays with
pinned module commits and the plan; Dagster's run storage is a record, not the
source of truth. And a Dagster deployment (daemon, run storage, workers) is a
service to operate, which the Snakemake backend is not; this is an alternative
backend for groups that want the UI, not a replacement.

## 5. First cut

`omni_dagster/` as a separate package: four factories plus `nextflow`, `plan()`
emitting `plan.yaml` only, Definitions with the default per-partition op calling
`run.sh`. Validate on the batch-integration plan: the emitted `plan.yaml` must
equal the one nf2ob produces, and a Dagster materialisation of one partition
must write the same file ob writes.
