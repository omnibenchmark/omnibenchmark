# 014 — One benchmark model across ob, Nextflow/nf-core, OpenEBench and Viash

Status: draft, 2026-10-24. Branch `feat/convergence` (on top of 013). Evidence:
two pipelines lowered to ob plans with scores matching the originals
(`~/phd/omni-label-projection`, OpenProblems via Viash; `~/phd/omni-batch-integration`,
nf-core/scdownstream via a Nextflow parser).

## 1. Claim

A benchmark is a plan-level DAG whose stages have roles. Everything the four
systems disagree on is either (a) how the DAG is written down, or (b) where
execution details live. Neither needs to be shared to share benchmarks: a
module written for one system runs in another when the stage contract matches,
and a benchmark wiring written in one system is compiled, not re-authored.

## 2. The shared model

| concept | ob | Nextflow / nf-core | OpenEBench | Viash / OpenProblems |
|---|---|---|---|---|
| unit of work | module (`run.sh` + env) | process (nf-core module) | step container | component (config.vsh.yaml) |
| contract | stage inputs/outputs by id | `input:`/`output:` declarations, `emit:` | fixed file layout + parameter names | `arguments` with `type: file`, `direction` |
| data role | stage `role: data` | upstream subworkflow output | validation step | `data_processor` |
| method role | stage `role: method`, N modules | dispatcher `if (tools.contains('x')) { X(...) }` | participant | `method` / `control_method` |
| metric role | stage `role: metric` | consumer process over the mixed channel | metrics step | `metric` |
| collect role | gather / collector | consolidation subworkflow | consolidation step | aggregation script |
| alternatives | modules in a stage; `exclude` | `contains()` blocks; `.filter` | challenges list | `info.variants` |
| lineage | path + labels (`provides`/`requires`) | `meta` map | participant id + challenge id | `dataset_id` / `method_id` in output |
| parameters | `parameters:` grid | `val` inputs, `ext.args` | per-step params | `arguments` |
| environment | `software_environments` (conda/apptainer/envmodule) | `container`/`conda` directive | container per step | `engines` (docker/native) |
| resources | `resources` on stage/module | `label 'process_high'` + `conf/*.config` | per-step | `resources` |

Rows 1 to 7 converge already: each system has the concept, the compilers
(013, nf2ob) translate it mechanically. Rows 8 to 10 are where ob should move.

## 3. What ob changes

**Roles on stages** (`role: data | method | metric | collect | report`). Today the
roles are implicit in position. Naming them lets any emitter (nf-core scaffold,
OpenEBench three-step, a report) read the plan without heuristics, and lets ob
validate shape ("a metric stage consumes a method stage").

**Execution profiles.** `Resources` exists on stage and module and is emitted
as Snakemake `threads`/`resources`. Move it, with `software_backend`, executor
and container arguments, into a companion file: `ob run plan.yaml --profile
omni.yaml`. A plan then runs unchanged on a laptop, omni and SLURM, and a
Nextflow `conf/*.config` lowers to a profile, which is what it is.

**Optional and multiple outputs.** `optional: true` (absent file is not a
failure; downstream sees no path) and `multiple: true` (a glob, downstream sees
a list). Nextflow has both (`[]`, `path(reads)`); nf-core modules assume them;
ob modules fake them with placeholder files today. Aligns with 012's typed
ports.

**Not changed: lineage.** `Stage.provides` + `Module.provides: {label: value}` +
`Module.requires` is `meta` declared in the plan, and `{label}` substitution
makes it readable downstream. The gap is run-time `meta` (values produced by
the data); that is the dynamic-filtering proposal, separate.

## 4. What the others would gain

- **nf-core**: a benchmark is one subworkflow shape (dispatcher + consumer).
  Stating the shape (roles, one tag key for the alternative, one shared channel)
  makes a pipeline benchmark-compatible and compilable to ob; ob modules compile
  to nf-core modules by the reverse mapping (row 2), which is the module-export
  path the group needs.
- **OpenEBench**: validation / metrics / consolidation is data / metric /
  collect with the participant as the method stage. An OEB-compatible workflow
  is an emitter target from the same plan; an ob metric module is an OEB metrics
  container with a fixed argv.
- **Viash/OpenProblems**: already mechanical (013).

## 5. Open

- Subworkflows as methods (sarek): recursive lowering vs subworkflow-as-module.
  First stop reason from the compiler bench, 2026-10-24.
- Run-time `meta` values (dynamic filtering).
- Where parameters live when a module is shared: ob grid vs nf `ext.args` string.
