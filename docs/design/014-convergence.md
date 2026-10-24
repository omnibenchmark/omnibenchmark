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

**Vectorized modules** (`Module.vectorized: true`, prototyped on this branch). One
invocation for all nodes of the module: the entrypoint gets `--nodes PATH`, a
JSON list of `{id, inputs, parameters, outputs}` records, and writes every
listed output. The Snakemake rule depends on the union of the nodes' inputs and
produces the union of their outputs, so lineage, `exclude` and the layout are
unchanged and downstream nodes are oblivious. This is how a Nextflow block runs
once with the fan-out inside (channels are the fan-out), how a method loads one
model for every dataset, and what a Dagster emitter maps to a partitioned asset
(partition key = lineage id).

**One rule per module** (default on this branch; `ob run --granularity node` keeps
the old layout for one release). The Snakefile is the plan-level DAG: a rule per
(stage, module) that depends on the union of its nodes' inputs and produces the
union of their outputs. Inside, ob's own driver runs each node as a bash function
with today's per-node CLI, in parallel up to the rule's threads, and skips a node
whose outputs exist and are newer than its inputs (outputs are flagged `update()`
so Snakemake does not delete them first; adding a dataset reruns only its nodes).
A module with `vectorized: true` is handed the records (`--nodes`) instead of the
driver. Gather and collector nodes keep their own rules. Lost on purpose:
Snakemake's per-node `benchmark:` file (performance.txt) and cross-module
pipelining; gained: a readable Snakefile, module-level resources and retries,
and a one-to-one mapping onto SLURM array jobs, nf-core subworkflows and Dagster
partitioned assets.

**A module kind that is a Nextflow workflow or process** (`kind: nextflow`). The
module names a pipeline, a workflow (or a guarded block of one) and the element
shape of each channel it reads; one generic runner builds the entry `main.nf`
from argv and publishes the emitted files. This is how sarek's variant callers
become methods without lowering their subworkflows, and how an nf-core module
runs in ob unchanged: no partial evaluation of `script:` blocks. The compiled
native module (template rendering) stays as an export path.

**The backend is a property of the environment** (on this branch, with #398's
`{type, uri}` environments cherry-picked). `software_backend` is the plan's
preference; an environment that does not declare it runs under one it does
declare, a module without an environment runs on the host, and `ob run` passes
one `--use-*` flag per backend actually in use. This is what lets one methods
stage hold a native conda module, an Apptainer image from OpenProblems and a
Nextflow block (host, bringing its own containers) side by side.

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

## 5. Emitters

Roles make the plan-level DAG readable without heuristics, so backends and
exports are the same kind of thing: Snakemake (today), an nf-core subworkflow
plus modules (the group's export need), an OpenEBench three-step workflow, and
Dagster (one asset per node, lineage as asset keys). None of them should read
ResolvedNodes; they read stages, modules and the lineage, which is why the
emitter lives at the plan level (worklog: plan-level DAG and stage roles).

## 6. Open

- Subworkflows as methods (sarek): recursive lowering vs subworkflow-as-module.
  First stop reason from the compiler bench, 2026-10-24.
- Run-time `meta` values (dynamic filtering).
- Where parameters live when a module is shared: ob grid vs nf `ext.args` string.
