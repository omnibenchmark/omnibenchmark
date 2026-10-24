# 013: Viash Component Adaptor

[![Status: Draft](https://img.shields.io/badge/Status-Draft-yellow.svg)](https://github.com/omnibenchmark/docs/design)
[![Version: 0.1](https://img.shields.io/badge/Version-0.1-blue.svg)](https://github.com/omnibenchmark/docs/design)

**Authors**: ben
**Date**: 2026-10-24
**Status**: Draft
**Version**: 0.1
**Supersedes**: N/A
**Reviewed-by**: TBD
**Related Issues**: TBD

## Changes

| Version | Date | Description | Author |
|---------|------|-------------|--------|
| 0.1 | 2026-10-24 | Initial draft, after the `task_label_projection` spike | ben |

## 1. Problem Statement

OpenProblems (OP) publishes benchmark components as [Viash](https://viash.io)
components: a script plus a `config.vsh.yaml` with typed arguments, a Docker
image and a Nextflow runner. We want to run them inside ob plans without
forking them, and eventually to run a whole OP task as an ob plan.

The two obvious bridges are expensive: an ob→Nextflow emitter (a second
backend, a second caching model, see §4), or an nf-core module shim (nf-core
modules are tool-shaped, Groovy-interpolated, and take parameters as one
`ext.args` string). A Viash component already is an ob module in all but
argv shape, so the adaptor is a thin argv rewrite plus metadata extraction.

## 2. Design Goals

- Run an unmodified Viash component as an ob module on api 0.8.0.
- Derive module metadata (inputs, outputs, typed parameters, image, resources)
  from `config.vsh.yaml`, not by hand.
- Derive the stages of an OP task's plan from its `src/api/` contracts.
- No changes to ob core: the output is an ordinary ob module.

### Non-Goals

- Executing Nextflow, or emitting Nextflow from ob.
- Publishing ob modules back to OP or nf-core (possible, see §6; not planned).
- Re-implementing OP's results site (`output_scores` etc.); we stop at an ob
  gather stage over the standard `file_score` outputs.

## 3. Proposed Solution

### 3.1 Component → module

`viash build --runner executable -c '.engines := [{"type":"native"}]'` turns a
component into a standalone CLI with the script inlined and no container
logic. The container becomes an ob `software_environment` (the image the
component's `engines[docker]` names, or the one OP publishes on ghcr). A
generated bash entrypoint rewrites ob's argv into the component's:

| ob argv (api 0.8) | component argv |
|---|---|
| `--output <id>=<path>` | `--<output arg> <path>` (map `id → arg`; `output` auto-mapped for single-output components) |
| `--<input id> <path>` | `--<input arg> <path>` (map `id → arg`) |
| `--output_dir`, `--name` | dropped |
| `--<param> <value>` | unchanged (ob params and Viash params share the `--key value` shape) |
| `resources.cores/memory` | `---cpus`, `---memory` (TODO, not in the spike) |

The mapping `ob id → Viash arg` is the only hand-written input; everything
else comes from the merged config (`viash config view --format json`).
Inputs the component does not declare are dropped, so a stage can offer a
wider input set than one of its modules consumes (§3.3, control methods).

The generated module is a git repo with `run.sh`, `omnibenchmark.yaml`,
`viash/<name>` (the built executable) and a `source:` record (OP repo,
commit, component path, viash version) so the build is reproducible.

### 3.2 Concept mapping

| OpenProblems / Viash | ob |
|---|---|
| task repo (`_viash.yaml`) | plan |
| `info.type` (`process_dataset`, `method`, `control_method`, `metric`) | stage id |
| `src/api/comp_<type>.yaml` (args with `direction`) | stage inputs / outputs |
| `src/api/file_<x>.yaml` (h5ad slot spec) | output id; later a stage validator (002) |
| component (`config.vsh.yaml` + script) | module |
| `engines[docker].image` + `setup`, built to ghcr | `software_environment` (apptainer/docker) |
| `runners[nextflow].directives.label` + `config_mods` labels | `resources` (cores/memory/time) |
| non-file arguments (typed, defaults) | module parameters |
| `info.variants` / `--paramsets` | parameter list of same-key dicts |
| `runEach(id: id + "." + comp.name + paramset)` | lineage path |
| `fromState`/`toState` wiring | stage `inputs:` by output id |
| `filter: preferred_normalization == dataset normalization` | `provides`/`requires` labels (needs the `provides` surface, see §5), or `exclude` |
| `--methods_include/exclude` | `-m`, `exclude`, obfilter |
| `extract_uns_metadata` + `joinStates` → `score_uns.yaml` | gather stage (010) with one shared collector module |
| `test_resources`, `check_config.py`, `run_and_check_output.py` | `ob validate module` (002) |
| `resources_test/` S3 subsample | small datasets for CI |

### 3.3 Task → plan

Three OP tasks were checked (`task_label_projection`, `task_batch_integration`,
`task_denoising`); their `run_benchmark/main.nf` share one skeleton:
`extract_uns_metadata` → `runEach(methods, filter by normalization)` →
`runEach(metrics)` → extract scores → `joinStates`. Everything but the
normalization filter and the collector is declared in `src/api/`, so a plan is
derivable:

```yaml
stages:
  - id: data        # one loader module per common dataset (S3 resources)
  - id: process     # comp_process_dataset: input → train, test, solution
    inputs: [dataset]
  - id: methods     # comp_method ∪ comp_control_method
    inputs: [train, test, solution]   # control methods consume solution; others drop it
  - id: metrics     # comp_metric
    inputs: [solution, prediction]
  - id: scores      # gather over metrics.score → score_uns.yaml
    gather: ...
```

Methods and metrics are the components under `src/methods`, `src/control_methods`,
`src/metrics`; parameters come from each component's arguments and
`info.variants`. Control methods and methods share a stage because they
produce the same `file_prediction` contract; two stages with one contract
collapse (010 §3.1).

### 3.4 First target

`task_label_projection` on its `resources_test` subsample
(`cxg_immune_cell_atlas`, 2.8 MB), minus the five foundation-model methods
(geneformer, scgpt ×2, scimilarity ×2, uce, scprint: GPU and multi-GB model
files). Leaves ~10 CPU methods, 3 control methods, 2 metrics, seconds each.

## 4. Alternatives Considered

- **Nextflow emitter**: a second backend. The resolver already yields
  backend-free `ResolvedNode`s, so it is feasible, but the path model (work
  dirs + `publishDir` vs ob's absolute lineage paths), `-resume` vs
  target-based reruns and `--until` all diverge. Only worth it to run whole
  OP/nf-core *workflows* as stages.
- **nf-core module shim**: run `nextflow run` per node. JVM per node, one
  `ext.args` string instead of typed parameters, Groovy before the `script:`
  block (36% of 2117 modules have ≥5 lines of it). Catalogue fit is also
  narrow: ~100–150 of the 2117 modules are competing methods on a shared
  output; the rest is plumbing or singletons.
- **Translate Viash scripts to plain ob modules**: loses the config-derived
  metadata and the published image; nothing gained over the wrapper.
- **A `kind: viash` module in ob core**: the argv rewrite in the Snakemake
  emitter instead of a wrapper. Saves a 25-line file per module, costs a
  core concept; rejected while the wrapper suffices.

## 5. Gaps

- `preferred_normalization` is a dataset attribute × method requirement: the
  `provides`/`requires` axis. `Stage.provides` is not in main (lives on
  `feat/dynamic-filtering`); until then it is an `exclude` list per method ×
  dataset.
- Static model files (`args: [model: file("s3://...")]`): no ob primitive for
  a constant file input. Options: bundle in the module repo (too big), a
  URL parameter the wrapper fetches, or a `models` data-like stage.
- Collector: OP's `score_uns.yaml` + `method_configs.yaml` + `dataset_uns.yaml`
  are task-independent (`file_score` is standard), so one shared gather module
  covers every task.
- Image pinning: ghcr tags are `build_main`; pin by digest in the plan.
- Resources: label → numbers via `_viash.yaml` `config_mods`; not forwarded yet.

## 6. Implementation Notes

- Spike: `scratch/viash_spike/viash2ob.py` (gitignored), runs knn + accuracy
  end to end on api 0.8.0 with the ghcr images; accuracy 0.497 on the subsample.
- Proposed home: `ob create module --from-viash <config.vsh.yaml> --input
  ob_id=arg ... --output ob_id=arg ...` (CLI layer; needs only the config
  parser), then `ob create benchmark --from-viash-task <dir>` for §3.3.
- Reverse direction (ob module → Viash component) is a config whose
  arguments are the stage contract and whose script calls the entrypoint;
  cheap if OP wants our modules, not planned.
- viash 0.9.7: `viash config view` exits 1 even on success; both commands must
  run from the directory holding `_viash.yaml` so `/src/api/...` merges resolve.
