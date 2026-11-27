# 015: Host Files — Per-Machine Budgets and Per-Stage Allocations

[![Status: Draft](https://img.shields.io/badge/Status-Draft-yellow.svg)](https://github.com/omnibenchmark/docs/design)
[![Version: 1](https://img.shields.io/badge/Version-1-blue.svg)](https://github.com/omnibenchmark/docs/design)

| | |
|---|---|
| **Authors** | btraven00 |
| **Date** | 2026-11-10 |
| **Status** | Draft |
| **Version** | 1 |
| **Supersedes** | N/A |
| **Reviewed-by** | TBD |
| **Related Issues** | [#418](https://github.com/omnibenchmark/omnibenchmark/pull/418) (default memory budget) |
| **Related designs** | [004](004-yaml-specification.md) §7 — `resources:` in the plan; [008](008-filtering.md) — `provides` labels and host capabilities |

## Changes

| Version | Date | Description | Author |
|---------|------|-------------|--------|
| 1 | 2026-11-10 | Initial draft | btraven00 |

## 1. Problem Statement

A plan's `resources:` block (004 §7) mixes two different things:

- **What a job needs.** A property of the method and the data: this stage peaks at 12 GB on a 28k-cell dataset.
- **What a machine can give.** A property of the host: 32 cores and 250 GB here, 8 cores and 16 GB on a laptop, a SLURM partition there.

Today both live in the plan, as one fixed number per stage. That has four costs.

1. **A plan is tuned for one machine.** `cores: 16` on a stage is right on a 192-core server and wrong on a laptop. People edit the plan per host, which changes a file that is supposed to describe the benchmark.
2. **One number per stage cannot follow the data.** Memory for an all-pairs metric grows with the square of the cell count. A stage that needs 12 GB on one dataset needs 150 GB on another, and the plan can only state one of them.
3. **Tuning changes the benchmark's identity.** `summary_hash` covers the `stages` section, `resources:` included. Raising a stage's `mem_mb` for a bigger machine yields a different hash for the same benchmark.
4. **The machine's total budget has no home.** It is passed by hand (`--cores`, `-- --resources mem_mb=…`) and forgotten. #418 adds a default memory budget, which closes the worst case (a declared `mem_mb` being ignored), but the number is still a guess at 90% of RAM.

The incident behind this document: a stage whose jobs each peak at 12 GB was started 32 times at once on a 251 GB machine. The machine was unreachable for about 70 minutes.

## 2. Design Goals

- **One file per machine** states that machine's budget and how it wants each stage sized.
- **The plan stays portable.** It may still carry `resources:` as a default, but nothing host-specific has to be written into it.
- **Allocations can follow the data**, through lineage labels (`provides`, 008) such as a dataset size tier.
- **Everything translates into Snakemake** directives and flags that exist today. ob adds no scheduler of its own.
- **Results do not move.** A host file never changes a node id, an output path or `summary_hash`.

### Non-Goals

- **Enforcing a job's limit.** A budget limits how many jobs start. Killing a job that exceeds its figure (cgroups, `ulimit`, SLURM) is a separate design; §6 says how it would attach.
- **Measuring needs automatically.** The numbers are written by a person, from `performance.txt` or a profiler.
- **Replacing Snakemake profiles.** A host file can name one (§3.6); it does not reimplement what profiles do.

## 3. Proposed Solution

### 3.1 The host file

A YAML file, one per machine:

```yaml
# hosts/omni.yaml
cores: 32                 # total budget -> snakemake --cores
mem_mb: 200000            # total budget -> --resources mem_mb=
gpu: 1                    # total budget -> --resources nvidia_gpu=
capabilities: [gpu]       # same as --with-capability gpu (008)

stages:
  PCA:      {cores: 4}
  PCA_M:    {cores: 16}
  EMBED_M:
    mem_mb: 16000                 # every EMBED_M node
    tier=l: {mem_mb: 64000}       # nodes whose lineage carries tier: l
  CLUST.cl_rapids: {gpu: 1}       # one module of a stage
```

Top-level keys are the machine's budget. `stages:` maps a stage id, or `STAGE.module_id`, to a resources block with the same fields as 004 §7.2 (`cores`, `mem_mb`, `disk_mb`, `runtime`, `gpu`).

Inside a stage entry, a key of the form `<label>=<value>` holds a resources block that applies only to nodes whose lineage carries that label value (§3.3).

### 3.2 Choosing the file

```
ob run plan.yaml --host hosts/omni.yaml     # a path
ob run plan.yaml --host omni                # a name
ob run plan.yaml                            # by hostname, if a file exists
```

A name is looked up as `<name>.yaml` in, in order:

1. `hosts/` next to the plan — allocations that belong to this benchmark.
2. `hosts/` in the omnibenchmark config directory — a machine's own defaults.

Without `--host`, ob looks for a file named after the short hostname in the same two places. If none exists, nothing changes: the run behaves as it does today, including the default memory budget of #418.

`--no-host` skips the lookup.

### 3.3 Tiers through labels

A stage advertises a label and each module gives it a value (008):

```yaml
- id: DATA
  provides: [tier]
  modules:
    - id: d_small
      provides: {tier: s}
    - id: d_atlas
      provides: {tier: l}
```

Every downstream node inherits the labels of its lineage, so when ob writes a rule for an `EMBED_M` node it knows that node's `tier`. The host file's `tier=l:` block then applies to exactly the nodes that descend from `d_atlas`.

The plan carries the tag and no numbers. The same tag can mean 64 GB on one host and "do not run here" on another (§3.5).

Any label works as a selector, not only one called `tier`.

### 3.4 Resolution order

For one node, each resource field is resolved separately; the first source that sets it wins:

1. Host file, `STAGE.module_id`, matching label selector
2. Host file, `STAGE.module_id`
3. Host file, `STAGE`, matching label selector
4. Host file, `STAGE`
5. Plan, module-level `resources:`
6. Plan, stage-level `resources:`
7. Default (`cores=2`, no memory figure)

If two label selectors in one entry match the same node and set the same field, the plan fails to load: there is no sensible order between labels.

For the budget: a `--cores` on the command line wins over the host file's `cores`; a `mem_mb` passed in `-- --resources` wins over the host file's `mem_mb`; the host file wins over the #418 default.

### 3.5 What the host file may not do

- It cannot add or remove nodes, with one exception: `capabilities:` takes part in the existing capability pruning (008), exactly as `--with-capability` does.
- A node whose resolved `mem_mb` or `cores` exceeds the host's budget is an error at load time, naming the node and both numbers. Today Snakemake silently clamps such a rule to the budget and runs it alone; for a job that needs 150 GB on a 64 GB laptop that is the wrong outcome. `--allow-oversized` restores the clamp.

### 3.6 Translation into Snakemake

Nothing below is new Snakemake behaviour.

| Host file | Becomes |
|---|---|
| `cores: N` | `--cores N` |
| `mem_mb: N` | `--resources mem_mb=N` |
| `gpu: N` | `--resources nvidia_gpu=N` |
| `capabilities: [...]` | ob's capability pruning; not passed to Snakemake |
| stage `cores: N` | the rule's `threads: N` and `resources: cores=N` |
| stage `mem_mb`, `disk_mb`, `runtime` | the rule's `resources:` entries of the same name |
| stage `gpu: N` | the rule's `resources: nvidia_gpu=N` |
| `snakemake_args: [...]` | appended to the Snakemake command line |

The per-stage values are applied when ob expands the plan: the resolved figure replaces the node's `resources` before the Snakefile is written. Per-rule overrides therefore do not go through `--set-resources` or `--set-threads`.

That choice is forced. Snakemake matches those overrides by exact rule name, and ob writes one rule per node: a 52-node stage is 52 differently named rules. A stage-level statement has to be resolved by ob.

`snakemake_args` is the escape hatch for what ob does not model: an executor, a profile, a default partition.

```yaml
# hosts/roland.yaml
capabilities: [gpu]
snakemake_args: ["--executor", "slurm", "--jobs", "40", "--default-resources", "slurm_partition=gpu"]
stages:
  EMBED_M: {mem_mb: 16000, runtime: 30}
```

With an executor, the top-level `cores` and `mem_mb` are not passed: a global `mem_mb` would cap the total submitted to the cluster, which is the scheduler's job.

### 3.7 Identity and provenance

- A host file is applied to the expanded nodes, after the plan is loaded. It therefore cannot change a node id, an output path or `summary_hash`, and Snakemake does not rerun a job because its resources changed. Switching host files leaves an existing output tree valid.
- The plan's own `resources:` blocks are a different matter: `summary_hash` hashes the whole `stages` section, so today editing a stage's `mem_mb` changes the benchmark's identity although no result can differ. This design proposes dropping `resources` from the hash, as was done for `validators` (004 §3.14).
- The run manifest records the host file's path and content hash, and the resolved Snakemake command, as it does today for the command.
- `ob run --dry` prints, per stage, the resolved resources and how many jobs fit the budget at once. This is the check to read before a large run.

## 4. Alternatives Considered

### Alternative 1: A hand-written Snakemake profile per host
- **Description**: `profiles/omni/config.yaml` with `cores`, `resources`, `set-resources`, `set-threads`.
- **Pros**: no ob change; a format Snakemake users know.
- **Cons**: works for the budget only. Per-rule overrides need exact rule names, and ob's rule names are per node and change with every parameter hash.
- **Reason for rejection**: it cannot express "stage EMBED_M".

### Alternative 2: Generate the profile from a mapping
- **Description**: a tool reads a stage mapping and the generated Snakefile, and writes a profile with one `set-resources` line per rule.
- **Pros**: no change to expansion; the output is a plain Snakemake artefact.
- **Cons**: a two-step run, and a generated file that is stale as soon as the plan or the filter changes. The Snakefile records stage and module per rule, but not labels, so tiers would need more plumbing.
- **Reason for rejection**: more moving parts than resolving the same mapping inside ob, for the same result.

### Alternative 3: Keep everything in the plan
- **Description**: status quo, plus the #418 budget.
- **Pros**: nothing to build.
- **Cons**: the four costs of §1 remain.
- **Reason for rejection**: it is the problem statement.

### Alternative 4: Ranges with retries
- **Description**: a stage gets `mem_mb: [8000, 32000]`; a job starts at the low figure and is retried with more when it fails.
- **Pros**: fits a tier that spans a wide range of sizes; no over-reservation at the small end.
- **Cons**: only works if a job is killed at its limit. On a plain local run nothing does that; the job keeps allocating.
- **Reason for rejection**: deferred, not rejected. It depends on the enforcement work in §6. The host-file shape allows it later: a scalar now, a pair then.

## 5. Implementation Plan

1. **Budget and stage map.** `--host`, lookup by name and hostname, top-level `cores` / `mem_mb` / `gpu` / `capabilities`, `stages:` by stage and by `STAGE.module_id`, resolution order, the oversize check, the manifest entry, and `resources` dropped from `summary_hash`.
2. **Label selectors.** `<label>=<value>` blocks inside a stage entry.
3. **`snakemake_args`** and the executor rule of §3.6.
4. **`ob run --dry` report** of resolved resources and concurrency per stage.

Phases 1 and 4 are enough to retire per-host edits of a plan. Phase 2 is what makes one host file serve datasets of different size.

### Testing Strategy

- **Unit**: resolution order field by field; selector match and conflict; oversize error; lookup order; a host file leaves node ids and `summary_hash` unchanged; a plan `resources:` edit leaves `summary_hash` unchanged.
- **Generated Snakefile**: a stage sized by the host file gets the expected `threads:` and `resources:`; the command line carries the budget.
- **End to end**: a small plan with a stage declared at 40% of the budget and `--cores 8` starts two jobs at a time.

## 6. Open Questions

1. **Labels per parameter set.** A module can produce several datasets as parameter sets (one module, two `dataset_name` values). `provides` is per module, so those datasets share a tier. Either such modules are split, which changes their ids and paths, or `provides` learns a per-parameter-set value. This decides whether tiers are usable on existing plans without moving outputs.
2. **Selector syntax.** `tier=l:` as a key is compact but is a small language inside YAML keys. The alternative is a list of `{when: {tier: l}, mem_mb: …}` entries: more verbose, no key parsing, and an explicit order.
3. **Enforcement.** Where a hard limit would attach: a wrapper in the rule's shell block (`ulimit -v`, coarse; or `systemd-run --scope -p MemoryMax=`, needs a user session), switched on by a host-file key such as `enforce: true`. SLURM already enforces.
4. **Should the plan keep `resources:` at all?** It is useful as a default for people without a host file. It could become advisory ("needs about this much") rather than a directive.
5. **Units.** 004 §7.2 documents `mem_mb` as MiB; Snakemake's `mem_mb` is megabytes and it derives `mem_mib` from it. The host file should say MB and 004 should be corrected.

## 7. References

1. [Snakemake: resources](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#resources)
2. [Snakemake: profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles)
3. [004: YAML specification](004-yaml-specification.md) §7
4. [008: Filtering](008-filtering.md)
