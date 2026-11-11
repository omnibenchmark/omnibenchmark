# 016: Run hooks

[![Status: Draft](https://img.shields.io/badge/Status-Draft-yellow.svg)](https://github.com/omnibenchmark/docs/design)
[![Version: 1](https://img.shields.io/badge/Version-1-blue.svg)](https://github.com/omnibenchmark/docs/design)

**Authors**: ben
**Date**: 2026-11-07
**Status**: Draft
**Version**: 1
**Supersedes**: N/A
**Reviewed-by**: TBD
**Related Issues**: [#397](https://github.com/omnibenchmark/omnibenchmark/issues/397) (this design), [#141](https://github.com/omnibenchmark/omnibenchmark/issues/141) (module validation hooks), [#407](https://github.com/omnibenchmark/omnibenchmark/pull/407) (run record, 007 §4.6)

## Changes

| Version | Date | Description | Author |
|---------|------|-------------|--------|
| 1 | 2026-11-07 | Initial draft | ben |

## 1. Problem Statement

After a run, information about it has to go somewhere omnibenchmark doesn't
own. The case in #397 is `obcost`: plan id, version and hash, dataset sizes,
and per-stage time and RAM go into a sqlite database or JSON file that feeds a
cost model. Others: a tracker that links a run to remote executions (a cloud
job, a SLURM job id, their costs), CI annotations, a notification. And before
a run, something may need to stop it (a budget or quota check).

The information is already on disk after a run: `manifest.json` (007 §4), the
run log `runs.jsonl` with the selection and node set (007 §4.6, #407), the
telemetry trace if enabled (006), and snakemake's per-rule `performance.txt`
(wall time, `max_rss`, I/O, CPU time) next to each node's outputs. What is
missing is a trigger.

## 2. Design Goals

- Consumers live outside omnibenchmark, in any language.
- Hooks observe a run and may stop it before it starts; they never change
  what it computes.
- Behave like git hooks, so there is nothing new to learn.
- No new data format: hooks read the files ob already writes.

### Non-Goals

- An ob-owned run database (4.1), an in-process plugin API (4.2), per-rule
  events (4.3).
- Marking which stages produce datasets (#397's "tricky part"): a plan
  annotation, e.g. `role: data` (014, on `feat/convergence`).
- Module validation hooks (#141): per module, inside the DAG.

## 3. Proposed Solution

As in git: an executable in the hooks directory, named after the event.

| Hook | Runs | Arguments | Exit status |
|---|---|---|---|
| `pre-run` | after the plan is loaded, before modules are fetched or anything is written or executed | `<plan> <out_dir>` | non-zero aborts: ob exits 3, nothing is written |
| `post-run` | after snakemake exits, whether it succeeded or not | `<plan> <out_dir> <exit_code>` | ignored |

- **Configuration**: in `omnibenchmark.cfg` (the file `ob config` reads and
  writes, `$XDG_CONFIG_HOME/omnibenchmark/`), section `[hooks]`:

  ```ini
  [hooks]
  ; a directory of executables named after events
  path = ~/.config/omnibenchmark/hooks
  ; or a command per event, run with the same arguments
  post-run = obcost ingest
  ; seconds, per hook
  timeout = 60
  ```

  Comments go on their own lines: the config file has no inline comments, so
  `post-run = obcost ingest ; x` would pass `;` and `x` to the command.
  Set with `ob config hooks.path ~/...` or `ob config hooks.post-run "obcost ingest"`.
  `OB_HOOKS_PATH` overrides `hooks.path` (git's `core.hooksPath`). For an
  event, `hooks.<event>` wins over a file in `hooks.path`; a missing or
  non-executable file is skipped.
- **No default directory**: with nothing configured, no hooks run. ob never
  runs hooks from a plan or module repository, for the reason git doesn't run
  hooks from a clone: whoever runs `ob run` chose them.
- **Timeout**: `hooks.timeout` (default 60 s). ob kills a hook that runs
  longer; for `pre-run` that counts as a veto. Unlike git, because a stuck
  hook would otherwise hold a long run's slot (or a CI job) forever.
- **Veto exit code**: 3, distinct from a failed run (snakemake's codes) and
  from ob's usage errors, so a script can tell "stopped by a hook" apart.
- **Visible**: when any hook is active, `ob run` logs it at info level before
  resolving, with where each one comes from, e.g.
  `Running with hooks: pre-run (~/.config/omnibenchmark/hooks/pre-run), post-run (hooks.post-run: obcost ingest)`.
  With `--no-hooks`, it logs that hooks are skipped. A run is never changed
  by hooks without saying so.
- **Skipping**: `ob run --no-hooks` (git's `--no-verify`). `--dry` runs no
  hooks, as it writes nothing under `.metadata/` (007 §4.6).
- **Environment**: the caller's, plus `OB_HOOK`,
  `OB_HOOKS_OUT` (`<out_dir>/.hooks/`, see 3.1), and the plan's identity
  (3.3): `OB_PLAN_URL` (`provenance.canonical_url`, empty if the plan has
  none), `OB_PLAN_ID`, `OB_PLAN_VERSION`; and the run's slice as
  `OB_SELECTION`, JSON shaped like the run record's `selection` (007 §4.6:
  `{"until": …}`, `{"module": …}`, `{"filter": <sha256>}`), empty for a full
  run. The run id is in
  `<out_dir>/.metadata/manifest.json`. The working directory is the one
  `ob run` was started in. Paths are absolute.
- **Output**: inherited. A hook's stdout goes to ob's stderr (git does the
  same), so it never mixes with `--telemetry` on stdout.

### 3.1 Hooks don't change results

A hook SHOULD NOT modify the plan, the module repositories, the outputs or
`.metadata/`: results must be the same with or without hooks, and the run
record must describe what actually ran. A hook that needs to write files uses
`<out_dir>/.hooks/` (`OB_HOOKS_OUT`), which ob never reads and creates before
`post-run` (a `pre-run` hook creates it if it needs it, so a stopped run
leaves no out dir); it
travels with the outputs, and `ob validate outputs` and `ob archive` treat it
as not part of the results (5, step 3).

ob checks this cheaply rather than trusting it:

- **pre-run**: ob hashes the plan file before the hook and again after it. If
  it changed, the run stops with exit 3 ("plan modified by pre-run hook"):
  the resolved DAG no longer matches the file.
- **post-run**: ob records size and mtime of every file in `<out_dir>` except
  `.hooks/` (outputs, `.metadata/`, logs) before the hook, and compares after
  it. Any difference
  is logged as a warning naming the files ("modified by post-run hook").
  It can't fail the run, which has already finished.

Size and mtime don't catch a same-size rewrite that restores the mtime; the
checks are a guard against mistakes, not against a hook set on tampering,
which runs with the user's rights anyway (see Configuration above).

### 3.2 What a hook reads

A `post-run` hook reads what it needs from `<out_dir>`: the last line of
`.metadata/runs.jsonl` is this run (007 §4.6), the node ids in it lead to the
output and `performance.txt` paths, and `.metadata/manifest.json` has the host.
That is the whole contract, so it changes only when 007 does.

`obcost` as a hook, either `ob config hooks.post-run "obcost ingest"` (it gets
`<plan> <out_dir> <exit_code>`), or as a file in `hooks.path`:

```sh
#!/bin/sh
# post-run
plan=$1 out=$2 rc=$3
[ "$rc" = 0 ] || exit 0
exec obcost ingest "$out"      # reads .metadata/runs.jsonl, */performance.txt
```

### 3.3 Plan identity

A consumer that keeps runs of many plans (a cost model, a tracker) needs a
key for "the same benchmark". `id` is not enough: it's unique only within
whoever wrote the plan, and two groups can both have a `clustering` plan.
The key is the plan's `provenance.canonical_url` (004 §9), the place its
authoritative version lives, plus `version`; `summary_hash` (007 §4.6) then
tells revisions of that version apart. Plans without a `canonical_url` fall
back to `id`.

This also groups variants of one benchmark the way their authors intend:
the conda, CI and reduced-CI plans of `omni-scrna/split-stages-plan` all set
the same `canonical_url`, so their runs land under one benchmark, told apart
by `summary_hash` and selection. The run record's `benchmark` block (007
§4.6) carries `id`, `version` and `summary_hash` today; it should carry
`canonical_url` as well (5, step 5).

### 3.4 Remote runs

Hooks run where `ob run` runs. A run on a remote runner (a cluster job, a
cloud job that calls `ob run`) runs that machine's hooks, or none; what a hook
writes to `<out_dir>/.hooks/` comes back with the outputs, and a tracker on the
user's machine can join it with what the remote's own API knows. ob knows
nothing about remotes.

## 4. Alternatives Considered

### 4.1 ob owns the run database
- **Pros**: no user code; one schema.
- **Cons**: ob maintains a database and its migrations, and each consumer
  still wants its own fields. The out dir already keeps an append-only run log.
- **Rejected**: the run log covers what ran; the rest is consumer-specific.

### 4.2 In-process plugins (pluggy, entry points)
- **Pros**: rich data and return values.
- **Cons**: Python only, installed into ob's environment, coupled to internal
  objects that change between releases, loaded on every run.
- **Rejected**: can be layered on later if a backend needs in-process access.

### 4.3 Telemetry consumers, snakemake logger plugins
- **Pros**: per-rule detail; both exist.
- **Cons**: telemetry needs `--telemetry` and has no veto; a logger plugin sees
  rules, not ob's plan or run record.
- **Rejected**: complementary, for per-rule needs.

### 4.4 A JSON payload on stdin
- **Pros**: one read for the hook.
- **Cons**: a second schema that restates 007 §4.6 and must track it.
- **Rejected**: the files are the contract.

## 5. Implementation Plan

1. `omnibenchmark/core/hooks.py`: `find(event)` looks up `hooks.<event>`,
   then `<hooks.path>/<event>` (`OB_HOOKS_PATH` first); `run(event, plan,
   out_dir, *args)` runs it with stdout to stderr and `hooks.timeout`, and
   does the 3.1 checks (plan hash; size and mtime of the out dir around
   `post-run`).
2. `cli/run.py`: `--no-hooks`. `pre-run` in `_run_benchmark` right after
   the plan loads (resolution already writes `out/Snakefile`), exit 3 on
   non-zero. `post-run` around `_run_snakemake`, which exits through
   `sys.exit(returncode)` in both outcomes: `except SystemExit` takes the
   code, `finally` runs the hook. Hooks set `OB_PLAN_URL`, `OB_PLAN_ID`,
   `OB_PLAN_VERSION`.
3. Later: `ob validate outputs` and `ob archive` skip `<out_dir>/.hooks/`.
4. With #407: `post-run` hooks read `runs.jsonl`, and its `benchmark` block
   gains `canonical_url` (3.3).

### Testing Strategy

Hooks written into a temporary `OB_HOOKS_PATH`, and a temporary
`omnibenchmark.cfg` (`XDG_CONFIG_HOME`):

- `pre-run` exiting non-zero makes ob exit 3 before `.metadata/` exists;
  exiting 0 runs. A `pre-run` that outlives `hooks.timeout` is killed and
  vetoes.
- `hooks.post-run` in the config wins over `<hooks.path>/post-run`.
- Nothing configured: no hooks run and nothing about hooks is logged; with
  hooks configured, the info line names each active hook and its source.
- A `pre-run` that edits the plan makes ob exit 3; a `post-run` that touches
  an output or `.metadata/` gets a warning naming the file; writing to
  `OB_HOOKS_OUT` is silent.
- `post-run` gets the snakemake exit code, also after a failed run, and its
  own exit status doesn't change ob's.
- `--no-hooks`, `--dry`, and a hook without the executable bit run nothing.
- A hook's stdout ends up on stderr (`--telemetry` output stays valid NDJSON).

## 6. Future Work

### 6.1 Purging aborted runs

An aborted run (`post-run` got a non-zero exit code, or no end was ever
recorded) leaves a line in `runs.jsonl`, logs, and possibly outputs. Purging
it is not a hook: hooks don't change results (3.1). It would be an ob command
(e.g. `ob prune`), and it has to respect incremental runs:

- Outputs of rules that finished are valid results, whether or not the run as
  a whole failed, and later runs reuse them (same node, same path). They stay.
- Outputs of rules that didn't finish are already removed by snakemake on
  failure, unless `--keep-incomplete` kept them; snakemake knows which ones
  are incomplete (its metadata, `--rerun-incomplete`). Those can go.
- A file is attributable to a run only through its node (the run's
  `nodes-*.json`) and its write time within the run's window; a later run of
  the same node may have rewritten it. When in doubt, keep.
- `runs.jsonl` is append-only (007 §4.6): a purge appends a status line for
  the run instead of deleting it.

This needs the run record to carry the end of a run (exit code, end time),
which 007 §4.6 doesn't have yet, and it overlaps with its open TODO on outputs
orphaned by earlier plan revisions (`snakemake --list-untracked`,
`--list-changes`).

### 6.2 Pulling remote records into local management

Runs on remote runners leave their records (`runs.jsonl`, `.hooks/`) in
remote out dirs, and the runners have their own APIs (job status, cost, data
sets). Gathering them into one local view belongs to a tool outside ob, fed
by `post-run` hooks and by pulling from each remote, so that ob stays unaware
of remotes. What makes this easy is already in place: run records are keyed
by `run_id` and append-only, so merging records pulled from several places is
a union by `run_id`, and pulling twice changes nothing. Runs are grouped by
plan identity (3.3), so the same benchmark run from a fork, a mirror or a CI
variant lands in one place.

## 7. References

1. [#397: design a hook system for ob run](https://github.com/omnibenchmark/omnibenchmark/issues/397)
2. [Git hooks](https://git-scm.com/docs/githooks) (`core.hooksPath`, stdout to stderr, `--no-verify`); `ob config` (`omnibenchmark/cli/config.py`)
3. 004-yaml-specification.md §9 (provenance: `canonical_url`, `derived_from`, `subset_of`); 006-telemetry.md; 007-output-layout.md §4, §4.6 (#407)
4. 014-convergence.md on `feat/convergence` (stage roles)
