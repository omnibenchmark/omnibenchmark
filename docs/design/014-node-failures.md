# 014: Node failures and timeouts

[![Status: Implemented](https://img.shields.io/badge/Status-Implemented-green.svg)](https://github.com/omnibenchmark/docs/design)
[![Version: 1](https://img.shields.io/badge/Version-1-blue.svg)](https://github.com/omnibenchmark/docs/design)

**Authors**: ben
**Date**: 2026-11-20
**Status**: Implemented
**Version**: 1
**Supersedes**: N/A
**Reviewed-by**: TBD
**Related Issues**: 008 §2 (no silent absence), 010 §3.5.1 (manifest), 013 (foreach), 015 (host files)

## Changes

| Version | Date | Description | Author |
|---------|------|-------------|--------|
| 1 | 2026-11-20 | Initial draft | ben |

## 1. Problem Statement

A module that exits non-zero fails its Snakemake job, and every node below it
is never run. Before `gather` that cost one branch. With a gather, one failed
member blocks the whole fan-in.

Observed on `clustering_example` (api 0.8.0, 1801 jobs): `FCPS_HDBSCAN_8`
crashes on 2 of 62 datasets. The `metrics` gather is grouped by dataset, so
those 2 datasets get no metrics for any of the other methods, and the report,
a gather over all metrics, is never rendered. `--continue-on-error` does not
help: it only keeps independent jobs going. The run of the same benchmark
archived with the paper lost only the two failing cells.

Three related gaps:

- **No timeout.** `ResolvedNode.timeout` and `LOCAL_TIMEOUT_VAR` exist but
  nothing reads them. Nine jobs of the same benchmark run for hours; the paper
  cut them at 7200 s.
- **No record.** What failed, and how, is only in the Snakemake log. Nothing
  in the output tree says "this cell is empty because the method timed out".
- **Workarounds live in modules.** The alternative today is for every method
  to catch its own errors and write an empty output, and for every consumer
  to skip empty files. An empty file can also be a legitimate result.

## 2. Design Goals

- A failed node is a **recorded result**, not a hole: one `Failure` per node,
  with kind, exit code and reason.
- Three questions answered for every failure: did it **fail to start**, did
  it **crash**, or did it **time out**.
- **Timeouts** enforced by ob.
- With `--continue-on-error`, the rest of the benchmark completes: a gather
  runs on the members that succeeded and is told which ones did not.
- **Modules do not change.** No error handling in methods, no skip logic in
  gather modules.
- **No silent absence** (008 §2): every failure is in the tree, in the gather
  manifest, in the run summary and in the exit status.
- Changing the timeout, or switching `--continue-on-error`, **reruns nothing
  that succeeded**.

### Non-Goals

- Automatic retries. `--retry-failed` reruns on request (§3.6).
- Failures outside the module process: environment creation, image pull,
  a scheduler killing the job. These stay Snakemake errors, reported as today.
- A per-stage or per-module timeout in the plan; it goes in host files (015).
- api < 0.8.0 plans and legacy metric collectors: their rules keep their
  text, so existing trees do not rerun.
- Validating that an output is well formed (002).

## 3. Proposed Solution

ob already generates the shell block of every rule. From api 0.8.0 that block
runs the module under a small wrapper which classifies the outcome and, on
failure, writes `failure.json` in the node directory. Everything else reads
that file.

### 3.1 The `Failure` model

`omnibenchmark/model/failure.py`. Plain data, immutable:

```python
class FailureKind(str, Enum):
    FAILED_TO_START = "failed_to_start"  # entrypoint, interpreter or package missing
    CRASHED = "crashed"                  # ran and exited non-zero, or was killed
    TIMED_OUT = "timed_out"              # killed by ob at the timeout
    MISSING_OUTPUT = "missing_output"    # exit 0, a declared output not written
    UPSTREAM = "upstream"                # not run: an input comes from a failed node


@dataclass(frozen=True)
class Failure:
    kind: FailureKind
    reason: str                      # one line, for people
    exit_code: Optional[int] = None  # None when the module was not run
    elapsed_s: Optional[int] = None  # wall time of the module process
    timeout_s: Optional[int] = None  # the limit in force, if any
    upstream: Tuple[str, ...] = ()   # node dirs of the failed producers
    log: Optional[str] = None        # the rule's log file, relative to out/
    node_dir: str = ""               # filled by the loader, not stored

    signal          # property: exit_code 137 -> 9
    load(node_dir, out_dir) -> Optional[Failure]
    to_dict() -> dict
```

The same module holds `collect` (every record of a plan), `clear` (for
`--retry-failed`), `summary` and `parse_duration`. `failure.json` is
`to_dict()`:

```json
{"kind": "timed_out", "reason": "timed out after 7200 s", "exit_code": 124,
 "elapsed_s": 7200, "timeout_s": 7200, "upstream": [],
 "log": ".logs/data_other_c0c1cea5_clustering_fcps_f20483e7.log"}
```

### 3.2 Classification

The wrapper decides, in this order:

| Condition | Kind | `exit_code` | `reason` |
|---|---|---|---|
| an input's producer has a `failure.json` | `upstream` | null | `input from a failed node` |
| a gather whose manifest has no member left | `upstream` | null | `every member failed` |
| timeout set and exit 124, or exit 137 at or past the limit | `timed_out` | as returned | `timed out after N s` |
| exit 126 / 127 | `failed_to_start` | 126 / 127 | `entrypoint not executable`, `entrypoint or interpreter not found` |
| non-zero exit and stderr matches a start pattern | `failed_to_start` | as returned | the matching line |
| exit above 128 | `crashed` | as returned | `killed by signal N`, `(out of memory?)` for 9 |
| any other non-zero exit | `crashed` | as returned | `exit code N: <last stderr line>` |
| exit 0 and a declared output does not exist | `missing_output` | 0 | `declared output <file> not written` |

Start patterns, a fixed list in the wrapper: `ModuleNotFoundError`,
`ImportError`, `can't open file` (Python), `there is no package called`,
`cannot open file` (R), `error while loading shared libraries`.

The reason carries one line of the module's stderr (the matching line, or the
last non-empty one; R's closing `Execution halted` is skipped), at most 300
characters. The full output stays in the rule log, which `log` points at.

### 3.3 The wrapper

`omnibenchmark/backend/ob_node.sh`, copied next to the Snakefile as
`.ob_node.sh`. Every rule sources it from the output root, then runs the
module through it:

```bash
source .ob_node.sh
…                                        # as before: paths, log, parameters.json
cd {params.module_dir}
OB_NODE=data.D/methods.M.ac0258f2 OB_UPSTREAM="data.D" OB_OUTPUTS="$OUTPUT_FILE_0" ob_run python3 {params.entrypoint} \
    --output_dir $OUTPUT_DIR …
```

A gather passes `OB_MANIFEST="$MANIFEST"` and only its own `inputs:`
producers as upstream.

On failure the wrapper writes `failure.json`, leaves every declared output as
an empty file, and:

- exits 0 when `OB_CONTINUE_ON_ERROR=1` (`ob run -k`): Snakemake sees a
  finished job and carries on;
- exits with the module's code otherwise: the job fails as it does today, and
  Snakemake removes the outputs.

The record is written in both modes. Downstream nodes and the manifest go by
`failure.json`, never by file size, so an empty output never passes as a
result.

- **Settings travel in the environment, not the rule text**: `OB_TIMEOUT` and
  `OB_CONTINUE_ON_ERROR`, set by `ob run`. Snakemake reruns a job when the
  rule's shell template changes, so changing either reruns nothing.
  Snakemake's apptainer call does not pass `--cleanenv`, so both reach the
  container.
- **Bash and coreutils only** (`timeout`, `grep`, `tr`, `mktemp`): the wrapper
  runs inside the module's environment, which may have no Python.
- **`timeout` signals its process group**, so children of the entrypoint die
  with it; `-k 30` sends KILL to whatever ignores TERM.
- **stderr is teed** to a temporary file in a synchronous pipeline, so the
  reason line is complete when the wrapper reads it; stdout and stderr still
  stream to the rule log.

Upgrading reruns existing api 0.8.0 trees once: every rule's text changes.

### 3.4 Timeouts

```
ob run plan.yml --timeout 2h      # 90s, 30m, 2h, 1d; the unit is required
```

One limit for every node of the run. The unit is required because
`resources.runtime` is in minutes, and a bare number would be read both ways.
The run manifest (`.metadata/manifest.json`) records it as `timeout_s`.

Per-stage timeouts belong in the host files of design 015; the command line
overrides them.

With `--timeout` and an executor or profile, `ob run`
adds `--default-resources runtime=<timeout + 5 min>`, so the scheduler's limit
outlasts ob's and the wrapper gets to write the record. A node whose plan
declares a shorter `runtime` keeps it, and Slurm may kill it first: it is then
an ordinary Snakemake error, without a record.

### 3.5 Downstream of a failure

Only with `--continue-on-error`; a strict run stops at the first failure.

- **A scatter node** whose input producer failed is not run. It gets its own
  `failure.json` (`upstream`) and empty outputs, so its children see it the
  same way.
- **A gather** runs on the members that succeeded. Its manifest rule now waits
  for the member files and, when it writes `lineage.json`, moves members that
  have a `failure.json` out of `members` into a new `failed` list, each with
  its record:

  ```json
  {"version": 2, "node_id": "…", "kind": "gather",
   "members": [ … ],
   "failed": [
     {"node_id": "…", "dir": "…", "files": { … }, "lineage": [ … ],
      "failure": {"kind": "crashed", "exit_code": 1, "reason": "exit code 1: boom", …}}
   ]}
  ```

  A gather module that ignores `failed` just sees fewer members. One that
  reads it can show the failure (a report: "HDBSCAN on iris5: crashed").
  The key is additive, so the manifest stays version 2.
- **A gather with no member left**, or whose own `inputs:` come from a failed
  node, is `upstream`.

### 3.6 Reruns

- With `-k`, `ob run` leaves recorded failures alone and lists them in its
  summary. (Without `-k` Snakemake removed the outputs, so the next run runs
  the node again anyway.)
- `ob run --retry-failed` deletes `failure.json` and the outputs of every
  failed node before starting Snakemake. `--retry-failed timed_out,crashed`
  limits it to those kinds.
- Nothing else is needed for the nodes below: when a retried node succeeds its
  outputs are newer, and Snakemake reruns its children and gathers.

### 3.7 Reporting

After Snakemake returns, in both modes, `ob run` loads `failure.json` for every
node of the plan and prints one line per failure, grouped by kind:

```
3 of 1801 nodes failed
  crashed (2)
    data.other.63b69cc8/clustering.fcps.76e2c75c   exit code 1: Error: Object was created without names.
    data.uci.63b4146f/clustering.fcps.76e2c75c   exit code 1: Error: Object was created without names.
  timed_out (1)
    data.other.c0c1cea5/clustering.fcps.f20483e7   timed out after 7200 s
Retry them with: ob run ... --retry-failed
```

`ob run` exits 1 when any node of the plan has a failure record, with or
without `-k`. A job Snakemake reports as failed without a record (the wrapper
never ran) is reported by Snakemake's own summary, as today.

## 4. Alternatives Considered

### Alternative 1: handle failures in each module
- **Description**: methods catch errors and write an empty output; gather
  modules skip empty files.
- **Pros**: no change to ob; about ten lines per module.
- **Cons**: repeated in every method and every consumer; empty is ambiguous;
  no kind, no exit code; does nothing for timeouts.
- **Reason for rejection**: the contract belongs in one place.

### Alternative 2: Snakemake checkpoints for every member
- **Description**: make each gathered stage a checkpoint and let the gather's
  input function list what exists.
- **Pros**: Snakemake-native; no tombstone files.
- **Cons**: the DAG is re-evaluated after every checkpoint job, 1,600 times
  on `clustering_example`; failed jobs still fail the run, so `--keep-going` semantics and the
  exit path get harder, not easier.
- **Reason for rejection**: cost and complexity for the same result.

### Alternative 3: tolerate failures only inside the gather
- **Description**: the gather rule lists members as optional inputs.
- **Pros**: scatter rules untouched.
- **Cons**: Snakemake has no optional inputs for a rule that must still wait
  for them; without a record there is nothing to tell "failed" from "not run
  yet"; no timeouts.
- **Reason for rejection**: needs the tombstone anyway.

### Alternative 4: always continue, no strict mode
- **Description**: the wrapper always exits 0.
- **Pros**: no environment variable.
- **Cons**: changes today's default; a typo in a new module would run the
  whole benchmark before anyone sees it.
- **Reason for rejection**: stop-at-first-failure is the right default while
  developing.

## 5. Implementation Plan

Done in one change:

1. `model/failure.py`: `FailureKind`, `Failure`, `collect`, `clear`,
   `summary`, `parse_duration`.
2. `backend/ob_node.sh`: the wrapper; `backend/snakemake.py` copies it, sources
   it and prefixes every api ≥ 0.8.0 module command with `ob_run`; the
   manifest rule waits for members and writes `failed`.
3. `cli/run.py`: `--timeout`, `--retry-failed`, the environment, the summary
   and exit status, `timeout_s` in the run manifest, `runtime` for executors.

Follow-ups, not in this change:

- Per-stage timeouts in host files (015).
- `ob collect`: a `failure` column in `performances.tsv`, so the two hours of
  a timed-out job are not read as a method's run time.
- Telemetry (006) marks the rule span failed from the record; `ob status`
  lists failures without running.

### Testing Strategy

- Unit (`tests/backend/test_node_failures.py`, about a second): the wrapper
  run with bash for every kind, the empty outputs, the JSON escaping, `-k`,
  the timeout; `Failure` round trip and `signal`; duration parsing; `collect`,
  `clear`, `summary`; the generated rules for 0.8.0 and 0.7.0.
- End to end (`tests/e2e/test_18_node_failures.py`, about 20 s): one module
  repo whose `--mode` picks the outcome. Five method nodes (ok, crash, sleep
  past a 3 s timeout, no output, failed import), a child stage under each and
  a gather over all five. With `-k`: the four kinds are recorded, the four
  children are `upstream`, the gather runs on one member and reads four
  `failed`, `ob run` exits 1. Without `-k`: the run stops and the record is
  there. Rerunning with another timeout runs nothing; `--retry-failed
  crashed` reruns that node and the gather above it.
- Acceptance: `clustering_example` at seed 20 with `-k --timeout 2h`: 62 of
  62 datasets in the report, and `tools/verify_zenodo.py` shows the same cells
  empty as the archive.

## 6. References

1. [010: gather](./010-gather.md), §3.5.1 the manifest
2. 015: host files (PR #419)
3. [008: filtering](./008-filtering.md), §2 no silent absence
4. GNU coreutils, `timeout`: exit status 124, 126, 127, 137
