# API versions

`api_version` in a benchmark [plan](glossary.md#plan) selects the contract the plan is written
against: which YAML fields are valid, what command line each module receives,
and how the output tree is laid out. It is independent of the installed
`omnibenchmark` version: a release keeps running plans written for older api
versions.

```yaml
api_version: "0.8.0"
```

- A plan without `api_version` is read as `0.4.0`.
- `ob create benchmark` writes the latest version.
- `ob describe plan` prints the version a plan declares.
- Raising `api_version` to 0.8.0 changes every output path, so the next run
  recomputes the whole benchmark.

## What each version changes

| api version | Plan | Modules | Output tree |
|---|---|---|---|
| 0.4.0 and below | Legacy contract. | `--name` is the dataset (first-stage module) id. | Nested, `stage/module/params/`. |
| 0.5.0 | `{dataset}` in an output path warns. | `--name` is the module's own id. | `performance.txt` instead of `{dataset}_performance.txt`. |
| 0.6.0 | No change. | No change. | No change. |
| 0.7.0 | `provides` / `requires` lineage labels; `gather:` stages; a stage may join inputs from divergent branches. | A gather module receives `--<id> p1 p2 …`. | Fan-in nodes write `lineage.json`. |
| 0.8.0 | Stage and module ids are letters, digits and underscores; output `validators`; `metric_collectors` rejected; `group_by` groups by ancestor node; a grouped gather may declare `inputs`. | `--output id=PATH` per stage output; a gather module receives `--from-files <lineage.json>`. | Flat, one directory `stage.module[.params]` per node; `out/by-params/` view. |

Template variables `{module.id}`, `{module.name}`, `{module.stage}` and
`{params.KEY}` work at every api version.

## Deprecated and removed

| Feature | Status | Use instead |
|---|---|---|
| `metric_collectors` | Rejected from 0.8.0. | A stage with `gather:` and no `group_by`. |
| `{dataset}` in output paths | Warns from 0.5.0. | `{module.id}` or `{params.KEY}`. |
| Stage or module ids with `-`, `.` or other punctuation | Warn below 0.8.0, rejected from 0.8.0. | Letters, digits, underscores. |
| Nested output layout | Plans below 0.8.0 warn at `ob run`; removed at 0.10.0. | `api_version: "0.8.0"`. |
| `--output_dir` and `--name` module flags | Still passed; dropped at 0.10.0. | `--output id=PATH`. |
| `benchmark_yaml_spec` | Deprecated. | `api_version`. |
| `inputs: [{entries: [...]}]` | Deprecated. | `inputs: [id, ...]`. |
| `parameters: [{values: [...]}]` | Deprecated. | `parameters: [{key: value}]`. |
| `storage_api`, `storage_bucket_name` | Deprecated. | The nested `storage:` block. |

To move a plan and its modules to 0.8.0, see
[Upgrade a benchmark to api 0.8.0](howto.md#upgrade-a-benchmark-to-api-080).
