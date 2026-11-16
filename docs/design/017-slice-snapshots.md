# 017: Slice snapshots

[![Status: Draft](https://img.shields.io/badge/Status-Draft-yellow.svg)](https://github.com/omnibenchmark/docs/design)
[![Version: 1](https://img.shields.io/badge/Version-1-blue.svg)](https://github.com/omnibenchmark/docs/design)

**Authors**: btraven
**Date**: 2026-11-16
**Status**: Draft
**Version**: 1
**Supersedes**: the mechanism of [003](./003-storage.md) §3 and §6 (its goals stand); the unmerged drafts "013 remote storage" (`scratch/`) and "012 incremental snapshots" (`feat/snapshots`)
**Reviewed-by**: TBD
**Related Issues**: #404

## Changes

| Version | Date | Description | Author |
|---------|------|-------------|--------|
| 1 | 2026-11-16 | Initial draft, written against `feat/slice-snapshots` | btraven |

## 1. Problem

- **A run recomputes everything.** The expensive, stable part of a plan (datasets, preprocessing) is identical across runs; only the cheap part (a new method, a parameter) changes.
- **Part of a run cannot be shared.** "The preprocessed datasets at v1.2" has no name and no file.
- **Publishing depends on S3 features.** 003 versions a live bucket with object versioning, locks and tags, and decides membership from predicted paths, which are wrong for joins and gathers.

## 2. Goals

- A **slice** names part of a run, by stage and by dataset.
- A slice is stored as **one file**, in a plain directory.
- A run **starts from** the stored slices that are valid for its plan and computes the rest.
- Validity is **exact**: decided by a hash, not by a version convention.
- No new dependency. No change to the plan YAML.

### Non-goals

- Running Snakemake against remote storage (`ob run --use-remote-storage`). Untouched here.
- Deduplication below the slice (per file, per block).
- Validating outputs (002).

## 3. Design

### 3.1 Node key

Every resolved node has a key:

```
key(n) = sha256(repository, commit, entrypoint, environment, parameters,
                node directory, output paths, sorted keys of n's parents)[:16]
```

- **Parents** are the node's lineage parent, its join partners and its gather members. The key therefore commits to everything upstream of the node.
- **Environment** is the content hash of the conda file, or the reference string for other backends.
- **A module run from a working copy** (`--dirty`) gets a random key: it is never reused.
- The key is **not** the hash in the directory name, which covers the node's own parameters only.

`ob run` writes all keys to `out/.metadata/nodes.json`, with each node's stage, directory, dataset label and parent ids.

### 3.2 Slice

A slice is a set of nodes closed under ancestors. It is selected by:

| Selector | Nodes |
|---|---|
| `--until STAGE` | nodes of stages up to `STAGE` in plan order |
| `--dataset D` | nodes whose `dataset` label is `D` |
| both | their intersection |
| neither | all nodes |

Ancestors of the selected nodes are always added. A gather over several datasets has no dataset label, so it belongs to no dataset slice.

```
slice hash = sha256(sorted node keys)[:8]
```

### 3.3 Snapshot file

A snapshot is one slice in one file:

```
<benchmark>-<version>-<plan hash:8>-<label>_<slice hash>[~<base hash>].<ext>
clustering-1.2-94686c86-D1.until-preprocess_deadbeef.tar.gz
```

- **Only the slice hash decides a match.** Benchmark, version, plan hash and label are for people.
- **Contents:** the files of the slice's nodes at their paths relative to `out/`, with no wrapping directory. A node's files exclude the directories of its child nodes, which the layout nests inside it.
- **Formats:** `.tar.gz` (default) and `.sqfs` (squashfs, zstd), chosen when packing.
- **Written once.** Packing onto an existing name fails.

### 3.4 Base

A snapshot packed with a base holds only the nodes its base slice lacks. Its hash is still the hash of the whole slice, and its name ends in `~<base hash>`. It is used only together with its base, which may itself have a base.

### 3.5 Store

Snapshots live in a local directory, resolved in this order:

1. `snapshots.path` in `<plan dir>/.ob/config` (not committed; `ob config --local`);
2. `snapshots.path` in the user's config (`ob config`);
3. `<plan dir>/snapshots`.

A relative path is relative to the plan. A remote (§6) mirrors this directory; nothing reads a remote directly.

### 3.6 Using snapshots

`ob run PLAN --use-snapshots`, after the plan is resolved and before Snakemake starts:

1. Compute the slice hash for every dataset (and none) × every stage cut (and none) of the resolved plan.
2. Keep the files in the store named by one of these hashes. A file with a base is kept only if its base is also there.
3. Apply them largest slice first, skipping a file that adds no node.
4. Apply a `.sqfs` by mounting it read-only under `out/.slices/<file>/` and symlinking each of its files into `out/`. Directories are created for real, so new child nodes are written next to cached files, and a write to a cached file fails. Without FUSE, or for a `.tar.gz`, the file is extracted.
5. Run Snakemake, which treats the nodes that are present as done.

- With `--dry`, steps 1–3 are reported and nothing is applied.
- Mounts stay after the run, since the results in `out/` point into them. `ob slice unmount` removes them; running again with `--use-snapshots` restores them.

### 3.7 Commands

| Command | Does |
|---|---|
| `ob archive PLAN --snapshot [--until S] [--dataset D] [--base-until S0] [--squashfs] [--out-dir DIR]` | pack a slice of `out/` into the store. `--snapshot` alone packs everything; the other options imply it |
| `ob run PLAN --use-snapshots [--dry]` | start from matching snapshots |
| `ob slice unmount` | unmount the squashfs snapshots of `out/` |
| `ob config --local KEY [VALUE]` | read or write `<cwd>/.ob/config` |

`ob archive` in snapshot mode shares no code with the zip archive or its S3 path: it returns before them, and rejects their options.

## 4. Alternatives

- **Version convention (012 draft).** Reuse when `(major, minor)` match, with a prefix hash as a warning. A promise where the node key gives a proof.
- **Content-addressed blobs and a manifest per version (013 draft).** Deduplicates per file and allows partial pulls, but needs a blob store and a manifest format, and a slice is no longer a file one can copy, mount or upload to an archive service.
- **One file per node.** Maximum reuse, thousands of files. Slices with a base get most of it.
- **One file per version only.** No reuse when anything changes. It is the slice with no selector.
- **Overlay filesystem over the squashfs image.** Makes the whole tree look writable and copies on write. The symlinks keep cached files read-only.
- **S3 object versioning (003).** See §1.

## 5. Implementation notes

State of `feat/slice-snapshots` (based on `integ/0.8`), verified on `tests/e2e/configs/15_flat_layout.yaml` with Snakemake 9.6 and the host backend:

- §3.1–§3.7 are implemented in `omnibenchmark/slices.py`, with hooks in `cli/run.py`, `cli/archive.py`, `cli/config.py` and `cli/slice.py`.
- Regression test: `tests/e2e/test_19_snapshots.py` (a base and a snapshot above it, tar.gz and squashfs, then a changed plan run from them).
- A full snapshot gives a run with 0 jobs. After a parameter change in a late stage, a snapshot up to an earlier stage still matches and only the later nodes run. The same holds for a mounted `.sqfs`, and for a base `.sqfs` under a `.tar.gz` that has it as base.
- Snakemake did not rerun present nodes although `out/.snakemake` was empty, for extracted files and for symlinks. The symlinks carry the mtimes of their targets.
- `mksquashfs` ≥ 4.5 is needed (it reads the tar stream). `squashfuse_ll` is looked up on `PATH` and in apptainer's `libexec`.
- Not exercised: `--dataset` end to end, conda and container backends, extraction of a `.sqfs` without FUSE, a base chain deeper than one.

## 6. Deferred

- **Remote.** `write = ssh://…` and `read = https://…` under `[remote "name"]` in the config: push copies files out of the store, fetch downloads the matching names into it. Listing over https uses the web server's JSON directory listing. The read URL may later also live in the plan.
- **S3 as a third transport** of the same four operations (list, exists, get, put once), replacing `remote/versioning.py`.
- **Absolute paths in `lineage.json`** of gathers and joins point into the `out/` that produced them. They must become relative before a snapshot leaves its machine.
- **Run record.** Which nodes came from which snapshot, in `runs.jsonl`.
- **Refuse to pack a node with missing outputs.** Today the files present are packed as they are.
- **Container jobs and mounts.** The mount under `out/` must be visible inside the job. Apptainer can also bind the image itself, without a host mount.
- **Slices by any lineage label**, not only `dataset`.
- **A local hardlink store** (012 draft) for reuse across out dirs on one volume without extracting twice.
- **zstd for tar.** Needs the `tar` binary or Python 3.14.

## 7. Open questions

- **Timings.** A node taken from a snapshot keeps the `performance.txt` of the machine that ran it. Is a run with mixed origins flagged, or are timings from snapshots dropped from `ob collect`?
- **Where the key lives.** In `nodes.json` only (today), or also in the directory name, which makes a stale directory impossible but changes the layout.
- **Stage order.** `--until` cuts by plan order, then closes under ancestors. Is a cut by graph depth needed?
- **Dangling links.** After a reboot the symlinks into an unmounted image dangle, and a run without `--use-snapshots` recomputes those nodes. Should `ob run` refuse, or remount on its own?
