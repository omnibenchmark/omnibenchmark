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

- A **slice** names part of a run: up to a stage, and below one node (one dataset, say).
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

`ob run` writes all keys to `out/.metadata/nodes.json`, with each node's stage, directory, outputs and parent ids.

From api 0.8.0 every rule also leaves its node's key in `<node dir>/.key`. The file says under which key a directory was actually produced; `nodes.json` only says what the plan resolves to now. They differ after a module reference moves, or after a changed plan is resolved into an old out dir.

### 3.2 Slice

A slice is a set of nodes closed under ancestors. It is selected by:

| Selector | Nodes |
|---|---|
| a cut, `until STAGE` | nodes of stages up to `STAGE` in plan order |
| a root node | the root's **group**: the nodes whose only ancestor at the root's stage is the root |
| both | their intersection |
| neither | all nodes |

Ancestors of the selected nodes are always added. A gather over several groups has several ancestors at the root's stage, so it is in no group.

```
slice hash = sha256(sorted node keys)[:8]
```

### 3.3 Snapshot file

A snapshot is one slice in one file:

```
<benchmark>-<version>-<plan hash:8>-<label>_<slice hash>[~<base hash>].<ext>
clustering-1.2-94686c86-data.D1.94686c86.until-preprocess_deadbeef.tar.gz
```

- **Only the slice hash decides a match.** Benchmark, version and plan hash are for people. The label is the root's directory segment, then the cut; it tells a reader of the store which nodes can be roots.
- **Contents:** the files of the slice's nodes at their paths relative to `out/`, with no wrapping directory. A node's files exclude the directories of its child nodes, which the layout nests inside it.
- **Formats:** `.tar.gz` (default) and `.sqfs` (squashfs, zstd), chosen when packing.
- **`lineage.json` is rewritten when packed:** the absolute paths of its member files become relative to `out/`. In a run they stay absolute, for the module that reads them.
- **Written once.** A name that exists is left alone: the same name is the same content.
- **Checked when packed.** Every node of the slice must hold its declared outputs and a `.key` equal to the key in `nodes.json`. Otherwise nothing is packed.

### 3.4 Base

A snapshot packed with a base holds only the nodes its base slice lacks. The base has the same root, or no root when the base cut lies below the root's stage. Its hash is still the hash of the whole slice, and its name ends in `~<base hash>`. It is used only together with its base, which may itself have a base.

A base is named by its slice hash, never by a location. Any file with that hash serves, wherever it was fetched from and whether or not it has a base of its own. Where to find it is the business of the remotes and the lock (§3.8), so a snapshot on one remote can build on a base published on another.

### 3.5 Store

Snapshots live in a local directory, resolved in this order:

1. `snapshots.path` in `<plan dir>/.ob/config` (not committed; `ob config --local`);
2. `snapshots.path` in the user's config (`ob config`);
3. `<plan dir>/snapshots`.

A relative path is relative to the plan. A data remote (§3.8) mirrors this directory; nothing else reads a remote.

Packing appends one row per snapshot to `README.md` in the store: the file, the benchmark and plan hash, what the slice covers, how many nodes and which stages the file holds, its size, and the base it needs. It is for people; matching never reads it.

### 3.6 Using snapshots

`ob run PLAN --use-snapshots`, after the plan is resolved and before Snakemake starts:

1. Compute the slice hash for every root named by a file in the store (and none) × every stage cut (and none) of the resolved plan.
2. Keep the files in the store named by one of these hashes. A file with a base is kept only if its base is also there.
3. Apply them largest slice first, skipping a file that adds no node.
4. Apply a `.sqfs` by mounting it read-only under `out/.slices/<file>/` and symlinking each of its files into `out/`. Directories are created for real, so new child nodes are written next to cached files, and a write to a cached file fails. Without FUSE, or for a `.tar.gz`, the file is extracted.
5. Run Snakemake, which treats the nodes that are present as done.

- With `--dry`, steps 1–3 are reported and nothing is applied.
- Mounts stay after the run, since the results in `out/` point into them. `ob slice unmount` removes them; running again with `--use-snapshots` restores them.

### 3.7 Commands

| Command | Does |
|---|---|
| `ob archive PLAN --snapshot [--until S] [--group-by G] [--base-until S0] [--squashfs] [--out-dir DIR]` | pack a slice of `out/` into the store. `--snapshot` alone packs everything; the other options imply it. `--group-by G` packs one snapshot per node of stage `G` |
| `ob run PLAN --use-snapshots [--dry]` | start from matching snapshots |
| `ob slice unmount` | unmount the squashfs snapshots of `out/` |
| `ob config --local KEY [VALUE]` | read or write `<cwd>/.ob/config` |

`ob archive` in snapshot mode shares no code with the zip archive or its S3 path: it returns before them, and rejects their options.

### 3.8 Data remotes (sketch, not implemented)

A data remote holds snapshot files somewhere else. It lists names and moves bytes; it never looks inside a snapshot. Running and packing only ever read the store (§3.5): a fetch fills it, a push copies out of it.

**Model**

```python
@dataclass(frozen=True)
class RemoteFile:            # one snapshot file, as a remote holds it
    name: str                # the snapshot file name (§3.3): its identity
    locator: str             # what this remote needs to fetch it: URL, path, backend id
    size: int | None
    checksum: str | None     # "sha256:..." or "md5:...", when the remote has one


class DataRemote(Protocol):
    name: str                                            # the config section's name
    def list(self) -> list[RemoteFile]: ...
    def fetch(self, file: RemoteFile, into: Path) -> Path: ...   # via <name>.part; checksum verified
    def push(self, path: Path) -> RemoteFile: ...        # optional; an existing name is left alone
```

**Config**, in `<plan dir>/.ob/config` or the user's config:

```ini
[remote "lab"]
kind  = http
read  = https://omni.example.org/snapshots/clustering
write = ssh://omni/srv/www/snapshots/clustering

[remote "paper"]
kind = zenodo
read = 10.5281/zenodo.1234567
```

Credentials come from the environment (`ZENODO_TOKEN`, `THOA_API_KEY`, the ssh agent), never from a config file.

**Kinds**

| kind | `list` | `fetch` | `push` | Constraints |
|---|---|---|---|---|
| `dir` | directory listing | copy | copy | another local folder or a mounted share |
| `http` | the server's JSON directory listing of `read` | GET | over ssh to `write`: stream to `<name>.part`, then a non-overwriting rename | the listing has no checksum: size only, unless a `.sha256` file is pushed alongside |
| `zenodo` | the record's file list | GET | deposit API; each push publishes a new version of the record | published files are immutable; md5 checksums; per-record size and file-count limits |
| `thoa` | the lock (below) | download by file id | upload, which creates a new dataset | datasets cannot be listed by snapshot name; a thoa job takes one dataset, so a snapshot and its bases must be in the same dataset to be staged for a job |
| `s3` | objects under a prefix | GET | conditional PUT | would replace `remote/versioning.py` |

**Lock.** A remote that cannot list by name needs an index, and a reader of a published plan should not need credentials or a listing. Both are served by `snapshots.lock`, committed next to the plan, one `RemoteFile` per line:

```
<remote>  <name>  <size>  <checksum>  <locator>
```

`list()` of a lock-backed remote is its lines. A push appends one. For a remote that can list, the lock is optional and pins which snapshots a plan was published with.

**Flow**

- `ob run PLAN --use-snapshots`: compute the wanted slice hashes (§3.6 step 1, with root candidates taken from remote names as well); for each remote in config order, fetch the wanted names the store lacks, with their bases; then continue as in §3.6. With `--dry`, report what would be fetched and its size.
- `ob archive PLAN --snapshot ... --push REMOTE`: pack into the store, then push what the remote lacks.

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
- Not exercised: extraction of a `.sqfs` without FUSE, a base chain deeper than two, container backends end to end.

## 6. Deferred

- **Data remotes (§3.8).** None is implemented. `dir` and `http` first; the `zenodo` and `thoa` rows are from their documentation and client, not from a test.
- **Run record.** Which nodes came from which snapshot, in `runs.jsonl`.
- **Check keys when a snapshot is used.** Each node directory in a snapshot carries its `.key` (16 hex digits); comparing it with the plan's key would catch a file that matches on the 8-digit slice hash by accident, or was packed wrongly. Today only packing checks.
- **Container jobs and mounts.** The mount under `out/` must be visible inside the job. Apptainer can also bind the image itself, without a host mount.
- **A local hardlink store** (012 draft) for reuse across out dirs on one volume without extracting twice.
- **zstd for tar.** Needs the `tar` binary or Python 3.14.

## 7. Open questions

- **Lock or listing.** Should every remote go through `snapshots.lock`, so that what a plan uses is always pinned and reviewable, or only those that cannot list?
- **Timings.** A node taken from a snapshot keeps the `performance.txt` of the machine that ran it. Is a run with mixed origins flagged, or are timings from snapshots dropped from `ob collect`?
- **The key in the directory name.** `.key` catches a stale directory when packing. Putting the key in the name would make one impossible, but changes the layout.
- **Stage order.** `--until` cuts by plan order, then closes under ancestors. Is a cut by graph depth needed?
- **Dangling links.** After a reboot the symlinks into an unmounted image dangle, and a run without `--use-snapshots` recomputes those nodes. Should `ob run` refuse, or remount on its own?
