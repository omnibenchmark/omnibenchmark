# Configuration

Omnibenchmark keeps per-user settings in one file. They are about your
machine (where things are cached, how conda is reached, who you are), not
about a benchmark: nothing here changes a benchmark's results. Every setting
is optional.

## The file

| Platform | Location |
|---|---|
| Linux | `~/.config/omnibenchmark/omnibenchmark.cfg` (or under `$XDG_CONFIG_HOME`) |
| macOS | `~/Library/Application Support/omnibenchmark/omnibenchmark.cfg` |

It is an INI file: `[section]` headers, then `key = value` lines. An example
with every setting:

```ini
[user]
# author prefilled by `ob create benchmark` and `ob create module`
name = Ada Lovelace
email = ada@example.org

[dirs]
# cache of fetched module repositories, shared by all your benchmarks
git_cache = ~/.cache/omnibenchmark/git
# where `ob archive` checks module code out, relative to where it runs
git_modules = .omnibenchmark/git

[conda]
# build each conda environment once, shared by all benchmarks
cache = ~/.cache/omnibenchmark/conda
# fetch conda-forge, bioconda, ... from a mirror
channel_alias = https://mirror.example.org/conda
# reach the internet through a forward proxy
proxy = http://proxy.example.org:3128
```

## Read and change settings

Edit the file, or use `ob config` with keys written as `section.key`:

```sh
ob config user.name "Ada Lovelace"   # set
ob config user.name                  # print
ob config --list                     # print everything
ob config --unset user.name          # remove
```

## Sections

### `[user]`: who you are

| Key | Meaning | Default |
|---|---|---|
| `name` | Author name for new benchmarks and modules | git's `user.name` |
| `email` | Author email for new benchmarks and modules | git's `user.email` |

Used by `ob create benchmark` and `ob create module`. `--author-name` and
`--author-email` on the command line override both.

### `[dirs]`: where repositories are kept

| Key | Meaning | Default |
|---|---|---|
| `git_cache` | Cache of fetched module repositories, shared by all your benchmarks | `~/.cache/omnibenchmark/git` (or under `$XDG_CACHE_HOME`) |
| `git_modules` | Directory where `ob archive` checks out module code | `.omnibenchmark/git`, relative to where `ob` runs |

`ob run` fetches each module repository into `git_cache` once, then copies the
pinned commit into the run's own `out/.modules/`. Point `git_cache` at a
larger disk if your home directory is small. It is safe to delete:
repositories are fetched again when needed.

### `[conda]`: conda-backed runs

Applies to `ob run` and `ob validate outputs` on a plan with
`software_backend: conda`. Each key sets an environment variable for
Snakemake; a variable you already exported in the shell wins over the file.

| Key | Meaning | Sets | Default |
|---|---|---|---|
| `cache` | Existing directory where conda environments are built | `SNAKEMAKE_CONDA_PREFIX` | unset: environments are built inside each output directory (`out/.snakemake/conda`) |
| `channel_alias` | Base URL of a conda channel mirror | `CONDA_CHANNEL_ALIAS` | unset |
| `proxy` | `http(s)://host:port` forward proxy | `HTTP_PROXY`, `HTTPS_PROXY` | unset |

**`cache`** lets several benchmarks (or several users on a shared filesystem)
build each environment once. Environments are keyed by the hash of their env
file, so the same file gives the same environment.

```sh
mkdir -p ~/.cache/omnibenchmark/conda
ob config conda.cache ~/.cache/omnibenchmark/conda
```

The cache only grows; remove environment directories by hand to reclaim
space. Do not pass `--conda-cleanup-envs` to Snakemake with a shared cache:
it deletes environments from the shared directory, including ones other
benchmarks still need.

**`channel_alias`** makes channels named in env files (`conda-forge`,
`bioconda`) come from `<channel_alias>/<channel>`. Channels written as full
URLs are not affected.

**`proxy`** is for compute nodes without direct internet access. Module jobs
inherit it too.

## For developers

The same settings are available from Python through `omnibenchmark.config`:

```python
from omnibenchmark.config import config

value = config.get("conda", "cache", default=None)   # read
config.set("user", "name", "Ada Lovelace")           # change in memory
config.save()                                        # write the file
config.sections()                                    # all sections
config.options("conda")                              # keys of a section
```

When adding a setting to omnibenchmark, document it on this page.
