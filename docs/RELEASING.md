# Releasing omnibenchmark

Internal dev guide for cutting a new release. Not part of the public docs site
(only `docs/src/` is published by mkdocs).

Versions are derived from git tags by `setuptools_scm` (see `pyproject.toml`:
`dynamic = ["version"]` and `[tool.setuptools_scm]`), so the **git tag is the
source of truth** — `pyproject.toml` is never bumped manually.

## Checklist

1. **Tests green.** Run both the fast suite and e2e:
   ```bash
   pixi run pytest -m "not e2e"
   pixi run pytest -m e2e tests/e2e/
   ```

2. **Finalize `CHANGELOG.md`.** Replace the `[X.Y.Z](UNRELEASED) ()` header
   with the release tag URL and date, e.g.:
   ```
   ## [0.5.1](https://github.com/omnibenchmark/omnibenchmark/releases/tag/v0.5.1) (May 8th 2026)
   ```

3. **Bump version strings in docs/README.** From the repo root:
   ```bash
   ./scripts/run_version_update.sh <OLD> <NEW>     # e.g. 0.5.0 0.5.1
   ```
   Updates `README.md`, `omni-environment.yml`, and `docs/src/*.md`.

4. **Refresh `ARCHITECTURE.md`.** Regenerate it from the package docstrings and
   imports so the layering diagram stays in sync (also fails on any import
   cycle):
   ```bash
   pixi run python scripts/gen_architecture.py
   pixi run python scripts/gen_architecture.py --check   # must exit 0
   ```
   Commit the result if it changed. (CI runs `--check` on every PR, so this is
   normally already current.)

5. **Commit and push to `main`.**
   ```bash
   git commit -am "chore: release vX.Y.Z"
   git push origin main
   ```

6. **Tag and push the tag.**
   ```bash
   git tag -a vX.Y.Z -m "Release X.Y.Z"
   git push origin vX.Y.Z
   ```

7. **Build the package.** Clean first to avoid stale artifacts polluting the
   upload:
   ```bash
   rm -rf dist/ build/ omnibenchmark.egg-info/
   pixi run python -m build
   pixi run twine check dist/*
   ```
   Verify the version embedded in the built artifacts matches the tag (no
   `+dirty` or `.devN` suffix — those mean the working tree wasn't clean or
   you forgot to tag).

8. **Upload to PyPI.** Optionally smoke-test on TestPyPI first:
   ```bash
   pixi run twine upload --repository testpypi dist/*
   ```
   Then the real upload:
   ```bash
   pixi run twine upload dist/*
   ```

9. **Create the GitHub release** for tag `vX.Y.Z` and paste the CHANGELOG entry
   as the body.

## Post-release

- Open a new `## [X.Y.(Z+1)](UNRELEASED) ()` section at the top of
  `CHANGELOG.md` so the next PR has somewhere to land its entry.
