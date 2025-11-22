# Release / packaging checklist (modern inmembrane)

This checklist is for **local / project releases** of the
`inmembrane` fork (SerraPHIM integration), not for the legacy
PyPI / DockerHub distribution.

---

## 1. Prepare a clean tree

- Use a fresh checkout or ensure your working tree is clean:

  ```bash
  git status
  ```

- Make sure `CHANGELOG` (if present) is updated with major changes
  since the last tag (scan the commit history if needed).

- Update the version string in:

  * `inmembrane/__init__.py`
  * `setup.py`

  They should both match (e.g. `0.96.0-dev`).

- Commit:

  ```bash
  git commit -a -m "Bump version to X.Y.Z and update changelog"
  ```

---

## 2. Tag the release

Create a tag and push it:

```bash
git tag inmembrane-X.Y.Z
git push
git push --tags
```

Use whatever tag naming scheme you prefer (`inmembrane-0.96.0-dev`,
`v0.96.0`, etc.), but keep it consistent.

---

## 3. Build and test from a clean virtualenv

Create a throwaway virtual environment and install from source:

```bash
python -m venv /tmp/inmembrane_venv
source /tmp/inmembrane_venv/bin/activate

# From the repo root:
pip install --upgrade pip
pip install .
```

Smoke test the CLI:

```bash
inmembrane_scan --help
inmembrane_scan --version
```

Run a **small Gram-negative test** (adjust paths as appropriate):

```bash
inmembrane_scan \
  --config ~/SerraPHIM_v2/tools/inmembrane/inmembrane.config \
  ~/SerraPHIM_v2/data/bakta_annotations/gram_neg_sample/proteome.faa
```

Then a **small Gram-positive test**:

```bash
inmembrane_scan \
  --config ~/SerraPHIM_v2/tools/inmembrane/inmembrane.config \
  ~/SerraPHIM_v2/data/bakta_annotations/gram_pos_sample/proteome.faa
```

Verify that:

- CSV and JSON files are created under the expected `out_dir`
- The `Category`, `LoopExtent`, and `Details` columns look reasonable
- `PhageReceptorCandidate=True/False` annotations appear where expected

If something fails here, fix it **before** pushing tags.

---

## 4. Optional: build source/wheel artifacts

If you want sdist/wheel artifacts (for internal distribution):

```bash
pip install build
python -m build
```

This creates:

- `dist/inmembrane-X.Y.Z.tar.gz`
- `dist/inmembrane-X.Y.Z-py3-none-any.whl`

You can share these internally or upload them to a private index if
needed.

---

## 5. Optional: publishing

At the moment this fork is intended primarily for **SerraPHIM internal
use**. If you later decide to push to PyPI or another public index:

- Create `~/.pypirc` and configure credentials
- Use `twine` to upload the files under `dist/`

This is **not** part of the standard workflow and should only be done
once the modern codebase is stable, documented, and decoupled from any
project-specific paths.
