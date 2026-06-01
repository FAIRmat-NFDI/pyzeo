# Making a release

pyzeo ships a Cython-based C++ extension. Following [Cython's
recommendation](https://cython.readthedocs.io/en/latest/src/userguide/source_files_and_compilation.html#distributing-cython-modules),
the **pre-generated `src/pyzeo/extension.cpp` is committed and shipped** so end
users do not need Cython to build from source (`setup.py` sets
`USE_CYTHON = False`).

Releases are published **manually with `twine`**: bump the version, regenerate
the C++ from the `.pyx`, build wheels via the `wheels` GitHub Actions workflow,
then upload the sdist + wheels to PyPI.

## Prerequisites

```sh
python -m pip install --upgrade pip
pip install --upgrade setuptools wheel cython build twine
```

`setuptools`, `wheel` and `cython` are needed to regenerate the C++ source and
compile (step 1); `build` and `twine` are used to package and upload (steps 3,
5). On Python 3.12+ `distutils` was removed from the standard library, so an
up-to-date `setuptools` (which vendors a `distutils` shim) is required.

- A C++ compiler and the Python development headers (step 1 compiles the
  extension after regenerating the C++). On Debian/Ubuntu, for the Python
  version you build with: `sudo apt install pythonX.Y-dev` (e.g.
  `python3.13-dev`). Without them the build fails with
  `fatal error: Python.h: No such file or directory`.
- A PyPI account with upload rights to `pyzeo`, configured via `~/.pypirc` or a
  token passed to `twine`.
- The [`gh`](https://cli.github.com/) CLI authenticated to this repo (used to
  trigger the workflow and download artifacts).

## 1. Regenerate the Cython C++ source

This keeps the shipped `extension.cpp` in sync with the current `extension.pyx`
and a current Cython.

```sh
# 1. In setup.py set: USE_CYTHON = True
python setup.py build_ext --inplace --force
# 2. Revert setup.py back to: USE_CYTHON = False
```

Commit the regenerated `src/pyzeo/extension.cpp`. The `.so` produced is
gitignored — do **not** commit it.

## 2. Bump the version

Edit `version` in `pyproject.toml`. Update the `requires-python` value and the
`Programming Language :: Python :: 3.x` classifiers if the supported Python
range changed (also update the matrix in `.github/workflows/test.yml` and
`CIBW_SKIP` in `.github/workflows/wheels.yml` accordingly).

## 3. Build & sanity-check the sdist locally

```sh
python -m build --sdist
tar tzf dist/pyzeo-<version>.tar.gz | grep extension.cpp   # must be present
```

Optional but recommended — prove a from-source install compiles without Cython:

```sh
python -m venv /tmp/pyzeo-test && source /tmp/pyzeo-test/bin/activate
pip install dist/pyzeo-<version>.tar.gz
cd tests && pytest
deactivate
```

> Note: only `extension.cpp` (not `.pyx`/`.pxd`) ships in the sdist. The bundled
> `.cc`/`.h`/`.hh`/Eigen files come from the `setup.py` Extension sources +
> `MANIFEST.in`.

## 4. Build the wheels (CI)

Push the version bump + regenerated `.cpp` to `main` (or a release branch and
merge it). Then trigger the **wheels** workflow — it is `workflow_dispatch`-only:

```sh
gh workflow run wheels.yml --ref main
```

(Or via the GitHub Actions UI → "wheels" → Run workflow.)

When the run finishes, download the built wheels into `dist/` (one artifact per
OS: ubuntu, macos-13, macos-14):

```sh
gh run download <run-id> -n cibw-wheels-* -D dist/
```

## 5. Upload to PyPI

```sh
twine check dist/pyzeo-<version>*           # validate sdist + wheels
twine upload dist/pyzeo-<version>.tar.gz dist/pyzeo-<version>-*.whl
```

Upload the sdist and all wheels together. To rehearse first, upload to TestPyPI:
`twine upload --repository testpypi dist/pyzeo-<version>*`.

## 6. Tag the release

```sh
git tag v<version>
git push origin v<version>
gh release create v<version> --generate-notes   # optional GitHub Release
```

## Checklist

- [ ] Regenerated `extension.cpp` (USE_CYTHON True → build → False) and committed it
- [ ] Bumped `version` in `pyproject.toml` (+ classifiers / Python range if changed)
- [ ] `python -m build --sdist`; confirmed `extension.cpp` is in the tarball
- [ ] (Optional) clean-venv install of the sdist + `pytest` passes
- [ ] Pushed to `main`, ran `gh workflow run wheels.yml`, downloaded artifacts
- [ ] `twine check` + `twine upload` of sdist and all wheels
- [ ] Tagged `v<version>` and pushed the tag
