## What changed

<!-- One or two sentences. -->

## Checklist

- [ ] `python -m pytest` passes locally (including the regression tests, which need the chain data)
- [ ] **Science:** does this change the processed profiles?
  - [ ] No: the regression test passes unchanged.
  - [ ] Yes, on purpose: the reference is regenerated (`python tests/regression/make_reference.py`) and the
        reason is in `CHANGELOG.md` under **Science:**.
- [ ] `CHANGELOG.md` has an entry under **Unreleased**
- [ ] New or changed public functions have NumPy-style docstrings (the tests check this)
- [ ] User-facing behavior or settings changed? The docs are updated (`README.md`, `docs/`, `src/ecco_insitu/defaults.yaml` comments)
- [ ] A known issue fixed? It moves to "Fixed" in `docs/known_issues.md`
