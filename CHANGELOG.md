# Changelog

All notable changes to `ecco_insitu`. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/). Every change to the code goes under **Unreleased**
until a version is tagged. Say whether a change alters the processed profiles ("**Science:**"),
since those changes need the regression reference regenerated on purpose
(`python tests/regression/make_reference.py`).

## [Unreleased]

### Added
- Production hardening, Phases 1–5 (details in `docs/production_hardening_report.md`):
  - regression and unit tests (`tests/`, `python -m pytest`);
  - settings files: `src/ecco_insitu/defaults.yaml` + `config/local.yaml` + `--config` / `--set`;
    every run saves `documentation/run_config.yaml`;
  - input checks against a written contract (`docs/input_contract.md`), `ecco-insitu validate`;
    per-run `documentation/run_summary.txt`; exit code 3 when input files failed;
  - climatology contract and `ecco-insitu climatology [--build | --prebake]`; pre-baked files are
    used only if made from the current climatology;
  - the `ecco_insitu` package (`src/` layout, `pip install -e .`) and the `ecco-insitu` command;
  - API documentation (Sphinx, `docs/sphinx/`).

### Added (work in progress)
- `src/ecco_insitu/model_comparison/` (`docs/model_comparison.md`): `ecco-insitu mds-to-nc` (model
  equivalents `.equi.data` → NetCDF, bit-identical to MITprof_gcm2nc.m on v4r6 iter47),
  `ecco-insitu model-cost` (offline cost as ECCO v4r6 `cost_profiles.F`) and `model-cost-bars`.
  Monthly means → `prof_Vestim` is planned.
- `ecco-insitu model-cost` also reports the anomaly cost with the spatial scaling factor and γ
  removed, and writes `cost_breakdown.csv` (by variable, file, year, latitude band, depth band);
  `ecco-insitu model-cost-versions` plots it. The model's own cost ("final") is unchanged.

### Changed
- Every NetCDF file the chain writes (yearly, step 11, step 13) has the global attribute
  `ecco_insitu_version` (e.g. `ecco_insitu 1.0.0.dev0`), and `run_config.yaml` names the version, so
  any output can be traced to the code that made it. Step 13 files' `created_by` is now
  `ecco_insitu.packing`. The regression reference was regenerated for this; **only these two global
  attributes changed**, and all data values are identical.
- Every script is now an `ecco-insitu` subcommand. `NCEI.py` → `ecco-insitu steps`;
  `run_ncei.sh` → `ecco-insitu run`, with options renamed: `-S` → `--no-scaling`, `-p` → `--plots`,
  `-c` → `--compact`, `-A` → `--no-pack`, `-C` → `--config`, `-s` → `--set`.
- Controller messages go through Python `logging`; the per-file step lines are unchanged.
- **Science:** none. Every data value is bit-identical to the outputs before the hardening work
  (commit 7fa2cb5); only the provenance attributes above were added.

### Known issues
- Recorded in `docs/known_issues.md` (KI-1 to KI-17), not yet fixed.

### Fixed
- `build_woa23_climatology` called `step06.calculate_potential_T`, which is in step05.
- Changing the climatology could silently keep using pre-baked files made from the old one.

### Removed
- The `pdb` imports and `set_trace` calls; the hard-coded `/Users/...` paths.

## Before 2026-10-23

See `PROJECT_LOG.md` for the history of the chain up to the hardening work (the double-conversion
fix of 2026-10-21, the high-latitude σ, step 13 packing, ...).
