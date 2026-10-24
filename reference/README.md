# Reference material (not part of the package)

- `mitgcm_cost_code/`: the model's profile-cost code that `ecco-insitu model-cost` follows.
  - `cost_profiles_ECCO_v4r6.F`: ECCO v4r6 (ECCO-GROUP/ECCO-v4-Configurations,
    `ECCOv4 Release 6/code/cost_profiles.F`); verified identical to the GitHub version on 2026-10-23.
  - `cost_profiles_MITgcm.F`: the current MITgcm `pkg/profiles` version, for comparison.
- `drafts/xmitgcm_gridded_mds_to_nc_draft.py`: an early draft (from another assistant) that converts
  *gridded* MDS output with `xmitgcm`. Not used: the task turned out to be the profiles package's
  `.equi.data` files (see `docs/model_comparison.md`).
