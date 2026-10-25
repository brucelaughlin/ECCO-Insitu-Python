"""
Comparing the processed profiles with the ECCO model. **Work in progress.**

These tools turn the model's values at the profiles ("model equivalents", ``prof_Vestim``) into
NetCDF and compute the profile cost offline, exactly as the model's ``cost_profiles.F`` does. They
work and were checked against ECCO v4r6, but their interfaces may still change and they don't yet
have regression tests of their own.

- :mod:`~ecco_insitu.model_comparison.mds_to_nc` (``ecco-insitu mds-to-nc``): MITgcm
  ``*.equi.data`` tile files to NetCDF ``prof_Vestim``; a Python replacement for MITprof_gcm2nc.m.
  *Working*; bit-identical to the MATLAB on all 21 v4r6 iter47 datasets.
- :mod:`~ecco_insitu.model_comparison.cost` (``ecco-insitu model-cost``): the offline profile cost,
  as ECCO v4r6 ``cost_profiles.F`` (SAMPLESPLIT + CLIMMASK; anomaly and mean terms). *Working*;
  reproduces the v4r6 per-type cost bars, not yet checked against a model run's own cost output.
- :mod:`~ecco_insitu.model_comparison.cost_bars` (``ecco-insitu model-cost-bars``): bar charts of
  cost per datum by dataset type. *Working.*
- :mod:`~ecco_insitu.model_comparison.cost_versions` (``ecco-insitu model-cost-versions``): the anomaly
  cost per datum with the spatial scaling factor and γ removed, by type, year, latitude band or depth
  band, from ``model-cost``'s ``cost_breakdown.csv``. *Working.*
- *Planned:* monthly-mean model fields to ``prof_Vestim`` (for cost checks without a model run).
  *Not started*; waiting for monthly-mean fields.

Reference material: ``reference/mitgcm_cost_code/`` (the v4r6 and current MITgcm ``cost_profiles.F``
the cost follows) and ``reference/drafts/`` (an early, unused draft for gridded MDS output).
Plan and status: ``docs/model_comparison.md``.
"""
