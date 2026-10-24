Outputs
=======

Each run writes one folder, ``<run.output_base>_<YYYYMMDD_HHMMSS>`` (or the ``-d`` folder):

``yearly/``
   One NetCDF file per input file, keeping the input's sub-folder layout (``yearly/CTD_WOD/...``),
   all on the input depth grid. Named ``<input name>__ncei_step_11.nc`` after step 11.

``ecco_arrays/``
   The model-ready profile files (step 13): one provenance per file, in time order, labeled in
   whole years, e.g. ``CTD_WOD_001_19920101_19921231.nc``.

``documentation/``
   - ``run_config.yaml``: the merged settings used;
   - ``run_summary.txt``: per input file, written / skipped (and why) / no data left / failed (at
     which step, with the error), plus its warnings;
   - ``PROFILES_SIZE_and_data.profiles.txt``: the ``NOBSGLOB`` / ``NFILESPROFMAX`` values and the
     ``data.profiles`` file list for the model;
   - ``ecco_arrays_manifest.json`` / ``.csv`` and ``optimize_sweep.csv``: what step 13 did.

Variables the chain adds
------------------------

- ``prof_bin_id_a``, ``prof_bin_id_b``: geodesic bins (10,242 and 2,562 bins)
- ``prof_Tclim``, ``prof_Sclim``: the climatology at each profile and depth
- ``prof_Tuncertainty``, ``prof_Suncertainty``: σ; ``prof_Vuncertainty_source``: 0 CTD σ,
  1 Arctic file, 2 southern constant
- ``prof_Tweight``, ``prof_Sweight``: the weights, γ/σ² × the step 11 factor
- ``prof_Tweight_code``, ``prof_Sweight_code``: why a weight is zero (step 06 QC codes, as bits)
- ``prof_area_gamma``: the LLC90 area factor γ (step 10)
- ``prof_Tspatial_scaling_factor``, ``prof_Sspatial_scaling_factor``: step 11 factors
- ``prof_Tcost``, ``prof_Scost``: climatology misfit cost, before γ and step 11
- ``prof_interp_*``, ``prof_point``: the model grid point of each profile (step 01)
- ``prof_T`` attribute ``T_conversion``: how step 05 decided whether to convert to potential T
