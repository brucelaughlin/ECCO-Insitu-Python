"""Profile analysis framework — three-step pipeline.

Step 1 (:mod:`~ecco_insitu.diagnostics.profile_analyses.compute`, ``ecco-insitu profile-analyze``):
    For each obs class (source), read all chain output files and compute the 10 standard
    analyses (T/S mean vs depth, cost vs depth/time/bin/latitude, histogram, depth-range
    summaries) for vs-climatology costs (and vs-model equivalents if supplied).  Saves
    per-source NetCDF files and figures.

Step 2 (:mod:`~ecco_insitu.diagnostics.profile_analyses.compare`, ``ecco-insitu profile-compare``):
    Load two or more Step-1 output directories and generate overlay comparison figures.
    No re-computation — reads the saved NetCDF directly.

Step 3 (:mod:`~ecco_insitu.diagnostics.profile_analyses.report`, ``ecco-insitu profile-report``):
    Assemble figures from Step 1 (and optionally Step 2) into a single HTML + PDF document,
    one section per obs class, numbered 1–10 by analysis.  No interpretation text.
"""
