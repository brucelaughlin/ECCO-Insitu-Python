"""The NCEI chain's steps.

Steps 01-10 run per input file (driven by :mod:`ecco_insitu.pipeline`); each has a ``main(ds, ...)``
that takes and returns an ``xarray.Dataset``:

====  ==========================================================================
01    model grid point and tile for each profile (LLC90)
02    geodesic bin ids
03    monthly climatology at each profile (WOA23)
04    uncertainties (sigma) and weights, including the high-latitude sigma
05    in-situ to potential temperature, if needed (metadata, then a climatology check)
06    quality-control tests that zero weights (dates, positions, climatology cost, ...)
07    drop profiles with no non-zero weight
08    optional removal of empty trailing depth levels
09    sub-daily profiles at one place reduced to once daily
10    LLC90 area factor gamma on the weights
====  ==========================================================================

Steps 11-13 run on a whole run folder after 01-10:

====  ==========================================================================
11    spatial scaling of the weights by profile density per bin and month
12    optional compaction of the files
13    packing into model-ready files (see :mod:`ecco_insitu.packing`)
====  ==========================================================================
"""
