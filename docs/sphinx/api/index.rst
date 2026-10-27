API reference
=============

Generated from the docstrings in ``src/ecco_insitu``.

.. automodule:: ecco_insitu
   :no-members:
   :no-index:

.. currentmodule:: ecco_insitu

Running the chain
-----------------

.. autosummary::
   :toctree: generated

   run
   pipeline
   cli

Settings, input checks and logging
----------------------------------

.. autosummary::
   :toctree: generated

   config
   validate
   log

The steps
---------

.. autosummary::
   :toctree: generated

   steps
   steps.step01
   steps.step02
   steps.step03
   steps.step04
   steps.step05
   steps.step06
   steps.step07
   steps.step08
   steps.step09
   steps.step10
   steps.step11_spatial_scaling
   steps.step12_profiles_compact
   steps.step13_pack_ecco_arrays

Model-ready files (step 13)
---------------------------

.. autosummary::
   :toctree: generated

   packing

Climatology
-----------

.. autosummary::
   :toctree: generated

   climatology
   climatology.setup
   climatology.build_woa23
   climatology.prebake
   clim_interp

WOD quality flags
-----------------

.. autosummary::
   :toctree: generated

   wod_flags

Uncertainties and utilities
---------------------------

.. autosummary::
   :toctree: generated

   high_latitude_sigma
   tools
   csv_to_nc

Diagnostics
-----------

.. autosummary::
   :toctree: generated

   diagnostics
   diagnostics.check_run_files
   diagnostics.check_ecco_arrays
   diagnostics.check_run_qc
   diagnostics.report_profile_counts
   diagnostics.report_costs
   diagnostics.report_spatial_scaling
   diagnostics.common
   diagnostics.costs_compute
   diagnostics.costs_report
   diagnostics.costs_scaling
   diagnostics.counts_index
   diagnostics.counts_maps
   diagnostics.counts_timeseries

Model comparison (work in progress)
-----------------------------------

.. autosummary::
   :toctree: generated

   model_comparison
   model_comparison.mds_to_nc
   model_comparison.cost
   model_comparison.cost_bars
   model_comparison.cost_versions
