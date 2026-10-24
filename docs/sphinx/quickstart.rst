Quick start
===========

1. **Install** (once, in the repository folder; see :doc:`installing`):

   .. code-block:: bash

      pip install -e .

2. **Point the chain at your data.** Copy ``config/local.example.yaml`` to ``config/local.yaml``
   and set ``data_root`` (and any paths that differ). See :doc:`settings`.

3. **Get a climatology: you don't need one already.** Either

   - **use your own:** point ``paths.climatology_file`` at it. Any NetCDF file with
     ``potential_T_monthly`` and ``S_monthly`` on (12 months, depth, lat, lon) and a regular lon/lat
     grid works (:doc:`input_contract`, section 2); or
   - **build the standard one** on the fly:

     .. code-block:: bash

        ecco-insitu climatology --build      # downloads WOA23 (~1.4 GB), builds it where paths.climatology_file points
        ecco-insitu climatology --prebake    # optional: faster step 03
        ecco-insitu climatology              # check the setup

4. **Check your input files** (optional; every run does this anyway):

   .. code-block:: bash

      ecco-insitu validate /path/to/inputs

5. **Run:**

   .. code-block:: bash

      ecco-insitu run

   The output goes to a new timestamped folder, with the settings used
   (``documentation/run_config.yaml``) and a summary of what happened to every input file
   (``documentation/run_summary.txt``). See :doc:`running` and :doc:`outputs`.
