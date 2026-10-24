Running the chain
=================

The whole chain
---------------

``ecco-insitu run`` (:mod:`ecco_insitu.run`) runs steps 01–10 on every input file, then steps 11–13
on the whole run folder:

.. code-block:: bash

   ecco-insitu run                              # input run.input_dir, output run.output_base_<timestamp>
   ecco-insitu run --config my_run.yaml         # settings for this run from a file
   ecco-insitu run --set chain.new_floor.prof_S=0.006
   ecco-insitu run -i /path/to/input -d /path/to/output -n 8
   ecco-insitu run --no-scaling                 # skip step 11
   ecco-insitu run -w 10day                     # step 11 with 10-day windows
   ecco-insitu run --plots                      # also the step 11 figures
   ecco-insitu run --compact                    # also step 12
   ecco-insitu run --no-pack                    # skip step 13

Before processing anything, the run checks the settings and data paths, the climatology, and every
input file (:doc:`input_contract`). A file that can't be processed is skipped with the reason; the
rest of the run carries on.

Everything printed goes to the screen and to ``z_logs/ncei_log_<YYYYMMDD_HHMMSS>.log``
(``run.log_dir``).

**Exit codes:** 0 all good; 3 finished, but some input files failed (see
``documentation/run_summary.txt``); 2 bad options or settings (nothing run); 1 a whole-run step
(11–13) failed.

All commands
------------

``ecco-insitu <command> --help`` shows each command's options.

=====================  =====================================================================
``run``                the whole chain (steps 01–13)
``steps``              steps 01–10 only (:mod:`ecco_insitu.pipeline`)
``scale``              step 11, spatial scaling, on a run folder
``compact``            step 12, compaction of one file
``pack``               step 13, the model-ready files
``optimize``           the step 13 packing program with all its options (:mod:`ecco_insitu.packing`)
``validate``           check input files, or a climatology, against the contract
``climatology``        check / ``--build`` / ``--prebake`` the climatology
``config``             print the merged settings
``check-files``        is a run complete and well-formed?
``check-arrays``       are the step 13 files a faithful repacking?
``check-qc``           slower scientific checks
``report-counts``      profile counts report
``report-costs``       cost report
``report-scaling``     step 11 figures
``csv-to-nc``          convert WOD CSV exports to NetCDF
=====================  =====================================================================
