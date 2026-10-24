"""
ecco_insitu: post-processing of in-situ ocean profiles for the ECCO state estimate.

The NCEI chain takes MITprof-format profile files (WOD, Argo, MEOP, ITP, moorings, ...) through
quality control, climatology comparison, uncertainty assignment and spatial scaling, and packs
the result into the files the ECCO model reads.

Main pieces:

- :mod:`ecco_insitu.pipeline` runs steps 01-10 on a folder of input files;
- :mod:`ecco_insitu.steps` holds the individual steps (01-13);
- :mod:`ecco_insitu.config` reads the settings (``defaults.yaml`` + ``config/local.yaml`` + run files);
- :mod:`ecco_insitu.validate` holds the input-file and climatology contracts;
- :mod:`ecco_insitu.climatology` checks, builds and pre-bakes the climatology;
- :mod:`ecco_insitu.packing` makes the model-ready files (step 13's engine);
- :mod:`ecco_insitu.diagnostics` checks and reports on a run.

Command line: ``ecco-insitu --help`` (see :mod:`ecco_insitu.cli`).
"""
__version__ = '1.0.0.dev0'


def provenance():
    """Text identifying the code that made a file: ``ecco_insitu <version>``.

    Written as the global attribute ``ecco_insitu_version`` in every NetCDF file the chain writes,
    and in ``run_config.yaml``, so any output can be traced to the code that produced it.
    """
    return f'ecco_insitu {__version__}'
