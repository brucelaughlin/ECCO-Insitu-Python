"""
The ``ecco-insitu`` command: one entry point for the chain and its tools.

Run ``ecco-insitu <command> --help`` for each command's options. Every command that uses the
settings accepts ``--config FILE`` and ``--set key=value``.

==========================  =============================================================
``run``                     the whole chain on a folder of inputs (steps 01-13)
``steps``                   steps 01-10 only
``scale``                   step 11, spatial scaling, on a run folder
``compact``                 step 12, optional compaction of one file
``pack``                    step 13, model-ready files
``optimize``                the step 13 packing program, with all its options
``validate``                check input files (or a climatology) against the contract
``climatology``             check / build (``--build``) / pre-bake (``--prebake``) the climatology
``wod-flags``               check / download (``--download``) / index (``--build``) the WOD quality flags
``config``                  print the merged settings
``check-files``             diagnostics: run output complete and consistent
``check-arrays``            diagnostics: step 13 files correct
``check-qc``                diagnostics: slower scientific checks
``report-counts``           diagnostics: profile counts report
``report-costs``            diagnostics: cost report
``report-scaling``          diagnostics: step 11 figures
``csv-to-nc``               convert WOD CSV exports to MITprof NetCDF
``mds-to-nc``               *work in progress:* model equivalents to NetCDF
``model-cost``              *work in progress:* offline profile cost (ECCO v4r6)
``model-cost-bars``         *work in progress:* cost bar charts
``model-cost-versions``     *work in progress:* cost with gamma / spatial scaling removed
``profile-analyze``         T/S profiles and cost distributions per obs class (NetCDF + figures)
``profile-compare``         overlay comparison of multiple profile-analyze outputs
``profile-report``          compile profile-analyze/compare figures into HTML + PDF
==========================  =============================================================
"""
import importlib
import sys

# command -> (module, function, one-line help)
COMMANDS = {
    'run': ('ecco_insitu.run', 'cli', 'the whole chain on a folder of inputs (steps 01-13)'),
    'steps': ('ecco_insitu.pipeline', 'cli', 'steps 01-10 on a folder of input files'),
    'scale': ('ecco_insitu.steps.step11_spatial_scaling', 'cli', 'step 11: spatial scaling of a run folder'),
    'compact': ('ecco_insitu.steps.step12_profiles_compact', 'cli', 'step 12: compact one output file'),
    'pack': ('ecco_insitu.steps.step13_pack_ecco_arrays', 'cli', 'step 13: model-ready files'),
    'optimize': ('ecco_insitu.packing', 'cli', 'the step 13 packing program, all options'),
    'validate': ('ecco_insitu.validate', 'cli', 'check input files or a climatology against the contract'),
    'climatology': ('ecco_insitu.climatology.setup', 'cli', 'check / --build / --prebake the climatology'),
    'wod-flags': ('ecco_insitu.wod_flags', 'cli', 'check / --download / --build the WOD quality-flag files'),
    'config': ('ecco_insitu.config', 'cli', 'print the merged settings'),
    'check-files': ('ecco_insitu.diagnostics.check_run_files', 'cli', 'diagnostics: run output complete'),
    'check-arrays': ('ecco_insitu.diagnostics.check_ecco_arrays', 'cli', 'diagnostics: step 13 files correct'),
    'check-qc': ('ecco_insitu.diagnostics.check_run_qc', 'cli', 'diagnostics: scientific checks (slow)'),
    'report-counts': ('ecco_insitu.diagnostics.report_profile_counts', 'cli', 'diagnostics: profile counts report'),
    'report-costs': ('ecco_insitu.diagnostics.report_costs', 'cli', 'diagnostics: cost report'),
    'report-scaling': ('ecco_insitu.diagnostics.report_spatial_scaling', 'cli', 'diagnostics: step 11 figures'),
    'csv-to-nc': ('ecco_insitu.csv_to_nc', 'cli', 'convert WOD CSV exports to NetCDF'),
    'mds-to-nc': ('ecco_insitu.model_comparison.mds_to_nc', 'cli',
                  '[work in progress] model equivalents (*.equi.data) to NetCDF'),
    'monthly-estim': ('ecco_insitu.model_comparison.monthly_estim', 'cli',
                      '[work in progress] approximate prof_Vestim from monthly-mean model fields'),
    'model-cost': ('ecco_insitu.model_comparison.cost', 'cli',
                   '[work in progress] offline profile cost as ECCO v4r6 cost_profiles.F'),
    'model-cost-bars': ('ecco_insitu.model_comparison.cost_bars', 'cli',
                        '[work in progress] bar charts of model-cost results'),
    'model-cost-versions': ('ecco_insitu.model_comparison.cost_versions', 'cli',
                            '[work in progress] cost per datum with weight factors removed, by type/year/...'),
    'profile-analyze': ('ecco_insitu.diagnostics.profile_analyses.compute', 'cli',
                        'profile analyses: T/S profiles, cost distributions, saved as NetCDF + figures'),
    'profile-compare': ('ecco_insitu.diagnostics.profile_analyses.compare', 'cli',
                        'profile analysis comparison: overlay figures from multiple analyze outputs'),
    'profile-report':  ('ecco_insitu.diagnostics.profile_analyses.report', 'cli',
                        'profile analysis report: compile figures into HTML + PDF'),
}



def _usage():
    from ecco_insitu import __version__
    lines = [f'ecco-insitu {__version__}: the ECCO in-situ profile chain', '',
             'usage: ecco-insitu <command> [options]     (ecco-insitu <command> --help for options)', '',
             'commands:']
    lines += [f'  {c:15s} {h}' for c, (_, _, h) in COMMANDS.items()]
    lines += ['', 'Settings: src/ecco_insitu/defaults.yaml, config/local.yaml (copy config/local.example.yaml), '
              '--config FILE, --set key=value.',
              "No climatology yet? 'ecco-insitu climatology --build' makes the standard one."]
    return '\n'.join(lines)


def main(argv=None):
    """Entry point of the ``ecco-insitu`` command: dispatch ``ecco-insitu <command> ...`` to that
    command's ``cli(argv)``.

    Parameters
    ----------
    argv : list of str, optional
        ``[command, *options]`` (default: ``sys.argv[1:]``). No command, ``-h``, ``--help`` or
        ``help`` prints the list of commands; ``--version`` prints the version.

    Returns
    -------
    int
        The command's exit code (2 for an unknown command).
    """
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ('-h', '--help', 'help'):
        print(_usage())
        return 0
    if argv[0] == '--version':
        from ecco_insitu import __version__
        print(__version__)
        return 0
    cmd, rest = argv[0], argv[1:]
    if cmd not in COMMANDS:
        print(f"ecco-insitu: unknown command '{cmd}'\n\n{_usage()}", file=sys.stderr)
        return 2
    mod, fn, _ = COMMANDS[cmd]
    sys.argv = [f'ecco-insitu {cmd}', *rest]          # so each command's --help shows the right name
    result = getattr(importlib.import_module(mod), fn)(rest)
    return int(result or 0)


if __name__ == '__main__':
    sys.exit(main())
