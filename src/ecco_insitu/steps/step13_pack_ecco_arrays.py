"""
Step 13 — pack the yearly chain output into ECCO-ready profile files.

A thin wrapper around :mod:`ecco_insitu.packing` (which explains the problem and does the
work), with the settings this chain uses. The settings are in the ``pack`` section of
``src/ecco_insitu/defaults.yaml`` (and ``paths.grid_dir``); override them in
``config/local.yaml``, a ``--config`` file or ``--set pack.<key>=<value>``. All of the packing
program's options are also available directly as ``ecco-insitu optimize``.

Run ('ecco-insitu run' does this after steps 11/12):
    ecco-insitu pack <yearly_dir> <out_dir> [-n 4] [--dry_run] [--config FILE]
"""

import argparse
import os
import sys

import ecco_insitu.config as config_loader  # noqa: E402
import ecco_insitu.packing as opt  # noqa: E402


def main(yearly_dir, out_dir, n_workers=None, dry_run=False, cfg=None):
    """Pack the files under ``yearly_dir`` into ``out_dir`` with the ``pack`` settings.

    Parameters
    ----------
    yearly_dir, out_dir : str or Path
        Input (chain output after step 11) and output folders.
    n_workers : int, optional
        Parallel workers; default ``run.n_workers``.
    dry_run : bool
        Plan and print, write nothing.
    cfg : dict, optional
        Settings from :func:`ecco_insitu.config.load`; default: the standard files.

    Returns
    -------
    list of dict
        The manifest of the files (see :func:`ecco_insitu.packing.run`).

    Raises
    ------
    ecco_insitu.config.ConfigError
        If ``paths.grid_dir`` is not set or doesn't exist.
    """
    if cfg is None:
        cfg = config_loader.load()
    config_loader.check_required_paths(cfg, ('paths.grid_dir',))
    pk = cfg['pack']
    if n_workers is None:
        n_workers = cfg['run']['n_workers']
    return opt.run(yearly_dir, None if dry_run else out_dir, llc=pk['llc'], tile=pk['tile'],
                   grid_dir=cfg['paths']['grid_dir'], position='interp', split_by_subdir=pk['split_by_subdir'],
                   year_rule=pk['year_rule'], nobsglob=pk['nobsglob'], tolerance=pk['tolerance'],
                   max_files=pk['max_files'], n_workers=n_workers)


def cli(argv=None):
    """Command-line entry point (``ecco-insitu pack``); see the module description."""
    ap = argparse.ArgumentParser(description='Pack yearly chain output into ECCO-ready profile files.')
    ap.add_argument('yearly_dir')
    ap.add_argument('out_dir')
    ap.add_argument('-n', '--n_workers', type=int, default=None)
    ap.add_argument('--dry_run', action='store_true', help='Plan and print, write nothing')
    config_loader.add_arguments(ap)
    a = ap.parse_args(argv)
    try:
        c = config_loader.from_args(a)
    except config_loader.ConfigError as e:
        sys.exit(f'config error: {e}')
    main(a.yearly_dir, a.out_dir, a.n_workers, a.dry_run, cfg=c)


if __name__ == '__main__':
    sys.exit(cli())
