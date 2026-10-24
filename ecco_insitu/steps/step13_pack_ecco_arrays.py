"""
Step 13 — pack the yearly chain output into ECCO-ready profile files.

A thin wrapper around optimize_profile_files.py (which explains the problem
and does the work), with the settings this chain uses. The settings are in the
``pack`` section of config/defaults.yaml (and ``paths.grid_dir``); override them in
config/local.yaml, a ``--config`` file or ``--set pack.<key>=<value>``. Any of
them can also be given to optimize_profile_files.py directly on the command line.

Run (run_ncei.sh does this after step 11/12):
    python step13_pack_ecco_arrays.py <yearly_dir> <out_dir> [-n 4] [--dry_run] [--config FILE]
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config_loader  # noqa: E402
import optimize_profile_files as opt  # noqa: E402


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
        Settings from :func:`config_loader.load`; default: the standard files.
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


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Pack yearly chain output into ECCO-ready profile files.')
    ap.add_argument('yearly_dir')
    ap.add_argument('out_dir')
    ap.add_argument('-n', '--n_workers', type=int, default=None)
    ap.add_argument('--dry_run', action='store_true', help='Plan and print, write nothing')
    config_loader.add_arguments(ap)
    a = ap.parse_args()
    try:
        c = config_loader.from_args(a)
    except config_loader.ConfigError as e:
        sys.exit(f'config error: {e}')
    main(a.yearly_dir, a.out_dir, a.n_workers, a.dry_run, cfg=c)
