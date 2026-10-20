"""
Step 13 — pack the yearly chain output into ECCO-ready profile files.

A thin wrapper around optimize_profile_files.py (which explains the problem
and does the work), with the settings this chain uses.  Edit them here; any of
them can also be given to optimize_profile_files.py directly on the command line.

Run (run_ncei.sh does this after step 11/12):
    python step13_pack_ecco_arrays.py <yearly_dir> <out_dir> [-n 4] [--dry_run]
"""

import argparse

import optimize_profile_files as opt

# ---------------------------------------------------------------------------
# Parameters (edit here)
# ---------------------------------------------------------------------------

LLC = 90                      # model grid: LLC90 (ECCO v4r6)
TILE = 30                     # per-rank tile size, 30x30 cells (v4r6 SIZE.h sNx = sNy = 30)
GRID_DIR = opt.DEFAULTS['grid_dir']
YEAR_RULE = 'whole_years'     # files labeled in whole years; a too-dense year is split into parts
NOBSGLOB = 'auto'             # or an integer
TOLERANCE = 0.05              # auto: fewest files within 5% of the smallest NFILESPROFMAX x NOBSGLOB
MAX_FILES = None              # auto: optional cap on the number of files
SPLIT_BY_SUBDIR = ['ITP']     # ITP_L2 / ITP_L3 separate; everything else one provenance per source


def main(yearly_dir, out_dir, n_workers=None, dry_run=False):
    return opt.run(yearly_dir, None if dry_run else out_dir, llc=LLC, tile=TILE, grid_dir=GRID_DIR,
                   position='interp', split_by_subdir=SPLIT_BY_SUBDIR, year_rule=YEAR_RULE,
                   nobsglob=NOBSGLOB, tolerance=TOLERANCE, max_files=MAX_FILES, n_workers=n_workers)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Pack yearly chain output into ECCO-ready profile files.')
    ap.add_argument('yearly_dir')
    ap.add_argument('out_dir')
    ap.add_argument('-n', '--n_workers', type=int, default=None)
    ap.add_argument('--dry_run', action='store_true', help='Plan and print, write nothing')
    a = ap.parse_args()
    main(a.yearly_dir, a.out_dir, a.n_workers, a.dry_run)
