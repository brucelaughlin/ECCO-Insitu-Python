"""
LEGACY — superseded by step11_spatial_scaling.py (per-variable counts).
Kept for reference / reproducing runs made before Oct 2026, e.g.
profile_files_NCEI_processed_20261007_144529. Not called by run_ncei.sh.

This version counts ALL profiles in a (bin, window) together, regardless of
whether they carry T, S, or both, and applies the same 1/n to T and S weights
(matches the original MATLAB). It also skips any file lacking prof_Sweight
(e.g. T-only Samoa moorings), leaving those files unscaled at step 10.

Step 11 — spatial density scaling.

Calculates prof_spatial_scaling_factor for every profile across all processed
NC files in a run directory, then writes the scaled weights back in-place.

Algorithm (matches Ian's MATLAB update_weights_based_on_spatial_density.m):
  1. Pass 1 — read prof_bin_id_a + prof_YYYYMMDD from every file; count
     how many profiles land in each (bin, window) cell.
  2. Compute scaling_factor(bin, window) = 1 / count  (or 1.0 if empty).
  3. Pass 2 — for each file, look up each profile's scaling factor, write
     prof_spatial_scaling_factor, and overwrite prof_Tweight / prof_Sweight
     with the scaled values.

Window modes
------------
WINDOW_MODE = 'monthly'   (default, matches MATLAB)
    One window per calendar month.  Key = (year, month).

WINDOW_MODE = '10day'     (set via --window 10day)
    10-day windows anchored to day 1 of each year:
      window_index = floor((day_of_year - 1) / 10)
    Key = (year, window_index).  Approximately matches the team's stated
    Argo reference of "one profile per 3°×3° box per 10-day window".

Run standalone:
    python step11_spatial_scaling_combined_LEGACY.py /path/to/run_dir
    python step11_spatial_scaling_combined_LEGACY.py /path/to/run_dir --window 10day
    python step11_spatial_scaling_combined_LEGACY.py /path/to/run_dir --dry_run

Called from run_ncei.sh via the -s flag.
"""

import argparse
import datetime
import sys
from pathlib import Path

import numpy as np
import xarray as xr

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

N_BINS = 10242          # fine geodesic resolution

# ---------------------------------------------------------------------------
# Window helpers
# ---------------------------------------------------------------------------

def _window_monthly(yyyymmdd):
    """Return (year, month) window key for a YYYYMMDD integer."""
    y = int(yyyymmdd) // 10000
    m = (int(yyyymmdd) % 10000) // 100
    return (y, m)


def _day_of_year(y, m, d):
    return (datetime.date(y, m, d) - datetime.date(y, 1, 1)).days  # 0-based


def _window_10day(yyyymmdd):
    """Return (year, 10day_index) window key for a YYYYMMDD integer.

    10-day windows: 0=Jan1-10, 1=Jan11-20, ..., 35=Dec-last
    (last window of the year absorbs any remainder).
    """
    v = int(yyyymmdd)
    y = v // 10000
    m = (v % 10000) // 100
    d = v % 100
    doy = _day_of_year(y, m, d)
    return (y, doy // 10)


_WINDOW_FN = {
    'monthly': _window_monthly,
    '10day':   _window_10day,
}

# ---------------------------------------------------------------------------
# Pass 1 — count profiles per (bin, window)
# ---------------------------------------------------------------------------

def _count_profiles(nc_files, window_fn):
    """Return dict {(bin_id, window_key): count} across all files."""
    counts = {}
    for f in nc_files:
        try:
            ds = xr.open_dataset(f, mask_and_scale=False)
        except Exception as e:
            print(f"  [warn] cannot open {f.name}: {e}")
            continue
        if 'prof_bin_id_a' not in ds or 'prof_YYYYMMDD' not in ds:
            print(f"  [warn] missing required variables in {f.name}, skipping")
            ds.close()
            continue
        bin_ids  = np.asarray(ds['prof_bin_id_a'],  dtype=np.int32)
        yyyymmdd = np.asarray(ds['prof_YYYYMMDD'],  dtype=np.int64)
        ds.close()
        for bid, date in zip(bin_ids, yyyymmdd):
            if bid < 1 or bid > N_BINS or date <= 0:
                continue
            key = (int(bid), window_fn(date))
            counts[key] = counts.get(key, 0) + 1
    return counts


# ---------------------------------------------------------------------------
# Pass 2 — apply scaling factors
# ---------------------------------------------------------------------------

def _apply_scaling(nc_files, counts, window_fn, dry_run=False):
    n_written = 0
    n_skipped = 0
    for f in nc_files:
        try:
            ds = xr.open_dataset(f, mask_and_scale=False)
        except Exception as e:
            print(f"  [warn] cannot open {f.name}: {e}")
            n_skipped += 1
            continue

        required = {'prof_bin_id_a', 'prof_YYYYMMDD', 'prof_Tweight', 'prof_Sweight'}
        if not required.issubset(ds.data_vars):
            missing = required - set(ds.data_vars)
            print(f"  [warn] {f.name} missing {missing}, skipping")
            ds.close()
            n_skipped += 1
            continue

        bin_ids  = np.asarray(ds['prof_bin_id_a'], dtype=np.int32)
        yyyymmdd = np.asarray(ds['prof_YYYYMMDD'], dtype=np.int64)
        n_prof   = len(bin_ids)

        sf = np.ones(n_prof, dtype=np.float64)
        for p in range(n_prof):
            bid  = int(bin_ids[p])
            date = int(yyyymmdd[p])
            if bid < 1 or bid > N_BINS or date <= 0:
                continue
            key = (bid, window_fn(date))
            c = counts.get(key, 0)
            sf[p] = 1.0 / c if c > 0 else 1.0

        Tw = np.asarray(ds['prof_Tweight'], dtype=np.float64)
        Sw = np.asarray(ds['prof_Sweight'], dtype=np.float64)

        # Multiply: sf is (iPROF,); weights are (iPROF, iDEPTH) — broadcast
        sf_col = sf[:, np.newaxis] if Tw.ndim == 2 else sf
        Tw_scaled = Tw * sf_col
        Sw_scaled = Sw * sf_col

        print(f"  {f.name}: {n_prof} profiles, "
              f"sf range [{sf.min():.4f}, {sf.max():.4f}]"
              + (" [dry run]" if dry_run else ""))

        if dry_run:
            ds.close()
            n_written += 1
            continue

        # Build updated dataset — keep all existing variables, overwrite weights,
        # add prof_spatial_scaling_factor
        ds_out = ds.copy()
        ds_out['prof_Tweight'] = xr.DataArray(
            Tw_scaled.astype(ds['prof_Tweight'].dtype),
            dims=ds['prof_Tweight'].dims,
            attrs=ds['prof_Tweight'].attrs,
        )
        ds_out['prof_Sweight'] = xr.DataArray(
            Sw_scaled.astype(ds['prof_Sweight'].dtype),
            dims=ds['prof_Sweight'].dims,
            attrs=ds['prof_Sweight'].attrs,
        )
        ds_out['prof_spatial_scaling_factor'] = xr.DataArray(
            sf.astype(np.float64),
            dims=('iPROF',),
            attrs={
                'long_name': 'spatial density scaling factor applied to T and S weights',
                'units':     ' ',
                '_FillValue': -9999.,
                'missing_value': -9999.,
            },
        )
        ds.close()

        # Derive output filename: replace __ncei_step_N with __ncei_step_11
        import re as _re
        new_name = _re.sub(r'__ncei_step_\d+', '__ncei_step_11', f.name)
        new_path = f.parent / new_name

        tmp = f.with_suffix('.spatial_scaling_tmp.nc')
        try:
            ds_out.to_netcdf(tmp)
            tmp.replace(new_path)
            # Remove old file if the name changed
            if new_path != f and f.exists():
                f.unlink()
            n_written += 1
        except Exception as e:
            print(f"  [error] failed to write {f.name}: {e}")
            if tmp.exists():
                tmp.unlink()
            n_skipped += 1
        finally:
            ds_out.close()

    return n_written, n_skipped


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(run_dir, window_mode='monthly', dry_run=False):
    run_dir = Path(run_dir)
    if not run_dir.is_dir():
        raise ValueError(f"run_dir does not exist: {run_dir}")

    window_fn = _WINDOW_FN.get(window_mode)
    if window_fn is None:
        raise ValueError(f"Unknown window_mode '{window_mode}'. "
                         f"Choose from: {list(_WINDOW_FN)}")

    nc_files = sorted(run_dir.rglob('*.nc'))
    if not nc_files:
        print(f"No .nc files found under {run_dir}")
        return

    print(f"\n{'='*60}")
    print(f"Step 11 — spatial density scaling")
    print(f"  run_dir    : {run_dir}")
    print(f"  window     : {window_mode}")
    print(f"  files      : {len(nc_files)}")
    print(f"  dry run    : {dry_run}")
    print(f"{'='*60}")

    print("\nPass 1: counting profiles per (bin, window)...")
    counts = _count_profiles(nc_files, window_fn)
    total_bin_windows = len(counts)
    max_count = max(counts.values()) if counts else 0
    print(f"  {total_bin_windows:,} occupied (bin, window) cells, "
          f"max count per cell = {max_count}")

    print("\nPass 2: applying scaling factors...")
    n_written, n_skipped = _apply_scaling(nc_files, counts, window_fn, dry_run)

    print(f"\n{'='*60}")
    status = "DRY RUN — no files modified" if dry_run else f"{n_written} files written"
    print(f"Step 11 complete: {status}, {n_skipped} skipped")
    print(f"{'='*60}\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Apply spatial density scaling to processed profile NC files.')
    parser.add_argument('run_dir',
                        help='Root directory of a completed NCEI run (searched recursively for *.nc)')
    parser.add_argument('--window', default='monthly',
                        choices=list(_WINDOW_FN),
                        help='Window mode for counting profile density '
                             '(default: monthly; 10day for 10-day windows)')
    parser.add_argument('--dry_run', action='store_true',
                        help='Print what would be done without modifying any files')
    args = parser.parse_args()
    main(args.run_dir, window_mode=args.window, dry_run=args.dry_run)
