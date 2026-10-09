"""
Step 11 — spatial density scaling, per observation type.

Calculates a spatial scaling factor for every profile, separately for each
observed variable (T, S), across all processed NC files in a run directory,
then writes the scaled weights back in-place.

Differs from the original MATLAB (update_weights_based_on_spatial_density.m,
reproduced in step11_spatial_scaling_combined_LEGACY.py), which counts every
profile together and applies one 1/n to both T and S. Here each variable is
counted on its own: a T-only profile (e.g. a Samoa mooring) raises the T count
in its bin but does not dilute the S weights of CTDs sharing that bin.

Algorithm:
  1. Pass 1 — for each variable V in SCALED_VARS, count how many profiles
     with at least one nonzero prof_Vweight land in each (bin, window) cell.
  2. Compute scaling_factor_V(bin, window) = 1 / count_V  (or 1.0 if empty).
  3. Pass 2 — for each file and each V it carries, look up each profile's
     factor, write prof_Vspatial_scaling_factor, and overwrite prof_Vweight
     with the scaled values. Files carrying only some of SCALED_VARS (e.g.
     T-only) are scaled for what they have.

A profile with no nonzero V weight is not counted for V and gets factor 1.0
(its V weights are all zero, so the factor has no effect).

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
    python step11_spatial_scaling.py /path/to/run_dir
    python step11_spatial_scaling.py /path/to/run_dir --window 10day
    python step11_spatial_scaling.py /path/to/run_dir --dry_run

Called from run_ncei.sh (on by default; -S to skip).
"""

import argparse
import datetime
import re
from pathlib import Path

import numpy as np
import xarray as xr

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

N_BINS = 10242          # fine geodesic resolution

# Observation types scaled independently. Each needs prof_<V>weight in a file
# to be counted/scaled there. (Samoa files also carry U/V; add them here if
# the team wants velocity weights scaled too.)
SCALED_VARS = ('T', 'S')

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
# Shared helpers
# ---------------------------------------------------------------------------

def _has_obs(weight):
    """(iPROF,) bool: profile has at least one nonzero weight for this var."""
    w = np.asarray(weight, dtype=np.float64)
    if w.ndim == 1:
        return w > 0
    return (w > 0).any(axis=1)


def _already_scaled(f):
    """True if f carries a scaling factor from this or the legacy step 11."""
    names = ['prof_spatial_scaling_factor'] + \
            [f'prof_{v}spatial_scaling_factor' for v in SCALED_VARS]
    try:
        with xr.open_dataset(f, mask_and_scale=False) as ds:
            return any(n in ds.variables for n in names)
    except Exception:
        return False


def _profile_keys(bin_ids, yyyymmdd, window_fn):
    """List of (bin_id, window_key) per profile, or None where invalid."""
    keys = []
    for bid, date in zip(bin_ids, yyyymmdd):
        if bid < 1 or bid > N_BINS or date <= 0:
            keys.append(None)
        else:
            keys.append((int(bid), window_fn(date)))
    return keys

# ---------------------------------------------------------------------------
# Pass 1 — count profiles per (bin, window), per variable
# ---------------------------------------------------------------------------

def _count_profiles(nc_files, window_fn):
    """Return {var: {(bin_id, window_key): count}} across all files."""
    counts = {v: {} for v in SCALED_VARS}
    for f in nc_files:
        try:
            ds = xr.open_dataset(f, mask_and_scale=False)
        except Exception as e:
            print(f"  [warn] cannot open {f.name}: {e}")
            continue
        if 'prof_bin_id_a' not in ds or 'prof_YYYYMMDD' not in ds:
            print(f"  [warn] missing bin id / date in {f.name}, skipping")
            ds.close()
            continue
        bin_ids  = np.asarray(ds['prof_bin_id_a'],  dtype=np.int32)
        yyyymmdd = np.asarray(ds['prof_YYYYMMDD'],  dtype=np.int64)
        has = {v: _has_obs(ds[f'prof_{v}weight'])
               for v in SCALED_VARS if f'prof_{v}weight' in ds}
        ds.close()

        keys = _profile_keys(bin_ids, yyyymmdd, window_fn)
        for v, mask in has.items():
            cv = counts[v]
            for key, ok in zip(keys, mask):
                if key is not None and ok:
                    cv[key] = cv.get(key, 0) + 1
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

        present = [v for v in SCALED_VARS if f'prof_{v}weight' in ds.data_vars]
        if 'prof_bin_id_a' not in ds or 'prof_YYYYMMDD' not in ds or not present:
            print(f"  [warn] {f.name} missing bin id / date / any of "
                  f"{[f'prof_{v}weight' for v in SCALED_VARS]}, skipping")
            ds.close()
            n_skipped += 1
            continue

        bin_ids  = np.asarray(ds['prof_bin_id_a'], dtype=np.int32)
        yyyymmdd = np.asarray(ds['prof_YYYYMMDD'], dtype=np.int64)
        n_prof   = len(bin_ids)
        keys     = _profile_keys(bin_ids, yyyymmdd, window_fn)

        ds_out = ds.copy()
        summary = []
        for v in present:
            wname = f'prof_{v}weight'
            has   = _has_obs(ds[wname])
            cv    = counts[v]

            sf = np.ones(n_prof, dtype=np.float64)
            for p in range(n_prof):
                if keys[p] is None or not has[p]:
                    continue
                c = cv.get(keys[p], 0)
                sf[p] = 1.0 / c if c > 0 else 1.0

            W = np.asarray(ds[wname], dtype=np.float64)
            # sf is (iPROF,); weights are (iPROF, iDEPTH) — broadcast
            sf_col = sf[:, np.newaxis] if W.ndim == 2 else sf
            ds_out[wname] = xr.DataArray(
                (W * sf_col).astype(ds[wname].dtype),
                dims=ds[wname].dims,
                attrs=ds[wname].attrs,
            )
            ds_out[f'prof_{v}spatial_scaling_factor'] = xr.DataArray(
                sf,
                dims=('iPROF',),
                attrs={
                    'long_name': f'spatial density scaling factor applied to {v} weights',
                    'units':     ' ',
                    '_FillValue': -9999.,
                    'missing_value': -9999.,
                },
            )
            sf_obs = sf[has]
            rng = (f"[{sf_obs.min():.4f}, {sf_obs.max():.4f}]"
                   if sf_obs.size else "n/a")
            summary.append(f"{v} sf {rng} ({int(has.sum())} w/ obs)")

        print(f"  {f.name}: {n_prof} profiles, " + "; ".join(summary)
              + (" [dry run]" if dry_run else ""))

        if dry_run:
            ds.close()
            ds_out.close()
            n_written += 1
            continue
        ds.close()

        # Derive output filename: replace __ncei_step_N with __ncei_step_11
        new_name = re.sub(r'__ncei_step_\d+', '__ncei_step_11', f.name)
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

    # Refuse to double-scale: step-11 files already carry scaled weights.
    already = [f for f in nc_files if _already_scaled(f)]
    if already:
        raise RuntimeError(
            f"{len(already)} file(s) already spatially scaled (e.g. {already[0].name}). "
            f"Run step 11 on an unscaled step-10 copy instead.")

    print(f"\n{'='*60}")
    print(f"Step 11 — spatial density scaling (per variable: {', '.join(SCALED_VARS)})")
    print(f"  run_dir    : {run_dir}")
    print(f"  window     : {window_mode}")
    print(f"  files      : {len(nc_files)}")
    print(f"  dry run    : {dry_run}")
    print(f"{'='*60}")

    print("\nPass 1: counting profiles per (bin, window), per variable...")
    counts = _count_profiles(nc_files, window_fn)
    for v in SCALED_VARS:
        cv = counts[v]
        print(f"  {v}: {len(cv):,} occupied (bin, window) cells, "
              f"{sum(cv.values()):,} profiles, "
              f"max count per cell = {max(cv.values()) if cv else 0}")

    print("\nPass 2: applying scaling factors...")
    n_written, n_skipped = _apply_scaling(nc_files, counts, window_fn, dry_run)

    print(f"\n{'='*60}")
    status = "DRY RUN — no files modified" if dry_run else f"{n_written} files written"
    print(f"Step 11 complete: {status}, {n_skipped} skipped")
    print(f"{'='*60}\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Apply per-variable spatial density scaling to processed profile NC files.')
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
