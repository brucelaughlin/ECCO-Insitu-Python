"""
Pre-bake WOA23 climatology onto the two fixed observation depth grids.

All 389 Scripps input files share exactly two prof_depth grids:
  GRID_97 : 97 levels, 2–6000 m  (CTD, GLD, MEOP, PFL, PFL_BGC, XBT, ITP, MRB/97_depths)
  GRID_36 : 36 levels, 1–750 m   (MRB/36_depths only)

This script runs the vertical interpolation (clim_interp._depth_interp + _fallback_fill)
once per grid cell, producing two output .nc files.  At runtime, step03 can then skip the
per-profile vertical interpolation entirely and do only time-blend + lat/lon lookup.

Usage:
    python prebake_woa23_climatology.py
    python prebake_woa23_climatology.py --source /path/to/fulldepth.nc --out_dir /path/to/out

Output files (written alongside the source by default):
    woa23_decav91C0_TS_clim_potential_T_1deg_97depths_prebaked.nc
    woa23_decav91C0_TS_clim_potential_T_1deg_36depths_prebaked.nc
"""

import argparse
import datetime
import sys
from pathlib import Path

import numpy as np
import xarray as xr

import clim_interp

# ==============================================================================
# Canonical observation depth grids (float, in metres, read from actual files).
# Exported so controller scripts can build the {tuple(GRID): path} lookup dict.
# ==============================================================================

GRID_97 = [
    2.0, 4.0, 7.0, 10.0, 13.0, 16.0, 20.0, 25.0, 30.0, 35.0,
    40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0,
    90.0, 95.0, 100.0, 100.0, 110.0, 120.0, 120.0, 120.0, 130.0, 140.0,
    140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0,
    240.0, 250.0, 260.0, 270.0, 280.0, 290.0, 300.0, 320.0, 350.0, 380.0,
    400.0, 420.0, 450.0, 480.0, 500.0, 520.0, 550.0, 600.0, 650.0, 700.0,
    750.0, 800.0, 850.0, 900.0, 950.0, 1000.0, 1000.0, 1100.0, 1200.0, 1300.0,
    1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2200.0, 2400.0, 2600.0,
    2800.0, 3000.0, 3200.0, 3400.0, 3600.0, 3800.0, 4000.0, 4200.0, 4400.0, 4600.0,
    4800.0, 5000.0, 5200.0, 5400.0, 5600.0, 5800.0, 6000.0,
]

GRID_36 = [
    1.0, 2.0, 5.0, 10.0, 13.0, 20.0, 25.0, 28.0, 30.0, 40.0,
    45.0, 48.0, 50.0, 53.0, 60.0, 75.0, 80.0, 83.0, 100.0, 103.0,
    120.0, 123.0, 125.0, 140.0, 150.0, 153.0, 175.0, 180.0, 200.0, 203.0,
    225.0, 250.0, 300.0, 400.0, 500.0, 750.0,
]

_DEFAULT_SOURCE = (
    "/Users/brucel/ecco/yip/woa23_climatology/"
    "woa23_decav91C0_TS_clim_potential_T_1deg_fulldepth.nc"
)


def prebake(source_path, out_dir):
    """Pre-bake WOA23 onto GRID_97 and GRID_36, writing one .nc per grid."""
    source_path = Path(source_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading {source_path} ...")
    ds = xr.open_dataset(source_path)
    T_all = ds['potential_T_monthly'].values.astype(float)   # (12, 102, 180, 360)
    S_all = ds['S_monthly'].values.astype(float)
    clim_depths = ds['depth'].values.astype(float)           # (102,)
    lat = ds['lat'].values.astype(np.float32)
    lon = ds['lon'].values.astype(np.float32)
    ds.close()

    nlat, nlon = lat.size, lon.size   # 180, 360
    ncell = nlat * nlon               # 64800

    grids = [
        ('97depths', np.array(GRID_97, dtype=float)),
        ('36depths', np.array(GRID_36, dtype=float)),
    ]

    stem = source_path.stem  # e.g. woa23_decav91C0_TS_clim_potential_T_1deg_fulldepth

    for tag, obs_depths in grids:
        ndepth_obs = obs_depths.size
        print(f"\nPre-baking onto {tag} ({ndepth_obs} levels, "
              f"{obs_depths.min():.0f}–{obs_depths.max():.0f} m) ...")

        T_out = np.full((12, ndepth_obs, nlat, nlon), np.nan, dtype=np.float32)
        S_out = np.full((12, ndepth_obs, nlat, nlon), np.nan, dtype=np.float32)

        for m in range(12):
            print(f"  month {m+1:2d}/12", end='\r', flush=True)
            for field_all, out_arr in [(T_all, T_out), (S_all, S_out)]:
                # Reshape (102, 180, 360) -> (ncell, 102): each cell is a "profile"
                clim_cols = field_all[m].reshape(102, ncell).T.copy()    # (ncell, 102)

                # Vertical interpolation onto obs grid (same ops as interpolate_climatology)
                interped = clim_interp._depth_interp(clim_cols, clim_depths, obs_depths)
                # _fallback_fill needs the original full-depth clim_cols as 2nd arg
                filled = clim_interp._fallback_fill(interped, clim_cols, clim_depths, obs_depths)

                # Reshape back: (ncell, ndepth_obs) -> (ndepth_obs, nlat, nlon)
                out_arr[m] = filled.T.reshape(ndepth_obs, nlat, nlon)

        print(f"  month 12/12 — done.")

        # Replace the 'fulldepth' tag in the stem with the obs-grid tag
        if 'fulldepth' in stem:
            out_stem = stem.replace('fulldepth', f'{tag}_prebaked')
        else:
            out_stem = f"{stem}_{tag}_prebaked"
        out_path = out_dir / f"{out_stem}.nc"

        out_ds = xr.Dataset(
            {
                'potential_T_monthly': (
                    ['month', 'obs_depth', 'lat', 'lon'],
                    T_out,
                ),
                'S_monthly': (
                    ['month', 'obs_depth', 'lat', 'lon'],
                    S_out,
                ),
            },
            coords={
                'month':     np.arange(1, 13, dtype=np.int32),
                'obs_depth': obs_depths.astype(np.float32),
                'lat':       lat,
                'lon':       lon,
            },
            attrs={
                'title': f'WOA23 pre-baked climatology on {ndepth_obs}-level observation depth grid',
                'source_file': str(source_path),
                'obs_depth_grid': tag,
                'build_date': datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ'),
                'note': (
                    'Depth interpolation and depth-fallback pre-applied. '
                    'Runtime needs only time-blend + bilinear space interpolation.'
                ),
            },
        )
        enc = {
            'potential_T_monthly': {'zlib': True, 'complevel': 4, 'dtype': 'float32'},
            'S_monthly':           {'zlib': True, 'complevel': 4, 'dtype': 'float32'},
        }
        out_ds.to_netcdf(out_path, encoding=enc)
        print(f"  Written: {out_path}")

    print("\nDone.")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--source', default=_DEFAULT_SOURCE,
                        help='Path to the WOA23 full-depth .nc climatology file')
    parser.add_argument('--out_dir', default=None,
                        help='Output directory (default: same as source file)')
    args = parser.parse_args()

    source = Path(args.source)
    if not source.exists():
        print(f"ERROR: source file not found: {source}", file=sys.stderr)
        sys.exit(1)
    out_dir = Path(args.out_dir) if args.out_dir else source.parent

    prebake(source, out_dir)


if __name__ == '__main__':
    main()
