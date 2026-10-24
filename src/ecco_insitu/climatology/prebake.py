"""
Pre-bake the climatology onto each observation depth grid found in an input folder (optional; makes
step 03 faster).

The input folder is scanned for all distinct ``prof_depth`` grids, and the climatology is
interpolated in depth onto each, once per climatology grid cell
(``clim_interp._depth_interp`` + ``clim_interp._fallback_fill``). That gives one
``<climatology stem>_<N>depths_prebaked.nc`` per grid. A new depth grid in new data just needs a
re-run.

At run time, step 03 then only blends in time and interpolates in lon/lat. The chain loads every
``*_prebaked.nc`` in ``paths.prebaked_clim_dir`` automatically, keyed by its depth grid, but only
those whose ``source_file`` attribute is the current ``paths.climatology_file`` (see
:func:`ecco_insitu.pipeline.check_climatologies`). Nothing needs to be registered by hand.

Usage::

    ecco-insitu climatology --prebake        # run.input_dir and paths.climatology_file from the settings

    # or directly, with explicit options:
    python -m ecco_insitu.climatology.prebake --input_dir /path/to/Interp_Profiles
    python -m ecco_insitu.climatology.prebake --source /path/to/fulldepth.nc --out_dir /path/to/out
    python -m ecco_insitu.climatology.prebake --min_levels 5     # skip grids with fewer levels
"""

import argparse
import collections
import datetime
import sys
from pathlib import Path

import numpy as np
import xarray as xr

import ecco_insitu.clim_interp as clim_interp

import ecco_insitu.config as config_loader

# Defaults from the chain settings (src/ecco_insitu/defaults.yaml + config/local.yaml)
_CFG = config_loader.load()
_DEFAULT_SOURCE = _CFG['paths']['climatology_file']
_DEFAULT_INPUT_DIR = _CFG['run']['input_dir']


def discover_grids(input_dir, min_levels=2):
    """Scan ``input_dir`` recursively for all distinct ``prof_depth`` grids.

    Parameters
    ----------
    input_dir : str or Path
        Folder of profile files.
    min_levels : int
        Grids with fewer levels are skipped (e.g. 1-level SOCAT files). Unreadable files are
        skipped too.

    Returns
    -------
    dict
        ``{depth tuple: [file names using it]}``.
    """
    grids = collections.defaultdict(list)
    for f in sorted(Path(input_dir).rglob("*.nc")):
        try:
            ds = xr.open_dataset(f)
            if 'prof_depth' in ds:
                depths = tuple(ds['prof_depth'].values.astype(float).tolist())
                if len(depths) >= min_levels:
                    grids[depths].append(str(f))
            ds.close()
        except Exception:
            pass
    return dict(grids)


def prebake_one_grid(obs_depths, T_all, S_all, clim_depths, lat, lon,
                     source_path, out_dir, stem):
    """Pre-bake T and S onto one observation depth grid and write one ``.nc`` file.

    Parameters
    ----------
    obs_depths : sequence of float
        The depth grid (meters).
    T_all, S_all : numpy.ndarray
        Climatology potential T and S, shape (12, depth, lat, lon).
    clim_depths, lat, lon : numpy.ndarray
        Climatology coordinates.
    source_path : Path
        The climatology file, recorded as the ``source_file`` attribute.
    out_dir : Path
        Output folder.
    stem : str
        Output name stem (``fulldepth`` in it is replaced by ``<N>depths_prebaked``).

    Returns
    -------
    Path
        The written file.
    """
    obs_depths = np.asarray(obs_depths, dtype=float)
    ndepth_obs = obs_depths.size
    nlat, nlon = lat.size, lon.size
    ncell = nlat * nlon

    tag = f"{ndepth_obs}depths"
    print(f"\nPre-baking onto {tag} ({ndepth_obs} levels, "
          f"{obs_depths[obs_depths > 0].min():.0f}–{obs_depths.max():.0f} m) ...")

    T_out = np.full((12, ndepth_obs, nlat, nlon), np.nan, dtype=np.float64)
    S_out = np.full((12, ndepth_obs, nlat, nlon), np.nan, dtype=np.float64)

    for m in range(12):
        print(f"  month {m+1:2d}/12", end='\r', flush=True)
        for field_all, out_arr in [(T_all, T_out), (S_all, S_out)]:
            # Reshape month slice (ndepth_clim, nlat, nlon) -> (ncell, ndepth_clim)
            clim_cols = field_all[m].reshape(field_all.shape[1], ncell).T.copy()

            # Same two operations that interpolate_climatology does at runtime
            interped = clim_interp._depth_interp(clim_cols, clim_depths, obs_depths)
            filled   = clim_interp._fallback_fill(interped, clim_cols, clim_depths, obs_depths)

            out_arr[m] = filled.T.reshape(ndepth_obs, nlat, nlon)

    print(f"  month 12/12 — done.")

    if 'fulldepth' in stem:
        out_stem = stem.replace('fulldepth', f'{tag}_prebaked')
    else:
        out_stem = f"{stem}_{tag}_prebaked"
    out_path = Path(out_dir) / f"{out_stem}.nc"

    out_ds = xr.Dataset(
        {
            'potential_T_monthly': (['month', 'obs_depth', 'lat', 'lon'], T_out),
            'S_monthly':           (['month', 'obs_depth', 'lat', 'lon'], S_out),
        },
        coords={
            'month':     np.arange(1, 13, dtype=np.int32),
            'obs_depth': obs_depths.astype(np.float64),
            'lat':       lat,
            'lon':       lon,
        },
        attrs={
            'title':           f'WOA23 pre-baked climatology on {ndepth_obs}-level obs depth grid',
            'source_file':     str(source_path),
            'obs_depth_grid':  tag,
            'build_date':      datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ'),
            'note': (
                'Depth interpolation and depth-fallback pre-applied. '
                'Runtime needs only time-blend + bilinear space interpolation.'
            ),
        },
    )
    enc = {
        'potential_T_monthly': {'zlib': True, 'complevel': 4, 'dtype': 'float64'},
        'S_monthly':           {'zlib': True, 'complevel': 4, 'dtype': 'float64'},
    }
    out_ds.to_netcdf(out_path, encoding=enc)
    print(f"  Written: {out_path}")
    return out_path


def prebake(source_path, input_dir, out_dir, min_levels=2):
    """Discover the depth grids in ``input_dir`` and pre-bake the climatology onto each.

    Parameters
    ----------
    source_path : str or Path
        The full-depth climatology (``potential_T_monthly``, ``S_monthly`` on month, depth,
        lat, lon).
    input_dir : str or Path
        Folder of profile files to take the depth grids from.
    out_dir : str or Path
        Output folder (created if needed).
    min_levels : int
        See :func:`discover_grids`.
    """
    source_path = Path(source_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Scanning {input_dir} for distinct depth grids (min_levels={min_levels}) ...")
    grids = discover_grids(input_dir, min_levels=min_levels)
    if not grids:
        print("No grids found — nothing to do.")
        return

    print(f"Found {len(grids)} distinct grid(s):")
    for depths, files in sorted(grids.items(), key=lambda x: -len(x[0])):
        arr = np.array(depths)
        print(f"  N={len(arr)}, {arr[arr>0].min():.0f}–{arr.max():.0f} m  "
              f"({len(files)} files, e.g. {Path(files[0]).name})")

    print(f"\nLoading {source_path} ...")
    ds = xr.open_dataset(source_path)
    T_all      = ds['potential_T_monthly'].values.astype(float)
    S_all      = ds['S_monthly'].values.astype(float)
    clim_depths = ds['depth'].values.astype(float)
    lat        = ds['lat'].values.astype(np.float32)
    lon        = ds['lon'].values.astype(np.float32)
    ds.close()

    stem = source_path.stem
    output_paths = {}

    for depths in grids:
        out_path = prebake_one_grid(
            depths, T_all, S_all, clim_depths, lat, lon,
            source_path, out_dir, stem)
        output_paths[depths] = out_path

    print("\n" + "="*60)
    print("Pre-bake complete.  To activate in ecco_insitu/pipeline.py, set:")
    print()
    print("    from prebake_woa23_climatology import discover_grids")
    print("    prebaked_clim_files = {")
    for depths, path in sorted(output_paths.items(), key=lambda x: -len(x[0])):
        arr = np.array(depths)
        comment = f"# {len(arr)}-level, {arr[arr>0].min():.0f}–{arr.max():.0f} m"
        print(f"        tuple({list(depths)}): '{path}',  {comment}")
    print("    }")
    print("="*60)


def cli(argv=None):
    """Command-line entry point (``python -m ecco_insitu.climatology.prebake``); ``ecco-insitu climatology --prebake`` runs :func:`prebake` with the settings' paths."""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--source', default=_DEFAULT_SOURCE,
                        help='Path to the WOA23 full-depth .nc climatology file')
    parser.add_argument('--input_dir', default=_DEFAULT_INPUT_DIR,
                        help='Input profile directory to scan for depth grids')
    parser.add_argument('--out_dir', default=None,
                        help='Output directory (default: same as source file)')
    parser.add_argument('--min_levels', type=int, default=2,
                        help='Skip grids with fewer than this many levels (default: 2)')
    args = parser.parse_args(argv)

    if not args.source or not args.input_dir:
        parser.error('give --source and --input_dir, or set paths.climatology_file and run.input_dir '
                     'in config/local.yaml')
    source = Path(args.source)
    if not source.exists():
        print(f"ERROR: source file not found: {source}", file=sys.stderr)
        sys.exit(1)
    out_dir = Path(args.out_dir) if args.out_dir else source.parent

    prebake(source, args.input_dir, out_dir, min_levels=args.min_levels)


if __name__ == '__main__':
    sys.exit(cli())
