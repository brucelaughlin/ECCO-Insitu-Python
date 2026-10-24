#!/usr/bin/env python3
"""
Validate the ECCO-ready files written by ecco_insitu/packing.py (step 13)
against the yearly files they were packed from.

Checks
  - one provenance per file; profiles in chronological order; one depth grid
  - no file puts more than the manifest's NOBSGLOB cap in any counted model
    tile (same grid, tile size and positions as the manifest)
  - year rule: whole_years files never mix part of a year with other years, and
    their labels are whole years; one_year files stay within one year
  - every yearly profile appears exactly once, and its data (every per-profile
    variable) is identical to the yearly file
  - the manifest's NOBSGLOB / NFILESPROFMAX match the files

Usage::

    ecco-insitu check-arrays <run_dir>     # <run_dir>/yearly, /ecco_arrays, /documentation
    ecco-insitu check-arrays --yearly <dir> --arrays <dir> [--docs <dir>]

The manifest is read from the documentation folder (default: beside the arrays
folder); older runs that kept it inside the arrays folder still work.
Exit code 0 if no errors.
"""

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import netCDF4
import numpy as np

import ecco_insitu.packing as opt  # noqa: E402


def cli(argv=None):
    """Command-line entry point (``ecco-insitu check-arrays``): validate a run's step 13 files.

    Reads the yearly files, the ECCO-ready files and ``ecco_arrays_manifest.json``, and checks that
    every yearly profile appears exactly once with identical data, that each file holds one
    provenance in time order with whole-year labels, that no tile exceeds ``NOBSGLOB``, and that
    the ``PROFILES_SIZE.h`` values match. Prints the errors and ``RESULT: PASS`` or ``FAIL``.

    The process exits with status 0 on PASS and 1 on FAIL (``sys.exit``).

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments (default: ``sys.argv[1:]``).
    """
    ap = argparse.ArgumentParser(description='Validate ECCO-ready profile files.')
    ap.add_argument('run_dir', nargs='?', default=None)
    ap.add_argument('--yearly', default=None)
    ap.add_argument('--arrays', default=None)
    ap.add_argument('--docs', default=None, help="folder with ecco_arrays_manifest.json (default: 'documentation' beside --arrays)")
    a = ap.parse_args(argv)
    yearly = Path(a.yearly) if a.yearly else Path(a.run_dir) / 'yearly'
    arrays = Path(a.arrays) if a.arrays else Path(a.run_dir) / 'ecco_arrays'
    errors = []
    err = errors.append

    files = sorted(arrays.glob('*.nc'))
    print(f'{arrays}: {len(files)} files')
    docs = Path(a.docs) if a.docs else arrays.parent / 'documentation'
    mpath = next((p for p in (docs / 'ecco_arrays_manifest.json', arrays / 'ecco_arrays_manifest.json') if p.exists()),
                 docs / 'ecco_arrays_manifest.json')
    others = [p.name for p in arrays.iterdir() if p.suffix != '.nc']
    if others:
        err(f'{arrays.name}/ holds non-profile files: {others} (summary files belong in {docs.name}/)')
    if not files or not mpath.exists():
        print('RESULT: FAIL (no files or no ecco_arrays_manifest.json)')
        sys.exit(1)
    meta = json.loads(mpath.read_text())
    cap = meta['nobsglob_cap']
    grid = opt.Grid(meta['grid_dir'], meta['llc'], meta['tile'])
    lon_v, lat_v = (('prof_interp_lon', 'prof_interp_lat') if meta['position'] == 'interp'
                    else ('prof_lon', 'prof_lat'))
    rule, split_sub = meta['year_rule'], set(meta['split_by_subdir'])
    years_used = {}                                    # (provenance, year) -> list of files holding it
    yearly_files = {f.name: f for f in yearly.rglob('*.nc')}
    seen = defaultdict(set)
    depths, busiest_all, n_packed = set(), 0, 0
    for f in files:
        with netCDF4.Dataset(f) as nc:
            nc.set_auto_mask(False)
            v = nc.variables
            n = nc.dimensions['iPROF'].size
            n_packed += n
            depths.add(tuple(np.round(v['prof_depth'][:], 3)))
            src_names = json.loads(nc.getncattr('yearly_files'))
            yf, yi = v['prof_yearly_file'][:], v['prof_yearly_index'][:]
            t = v['prof_YYYYMMDD'][:].astype(np.int64) * 1_000_000 + v['prof_HHMMSS'][:].astype(np.int64)
            if (np.diff(t) < 0).any():
                err(f'{f.name}: profiles not in chronological order')
            provs = {opt.provenance(yearly_files[nm], yearly, split_sub) for nm in src_names}
            if len(provs) != 1:
                err(f'{f.name}: mixes provenances {sorted(provs)}')
            tl = grid.tiles_of(v[lon_v][:].ravel(), v[lat_v][:].ravel())
            yrs = sorted(set((t // 10_000_000_000).tolist()))
            prov = sorted(provs)[0]
            for y in yrs:
                years_used.setdefault((prov, y), []).append((f.name, len(yrs)))
            t0, t1 = nc.getncattr('time_start'), nc.getncattr('time_end')
            if rule == 'whole_years' and not (t0.endswith('0101') and t1.endswith('1231')
                                              and int(t0[:4]) <= yrs[0] and int(t1[:4]) >= yrs[-1]):
                err(f'{f.name}: label {t0}-{t1} is not whole years covering its data {yrs[0]}-{yrs[-1]}')
            if rule == 'one_year' and len(yrs) > 1:
                err(f'{f.name}: spans years {yrs[0]}-{yrs[-1]}')
            busiest = int(np.bincount(tl).max())
            busiest_all = max(busiest_all, busiest)
            if busiest > cap:
                err(f'{f.name}: a tile holds {busiest:,} profiles > NOBSGLOB cap {cap:,}')
            sources = {}
            for j, name in enumerate(src_names):
                sel = np.flatnonzero(yf == j)
                idx = yi[sel]
                dup = seen[name] & set(idx.tolist())
                if dup:
                    err(f'{f.name}: {len(dup)} profiles of {name} already packed elsewhere')
                seen[name] |= set(idx.tolist())
                src = netCDF4.Dataset(yearly_files[name]); src.set_auto_mask(False)
                sources[name] = (src, sel, idx)
            # each packed variable is read once, then compared with every contributing yearly file
            for var in v:
                if 'iPROF' not in v[var].dimensions or var.startswith('prof_yearly_'):
                    continue
                ax = v[var].dimensions.index('iPROF')
                packed = v[var][:]
                for name, (src, sel, idx) in sources.items():
                    if var not in src.variables:
                        continue
                    want = np.take(src.variables[var][:], idx, axis=ax)
                    got = np.take(packed, sel, axis=ax)
                    got = got[tuple(slice(0, L) for L in want.shape)]
                    same = (np.array_equal(want, got, equal_nan=True) if want.dtype.kind == 'f'
                            else np.array_equal(want, got))
                    if not same:
                        err(f'{f.name}: {var} differs from {name}')
                del packed
            for name, (src, _, _) in sources.items():
                for var in src.variables:            # every yearly variable must be in the packed file
                    if 'iPROF' in src.variables[var].dimensions and var not in v:
                        err(f'{f.name}: missing variable {var} (present in {name})')
                src.close()
        print(f'  {f.name}: {n:,} profiles, busiest tile {busiest:,}', flush=True)

    if len(depths) != 1:
        err(f'files have {len(depths)} different depth grids')
    if rule == 'whole_years':                          # a split year must not share a file with other years
        for (p, y), holders in years_used.items():
            if len(holders) > 1 and any(n_years > 1 for _, n_years in holders):
                err(f'{p} {y} is split, but part of it shares a file with other years: {[h for h, _ in holders]}')
    n_yearly = 0
    for name, f in yearly_files.items():
        with netCDF4.Dataset(f) as nc:
            m = nc.dimensions['iPROF'].size
        n_yearly += m
        if len(seen[name]) != m:
            err(f'{name}: {m - len(seen[name]):,} of {m:,} profiles not packed')
    if meta['nfilesprofmax_needed'] != len(files):
        err(f"manifest NFILESPROFMAX {meta['nfilesprofmax_needed']} != {len(files)} files")
    if meta['nobsglob_needed'] != busiest_all:
        err(f"manifest NOBSGLOB {meta['nobsglob_needed']:,} != busiest tile found {busiest_all:,}")
    print(f'\nyearly profiles {n_yearly:,}; packed {n_packed:,}; NFILESPROFMAX {len(files)} x NOBSGLOB '
          f'{busiest_all:,} = {len(files) * busiest_all:,}; depth grids {len(depths)}')
    print(f'ERRORS ({len(errors)}):')
    for e in errors[:50]:
        print('  ' + e)
    print('RESULT:', 'PASS' if not errors else 'FAIL')
    sys.exit(1 if errors else 0)


if __name__ == '__main__':
    sys.exit(cli())
