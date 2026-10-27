"""
**Work in progress** (see :mod:`ecco_insitu.model_comparison`). Approximate model equivalents from
monthly-mean fields: ``prof_Testim`` / ``prof_Sestim`` without a model run on the profiles.

For every profile file, the model's monthly-mean THETA and SALT (MITgcm MDS ``.data``/``.meta``, LLC90,
compact layout) are sampled at the profile's grid cell, month and depths, and written to
``<name>_model.nc`` in the same layout ``ecco-insitu mds-to-nc`` writes (so ``model-cost`` and
``profile-analyze --model-equivalents`` take either).

How each value is made:

- **Cell:** the LLC90 cell the chain put the profile in (``prof_interp_lon`` / ``prof_interp_lat``, the
  cell center; one cell per profile, as ``prof_interp_weights`` = 1).
- **Month:** the monthly mean of the profile's calendar month (``prof_YYYYMMDD``). A record whose
  time step is ``s`` (``deltaT`` = 3600 s, model start ``--start``, default 1992-01-01) ends at
  ``start + s hours``; it is the mean of the month before that.
- **Depth:** as the model's ``profiles_interp.F``. Linear between the two model levels around the depth
  (level centers ``-RC``), when the deeper one is wet; the shallowest level's value above it; the deepest
  level's value down to half a level below it; otherwise the upper level's value if the depth is above
  that cell's bottom (``-RF``), else none. A dry cell (``hFacC`` = 0) has no value.
- **Missing:** -9999 (no model value), as in ``mds-to-nc``'s output.

**It is approximate.** The model samples its own state at the profile's time (daily), not a monthly mean,
so these values differ from a run's ``.equi`` output; the output says so in its attributes. Use it for
cost analyses (``profile-analyze``, ``model-cost``), not where exact model equivalents are needed.

Usage::

    ecco-insitu monthly-estim --theta <dir of THETA_mon_mean.*> --salt <dir of SALT_mon_mean.*>
                              --data <dir of profile .nc> --out <dir> [--grid <LLC90 grid dir>] [-n 6]

``--grid`` defaults to ``paths.grid_dir`` in the settings.
"""
import argparse
import datetime
import glob
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor

import netCDF4
import numpy as np

FILL = -9999.0
NX = 90
NR = 50
_VARS = {'T': 'theta', 'S': 'salt'}


# --------------------------------------------------------------------------------------------------
def read_grid(grid_dir):
    """LLC90 cell centers (XC, YC; compact, flat), level centers and faces (RC, RF, positive down), and the
    wet mask (hFacC > 0) per level, flat over the compact layout."""
    g = lambda n: np.fromfile(os.path.join(grid_dir, n + '.data'), '>f4')
    xc, yc, rc, rf, h = g('XC'), g('YC'), -g('RC'), -g('RF'), g('hFacC')
    n2 = xc.size
    return {'xc': xc, 'yc': yc, 'zc': rc.astype(float), 'zf': rf.astype(float),
            'wet': (h.reshape(NR, n2) > 0)}


def month_records(field_dir, prefix, start):
    """{yyyymm: path} for the monthly-mean files ``<prefix>.<step>.data`` in ``field_dir``."""
    out = {}
    for p in glob.glob(os.path.join(field_dir, f'{prefix}.*.data')):
        step = int(re.search(r'\.(\d+)\.data$', p).group(1))
        end = start + datetime.timedelta(hours=step)
        mid = end - datetime.timedelta(days=14)          # inside the averaged month
        out[mid.year * 100 + mid.month] = p
    return out


def cell_index(grid, lon, lat):
    """Flat compact-grid index of each profile's cell, from the cell center the chain stored."""
    from scipy.spatial import cKDTree
    def xyz(lo, la):
        lo, la = np.radians(lo), np.radians(la)
        return np.c_[np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
    d, i = cKDTree(xyz(grid['xc'], grid['yc'])).query(xyz(lon, lat))
    return i, d


def vertical_weights(zc, zf, depth):
    """For each profile depth: (k0, k1, a, rule) so that value = (1 - a) * F[k0] + a * F[k1], following
    profiles_interp.F. rule 0 = above the top center, 1 = between centers, 2 = below the deepest center
    (within half a level), -1 = deeper than that (no value)."""
    k0 = np.zeros(len(depth), int); k1 = np.zeros(len(depth), int); a = np.zeros(len(depth)); rule = np.full(len(depth), -1)
    nr = len(zc)
    half = 0.5 * (zf[nr] - zf[nr - 1])
    for j, z in enumerate(depth):
        if z < zc[0]:
            rule[j] = 0
        elif z >= zc[-1]:
            if z < zc[-1] + half:
                k0[j] = k1[j] = nr - 1; rule[j] = 0
        else:
            k = int(np.searchsorted(zc, z, side='right') - 1)       # zc[k] <= z < zc[k+1]
            k0[j], k1[j] = k, k + 1; a[j] = (z - zc[k]) / (zc[k + 1] - zc[k]); rule[j] = 1
    return k0, k1, a, rule


def _estimate_file(args):
    path, out_dir, theta, salt, grid_dir, start = args
    grid = read_grid(grid_dir)
    T_rec, S_rec = month_records(theta, 'THETA_mon_mean', start), month_records(salt, 'SALT_mon_mean', start)
    name = os.path.basename(path)[:-3]
    with netCDF4.Dataset(path) as ds:
        ds.set_auto_mask(False)
        depth = np.asarray(ds['prof_depth'][:], float)
        lon = np.asarray(ds['prof_interp_lon'][:], float).reshape(-1, ds['prof_interp_lon'].shape[-1] if ds['prof_interp_lon'].ndim > 1 else 1)[:, 0]
        lat = np.asarray(ds['prof_interp_lat'][:], float).reshape(len(lon), -1)[:, 0]
        ym = (np.asarray(ds['prof_YYYYMMDD'][:], float) // 100).astype(int)
        has = {V: f'prof_{V}' in ds.variables for V in _VARS}
    cells, dist = cell_index(grid, lon, lat)
    k0, k1, a, rule = vertical_weights(grid['zc'], grid['zf'], depth)
    # profiles_interp.F: between centers with the deeper cell dry, use the upper value if above its bottom
    upper_ok = (rule == 1) & (depth < grid['zf'][k1])
    est = {V: np.full((len(lon), len(depth)), FILL) for V in _VARS if has[V]}
    n2 = grid['xc'].size
    missing_month = 0
    for m in np.unique(ym):
        sel = np.nonzero(ym == m)[0]
        for V in est:
            rec = (T_rec if V == 'T' else S_rec).get(int(m))
            if rec is None:
                missing_month += len(sel) if V == 'T' else 0
                continue
            F = np.fromfile(rec, '>f4').reshape(NR, n2).astype(float)
            c = cells[sel]
            f0, f1 = F[k0][:, c].T, F[k1][:, c].T
            w0, w1 = grid['wet'][k0][:, c].T, grid['wet'][k1][:, c].T
            v = np.full(f0.shape, np.nan)
            t = rule == 0
            v[:, t] = np.where(w0[:, t], f0[:, t], np.nan)
            b = rule == 1
            v[:, b] = np.where(w1[:, b], (1 - a[b]) * f0[:, b] + a[b] * f1[:, b],
                               np.where(w0[:, b] & upper_ok[b], f0[:, b], np.nan))
            est[V][sel] = np.where(np.isfinite(v), v, FILL)
    dst = os.path.join(out_dir, name + '_model.nc'); tmp = dst + '.part'
    with netCDF4.Dataset(path) as src, netCDF4.Dataset(tmp, 'w') as o:
        o.createDimension('iPROF', len(lon)); o.createDimension('iDEPTH', len(depth))
        for v in ('prof_depth', 'prof_YYYYMMDD', 'prof_HHMMSS', 'prof_lon', 'prof_lat', 'prof_descr', 'prof_bin_id_a'):
            if v in src.variables:
                sv = src[v]
                for d in sv.dimensions:
                    if d not in o.dimensions:
                        o.createDimension(d, len(src.dimensions[d]))
                ov = o.createVariable(v, sv.dtype, sv.dimensions, zlib=True)
                ov[:] = sv[:]
        for V, arr in est.items():
            ov = o.createVariable(f'prof_{V}estim', 'f8', ('iPROF', 'iDEPTH'), fill_value=FILL, zlib=True, complevel=4)
            ov[:] = arr
            ov.long_name = f'approximate model {_VARS[V]} at the profile (monthly mean; see global attributes)'
        o.setncattr('model_equivalents', 'APPROXIMATE: sampled from monthly-mean fields, not a model run on these '
                    'profiles (ecco_insitu.model_comparison.monthly_estim)')
        o.setncattr('monthly_mean_theta_dir', os.path.abspath(theta)); o.setncattr('monthly_mean_salt_dir', os.path.abspath(salt))
        o.setncattr('source_file', os.path.abspath(path))
        o.setncattr('horizontal', 'the LLC90 cell of prof_interp_lon/lat (nearest neighbor, as the chain)')
        o.setncattr('vertical', "as MITgcm profiles_interp.F (linear between level centers; dry cells have no value)")
        o.setncattr('time', f"the profile's calendar month; model start {start:%Y-%m-%d}, deltaT 3600 s")
    os.replace(tmp, dst)
    n = {V: int((arr != FILL).sum()) for V, arr in est.items()}
    far = int((np.degrees(dist) > 0.01).sum())
    return (f'{name}: {len(lon):,} profiles; model values T {n.get("T", 0):,}, S {n.get("S", 0):,}; '
            f'profiles in months with no field: {missing_month:,}; cell mismatch > 0.01 deg: {far:,}')


def run(data_dir, out_dir, theta, salt, grid_dir, start=datetime.datetime(1992, 1, 1), n_workers=4):
    """Write ``<out_dir>/<name>_model.nc`` for every ``*.nc`` in ``data_dir``; returns the summary lines."""
    os.makedirs(out_dir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(data_dir, '*.nc')))
    with ProcessPoolExecutor(n_workers) as ex:
        out = []
        for line in ex.map(_estimate_file, [(f, out_dir, theta, salt, grid_dir, start) for f in files]):
            print(line, flush=True); out.append(line)
    return out


def cli(argv=None):
    """Command-line entry point (``ecco-insitu monthly-estim``); see the module description."""
    import ecco_insitu.config as config_loader
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[1],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument('--theta', required=True, help='folder of THETA_mon_mean.<step>.data/.meta')
    ap.add_argument('--salt', required=True, help='folder of SALT_mon_mean.<step>.data/.meta')
    ap.add_argument('--data', required=True, help='folder of profile files (e.g. a run\'s ecco_arrays/)')
    ap.add_argument('--out', required=True, help='output folder for <name>_model.nc')
    ap.add_argument('--grid', help='LLC90 grid folder (XC, YC, RC, RF, hFacC); default paths.grid_dir')
    ap.add_argument('--start', default='1992-01-01', help='model start date (default 1992-01-01, ECCO v4)')
    ap.add_argument('-n', type=int, default=4, help='parallel workers (default 4)')
    config_loader.add_arguments(ap)
    a = ap.parse_args(argv)
    grid = a.grid or config_loader.from_args(a)['paths']['grid_dir']
    run(a.data, a.out, a.theta, a.salt, grid, datetime.datetime.strptime(a.start, '%Y-%m-%d'), a.n)
    return 0


if __name__ == '__main__':
    sys.exit(cli())
