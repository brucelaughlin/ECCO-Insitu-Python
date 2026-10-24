"""
**Work in progress** (see :mod:`ecco_insitu.model_comparison`). Reproduces the v4r6 iter47
per-type cost bars; not yet checked against a model run's own cost output.

Offline profile cost, computed as ECCO v4r6's cost_profiles.F does it (ECCO-GROUP/ECCO-v4-Configurations,
"ECCOv4 Release 6/code/cost_profiles.F", with ALLOW_PROFILES_SAMPLESPLIT_COST and
ALLOW_PROFILES_CLIMMASK defined, as in that release's PROFILES_OPTIONS.h).

Per variable V (T, S, ...), with  d = prof_V,  m = model equivalent (prof_Vestim),
w = prof_Vweight, and the averaging bin b = prof_bin_id_a (truncated to an integer)::

  CLIMMASK      w = 0 where prof_Vclim <= -990                       (profiles_readvector.F)
  model mask    mask = 1 where m != 0, else 0                         (no model value is 0)
  bin means     over all files, per (bin, depth value), using points with w > 0 and mask = 1:
                  m_mean = mean(m),  d_mean = mean(d),  s = mean(1/w),  n = count
                  w_mean = 1 / max(s / n, s / profiles_mean_indsamples)        (default 9)
                Depths are matched across files by exact value (prof_depth_comb in the model).
  anomaly cost  sum over points with w > 0 and d_mean != 0 of
                  w * mask * (m - d - m_mean + d_mean)^2,     count = sum of mask
                (reported per file and variable, as 'cost_profiles(file,var)')
  mean cost     sum over (bin, depth) with w_mean > 0, d_mean != 0 and m_mean != 0 of
                  w_mean * (m_mean - d_mean)^2,               count = number of terms
                (reported per variable, as 'cost_profiles_mean(var)')
  total         sum of mult_profiles * anomaly cost + mult_profiles_mean * mean cost
                (all multipliers are 1 in v4r6's data.profiles)

Model values are read from the data file if it has prof_Vestim, else from <name>_model.nc in
--estim (e.g. made by ``ecco-insitu mds-to-nc``, or MITprof_gcm2nc.m).
-9999 in prof_Vestim means "no model value" and is treated as 0, as in the model's binary output.

Usage::

    ecco-insitu model-cost --data <dir with profile .nc> [--estim <dir with *_model.nc>]
                           [--vars T S] [--indsamples 9] [--out <dir>] [-n 8]

"""
import argparse
import sys
import glob
import json
import os
from concurrent.futures import ProcessPoolExecutor

import netCDF4
import numpy as np

CHUNK = 20000            # profiles read at a time
NAVGBIN = 10246          # bin ids 1..10245 (NAVGBINMAX in v4r6's PROFILES_SIZE.h); 0 unused


def estim_path(path, estim_dir):
    """File holding the model values for data file ``path``: ``path`` itself if ``estim_dir`` is
    None, else ``<estim_dir>/<name>_model.nc``."""
    if estim_dir is None:
        return path
    name = os.path.basename(path)[:-3]
    return os.path.join(estim_dir, name + '_model.nc')


def read_chunk(path, epath, V, i0, i1):
    """d, m, w (after CLIMMASK) for profiles i0:i1, and bin ids (int, -1 where unusable)."""
    with netCDF4.Dataset(path) as ds:
        ds.set_auto_mask(False)
        d = np.asarray(ds[f'prof_{V}'][i0:i1], float)
        w = np.asarray(ds[f'prof_{V}weight'][i0:i1], float)
        if f'prof_{V}clim' in ds.variables:
            c = np.asarray(ds[f'prof_{V}clim'][i0:i1], float)
            w = np.where(c <= -990, 0.0, w)
        b = np.asarray(ds['prof_bin_id_a'][i0:i1], float)
    with netCDF4.Dataset(epath) as ds:
        ds.set_auto_mask(False)
        m = np.asarray(ds[f'prof_{V}estim'][i0:i1], float)
    m = np.where(m == -9999, 0.0, m)
    with np.errstate(invalid='ignore'):
        bi = np.where(np.isfinite(b) & (b >= 1) & (b < NAVGBIN), np.trunc(np.nan_to_num(b)), -1).astype(np.int64)
    return d, m, w, bi


def pass1(args):
    """Bin sums for one chunk: sum m, sum d, sum 1/w, count, per flat (bin, comb level)."""
    path, epath, V, i0, i1, lev, nlev = args
    d, m, w, bi = read_chunk(path, epath, V, i0, i1)
    with np.errstate(invalid='ignore'):
        use = (w > 0) & (m != 0) & (bi[:, None] >= 1)
    fin = np.isfinite(d) & np.isfinite(m) & np.isfinite(w)
    ok = use & fin
    flat = (bi[:, None] * nlev + lev[None, :])[ok]
    size = NAVGBIN * nlev
    bad = int((use & ~fin).sum())
    return (np.bincount(flat, weights=m[ok], minlength=size),
            np.bincount(flat, weights=d[ok], minlength=size),
            np.bincount(flat, weights=1.0 / w[ok], minlength=size),
            np.bincount(flat, minlength=size).astype(float), bad)


def pass2(args):
    """Anomaly cost and count for one chunk, given the bin means."""
    path, epath, V, i0, i1, lev, nlev, mmean, dmean = args
    d, m, w, bi = read_chunk(path, epath, V, i0, i1)
    flat = np.where(bi[:, None] >= 1, bi[:, None] * nlev + lev[None, :], 0)
    mm, dm = mmean[flat], dmean[flat]            # flat 0 is bin 0, whose means are 0
    mask = (m != 0).astype(float)
    with np.errstate(invalid='ignore'):
        sel = (w > 0) & (dm != 0) & np.isfinite(d) & np.isfinite(m) & np.isfinite(w)
    r = (m - d - mm + dm)[sel]
    return float(np.sum(w[sel] * mask[sel] * r * r)), float(mask[sel].sum())


def cli(argv=None):
    """Command-line entry point (``ecco-insitu`` subcommand)."""
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--data', required=True, help='folder with the profile files (searched recursively)')
    p.add_argument('--estim', help='folder with <name>_model.nc files holding prof_Vestim (if the data files lack it)')
    p.add_argument('--vars', nargs='+', default=['T', 'S'])
    p.add_argument('--indsamples', type=float, default=9.0, help='profiles_mean_indsamples (v4r6: 9)')
    p.add_argument('--out', help='folder for cost_summary.json (default: print only)')
    p.add_argument('-n', type=int, default=4, help='parallel workers (default 4)')
    a = p.parse_args(argv)

    files = sorted(f for f in glob.glob(os.path.join(a.data, '**', '*.nc'), recursive=True)
                   if not f.endswith('_model.nc'))
    # per file: profile count, depth levels, and which variables it has
    info = {}
    for f in files:
        ep = estim_path(f, a.estim)
        if not os.path.exists(ep):
            print(f'{os.path.basename(f)}: no model file {ep}, skipped'); continue
        with netCDF4.Dataset(f) as ds, netCDF4.Dataset(ep) as es:
            vs = [V for V in a.vars if f'prof_{V}' in ds.variables and f'prof_{V}weight' in ds.variables
                  and f'prof_{V}estim' in es.variables]
            if 'prof_bin_id_a' not in ds.variables:
                print(f'{os.path.basename(f)}: no prof_bin_id_a, skipped'); continue
            info[f] = (ep, len(ds.dimensions['iPROF']), np.asarray(ds['prof_depth'][:], float), vs)
    comb = np.unique(np.concatenate([v[2] for v in info.values()]))      # combined depths, exact values
    comb = comb[comb >= 0]
    nlev = len(comb)
    lev_of = {f: np.searchsorted(comb, v[2]) for f, v in info.items()}

    result = {'settings': {'indsamples': a.indsamples, 'files': len(info), 'combined_levels': nlev}}
    total = 0.0
    with ProcessPoolExecutor(a.n) as ex:
        for V in a.vars:
            fl = [f for f in info if V in info[f][3]]
            jobs = [(f, info[f][0], V, i0, min(i0 + CHUNK, info[f][1]), lev_of[f], nlev)
                    for f in fl for i0 in range(0, info[f][1], CHUNK)]
            if not jobs:
                continue
            Sm = Sd = Ss = N = 0; bad = 0
            for sm, sd, ss, n, bd in ex.map(pass1, jobs):
                Sm, Sd, Ss, N, bad = Sm + sm, Sd + sd, Ss + ss, N + n, bad + bd
            with np.errstate(invalid='ignore', divide='ignore'):
                mmean = np.where(N > 0, Sm / N, 0.0); dmean = np.where(N > 0, Sd / N, 0.0)
                s = np.where(N > 0, Ss / N, 0.0)
                sig2 = np.maximum(s / np.where(N > 0, N, 1), s / a.indsamples)
                wmean = np.where((N > 0) & (sig2 != 0), 1.0 / np.where(sig2 != 0, sig2, 1), 0.0)
            sel = (wmean > 0) & (dmean != 0) & (mmean != 0)
            jmean = float(np.sum(wmean[sel] * (mmean[sel] - dmean[sel]) ** 2)); nmean = int(sel.sum())

            per_file = {}
            jobs2 = [j + (mmean, dmean) for j in jobs]
            for j, (c, n) in zip(jobs2, ex.map(pass2, jobs2)):
                k = os.path.relpath(j[0], a.data)
                pc, pn = per_file.get(k, (0.0, 0.0)); per_file[k] = (pc + c, pn + n)
            janom = sum(c for c, _ in per_file.values()); nanom = sum(n for _, n in per_file.values())
            total += janom + jmean
            result[V] = {'anomaly_cost': janom, 'anomaly_count': nanom, 'mean_cost': jmean, 'mean_count': nmean,
                         'per_file': {k: {'cost': c, 'count': n} for k, (c, n) in per_file.items()},
                         'nonfinite_points_used': bad}

            print(f'\n{V}:')
            for i, (k, (c, n)) in enumerate(sorted(per_file.items()), 1):
                print(f'  cost_profiles({i:2d},{V})= {c:12.5E} {n:12.5E}   {k}')
            print(f'  cost_profiles_mean({V})= {jmean:12.5E} {float(nmean):12.5E}')
            print(f'  anomaly total {janom:.6E} (count {nanom:.0f}), mean {jmean:.6E}; cost per anomaly point '
                  f'{janom / nanom if nanom else float("nan"):.4f}')
            if bad:
                print(f'  WARNING: {bad} points with weight > 0 have non-finite data, model or weight values; left out')
    result['total'] = total
    print(f'\nTotal profiles cost (all multipliers 1): {total:.6E}')
    if a.out:
        os.makedirs(a.out, exist_ok=True)
        with open(os.path.join(a.out, 'cost_summary.json'), 'w') as fh:
            json.dump(result, fh, indent=1)
        print('wrote', os.path.join(a.out, 'cost_summary.json'))


if __name__ == '__main__':
    sys.exit(cli())
