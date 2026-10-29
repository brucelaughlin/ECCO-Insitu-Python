"""
Profile analysis framework — compute step.

For each obs class (source), read all chain output files and compute the 10
standard analyses.  Save as per-source NetCDF files plus figures.

    ecco-insitu profile-analyze --run <run_dir> --out <out_dir>
                                [--sources XBT CTD_WOD ...]
                                [--model-equivalents <dir>]
                                [--run-b <dir>]
                                [-n N]

    ecco-insitu profile-analyze --replot --out <out_dir>    # redraw the figures only, from the
                                                            # saved NetCDF (about a minute)

The 10 analyses (per obs class, for each cost type and split/unsplit mode):

  1.  Global mean T and S vs depth
  2.  Global mean T and S cost vs depth
  3.  Global mean cost per datum vs time  (factors removed = 'true' weight version)
  4.  Global mean total cost vs time      (final weights)
  5.  Vertically integrated cost per geodesic bin
  6.  Vertically integrated cost per datum per geodesic bin  (factors removed)
  7.  Zonal mean cost per datum vs depth
  8.  Histogram of per-datum cost
  9.  Global mean cost per datum vs depth  (analysis 2 with 'true' weights, separate figure)
  10. Global mean cost per datum summed into depth ranges
      (0–100, 0–250, 0–500, 0–700, 0–2000 m)

Cost types:
  vs_clim    (prof_V - prof_Vclim)²    always
  vs_model   (prof_V - prof_Vestim)²   if --model-equivalents given
  vs_run_b   same but --run-b          if --run-b given

Model-equivalent folders (--model-equivalents, --run-b) may hold either
  - ``<yearly file name>_model.nc``, one per yearly file; or
  - ``<packed file name>_model.nc``, one per model-ready (step 13) file in ``<run>/ecco_arrays/``, as
    ``ecco-insitu mds-to-nc`` (a model run) and ``ecco-insitu monthly-estim`` (approximate, from monthly
    means) write them. Their values are mapped back to the yearly files through the packed files'
    ``prof_yearly_file`` / ``prof_yearly_index``.
A yearly file with no model values for some profiles leaves those out of the model costs, and the count
is printed.

Split vs unsplit (per point; w = weight, V = observation, Vref = reference):
  unsplit    w (V - Vref)²                       the mean misfit is not removed
  split      w (V - Vref - b)²,  b = mean(V) - mean(Vref) per (geodesic bin, depth), pooled
             over all obs classes, points with w > 0 (and, for model costs, Vref != 0):
             ECCO v4r6 cost_profiles.F's anomaly term. Costs are still tallied per obs class.
             ECCO's mean term, w_mean (mean(Vref) - mean(V))², w_mean = 1/max(s/n, s/9),
             s = mean(1/w), is pooled, so it is reported once per cost type (meta.json and
             the NetCDF attributes), not per obs class.
  split_own  as split, but b from this obs class only.
  Bias arrays (pooled and own) are written to the NetCDF.

Weight versions:
  final   γ × sf × (ΔV/σ)²  what ECCO sees
  nosf    γ × (ΔV/σ)²       spatial scaling removed
  true    (ΔV/σ)²           both factors removed ('per datum')
"""
import argparse
import datetime
import json
import os
import time as _time
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'axes.grid': True, 'grid.alpha': 0.5, 'grid.linewidth': 0.6})
import matplotlib.colors as mcolors
import netCDF4
import numpy as np

import ecco_insitu.diagnostics.common as common

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
VARS = ('T', 'S')
WEIGHT_VERS = ('final', 'nosf', 'true')   # matches costs_compute.py KINDS
SPLIT_MODES = ('unsplit', 'split', 'split_own')
PROFILES_MEAN_INDSAMPLES = 9   # v4r6 data.profiles
# Analysis 11: data inside vs after the climatology's years; read from the run's climatology file
# (its 'title', e.g. 'WOA23 1991-2020 ...'); this default is WOA23 decav91C0's
CLIM_YEARS = (1991, 2020)
PERIODS = ('in_clim_years', 'after_clim_years')
COST_TYPES  = ('vs_clim', 'vs_model', 'vs_run_b')

N_BINS   = 10242          # geodesic bin ids 1..10242; index 0 unused
N_DEPTHS = 97             # expected; actual depth array read from first file

# Analysis 10 depth ranges: upper bound in metres
DEPTH_RANGE_BOUNDS  = [100, 250, 500, 700, 2000]
DEPTH_RANGE_LABELS  = ['0–100 m', '0–250 m', '0–500 m', '0–700 m', '0–2000 m']
N_DEPTH_RANGES      = len(DEPTH_RANGE_BOUNDS)

# Analysis 7 latitude bands: 5° bands
LAT_EDGES    = np.linspace(-90, 90, 37).astype(np.float32)   # 36 bands
N_LAT_BANDS  = 36

# Analysis 8 histogram: log10(cost) bins
HIST_EDGES   = np.linspace(-4, 5, 181)   # 180 bins covering 10^-4 to 10^5
N_HIST_BINS  = len(HIST_EDGES) - 1

_FILL_9S = np.array([999, 9999, 99999, 999999, 9999999], dtype=np.float64)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _is_fill_value(x):
    return ~np.isfinite(x) | (np.abs(x) >= 900)


def _is_fill_weight(w):
    a = np.abs(w)
    return ~np.isfinite(w) | (a >= 1e30) | np.isin(a, _FILL_9S)


def _read(v, name, dtype=np.float64):
    if name not in v:
        return None
    x = v[name][:]
    x = np.ma.filled(np.ma.asarray(x, dtype=dtype), np.nan) if dtype == np.float64 else np.asarray(x, dtype=dtype)
    return x


def _lat_band_idx(lat):
    """Map latitudes to band index 0..N_LAT_BANDS-1."""
    return np.clip(np.searchsorted(LAT_EDGES[1:-1], lat), 0, N_LAT_BANDS - 1).astype(np.int32)


def _yyyymm(ymd):
    """YYYYMMDD integer array → YYYYMM integer array."""
    return (np.asarray(ymd, dtype=np.int64) // 100).astype(np.int32)

# ---------------------------------------------------------------------------
# File loading
# ---------------------------------------------------------------------------

class EstimIndex:
    """Model equivalents for the yearly files of a run, from a folder of ``*_model.nc`` files (see the module
    description): per-yearly files are used directly; per-packed files are mapped back to the yearly files
    through ``<run>/ecco_arrays``."""

    def __init__(self, estim_dir, run_dir):
        self.dir = Path(estim_dir)
        self.map = {}                    # yearly file name -> list of (model file, packed rows, yearly rows)
        arrays = common.run_root(run_dir) / 'ecco_arrays'
        if not arrays.is_dir():
            return
        import json
        for packed in sorted(arrays.glob('*.nc')):
            mf = self.dir / (packed.stem + '_model.nc')
            if not mf.exists():
                continue
            with netCDF4.Dataset(packed) as ds:
                ds.set_auto_mask(False)
                names = json.loads(ds.getncattr('yearly_files'))
                fi = np.asarray(ds['prof_yearly_file'][:], np.int64)
                ix = np.asarray(ds['prof_yearly_index'][:], np.int64)
            for k, name in enumerate(names):
                rows = np.nonzero(fi == k)[0]
                self.map.setdefault(name, []).append((mf, rows, ix[rows]))

    def get(self, yearly_path, n, n_depth):
        """{V: (n, n_depth) model values, NaN where none} for one yearly file, or None if there are none."""
        own = self.dir / (Path(yearly_path).stem + '_model.nc')
        if own.exists():
            out = {}
            with netCDF4.Dataset(own) as ds:
                ds.set_auto_mask(False)
                for V in VARS:
                    out[V] = _read(ds.variables, f'prof_{V}estim')
            return out
        parts = self.map.get(Path(yearly_path).name)
        if not parts:
            return None
        out = {V: np.full((n, n_depth), np.nan) for V in VARS}
        for mf, rows, yrows in parts:
            with netCDF4.Dataset(mf) as ds:
                ds.set_auto_mask(False)
                for V in VARS:
                    if f'prof_{V}estim' in ds.variables:
                        lo, hi = rows.min(), rows.max() + 1          # read one block, then pick the rows
                        out[V][yrows] = np.asarray(ds[f'prof_{V}estim'][lo:hi], float)[rows - lo]
        return out


def _load_file(path, estim_path=None, estim_index=None):
    """Load arrays from one profile NetCDF.  Returns a dict or None on error.

    Model equivalents come from ``estim_index`` (an :class:`EstimIndex`) if given, else from
    ``estim_path``, else from ``prof_Vestim`` in the file itself."""
    out = {}
    try:
        with netCDF4.Dataset(path) as ds:
            ds.set_auto_mask(False)
            v = ds.variables
            n = ds.dimensions['iPROF'].size
            if n == 0:
                return None
            out['depth']     = _read(v, 'prof_depth', np.float32)
            out['lat']       = _read(v, 'prof_lat')
            out['lon']       = _read(v, 'prof_lon')
            out['ymd']       = _read(v, 'prof_YYYYMMDD', np.int32)
            out['bin_a']     = _read(v, 'prof_bin_id_a', np.int32)
            out['gamma']     = _read(v, 'prof_area_gamma')
            for V in VARS:
                out[f'val_{V}']   = _read(v, f'prof_{V}')
                out[f'clim_{V}']  = _read(v, f'prof_{V}clim')
                out[f'w_{V}']     = _read(v, f'prof_{V}weight')
                out[f'sf_{V}']    = _read(v, f'prof_{V}spatial_scaling_factor')
                # model equivalents embedded in the data file
                out[f'estim_{V}'] = _read(v, f'prof_{V}estim')
        if estim_index is not None:
            nd = out['depth'].size if out['depth'] is not None else 0
            e = estim_index.get(path, n, nd)
            if e is not None:
                for V in VARS:
                    if e[V] is not None:
                        out[f'estim_{V}'] = e[V]
        elif estim_path and Path(estim_path).exists():
            with netCDF4.Dataset(estim_path) as ds:
                ds.set_auto_mask(False)
                for V in VARS:
                    e = _read(ds.variables, f'prof_{V}estim')
                    if e is not None:
                        out[f'estim_{V}'] = e
    except Exception as exc:
        print(f'  [warn] {Path(path).name}: {exc}')
        return None
    # defaults for missing optional arrays
    if out.get('gamma') is None:
        out['gamma'] = np.ones(n)
    for V in VARS:
        if out[f'sf_{V}'] is None:
            out[f'sf_{V}'] = np.ones(n)
        # treat -9999 sentinel as missing model value
        if out[f'estim_{V}'] is not None:
            out[f'estim_{V}'] = np.where(out[f'estim_{V}'] == -9999, np.nan, out[f'estim_{V}'])
    out['n'] = n
    return out


def _weight_versions(w, sf, gamma):
    """Return dict of three weight-version arrays from the final weight w."""
    sf_ok    = np.where(np.isfinite(sf) & (sf > 0), sf, 1.0)
    gamma_ok = np.where(np.isfinite(gamma) & (gamma > 0), gamma, 1.0)
    return {
        'final': w,
        'nosf':  w / sf_ok[:, None],
        'true':  w / (gamma_ok[:, None] * sf_ok[:, None]),
    }

# ---------------------------------------------------------------------------
# Accumulators
# ---------------------------------------------------------------------------

def _new_acc(n_depths):
    """Fresh accumulator dict for one (source, cost_type) pass."""
    acc = {
        # --- analysis 1: T/S vs depth (no cost) ---
        'a01_T_sum':   np.zeros(n_depths),
        'a01_T_count': np.zeros(n_depths, dtype=np.int64),
        'a01_S_sum':   np.zeros(n_depths),
        'a01_S_count': np.zeros(n_depths, dtype=np.int64),
        # --- bias accumulators for split cost, shape (N_BINS+1, n_depths) ---
        'bias_sum_T':   np.zeros((N_BINS + 1, n_depths)),
        'bias_sum_S':   np.zeros((N_BINS + 1, n_depths)),
        'bias_count_T': np.zeros((N_BINS + 1, n_depths), dtype=np.int64),
        'bias_count_S': np.zeros((N_BINS + 1, n_depths), dtype=np.int64),
    }
    # cost analyses: a02..a10, per variable, per weight version, per split mode
    for V in VARS:
        for sp in SPLIT_MODES:
            for wv in WEIGHT_VERS:
                k = (V, sp, wv)
                # a02: mean cost vs depth   (n_depths,)
                acc[('a02', *k, 'sum')]   = np.zeros(n_depths)
                acc[('a02', *k, 'count')] = np.zeros(n_depths, dtype=np.int64)
                # a03: mean cost per datum vs time  (YYYYMM → running total)
                acc[('a03', *k, 'sum')]   = defaultdict(float)
                acc[('a03', *k, 'count')] = defaultdict(int)
                # a05: integrated cost per bin  (N_BINS+1,)
                acc[('a05', *k, 'sum')]   = np.zeros(N_BINS + 1)
                acc[('a05', *k, 'count')] = np.zeros(N_BINS + 1, dtype=np.int64)
                # a05 thickness-weighted: per profile sum_k c_k dz_k / sum_k dz_k (non-empty k), summed per bin
                acc[('a05tw', *k, 'sum')]   = np.zeros(N_BINS + 1)
                acc[('a05tw', *k, 'count')] = np.zeros(N_BINS + 1, dtype=np.int64)
                # a07: zonal mean cost vs depth  (N_LAT_BANDS, n_depths)
                acc[('a07', *k, 'sum')]   = np.zeros((N_LAT_BANDS, n_depths))
                acc[('a07', *k, 'count')] = np.zeros((N_LAT_BANDS, n_depths), dtype=np.int64)
                # a08: histogram of log10 cost  (N_HIST_BINS,)
                acc[('a08', *k, 'hist')]  = np.zeros(N_HIST_BINS)
                # a11: cost vs depth, by period relative to the climatology's years
                for pe in PERIODS:
                    acc[('a11', *k, pe, 'sum')]   = np.zeros(n_depths)
                    acc[('a11', *k, pe, 'count')] = np.zeros(n_depths, dtype=np.int64)
    return acc


def _new_pooled(n_depths):
    """Sums per (bin, depth) over every obs class, for the Fortran-style split (cost_profiles.F)."""
    def z():
        return np.zeros((N_BINS + 1, n_depths))
    return {V: {'d': z(), 'm': z(), 'invw': z(), 'n': np.zeros((N_BINS + 1, n_depths), np.int64)}
            for V in VARS}


def _pooled_stats(pooled, indsamples=PROFILES_MEAN_INDSAMPLES):
    """Bin means and ECCO's mean term from pooled sums, as cost_profiles.F (and ``model-cost``).

    b = d_mean - m_mean (the time-mean misfit, NaN where undefined); mean term
    sum of w_mean * (m_mean - d_mean)^2 with w_mean = 1 / max(s/n, s/indsamples), s = mean(1/w).
    """
    out = {}
    for V, p in pooled.items():
        n = p['n']
        with np.errstate(invalid='ignore', divide='ignore'):
            dmean = np.where(n > 0, p['d'] / n, 0.0)
            mmean = np.where(n > 0, p['m'] / n, 0.0)
            s = np.where(n > 0, p['invw'] / n, 0.0)
            sig2 = np.maximum(s / np.where(n > 0, n, 1), s / indsamples)
            wmean = np.where((n > 0) & (sig2 != 0), 1.0 / np.where(sig2 != 0, sig2, 1), 0.0)
        sel = (wmean > 0) & (dmean != 0) & (mmean != 0)
        out[V] = {
            'bias': np.where((n > 0) & (dmean != 0), dmean - mmean, np.nan),
            'count': n, 'wmean': wmean,
            'mean_cost': float(np.sum(wmean[sel] * (mmean[sel] - dmean[sel]) ** 2)),
            'mean_count': int(sel.sum()),
        }
    return out


def _layer_thickness(depth):
    """Thickness (m) of the layer around each depth level: edges at the surface, the midpoints between
    levels, and half a spacing below the deepest level."""
    d = np.asarray(depth, float)
    edges = np.concatenate([[0.0], (d[1:] + d[:-1]) / 2, [d[-1] + (d[-1] - d[-2]) / 2]])
    return np.diff(edges)


def clim_years(run_dir):
    """(first, last) year of the climatology the run used, from the first 'YYYY-YYYY' in its file's
    title; ``CLIM_YEARS`` if the file or the years can't be found."""
    import re
    import yaml
    cfg = common.run_root(run_dir) / 'documentation' / 'run_config.yaml'
    try:
        path = (yaml.safe_load(cfg.read_text()) or {})['paths']['climatology_file']
        with netCDF4.Dataset(path) as ds:
            m = re.search(r'(\d{4})\s*[-–]\s*(\d{4})', getattr(ds, 'title', ''))
        return (int(m.group(1)), int(m.group(2))) if m else CLIM_YEARS
    except (OSError, KeyError, TypeError):
        return CLIM_YEARS


def _add_costs(acc, V, mode, mf, ok, w, sf, gamma, ym, bi_safe, lb, dz):
    """Add one file's costs (misfit ``mf`` at points ``ok``) to every analysis, for each weight version."""
    weights = _weight_versions(w, sf, gamma)
    yr = ym // 100
    in_clim = (yr >= CLIM_YEARS[0]) & (yr <= CLIM_YEARS[1])
    period_ok = {'in_clim_years': ok & in_clim[:, None], 'after_clim_years': ok & (yr > CLIM_YEARS[1])[:, None]}
    dz_ok = np.where(ok, dz[None, :], 0.0)
    dz_tot = dz_ok.sum(axis=1)
    has = dz_tot > 0
    rows_ok, cols_ok = np.where(ok)
    n_d = mf.shape[1]
    for wv in WEIGHT_VERS:
        c2d = np.where(ok, mf, 0.0) ** 2 * weights[wv]
        k = (V, mode, wv)
        np.add.at(acc[('a02', *k, 'sum')],   cols_ok, c2d[ok])
        np.add.at(acc[('a02', *k, 'count')], cols_ok, 1)
        for pe, okp in period_ok.items():
            cols_p = np.where(okp)[1]
            np.add.at(acc[('a11', *k, pe, 'sum')],   cols_p, c2d[okp])
            np.add.at(acc[('a11', *k, pe, 'count')], cols_p, 1)
        c_per_prof = np.where(ok, c2d, 0.0).sum(axis=1)
        n_per_prof = ok.sum(axis=1).astype(np.int64)
        for i in np.nonzero(n_per_prof)[0]:
            ymv = int(ym[i])
            acc[('a03', *k, 'sum')][ymv]   += float(c_per_prof[i])
            acc[('a03', *k, 'count')][ymv] += int(n_per_prof[i])
        np.add.at(acc[('a05', *k, 'sum')],   bi_safe, c_per_prof)
        np.add.at(acc[('a05', *k, 'count')], bi_safe, n_per_prof)
        tw = (np.where(ok, c2d, 0.0) * dz_ok).sum(axis=1)[has] / dz_tot[has]
        np.add.at(acc[('a05tw', *k, 'sum')],   bi_safe[has], tw)
        np.add.at(acc[('a05tw', *k, 'count')], bi_safe[has], 1)
        if rows_ok.size:
            flat_idx = lb[rows_ok] * n_d + cols_ok
            np.add.at(acc[('a07', *k, 'sum')].ravel(),   flat_idx, c2d[ok])
            np.add.at(acc[('a07', *k, 'count')].ravel(), flat_idx, 1)
        pos = c2d[ok & (c2d > 0)]
        if pos.size:
            hist, _ = np.histogram(np.log10(pos), bins=HIST_EDGES)
            acc[('a08', *k, 'hist')] += hist


def _accumulate_file(acc, fd, cost_type, pass_, biases=None, pooled=None):
    """Add one file to the accumulators of one cost type.

    ``pass_ = 'unsplit'``: unsplit costs, analysis 1, the class's own bias sums, and (if given) the
    ``pooled`` sums over all obs classes. ``pass_ = 'split'``: the split costs, one per entry of
    ``biases`` ({split mode: {V: (N_BINS+1, n_depths) bias, NaN where undefined}}).
    """
    n = fd['n']
    bin_a = fd['bin_a'] if fd.get('bin_a') is not None else np.zeros(n, dtype=np.int32)
    valid_bin = (bin_a >= 1) & (bin_a <= N_BINS)
    bi_safe = np.clip(bin_a, 0, N_BINS).astype(np.int32)
    ym = _yyyymm(fd['ymd'])
    lb = _lat_band_idx(fd['lat'])
    dz = _layer_thickness(fd['depth'])

    for V in VARS:
        val, w, sf = fd[f'val_{V}'], fd[f'w_{V}'], fd[f'sf_{V}']
        if val is None or fd[f'clim_{V}'] is None or w is None:
            continue
        if pass_ == 'unsplit':
            # analysis 1: observations only, the same for every cost type
            obs = ~_is_fill_value(val) & ~_is_fill_weight(w) & (w > 0)
            cols = np.where(obs)[1]
            np.add.at(acc[f'a01_{V}_sum'],   cols, val[obs])
            np.add.at(acc[f'a01_{V}_count'], cols, 1)
        if cost_type == 'vs_clim':
            ref = fd[f'clim_{V}']
        else:
            ref = fd.get(f'estim_{V}' if cost_type == 'vs_model' else f'estimb_{V}')
            if ref is None:
                continue
        ok = ~_is_fill_value(val) & ~_is_fill_value(ref) & ~_is_fill_weight(w) & (w > 0)
        if cost_type != 'vs_clim':
            ok &= ref != 0          # the model's "no value" (cost_profiles.F mask)
        mf = val - ref

        if pass_ == 'unsplit':
            okb = ok & valid_bin[:, None]
            rows, cols = np.where(okb)
            flat = bi_safe[rows] * mf.shape[1] + cols
            np.add.at(acc[f'bias_sum_{V}'].ravel(),   flat, mf[okb])
            np.add.at(acc[f'bias_count_{V}'].ravel(), flat, 1)
            if pooled is not None:
                p = pooled[V]
                np.add.at(p['d'].ravel(),    flat, val[okb])
                np.add.at(p['m'].ravel(),    flat, ref[okb])
                np.add.at(p['invw'].ravel(), flat, 1.0 / w[okb])
                np.add.at(p['n'].ravel(),    flat, 1)
            _add_costs(acc, V, 'unsplit', mf, ok, w, sf, fd['gamma'], ym, bi_safe, lb, dz)
        else:
            for mode, bias in biases.items():
                b = bias[V][bi_safe]                     # (n_prof, n_depths); bin 0 has no mean
                okm = ok & valid_bin[:, None] & np.isfinite(b)
                _add_costs(acc, V, mode, mf - np.where(okm, b, 0.0), okm, w, sf, fd['gamma'],
                           ym, bi_safe, lb, dz)


# ---------------------------------------------------------------------------
# Finalization
# ---------------------------------------------------------------------------

def _finalize(acc, depth):
    """Convert running sums to means / normalized quantities.  Returns results dict."""
    n_d = len(depth)
    res = {'depth': depth}

    # analysis 1
    res['a01'] = {}
    for V in VARS:
        s = acc[f'a01_{V}_sum']
        c = acc[f'a01_{V}_count']
        with np.errstate(invalid='ignore'):
            res['a01'][V] = {'mean': np.where(c > 0, s / c, np.nan), 'count': c.copy()}

    # biases (for documentation)
    res['bias'] = {}
    for V in VARS:
        s = acc[f'bias_sum_{V}']
        c = acc[f'bias_count_{V}']
        with np.errstate(invalid='ignore'):
            res['bias'][V] = np.where(c > 0, s / c, np.nan)   # (N_BINS+1, n_depths)
        res[f'bias_count_{V}'] = c

    # analyses 2–10
    for V in VARS:
        for sp in SPLIT_MODES:
            for wv in WEIGHT_VERS:
                k = (V, sp, wv)
                # a02: mean cost vs depth
                s = acc[('a02', *k, 'sum')]
                c = acc[('a02', *k, 'count')]
                with np.errstate(invalid='ignore'):
                    res[('a02', *k)] = {'mean': np.where(c > 0, s / c, np.nan), 'count': c.copy()}

                # a03: mean cost per datum vs time
                ym_d = acc[('a03', *k, 'sum')]
                cn_d = acc[('a03', *k, 'count')]
                yms = sorted(ym_d.keys())
                with np.errstate(invalid='ignore'):
                    a03_vals = np.array([ym_d[y] / cn_d[y] if cn_d[y] > 0 else np.nan for y in yms])
                res[('a03', *k)] = {'yyyymm': np.array(yms, dtype=np.int32),
                                    'mean_per_datum': a03_vals,
                                    'count': np.array([cn_d[y] for y in yms], dtype=np.int64)}

                # a04: total cost vs time (same accumulators as a03 — see comment below)
                # a04 accumulates the same sum and count, but is reported as total cost
                # (not normalised by count); same dict for sum
                res[('a04', *k)] = {'yyyymm': np.array(yms, dtype=np.int32),
                                    'total': np.array([ym_d[y] for y in yms]),
                                    'count': np.array([cn_d[y] for y in yms], dtype=np.int64)}

                # a05: integrated cost per bin
                s5 = acc[('a05', *k, 'sum')]
                c5 = acc[('a05', *k, 'count')]
                with np.errstate(invalid='ignore'):
                    res[('a05', *k)] = {
                        'total':     np.where(c5[1:] > 0, s5[1:], np.nan),    # bins 1..N_BINS; NaN if empty
                        'per_datum': np.where(c5[1:] > 0, s5[1:] / c5[1:], np.nan),
                        'count':     c5[1:].copy().astype(np.int64),
                    }
                st = acc[('a05tw', *k, 'sum')]
                ct_ = acc[('a05tw', *k, 'count')]
                res[('a05tw', *k)] = {'total': np.where(ct_[1:] > 0, st[1:], np.nan),
                                      'count': ct_[1:].copy().astype(np.int64)}

                # a07: zonal mean cost vs depth
                s7 = acc[('a07', *k, 'sum')]
                c7 = acc[('a07', *k, 'count')]
                with np.errstate(invalid='ignore'):
                    res[('a07', *k)] = {'mean': np.where(c7 > 0, s7 / c7, np.nan), 'count': c7.copy()}

                # a08: histogram
                res[('a08', *k)] = acc[('a08', *k, 'hist')].copy()

                # a09: same as a02 (cost per datum vs depth) — referenced from a02 with 'true'
                # (no separate accumulator; figure will use a02 'true')

                # a10: depth-range summaries from a02
                a02_mean = res[('a02', *k)]['mean']
                a02_cnt  = res[('a02', *k)]['count']
                depth_ranges = []
                for bound in DEPTH_RANGE_BOUNDS:
                    mask = depth <= bound
                    n_pts = a02_cnt[mask].sum()
                    if n_pts > 0:
                        val = np.nansum(a02_mean[mask] * a02_cnt[mask]) / n_pts
                    else:
                        val = np.nan
                    depth_ranges.append(val)
                res[('a10', *k)] = np.array(depth_ranges)

                # a11: cost vs depth by period, and its 0-2000 m mean
                for pe in PERIODS:
                    s11 = acc[('a11', *k, pe, 'sum')]
                    c11 = acc[('a11', *k, pe, 'count')]
                    top = depth <= DEPTH_RANGE_BOUNDS[-1]
                    with np.errstate(invalid='ignore'):
                        res[('a11', *k, pe)] = {
                            'mean': np.where(c11 > 0, s11 / c11, np.nan), 'count': c11.copy(),
                            'mean_0_2000': s11[top].sum() / c11[top].sum() if c11[top].sum() else np.nan}

    return res

# ---------------------------------------------------------------------------
# NetCDF writer
# ---------------------------------------------------------------------------

def _write_nc(res, source, run_ts, path, cost_types_used, pooled_sources=None):
    """Write all analysis results to a NetCDF file."""
    depth   = res['depth']
    n_d     = len(depth)
    # collect all yyyymm keys across all (V, sp, wv, cost_type)
    all_ym = set()
    for k in res:
        if isinstance(k, tuple) and k[0] == 'a03':
            all_ym.update(res[k]['yyyymm'].tolist())
    months = np.array(sorted(all_ym), dtype=np.int32)

    path.parent.mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset(path, 'w') as nc:
        nc.setncatts({
            'source': source, 'run_ts': run_ts,
            'split_method': 'cost_profiles_F_pooled_bin_mean',
            'split_note': (
                'unsplit = w*(V-Vref)^2. split = w*(V-Vref-b)^2 with b = mean(V) - mean(Vref) per '
                '(geodesic bin, depth) pooled over all obs classes in split_pooled_over, as ECCO '
                'cost_profiles.F (the anomaly term). split_own = the same with b from this obs class '
                'only. ECCO\'s mean term w_mean*(mean(Vref)-mean(V))^2 is in the attributes '
                'ecco_mean_term_<V> (pooled, not per obs class). Bias arrays are in this file.'
            ),
            'split_pooled_over': ', '.join(pooled_sources or [source]),
            'profiles_mean_indsamples': PROFILES_MEAN_INDSAMPLES,
            'a11_clim_years': f'{CLIM_YEARS[0]}-{CLIM_YEARS[1]}',
            'weight_versions': 'final = gamma*sf*(dV/sigma)^2; nosf = gamma*(dV/sigma)^2; true = (dV/sigma)^2',
            'cost_types_used': ', '.join(cost_types_used),
            'created': datetime.datetime.now().isoformat(timespec='seconds'),
        })
        nc.createDimension('depth',        n_d)
        nc.createDimension('month',        len(months))
        nc.createDimension('bin_a',        N_BINS)
        nc.createDimension('lat_band',     N_LAT_BANDS)
        nc.createDimension('depth_range',  N_DEPTH_RANGES)
        nc.createDimension('hist_edge',    len(HIST_EDGES))
        nc.createDimension('hist_bin',     N_HIST_BINS)
        nc.createDimension('bias_bin_a',   N_BINS + 1)

        def cv(name, data, dims, dtype='f4', **attrs):
            kw = dict(zlib=True, complevel=4)
            if dtype in ('f4', 'f8'):
                kw['fill_value'] = np.nan
            x = nc.createVariable(name, dtype, dims, **kw)
            x[:] = data
            x.setncatts(attrs)

        cv('depth_m',        depth, ('depth',),       'f4', units='m', long_name='depth level')
        cv('yyyymm',         months, ('month',),      'i4', long_name='year-month (YYYYMM)')
        cv('hist_edges',     HIST_EDGES, ('hist_edge',), 'f4',
           long_name='log10 histogram bin edges for cost per datum')
        cv('lat_band_center',
           (LAT_EDGES[:-1] + LAT_EDGES[1:]) / 2, ('lat_band',), 'f4', units='degrees_north')
        depth_range_upper = np.array(DEPTH_RANGE_BOUNDS, dtype=np.int32)
        nc.createVariable('depth_range_upper', 'i4', ('depth_range',))[:] = depth_range_upper
        labels_var = nc.createVariable('depth_range_label', str, ('depth_range',))
        for i, lbl in enumerate(DEPTH_RANGE_LABELS):
            labels_var[i] = lbl

        # analysis 1
        cv('a01_T_mean_depth', res['a01']['T']['mean'], ('depth',), 'f4',
           long_name='Global mean T vs depth (all weighted observations; the same for every cost type)')
        cv('a01_S_mean_depth', res['a01']['S']['mean'], ('depth',), 'f4',
           long_name='Global mean S vs depth (all weighted observations; the same for every cost type)')
        nc.createVariable('a01_T_count_depth', 'i8', ('depth',))[:] = res['a01']['T']['count']
        nc.createVariable('a01_S_count_depth', 'i8', ('depth',))[:] = res['a01']['S']['count']

        # biases
        for V in VARS:
            cv(f'bias_own_{V}', res['bias'][V][1:], ('bin_a', 'depth'), 'f4',
               long_name=f'Time-mean (prof_{V} - ref) per (bin, depth), this obs class only (split_own)')
            nc.createVariable(f'bias_own_count_{V}', 'i8', ('bin_a', 'depth'))[:] = \
                res[f'bias_count_{V}'][1:]
            pl = res.get('pooled', {}).get(V)
            if pl is not None:
                cv(f'bias_pooled_{V}', pl['bias'][1:], ('bin_a', 'depth'), 'f4',
                   long_name=f'mean(prof_{V}) - mean(ref) per (bin, depth), all obs classes pooled (split)')
                nc.createVariable(f'bias_pooled_count_{V}', 'i8', ('bin_a', 'depth'))[:] = pl['count'][1:]
                nc.setncattr(f'ecco_mean_term_{V}', pl['mean_cost'])
                nc.setncattr(f'ecco_mean_term_count_{V}', pl['mean_count'])

        # analyses 2–10
        for V in VARS:
            for sp in SPLIT_MODES:
                for wv in WEIGHT_VERS:
                    k = (V, sp, wv)
                    pfx = f'a__{V}__{sp}__{wv}'

                    # a02
                    r2 = res[('a02', *k)]
                    cv(f'a02__{pfx}__mean_depth', r2['mean'], ('depth',), 'f4',
                       long_name=f'[a02] {V} mean cost vs depth ({sp}, {wv})')
                    nc.createVariable(f'a02__{pfx}__count_depth', 'i8', ('depth',))[:] = r2['count']

                    # a03
                    r3 = res[('a03', *k)]
                    ym_idx = np.searchsorted(months, r3['yyyymm'])
                    tmp3 = np.full(len(months), np.nan, dtype=np.float32)
                    tmp3[ym_idx] = r3['mean_per_datum']
                    cv(f'a03__{pfx}__mean_per_datum_time', tmp3, ('month',), 'f4',
                       long_name=f'[a03] {V} mean cost per datum vs time ({sp}, {wv})')
                    cnt3 = np.zeros(len(months), dtype=np.int64)
                    cnt3[ym_idx] = r3['count']
                    nc.createVariable(f'a03__{pfx}__count_time', 'i8', ('month',))[:] = cnt3

                    # a04
                    r4 = res[('a04', *k)]
                    tmp4 = np.full(len(months), np.nan, dtype=np.float32)
                    tmp4[ym_idx] = r4['total']
                    cv(f'a04__{pfx}__total_cost_time', tmp4, ('month',), 'f4',
                       long_name=f'[a04] {V} total cost vs time ({sp}, {wv})')

                    # a05
                    r5 = res[('a05', *k)]
                    cv(f'a05__{pfx}__integrated_by_bin', r5['total'], ('bin_a',), 'f4',
                       long_name=f'[a05] {V} vertically integrated cost per geodesic bin ({sp}, {wv})')
                    cv(f'a06__{pfx}__per_datum_by_bin', r5['per_datum'], ('bin_a',), 'f4',
                       long_name=f'[a06] {V} integrated cost per datum per geodesic bin ({sp}, {wv})')
                    nc.createVariable(f'a05__{pfx}__count_bin', 'i8', ('bin_a',))[:] = r5['count']
                    rt = res[('a05tw', *k)]
                    cv(f'a05tw__{pfx}__thickness_weighted_by_bin', rt['total'], ('bin_a',), 'f4',
                       long_name=f'[a05] {V} cost per geodesic bin, thickness weighted: sum over profiles of '
                                 f'sum_k c_k dz_k / sum_k dz_k (non-empty levels) ({sp}, {wv})')
                    nc.createVariable(f'a05tw__{pfx}__profiles_bin', 'i8', ('bin_a',))[:] = rt['count']

                    # a07
                    r7 = res[('a07', *k)]
                    cv(f'a07__{pfx}__zonal_mean_depth', r7['mean'], ('lat_band', 'depth'), 'f4',
                       long_name=f'[a07] {V} zonal mean cost per datum vs depth ({sp}, {wv})')

                    # a08
                    nc.createVariable(f'a08__{pfx}__histogram', 'f4', ('hist_bin',)
                                      )[:] = res[('a08', *k)]

                    # a10
                    cv(f'a10__{pfx}__by_depth_range', res[('a10', *k)], ('depth_range',), 'f4',
                       long_name=f'[a10] {V} mean cost per datum by depth range ({sp}, {wv})')

                    # a11
                    for pe in PERIODS:
                        r11 = res[('a11', *k, pe)]
                        cv(f'a11__{pfx}__{pe}__mean_depth', r11['mean'], ('depth',), 'f4',
                           long_name=f'[a11] {V} mean cost per datum vs depth, {pe} ({sp}, {wv})')
                        nc.createVariable(f'a11__{pfx}__{pe}__count_depth', 'i8', ('depth',))[:] = r11['count']

    print(f'  wrote {path}')

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

_V_COLOR = {'T': 'tab:red', 'S': 'tab:blue'}
_SP_LS   = {'unsplit': '--', 'split': '-', 'split_own': ':'}
_WV_COLOR = {'final': 'tab:blue', 'nosf': 'tab:orange', 'true': 'tab:red'}
# Time-series split/unsplit: slightly darker shade + circle markers for split
_TS_SP_COLOR = {
    ('T', 'unsplit'): 'tab:red',  ('T', 'split'): '#8B0000', ('T', 'split_own'): '#E8865A',
    ('S', 'unsplit'): 'tab:blue', ('S', 'split'): '#00008B', ('S', 'split_own'): '#5AA8E8',
}
_SP_MARKER = {'unsplit': '', 'split': ''}
# Only suppress depth levels with zero observations. The Mediterranean spike in XBT below
# ~2200 m (T ≈ 13 °C) is real data and is explained in the report's Analysis 1 note.
_MIN_DEPTH_COUNT = 1


# Fixed margins (fractions of the figure): every cost type of a figure is then saved at the same size
# with its axes at the same pixel positions, so switching tabs changes only the plotted data.
# (bbox_inches='tight' would crop each image to its own titles and legends, which shifts the frame.)
_MARGINS = dict(left=0.09, right=0.97, bottom=0.11, top=0.88)


def _savefig(fig, path, margins=True):
    path.parent.mkdir(parents=True, exist_ok=True)
    if margins:
        fig.subplots_adjust(**_MARGINS)
    fig.savefig(path, dpi=110)
    plt.close(fig)


def _decimal_year(ym):
    ym = np.asarray(ym, dtype=int)
    return (ym // 100) + (ym % 100 - 0.5) / 12.0


def _span(arrays, positive=False):
    vals = [np.asarray(a, dtype=float).ravel() for a in arrays]
    v = np.concatenate(vals) if vals else np.array([])
    v = v[np.isfinite(v) & ((v > 0) if positive else True)]
    return (float(v.min()), float(v.max())) if v.size else None


def _pad(span, frac=0.05, zero=False):
    if span is None:
        return None
    lo, hi = span
    if zero:
        lo = min(lo, 0.0)
    d = (hi - lo) * frac or abs(hi) * frac or 1.0
    return (0.0 if zero and lo == 0 else lo - d, hi + d)


def _p98(arrays):
    vals = [p for a in arrays
            for p in [np.nanpercentile(a[np.isfinite(a)], 98) if np.isfinite(a).any() else np.nan]]
    vmax = np.nanmax(vals) if np.isfinite(vals).any() else np.nan
    return float(vmax) if np.isfinite(vmax) and vmax > 0 else 1.0


def _binned(res, V, wv, field, key='a05'):
    """Per-bin map values (analyses 5/6), NaN in bins with no data."""
    r = res[(key, V, 'unsplit', wv)]
    return np.where(r['count'] > 0, r[field], np.nan)


def _shared_limits(results):
    """Axis limits per (analysis, variable) spanning every cost type's results."""
    rs = list(results.values())
    lims = {}
    for V in VARS:
        lims[('a01', V)] = _pad(_span(
            [np.where(r['a01'][V]['count'] >= _MIN_DEPTH_COUNT, r['a01'][V]['mean'], np.nan) for r in rs]))
        for an, wvs in (('02', WEIGHT_VERS), ('09', ('true',))):
            lims[(f'a{an}', V)] = _pad(_span(
                [np.where(r[('a02', V, sp, wv)]['count'] >= _MIN_DEPTH_COUNT, r[('a02', V, sp, wv)]['mean'], np.nan)
                 for r in rs for sp in SPLIT_MODES for wv in wvs]), zero=True)
        lims[('a03', V)] = _pad(_span([r[('a03', V, sp, 'true')]['mean_per_datum']
                                       for r in rs for sp in SPLIT_MODES]), zero=True)
        s4 = _span([r[('a04', V, 'unsplit', wv)]['total'] for r in rs for wv in WEIGHT_VERS], positive=True)
        lims[('a04', V)] = None if s4 is None else (s4[0] / 1.5, s4[1] * 1.5)
        lims[('a05', V)] = _p98([_binned(r, V, 'final', 'total') for r in rs])
        lims[('a05tw', V)] = _p98([_binned(r, V, 'final', 'total', 'a05tw') for r in rs])
        lims[('a06', V)] = _p98([_binned(r, V, 'true', 'per_datum') for r in rs])
        lims[('a07', V)] = _p98([r[('a07', V, 'unsplit', 'true')]['mean'] for r in rs])
        lims[('a08', V)] = _pad(_span([h / h.sum() for r in rs for wv in WEIGHT_VERS
                                       for h in [r[('a08', V, 'unsplit', wv)]] if h.sum() > 0]), zero=True)
        lims[('a11', V)] = _pad(_span([np.where(r[('a11', V, sp, 'true', pe)]['count'] > 0,
                                                r[('a11', V, sp, 'true', pe)]['mean'], np.nan)
                                       for r in rs for sp in ('unsplit', 'split') for pe in PERIODS]), zero=True)
        lims[('a10', V)] = _pad(_span([r[('a10', V, sp, 'true')] for r in rs for sp in SPLIT_MODES]),
                                frac=0.1, zero=True)
    for an, key in (('a01', None), ('a02', ('a02', 'unsplit', 'final'))):
        lims[(an, 'n')] = _log_lim([r['a01'][V]['count'] if key is None else r[(key[0], V, *key[1:])]['count']
                                    for r in rs for V in VARS])
    lims[('a11', 'n')] = _log_lim([r[('a11', V, 'unsplit', 'true', pe)]['count'] for r in rs for V in VARS
                                   for pe in PERIODS])
    lims[('a04', 'n')] = _log_lim([r[('a03', V, 'unsplit', 'true')]['count'] for r in rs for V in VARS])
    lims['year'] = _pad(_span([_decimal_year(r[('a03', V, sp, 'true')]['yyyymm'])
                               for r in rs for V in VARS for sp in SPLIT_MODES]), frac=0.01)
    return lims


_CT_LABEL = {'vs_clim': 'vs climatology (WOA23)', 'vs_model': 'vs model (v4r6 monthly means)',
             'vs_run_b': 'vs model run B (monthly means)'}


def _log_lim(arrays):
    s = _span(arrays, positive=True)
    return None if s is None else (max(s[0] / 1.5, 0.8), s[1] * 1.5)


def _depth_figure(n_main, width_main=4.5, height=6):
    """Figure with ``n_main`` depth panels and a narrow panel on the right for the number of points."""
    fig = plt.figure(figsize=(width_main * n_main + 1.8, height))
    gs = fig.add_gridspec(1, n_main + 1, width_ratios=[4] * n_main + [1.2], wspace=0.12)
    axes = [fig.add_subplot(gs[0, 0])]
    axes += [fig.add_subplot(gs[0, i], sharey=axes[0]) for i in range(1, n_main + 1)]
    return fig, axes[:-1], axes[-1]


def _count_ticks(axis, lim, max_labels=5):
    """Readable log ticks for a count axis: at most ``max_labels`` labels (1, 1-3, 1-2-5 or every few
    powers of ten, depending on the range), short labels (2k, 50k, 1M), and a faint minor grid."""
    from matplotlib.ticker import FixedLocator, FuncFormatter, LogLocator, NullFormatter

    def fmt(v, _):
        for div, suf in ((1e9, 'G'), (1e6, 'M'), (1e3, 'k')):
            if v >= div:
                return f'{v / div:g}{suf}'
        return f'{v:g}'
    lo, hi = lim if lim else axis.get_view_interval()
    lo, hi = max(lo, 1e-9), max(hi, lo * 10)
    e0, e1 = int(np.floor(np.log10(lo))), int(np.ceil(np.log10(hi)))
    for subs in ((1, 2, 5), (1, 3), (1,)):
        ticks = [m * 10.0 ** e for e in range(e0, e1 + 1) for m in subs if lo <= m * 10.0 ** e <= hi]
        if len(ticks) <= max_labels:
            break
    else:
        step = int(np.ceil(len(ticks) / max_labels))
        ticks = ticks[::step]
    axis.set_major_locator(FixedLocator(ticks))
    axis.set_major_formatter(FuncFormatter(fmt))
    axis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10), numticks=30))
    axis.set_minor_formatter(NullFormatter())


def _count_panel(ax, counts, depth, lim):
    """Number of points per depth level (log scale), one line per variable."""
    for V, c in counts.items():
        ax.semilogx(np.where(c > 0, c, np.nan), depth, color=_V_COLOR[V], label=V, lw=1.1)
    _apply(ax.set_xlim, lim)
    _count_ticks(ax.xaxis, lim, 4)
    ax.tick_params(axis='x', labelsize=7, labelrotation=90)
    ax.grid(True, which='minor', axis='x', alpha=0.2, lw=0.4)
    ax.set_xlabel('# points')
    ax.set_title('n')
    ax.tick_params(labelleft=False)
    if len(counts) > 1:
        ax.legend(fontsize=7)


def _apply(set_lim, lim):
    if lim is not None:
        set_lim(*lim)


def _make_figures(res, source, fig_dir, lims=None, cost_type=None):
    """Generate all 10 analysis figures; ``lims`` (from _shared_limits) fixes the axis ranges."""
    lims = lims if lims is not None else _shared_limits({'_': res})
    ctl = f' — {_CT_LABEL.get(cost_type, cost_type)}' if cost_type else ''
    depth = res['depth']
    _dmax = float(np.nanmax(depth)) * 1.05  # bottom of every depth axis

    # Only plot variables that have valid data in this obs class (e.g. XBT has no S)
    active_vars = [V for V in VARS if np.isfinite(res['a01'][V]['mean']).any()]
    if not active_vars:
        return
    nv = len(active_vars)

    # ---- Analysis 1: T and S mean profile ----
    fig, axes, cax = _depth_figure(nv, 4)
    for ax, V in zip(axes, active_vars):
        m = res['a01'][V]['mean']
        c = res['a01'][V]['count']
        m_plot = np.where(c >= _MIN_DEPTH_COUNT, m, np.nan)
        ax.plot(m_plot, depth, color=_V_COLOR[V])
        ax.set_xlabel(f'Mean {V} (°C)' if V == 'T' else f'Mean {V} (psu)')
        ax.set_ylabel('Depth (m)') if ax is axes[0] else None
        ax.set_title(f'{V} (n_max={c.max():,.0f})')
        _apply(ax.set_xlim, lims[('a01', V)])
    _count_panel(cax, {V: res['a01'][V]['count'] for V in active_vars}, depth, lims[('a01', 'n')])
    axes[0].set_ylim(_dmax, 0)
    fig.suptitle(f'{source}: Analysis 1 — Global mean T & S vs depth (observations only)')
    _savefig(fig, fig_dir / '01_mean_T_S_depth.png')

    # ---- Analysis 2/9: cost vs depth ----
    for an, wvs, title_sfx in [
        ('02', WEIGHT_VERS, 'mean cost'),
        ('09', ('true',),   'mean cost per datum'),
    ]:
        fig, axes, cax = _depth_figure(nv, 5)
        for ax, V in zip(axes, active_vars):
            c = res[('a02', V, 'unsplit', 'final')]['count']  # count is same across sp/wv
            for sp in SPLIT_MODES:
                for wv in wvs:
                    if an == '09' and wv != 'true':
                        continue
                    m = res[('a02', V, sp, wv)]['mean']
                    m_plot = np.where(c >= _MIN_DEPTH_COUNT, m, np.nan)
                    ax.plot(m_plot, depth, ls=_SP_LS[sp], color=_WV_COLOR[wv],
                            label=f'{wv} {sp}', lw=1.2)
            ax.set_xlabel(f'{V} {title_sfx}')
            ax.set_ylabel('Depth (m)') if ax is axes[0] else None
            ax.set_title(V)
            ax.legend(fontsize=7)
            _apply(ax.set_xlim, lims[(f'a{an}', V)])
        _count_panel(cax, {V: res[('a02', V, 'unsplit', 'final')]['count'] for V in active_vars},
                     depth, lims[('a02', 'n')])
        axes[0].set_ylim(_dmax, 0)
        fig.suptitle(f'{source}: Analysis {an} — {title_sfx} vs depth{ctl}')
        _savefig(fig, fig_dir / f'{an}_cost_depth_{source}.png')

    # ---- Analysis 3: mean cost per datum vs time ----
    fig, axes = plt.subplots(nv, 1, figsize=(10, 3.5*nv + 1), sharex=True, squeeze=False)
    axes = axes.flatten()
    for ax, V in zip(axes, active_vars):
        any_data = False
        for sp in SPLIT_MODES:
            r = res[('a03', V, sp, 'true')]
            if len(r['yyyymm']) == 0:
                continue
            any_data = True
            x = _decimal_year(r['yyyymm'])
            ax.plot(x, r['mean_per_datum'],
                    color=_TS_SP_COLOR[(V, sp)], ls=_SP_LS[sp],
                    label=sp, lw=1.1)
        ax.set_ylabel(f'{V} mean cost/datum')
        if any_data:
            ax.legend(fontsize=7)
        _apply(ax.set_ylim, lims[('a03', V)])
        _apply(ax.set_xlim, lims['year'])
    axes[-1].set_xlabel('Year')
    fig.suptitle(f'{source}: Analysis 3 — Mean cost per datum vs time (true weights){ctl}')
    _savefig(fig, fig_dir / f'03_cost_per_datum_time_{source}.png')

    # ---- Analysis 4: total cost vs time ----
    fig = plt.figure(figsize=(10, 3.5 * nv + 2.2))
    gs = fig.add_gridspec(nv + 1, 1, height_ratios=[1.1] + [3] * nv, hspace=0.12)
    nax = fig.add_subplot(gs[0])
    axes = [fig.add_subplot(gs[i + 1], sharex=nax) for i in range(nv)]
    for V in active_vars:
        r = res[('a03', V, 'unsplit', 'true')]
        if len(r['yyyymm']):
            nax.semilogy(_decimal_year(r['yyyymm']), np.where(r['count'] > 0, r['count'], np.nan),
                         color=_V_COLOR[V], label=V, lw=1.0)
    nax.set_ylabel('# points\nper month')
    _count_ticks(nax.yaxis, lims[('a04', 'n')], 4)
    nax.tick_params(axis='y', labelsize=7)
    nax.grid(True, which='minor', axis='y', alpha=0.2, lw=0.4)
    _apply(nax.set_ylim, lims[('a04', 'n')])
    nax.tick_params(labelbottom=False)
    if nv > 1:
        nax.legend(fontsize=7)
    for ax, V in zip(axes, active_vars):
        if ax is not axes[-1]:
            ax.tick_params(labelbottom=False)
        for wv in WEIGHT_VERS:
            r = res[('a04', V, 'unsplit', wv)]
            if len(r['yyyymm']) == 0:
                continue
            x = _decimal_year(r['yyyymm'])
            ax.semilogy(x, np.where(np.array(r['total']) > 0, r['total'], np.nan),
                        color=_WV_COLOR[wv], label=wv, lw=1.1)
        ax.set_ylabel(f'{V} total cost\n(sum over the month)')
        ax.legend(fontsize=7)
        _apply(ax.set_ylim, lims[('a04', V)])
        _apply(ax.set_xlim, lims['year'])
    axes[-1].set_xlabel('Year')
    fig.suptitle(f'{source}: Analysis 4 — Total cost per month (sum over all points), unsplit{ctl}')
    _savefig(fig, fig_dir / f'04_total_cost_time_{source}.png')

    # ---- Analysis 5/6: spatial map ----
    try:
        _make_spatial_map(res, source, fig_dir, active_vars, lims, ctl)
    except Exception as exc:
        print(f'  [warn] spatial map failed: {exc}')

    # ---- Analysis 7: zonal mean cost per datum vs depth ----
    for V in active_vars:
        k = (V, 'unsplit', 'true')
        r = res[('a07', *k)]['mean']   # (N_LAT_BANDS, n_depths)
        lat_centers = (LAT_EDGES[:-1] + LAT_EDGES[1:]) / 2
        valid = np.isfinite(r).any(axis=1)
        if not valid.any():
            continue
        fig = plt.figure(figsize=(11, 6))
        gs = fig.add_gridspec(1, 2, width_ratios=[5, 1], wspace=0.08)
        ax = fig.add_subplot(gs[0])
        cax = fig.add_subplot(gs[1], sharey=ax)
        _count_panel(cax, {V: res[('a02', V, 'unsplit', 'true')]['count']}, depth, lims[('a02', 'n')])
        vmax = lims[('a07', V)]
        im = ax.contourf(lat_centers, depth, r.T, levels=np.linspace(0, vmax, 21),
                         cmap='YlOrRd', extend='max')
        ax.set_ylim(_dmax, 0)
        ax.set_xlim(LAT_EDGES[0], LAT_EDGES[-1])
        plt.colorbar(im, ax=ax, label=f'{V} cost per datum (ΔV/σ)²', pad=0.02)
        ax.set_xlabel('Latitude (°N)')
        ax.set_ylabel('Depth (m)')
        fig.suptitle(f'{source}: Analysis 7 — Zonal mean {V} cost per datum (unsplit, true){ctl}')
        _savefig(fig, fig_dir / f'07_zonal_{V}_{source}.png')

    # ---- Analysis 8: histogram ----
    fig, axes = plt.subplots(1, nv, figsize=(5*nv + 1, 5), squeeze=False)
    axes = axes.flatten()
    bin_centers = (HIST_EDGES[:-1] + HIST_EDGES[1:]) / 2
    for ax, V in zip(axes, active_vars):
        for wv in WEIGHT_VERS:
            h = res[('a08', V, 'unsplit', wv)]
            if h.sum() > 0:
                ax.plot(bin_centers, h / h.sum(), color=_WV_COLOR[wv], label=wv, lw=1.2)
        ax.axvline(np.log10(100), color='k', ls=':', lw=1, label='QC limit=100')
        ax.set_xlabel(f'log₁₀({V} cost per datum)')
        ax.set_ylabel('Density')
        ax.set_title(V)
        ax.legend(fontsize=7)
        ax.set_xlim(HIST_EDGES[0], HIST_EDGES[-1])
        _apply(ax.set_ylim, lims[('a08', V)])
    fig.suptitle(f'{source}: Analysis 8 — Histogram of cost per datum{ctl}')
    _savefig(fig, fig_dir / f'08_histogram_{source}.png')

    # ---- Analysis 10: depth-range summary ----
    fig, axes = plt.subplots(1, nv, figsize=(5*nv + 1, 5), squeeze=False)
    axes = axes.flatten()
    x = np.arange(N_DEPTH_RANGES)
    width = 0.8 / len(SPLIT_MODES)
    for ax, V in zip(axes, active_vars):
        for i, sp in enumerate(SPLIT_MODES):
            vals = res[('a10', V, sp, 'true')]
            ax.bar(x + i * width, vals, width, label=sp, alpha=0.8)
        ax.set_xticks(x + width * (len(SPLIT_MODES) - 1) / 2)
        ax.set_xticklabels(DEPTH_RANGE_LABELS, rotation=30, ha='right', fontsize=8)
        ax.set_ylabel(f'{V} mean cost per datum (ΔV/σ)²')
        ax.set_title(V)
        ax.legend(fontsize=7)
        _apply(ax.set_ylim, lims[('a10', V)])
    fig.suptitle(f'{source}: Analysis 10 — Mean cost per datum by depth range{ctl}')
    _savefig(fig, fig_dir / f'10_depth_ranges_{source}.png')

    # ---- Analysis 11 (extra): cost per datum vs depth, inside vs after the climatology's years ----
    has_after = any(res[('a11', V, 'unsplit', 'true', 'after_clim_years')]['count'].sum() for V in active_vars)
    if has_after:
        last = max(int(np.max(res[('a03', V, 'unsplit', 'true')]['yyyymm'])) // 100
                   for V in active_vars if len(res[('a03', V, 'unsplit', 'true')]['yyyymm']))
        after = f'{CLIM_YEARS[1] + 1}' + (f'–{last}' if last > CLIM_YEARS[1] + 1 else '')
        pe_style = {'in_clim_years': ('tab:blue', f'{CLIM_YEARS[0]}–{CLIM_YEARS[1]} (in climatology)'),
                    'after_clim_years': ('tab:red', f'{after} (after climatology)')}
        fig, axes, cax = _depth_figure(nv, 5)
        for ax, V in zip(axes, active_vars):
            for pe, (col, lab) in pe_style.items():
                for sp in ('unsplit', 'split'):
                    r = res[('a11', V, sp, 'true', pe)]
                    m = np.where(r['count'] > 0, r['mean'], np.nan)
                    ax.plot(m, depth, color=col, ls=_SP_LS[sp], lw=1.2,
                            label=f'{lab}, {sp}: {r["mean_0_2000"]:.2f}')
            ax.set_xlabel(f'{V} mean cost per datum (true)')
            ax.set_ylabel('Depth (m)') if ax is axes[0] else None
            ax.set_title(V)
            ax.legend(fontsize=7, title='period, split mode: 0–2000 m mean', title_fontsize=7)
            _apply(ax.set_xlim, lims[('a11', V)])
        for pe, (col, lab) in pe_style.items():
            c = sum(res[('a11', V, 'unsplit', 'true', pe)]['count'] for V in active_vars[:1])
            cax.semilogx(np.where(c > 0, c, np.nan), depth, color=col, lw=1.1, label=lab.split(' (')[0])
        _apply(cax.set_xlim, lims[('a11', 'n')])
        _count_ticks(cax.xaxis, lims[('a11', 'n')], 4)
        cax.tick_params(axis='x', labelsize=7, labelrotation=90)
        cax.grid(True, which='minor', axis='x', alpha=0.2, lw=0.4)
        cax.set_xlabel('# points')
        cax.set_title(f'n ({active_vars[0]})')
        cax.tick_params(labelleft=False)
        cax.legend(fontsize=7)
        axes[0].set_ylim(_dmax, 0)
        fig.suptitle(f'{source}: Analysis 11 — cost per datum vs depth, data inside vs after the '
                     f'climatology\'s years ({CLIM_YEARS[0]}–{CLIM_YEARS[1]}){ctl}')
        _savefig(fig, fig_dir / f'11_clim_years_{source}.png')


_RASTERS = {}


def _rasters():
    """Bin-id rasters for the three map views (global Robinson, north and south polar stereographic, to 50°)."""
    if not _RASTERS:
        bl = common.BinLookup()
        for view, npix in (('global', 700), ('NH', 520), ('SH', 520)):
            _RASTERS[view] = bl.projected_raster(view, npix=npix)
    return _RASTERS


def _make_spatial_map(res, source, fig_dir, active_vars, lims, ctl=''):
    """Analyses 5 and 6: per-bin maps, global plus north and south polar stereographic panels."""
    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
    except ImportError:
        print('  [warn] cartopy not available; skipping spatial maps')
        return

    views = [('global', ccrs.Robinson(), 'Global'),
             ('NH', ccrs.NorthPolarStereo(), 'North of 50°N'),
             ('SH', ccrs.SouthPolarStereo(), 'South of 50°S')]
    rasters = _rasters()

    for V in active_vars:
        for an, key, field, sp, wv, title, sfx in [
            ('05', 'a05',   'total',     'unsplit', 'final', 'total cost (sum over depth levels and profiles)', ''),
            ('05', 'a05tw', 'total',     'unsplit', 'final',
             'total cost, thickness weighted (sum over profiles of each profile\'s thickness-weighted mean)',
             '_thickness_weighted'),
            ('06', 'a05',   'per_datum', 'unsplit', 'true',  'cost per datum', ''),
        ]:
            vals = _binned(res, V, wv, field, key)            # (N_BINS,) bins 1..N_BINS, NaN if empty
            per_bin_full = np.concatenate([[np.nan], vals])   # index 0 = no bin
            vmax = lims[('a06' if an == '06' else key, V)]
            norm = mcolors.Normalize(0, vmax)

            # fixed axes rectangles [left, bottom, width, height] (figure fractions), so every cost type
            # puts the maps and color bar at the same pixels
            fig = plt.figure(figsize=(16, 5.6))
            rects = [[0.02, 0.22, 0.50, 0.66], [0.54, 0.22, 0.21, 0.66], [0.77, 0.22, 0.21, 0.66]]
            for i, (view, proj, label) in enumerate(views):
                raster, land, outside, extent = rasters[view]
                grid_a = raster['a']
                img = np.where(land | outside | (grid_a == 0), np.nan, per_bin_full[grid_a])
                ax = fig.add_axes(rects[i], projection=proj)
                ax.imshow(img, origin='lower', extent=extent, transform=proj, cmap='YlOrRd',
                          norm=norm, zorder=0, interpolation='nearest')
                ax.add_feature(cfeature.LAND, color='lightgray', zorder=1)
                if view == 'global':
                    ax.set_global()
                else:
                    ax.set_extent(extent, crs=proj)
                ax.gridlines(linewidth=0.5, alpha=0.5, color='gray')
                ax.set_title(label, fontsize=9)
            fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='YlOrRd'), cax=fig.add_axes([0.3, 0.1, 0.4, 0.035]),
                         orientation='horizontal', label=f'{V} {title} ({sp}, {wv})')
            fig.suptitle(f'{source}: Analysis {an} — {V} {title} per geodesic bin{ctl}', fontsize=11, y=0.97)
            _savefig(fig, fig_dir / f'{an}_spatial_{V}_{source}{sfx}.png', margins=False)

# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def _source_files(source, files, run_dir, estim_dir, estim_dir_b):
    """Model-equivalent indexes for one obs class, and a loader that attaches them to each file."""
    est_index = EstimIndex(estim_dir, run_dir) if estim_dir else None
    est_index_b = EstimIndex(estim_dir_b, run_dir) if estim_dir_b else None
    if est_index is not None:
        found = sum(1 for f in files if (Path(estim_dir) / (Path(f).stem + '_model.nc')).exists()
                    or Path(f).name in est_index.map)
        print(f'    {source}: model equivalents for {found} of {len(files)} files')

    def load(f):
        fd = _load_file(f, estim_index=est_index)
        if fd is not None and est_index_b is not None:
            fb = est_index_b.get(f, fd['n'], fd['depth'].size)
            for V in VARS:
                e = None if fb is None else fb[V]
                fd[f'estimb_{V}'] = None if e is None else np.where(e == -9999, np.nan, e)
        return fd
    return load


def _progress(i):
    if (i + 1) % 50 == 0:
        print(f' {i+1}', end='', flush=True)


def run(run_dir, out_dir, sources=None, estim_dir=None, estim_dir_b=None, n_workers=None):
    """Compute all standard analyses for every obs class (source) in ``run_dir``.

    Pass 1 reads every file once: unsplit costs, each class's own bias sums, and sums per
    (geodesic bin, depth) pooled over **all** obs classes (as ECCO's ``cost_profiles.F``).
    Pass 2 reads each file again for the two split costs: ``split`` (the pooled, Fortran-style
    mean removed) and ``split_own`` (the class's own mean removed). Results are per obs class.
    The pooled means come from the obs classes processed, so ``sources`` changes ``split``.

    Parameters
    ----------
    run_dir : str or Path
        Processed run directory (or its yearly/ sub-folder).
    out_dir : str or Path
        Output directory.
    sources : list of str, optional
        Only these sources.  Default: all.
    estim_dir : str or Path, optional
        Directory with ``<name>_model.nc`` model equivalents (enables vs_model cost).
    estim_dir_b : str or Path, optional
        Second model-equivalent directory (enables vs_run_b cost).
    n_workers : int, optional
        Unused (kept for the command line); the passes run serially.

    Returns
    -------
    Path
        The output directory.
    """
    global CLIM_YEARS
    run_dir = common.yearly_dir(run_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    CLIM_YEARS = clim_years(run_dir)

    cost_types = ['vs_clim'] + (['vs_model'] if estim_dir else []) + (['vs_run_b'] if estim_dir_b else [])

    groups = {}
    for f in sorted(run_dir.rglob('*.nc')):
        src, _ = common.classify(f, run_dir)
        if sources and src not in sources:
            continue
        groups.setdefault(src, []).append(f)
    if not groups:
        print('No files found.')
        return out_dir

    t0 = _time.time()
    print(f'Run: {run_dir}\nOut: {out_dir}')
    print(f'{sum(map(len, groups.values()))} files across {len(groups)} sources')
    print(f'Cost types: {cost_types}')

    depth = None
    for files in groups.values():
        for f in files:
            fd = _load_file(f)
            if fd is not None and fd.get('depth') is not None:
                depth = fd['depth'].astype(np.float32)
                break
        if depth is not None:
            break
    n_d = len(depth)
    loaders = {src: _source_files(src, files, run_dir, estim_dir, estim_dir_b) for src, files in groups.items()}

    accs = {src: {ct: _new_acc(n_d) for ct in cost_types} for src in groups}
    pooled = {ct: _new_pooled(n_d) for ct in cost_types}

    print('Pass 1 (unsplit, bin means over all obs classes)')
    for src, files in sorted(groups.items()):
        print(f'  {src}:', end='', flush=True)
        for i, f in enumerate(files):
            fd = loaders[src](f)
            if fd is None:
                continue
            if fd['depth'].size != n_d or not np.array_equal(fd['depth'], depth):
                raise ValueError(f'{f}: depth levels differ from the other files')
            for ct in cost_types:
                _accumulate_file(accs[src][ct], fd, ct, 'unsplit', pooled=pooled[ct])
            _progress(i)
        print()

    stats = {ct: _pooled_stats(pooled[ct]) for ct in cost_types}
    for ct in cost_types:
        print(f'  {ct}: ECCO mean term T {stats[ct]["T"]["mean_cost"]:.6E} ({stats[ct]["T"]["mean_count"]} terms), '
              f'S {stats[ct]["S"]["mean_cost"]:.6E} ({stats[ct]["S"]["mean_count"]} terms)')

    print('Pass 2 (split: pooled mean and own-class mean removed)')
    for src, files in sorted(groups.items()):
        own = {}
        for ct in cost_types:
            own[ct] = {}
            for V in VARS:
                s_, c_ = accs[src][ct][f'bias_sum_{V}'], accs[src][ct][f'bias_count_{V}']
                with np.errstate(invalid='ignore'):
                    own[ct][V] = np.where(c_ > 0, s_ / c_, np.nan)
        print(f'  {src}:', end='', flush=True)
        for i, f in enumerate(files):
            fd = loaders[src](f)
            if fd is None:
                continue
            for ct in cost_types:
                biases = {'split': {V: stats[ct][V]['bias'] for V in VARS}, 'split_own': own[ct]}
                _accumulate_file(accs[src][ct], fd, ct, 'split', biases=biases)
            _progress(i)
        print()

    run_ts_name = common.run_timestamp(run_dir)
    for src in sorted(groups):
        results = {ct: _finalize(accs[src][ct], depth) for ct in cost_types}
        for ct in cost_types:
            results[ct]['pooled'] = stats[ct]
        lims = _shared_limits(results)
        for ct, res in results.items():
            out_nc = out_dir / src / f'analysis_{src}_{run_ts_name}_{ct}.nc'
            _write_nc(res, src, run_ts_name, out_nc, [ct], pooled_sources=sorted(groups))
            _make_figures(res, src, out_dir / src / 'figures' / ct, lims, ct)
        print(f'  {src}: figures done ({_time.time() - t0:.0f}s)')

    meta = {
        'run': str(run_dir), 'out': str(out_dir),
        'sources': sorted(groups.keys()),
        'cost_types': cost_types,
        'estim_dir': str(estim_dir) if estim_dir else None,
        'estim_dir_b': str(estim_dir_b) if estim_dir_b else None,
        'split_pooled_over': sorted(groups.keys()),
        'clim_years': list(CLIM_YEARS),
        'ecco_mean_term': {ct: {V: {'cost': stats[ct][V]['mean_cost'], 'count': stats[ct][V]['mean_count']}
                                for V in VARS} for ct in cost_types},
        'created': datetime.datetime.now().isoformat(timespec='seconds'),
    }
    (out_dir / 'meta.json').write_text(json.dumps(meta, indent=2))
    print(f'All done in {_time.time() - t0:.0f}s → {out_dir}')
    return out_dir


def load_results(path):
    """Rebuild the in-memory results of one obs class and cost type from its saved NetCDF
    (what ``_make_figures`` needs), so figures can be redrawn without recomputing."""
    with netCDF4.Dataset(path) as nc:
        nc.set_auto_mask(False)
        g = lambda name: np.asarray(nc[name][:])
        depth = g('depth_m').astype(np.float32)
        months = g('yyyymm').astype(np.int32)
        res = {'depth': depth,
               'a01': {V: {'mean': g(f'a01_{V}_mean_depth').astype(float),
                           'count': g(f'a01_{V}_count_depth')} for V in VARS}}
        top = depth <= DEPTH_RANGE_BOUNDS[-1]
        for V in VARS:
            for sp in SPLIT_MODES:
                for wv in WEIGHT_VERS:
                    k = (V, sp, wv)
                    pfx = f'a__{V}__{sp}__{wv}'
                    res[('a02', *k)] = {'mean': g(f'a02__{pfx}__mean_depth').astype(float),
                                        'count': g(f'a02__{pfx}__count_depth')}
                    m3 = g(f'a03__{pfx}__mean_per_datum_time').astype(float)
                    t4 = g(f'a04__{pfx}__total_cost_time').astype(float)
                    if f'a03__{pfx}__count_time' in nc.variables:
                        c3 = g(f'a03__{pfx}__count_time')
                    else:   # files written before the counts were saved: total / mean
                        with np.errstate(invalid='ignore', divide='ignore'):
                            c3 = np.where(m3 > 0, np.rint(t4 / m3), 0).astype(np.int64)
                    have = np.isfinite(m3)
                    res[('a03', *k)] = {'yyyymm': months[have], 'mean_per_datum': m3[have], 'count': c3[have]}
                    res[('a04', *k)] = {'yyyymm': months[have], 'total': t4[have], 'count': c3[have]}
                    res[('a05', *k)] = {'total': g(f'a05__{pfx}__integrated_by_bin').astype(float),
                                        'per_datum': g(f'a06__{pfx}__per_datum_by_bin').astype(float),
                                        'count': g(f'a05__{pfx}__count_bin')}
                    if f'a05tw__{pfx}__thickness_weighted_by_bin' in nc.variables:
                        res[('a05tw', *k)] = {'total': g(f'a05tw__{pfx}__thickness_weighted_by_bin').astype(float),
                                              'count': g(f'a05tw__{pfx}__profiles_bin')}
                    res[('a07', *k)] = {'mean': g(f'a07__{pfx}__zonal_mean_depth').astype(float)}
                    res[('a08', *k)] = g(f'a08__{pfx}__histogram').astype(float)
                    res[('a10', *k)] = g(f'a10__{pfx}__by_depth_range').astype(float)
                    for pe in PERIODS:
                        name = f'a11__{pfx}__{pe}__mean_depth'
                        if name not in nc.variables:
                            continue
                        m = g(name).astype(float)
                        c = g(f'a11__{pfx}__{pe}__count_depth')
                        n = c[top].sum()
                        res[('a11', *k, pe)] = {'mean': m, 'count': c,
                                                'mean_0_2000': float(np.nansum(m[top] * c[top]) / n) if n else np.nan}
    return res


def replot(out_dir, sources=None):
    """Redraw every figure of a ``profile-analyze`` output folder from its saved NetCDF files
    (minutes, instead of re-reading the run). Uses the cost types and climatology years in meta.json."""
    global CLIM_YEARS
    out_dir = Path(out_dir)
    meta = json.loads((out_dir / 'meta.json').read_text())
    CLIM_YEARS = tuple(meta.get('clim_years', CLIM_YEARS))
    t0 = _time.time()
    for src in meta['sources']:
        if sources and src not in sources:
            continue
        results = {}
        for ct in meta['cost_types']:
            f = sorted((out_dir / src).glob(f'analysis_{src}_*_{ct}.nc'))
            if f:
                results[ct] = load_results(f[-1])
        if not results:
            continue
        lims = _shared_limits(results)
        for ct, res in results.items():
            _make_figures(res, src, out_dir / src / 'figures' / ct, lims, ct)
        print(f'  {src}: figures redrawn ({_time.time() - t0:.0f}s)')
    print(f'Replotted {out_dir} in {_time.time() - t0:.0f}s')
    return out_dir


def cli(argv=None):
    """CLI entry point (``ecco-insitu profile-analyze``)."""
    p = argparse.ArgumentParser(
        description='Compute profile analyses (T/S profiles, cost distributions) '
                    'for each obs class in a processed run.')
    p.add_argument('--run', help='Processed run directory (or its yearly/ sub-folder); not needed with --replot')
    p.add_argument('--out', required=True, help='Output directory')
    p.add_argument('--sources', nargs='*', default=None,
                   help='Obs classes to process (default: all)')
    p.add_argument('--model-equivalents', default=None, dest='estim_dir',
                   help='Dir with <name>_model.nc files for vs-model cost')
    p.add_argument('--run-b', default=None, dest='estim_dir_b',
                   help='Second model-equivalent dir for vs-run-B cost')
    p.add_argument('-n', type=int, default=None,
                   help='Parallel workers (default: 1)')
    p.add_argument('--replot', action='store_true',
                   help='Only redraw the figures in --out from its saved NetCDF files (no recomputing)')
    a = p.parse_args(argv)
    if a.replot:
        replot(a.out, a.sources)
        return
    if not a.run:
        p.error('--run is required (unless --replot)')
    run(a.run, a.out, a.sources, a.estim_dir, a.estim_dir_b, a.n)
