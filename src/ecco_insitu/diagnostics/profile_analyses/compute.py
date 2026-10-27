"""
Profile analysis framework — compute step.

For each obs class (source), read all chain output files and compute the 10
standard analyses.  Save as per-source NetCDF files plus figures.

    ecco-insitu profile-analyze --run <run_dir> --out <out_dir>
                                [--sources XBT CTD_WOD ...]
                                [--model-equivalents <dir>]
                                [--run-b <dir>]
                                [-n N]

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

Split vs unsplit:
  unsplit  raw misfit  (V - Vref)² × w
  split    Fortran-style: (V - Vref - bias)² × w
           bias[bin,depth] = time-mean (V - Vref) across all obs at that (bin,depth)
  Both saved; bias arrays are also written to NetCDF for reproducibility.
  Global attribute split_method = 'mean_bias_per_bin_depth' documents the formula.

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
import matplotlib.colors as mcolors
import netCDF4
import numpy as np

import ecco_insitu.diagnostics.common as common

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
VARS = ('T', 'S')
WEIGHT_VERS = ('final', 'nosf', 'true')   # matches costs_compute.py KINDS
SPLIT_MODES = ('unsplit', 'split')
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
                # a04: total cost vs time
                acc[('a04', *k, 'sum')]   = defaultdict(float)
                acc[('a04', *k, 'count')] = defaultdict(int)
                # a05: integrated cost per bin  (N_BINS+1,)
                acc[('a05', *k, 'sum')]   = np.zeros(N_BINS + 1)
                acc[('a05', *k, 'count')] = np.zeros(N_BINS + 1, dtype=np.int64)
                # a07: zonal mean cost vs depth  (N_LAT_BANDS, n_depths)
                acc[('a07', *k, 'sum')]   = np.zeros((N_LAT_BANDS, n_depths))
                acc[('a07', *k, 'count')] = np.zeros((N_LAT_BANDS, n_depths), dtype=np.int64)
                # a08: histogram of log10 cost  (N_HIST_BINS,)
                acc[('a08', *k, 'hist')]  = np.zeros(N_HIST_BINS)
    return acc


def _accumulate_cost(acc, V, sp, cost2d, ok, weights, lat, bin_a, yyyymm_arr):
    """Add one file's cost values to all cost-analysis accumulators.

    cost2d  : (n_prof, n_depths) cost array for one weight version → supplied for all 3
    ok      : (n_prof, n_depths) bool mask of valid, nonzero-weight points
    weights : dict of {wver: (n_prof, n_depths) weight array}
    lat     : (n_prof,)
    bin_a   : (n_prof,) int32 bin ids (1..10242; clipped to valid range)
    yyyymm_arr: (n_prof,) YYYYMM int
    """
    n_prof, n_depths = cost2d.shape
    # safe bin index (clip to 0..N_BINS; 0 = unused sentinel)
    bi = np.clip(bin_a, 0, N_BINS).astype(np.int32)
    lb = _lat_band_idx(lat)

    for wv in WEIGHT_VERS:
        w2d = weights[wv]  # (n_prof, n_depths)
        # compute cost for this weight version: base misfit^2 × w_ver
        # 'final' weight already contains misfit^2 in cost2d; other versions scale it
        # → cost2d IS (misfit^2 * w_final); divide by appropriate factor
        if wv == 'final':
            c2d = cost2d
        elif wv == 'nosf':
            # nosf = final / sf  →  cost2d_nosf = cost2d_final / sf
            # But we passed the per-weight-version costs already. Actually cost2d is
            # (misfit^2 * w_final), and w_nosf = w_final/sf. So cost_nosf = misfit^2 * w_nosf.
            # We need misfit^2 from cost2d_final / w_final, then * w_nosf.
            # Simpler: pass cost arrays for all versions together. We'll recompute below.
            c2d = cost2d   # placeholder; replaced by caller that supplies per-version array
        else:
            c2d = cost2d

        k = (V, sp, wv)
        c_ok = np.where(ok, c2d, np.nan)

        # a02: mean cost vs depth
        np.add.at(acc[('a02', *k, 'sum')],   np.where(ok)[1], c2d[ok])
        np.add.at(acc[('a02', *k, 'count')], np.where(ok)[1], 1)

        # a03, a04: per-time aggregates
        for i in range(n_prof):
            ym = int(yyyymm_arr[i])
            ok_i = ok[i]
            c_i  = c2d[i, ok_i]
            n_i  = ok_i.sum()
            if n_i == 0:
                continue
            s = float(c_i.sum())
            acc[('a03', *k, 'sum')][ym]   += s
            acc[('a03', *k, 'count')][ym] += int(n_i)
            acc[('a04', *k, 'sum')][ym]   += s
            acc[('a04', *k, 'count')][ym] += int(n_i)

        # a05: integrated cost per bin
        #  sum all depth levels for each profile, then scatter by bin
        c_per_prof = np.where(ok, c2d, 0.0).sum(axis=1)
        count_per_prof = ok.sum(axis=1).astype(np.int64)
        np.add.at(acc[('a05', *k, 'sum')],   bi, c_per_prof)
        np.add.at(acc[('a05', *k, 'count')], bi, count_per_prof)

        # a07: zonal mean cost per datum vs depth
        rows_ok, cols_ok = np.where(ok)
        if rows_ok.size:
            flat_idx = lb[rows_ok] * n_depths + cols_ok
            np.add.at(acc[('a07', *k, 'sum')].ravel(),   flat_idx, c2d[ok])
            np.add.at(acc[('a07', *k, 'count')].ravel(), flat_idx, 1)

        # a08: histogram of log10 cost (only where c2d > 0)
        pos = c2d[ok & (c2d > 0)]
        if pos.size:
            hist, _ = np.histogram(np.log10(pos), bins=HIST_EDGES)
            acc[('a08', *k, 'hist')] += hist


def _accumulate_file(acc, fd, bias_T, bias_S, cost_type, split_mode):
    """Update accumulators for one file, one (cost_type, split_mode) combo."""
    n      = fd['n']
    n_d    = len(fd['depth'])
    ymd    = fd['ymd']
    lat    = fd['lat']
    gamma  = fd['gamma']
    bin_a  = fd['bin_a'] if fd.get('bin_a') is not None else np.zeros(n, dtype=np.int32)
    ym     = _yyyymm(ymd)

    # Only valid bin ids contribute to spatial analyses
    valid_bin = (bin_a >= 1) & (bin_a <= N_BINS)

    for V in VARS:
        val   = fd[f'val_{V}']
        clim  = fd[f'clim_{V}']
        w     = fd[f'w_{V}']
        sf    = fd[f'sf_{V}']
        if val is None or clim is None or w is None:
            continue
        # reference for this cost type
        if cost_type == 'vs_clim':
            ref = clim
        elif cost_type in ('vs_model', 'vs_run_b'):
            estim = fd.get(f'estim_{V}' if cost_type == 'vs_model' else f'estimb_{V}')
            if estim is None:
                continue
            ref = estim
        else:
            continue

        # valid mask
        fv  = _is_fill_value(val)
        fref = _is_fill_value(ref)
        fw   = _is_fill_weight(w)
        ok   = ~fv & ~fref & ~fw & (w > 0)

        # analysis 1: T/S mean vs depth (only once, during unsplit pass, independent of split)
        if split_mode == 'unsplit':
            if V == 'T':
                np.add.at(acc['a01_T_sum'],   np.where(ok)[1], val[ok])
                np.add.at(acc['a01_T_count'], np.where(ok)[1], 1)
            else:
                np.add.at(acc['a01_S_sum'],   np.where(ok)[1], val[ok])
                np.add.at(acc['a01_S_count'], np.where(ok)[1], 1)

        # bias accumulation (only during unsplit pass, for computing split biases)
        if split_mode == 'unsplit' and valid_bin.any():
            misfit = val - ref   # (n_prof, n_depths)
            bi_c = np.clip(bin_a, 0, N_BINS).astype(np.int32)
            rows_ok, cols_ok = np.where(ok & valid_bin[:, None])
            flat = bi_c[rows_ok] * n_d + cols_ok
            mf_flat = misfit[rows_ok, cols_ok]
            bs_key = f'bias_sum_{V}'
            bc_key = f'bias_count_{V}'
            np.add.at(acc[bs_key].ravel(), flat, mf_flat)
            np.add.at(acc[bc_key].ravel(), flat, 1)

        # misfit
        if split_mode == 'unsplit':
            mf = val - ref
        else:
            # split: subtract bias at each (bin, depth)
            bi_c = np.clip(bin_a, 0, N_BINS).astype(np.int32)
            mf = val - ref
            if V == 'T':
                b_arr = bias_T
            else:
                b_arr = bias_S
            if b_arr is not None:
                # broadcast bias for each profile's bin
                mf_adj = mf.copy()
                for i in range(n):
                    b = bi_c[i]
                    mf_adj[i] -= b_arr[b]   # shape (n_depths,)
                mf = mf_adj

        # three weight-version costs
        weights = _weight_versions(w, sf, gamma)
        costs_by_ver = {wv: mf ** 2 * weights[wv] for wv in WEIGHT_VERS}

        # for _accumulate_cost we pass the 'final' cost and let it use weight-version arrays
        # Actually, we need to call _accumulate_cost per weight version.  Do it here directly
        # to avoid the inner loop overhead of the helper function approach.
        bi_safe = np.clip(bin_a, 0, N_BINS).astype(np.int32)
        lb = _lat_band_idx(lat)

        for wv in WEIGHT_VERS:
            c2d = costs_by_ver[wv]
            k = (V, split_mode, wv)

            rows_ok, cols_ok = np.where(ok)

            # a02
            np.add.at(acc[('a02', *k, 'sum')],   cols_ok, c2d[ok])
            np.add.at(acc[('a02', *k, 'count')], cols_ok, 1)

            # a03, a04
            c_per_prof = np.where(ok, c2d, 0.0).sum(axis=1)
            n_per_prof = ok.sum(axis=1).astype(np.int64)
            for i in range(n):
                ymv = int(ym[i])
                ni  = int(n_per_prof[i])
                if ni == 0:
                    continue
                s = float(c_per_prof[i])
                acc[('a03', *k, 'sum')][ymv]   += s
                acc[('a03', *k, 'count')][ymv] += ni
                acc[('a04', *k, 'sum')][ymv]   += s
                acc[('a04', *k, 'count')][ymv] += ni

            # a05
            np.add.at(acc[('a05', *k, 'sum')],   bi_safe, c_per_prof)
            np.add.at(acc[('a05', *k, 'count')], bi_safe, n_per_prof)

            # a07
            if rows_ok.size:
                flat_idx = lb[rows_ok] * n_d + cols_ok
                np.add.at(acc[('a07', *k, 'sum')].ravel(),   flat_idx, c2d[ok])
                np.add.at(acc[('a07', *k, 'count')].ravel(), flat_idx, 1)

            # a08
            pos = c2d[ok & (c2d > 0)]
            if pos.size:
                hist, _ = np.histogram(np.log10(pos), bins=HIST_EDGES)
                acc[('a08', *k, 'hist')] += hist

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
                        'total':     s5[1:].copy(),    # bins 1..N_BINS
                        'per_datum': np.where(c5[1:] > 0, s5[1:] / c5[1:], np.nan),
                        'count':     c5[1:].copy().astype(np.int64),
                    }

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

    return res

# ---------------------------------------------------------------------------
# NetCDF writer
# ---------------------------------------------------------------------------

def _write_nc(res, source, run_ts, path, cost_types_used):
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
            'split_method': 'mean_bias_per_bin_depth',
            'split_note': (
                'split = Fortran-style anomaly cost: misfit minus its time-mean per (bin,depth). '
                'unsplit = raw (V-Vref)^2 * w.  Both are saved; bias arrays are in this file '
                'for reproducibility.  To recover the raw-misfit view use the unsplit variables.'
            ),
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
           long_name='Global mean T vs depth (all profiles; not cost-type dependent)')
        cv('a01_S_mean_depth', res['a01']['S']['mean'], ('depth',), 'f4',
           long_name='Global mean S vs depth')
        nc.createVariable('a01_T_count_depth', 'i8', ('depth',))[:] = res['a01']['T']['count']
        nc.createVariable('a01_S_count_depth', 'i8', ('depth',))[:] = res['a01']['S']['count']

        # biases
        for V in VARS:
            b_arr = res['bias'][V][1:]   # bins 1..N_BINS, shape (N_BINS, n_depths)
            cv(f'bias_{V}', b_arr, ('bin_a', 'depth'), 'f4',
               long_name=f'Time-mean (prof_{V} - ref) per (bin, depth) used for split cost')
            nc.createVariable(f'bias_count_{V}', 'i8', ('bin_a', 'depth'))[:]  = \
                res[f'bias_count_{V}'][1:]

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

    print(f'  wrote {path}')

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

_V_COLOR = {'T': 'tab:red', 'S': 'tab:blue'}
_SP_LS   = {'unsplit': '--', 'split': '-'}
_WV_COLOR = {'final': 'tab:blue', 'nosf': 'tab:orange', 'true': 'tab:red'}
# Time-series split/unsplit: slightly darker shade + circle markers for split
_TS_SP_COLOR = {
    ('T', 'unsplit'): 'tab:red',  ('T', 'split'): '#8B0000',
    ('S', 'unsplit'): 'tab:blue', ('S', 'split'): '#00008B',
}
_SP_MARKER = {'unsplit': '', 'split': ''}
# Only suppress depth levels with zero observations. The Mediterranean spike in XBT below
# ~2200 m (T ≈ 13 °C) is real data and is explained in the report's Analysis 1 note.
_MIN_DEPTH_COUNT = 1


def _savefig(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=110, bbox_inches='tight')
    plt.close(fig)


def _make_figures(res, source, fig_dir):
    """Generate all 10 analysis figures."""
    depth = res['depth']
    _dmax = float(np.nanmax(depth)) * 1.05  # bottom of every depth axis

    # Only plot variables that have valid data in this obs class (e.g. XBT has no S)
    active_vars = [V for V in VARS if np.isfinite(res['a01'][V]['mean']).any()]
    if not active_vars:
        return
    nv = len(active_vars)

    # ---- Analysis 1: T and S mean profile ----
    fig, axes = plt.subplots(1, nv, figsize=(4*nv + 2, 6), sharey=True, squeeze=False)
    axes = axes.flatten()
    for ax, V in zip(axes, active_vars):
        m = res['a01'][V]['mean']
        c = res['a01'][V]['count']
        m_plot = np.where(c >= _MIN_DEPTH_COUNT, m, np.nan)
        ax.plot(m_plot, depth, color=_V_COLOR[V])
        ax.set_xlabel(f'Mean {V} (°C)' if V == 'T' else f'Mean {V} (psu)')
        ax.set_ylabel('Depth (m)') if ax is axes[0] else None
        ax.set_title(f'{V} (n_max={c.max():,.0f})')
    axes[0].set_ylim(_dmax, 0)
    fig.suptitle(f'{source}: Analysis 1 — Global mean T & S vs depth')
    _savefig(fig, fig_dir / '01_mean_T_S_depth.png')

    # ---- Analysis 2/9: cost vs depth ----
    for an, wvs, title_sfx in [
        ('02', WEIGHT_VERS, 'mean cost'),
        ('09', ('true',),   'mean cost per datum'),
    ]:
        fig, axes = plt.subplots(1, nv, figsize=(5*nv + 1, 6), sharey=True, squeeze=False)
        axes = axes.flatten()
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
        axes[0].set_ylim(_dmax, 0)
        fig.suptitle(f'{source}: Analysis {an} — {title_sfx} vs depth')
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
            ym = np.asarray(r['yyyymm'], dtype=int)
            x = (ym // 100) + (ym % 100 - 0.5) / 12.0  # evenly-spaced decimal year
            ax.plot(x, r['mean_per_datum'],
                    color=_TS_SP_COLOR[(V, sp)], ls=_SP_LS[sp],
                    label=sp, lw=1.1)
        ax.set_ylabel(f'{V} mean cost/datum')
        if any_data:
            ax.legend(fontsize=7)
    axes[-1].set_xlabel('Year')
    fig.suptitle(f'{source}: Analysis 3 — Mean cost per datum vs time (true weights)')
    _savefig(fig, fig_dir / f'03_cost_per_datum_time_{source}.png')

    # ---- Analysis 4: total cost vs time ----
    fig, axes = plt.subplots(nv, 1, figsize=(10, 3.5*nv + 1), sharex=True, squeeze=False)
    axes = axes.flatten()
    for ax, V in zip(axes, active_vars):
        for wv in WEIGHT_VERS:
            r = res[('a04', V, 'unsplit', wv)]
            if len(r['yyyymm']) == 0:
                continue
            ym = np.asarray(r['yyyymm'], dtype=int)
            x = (ym // 100) + (ym % 100 - 0.5) / 12.0  # evenly-spaced decimal year
            ax.semilogy(x, np.where(np.array(r['total']) > 0, r['total'], np.nan),
                        color=_WV_COLOR[wv], label=wv, lw=1.1)
        ax.set_ylabel(f'{V} total cost')
        ax.legend(fontsize=7)
    axes[-1].set_xlabel('Year')
    fig.suptitle(f'{source}: Analysis 4 — Total cost vs time')
    _savefig(fig, fig_dir / f'04_total_cost_time_{source}.png')

    # ---- Analysis 5/6: spatial map ----
    try:
        _make_spatial_map(res, source, fig_dir, active_vars)
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
        fig, ax = plt.subplots(figsize=(9, 6))
        with np.errstate(invalid='ignore'):
            vmax = np.nanpercentile(r[valid], 98)
            vmax = vmax if np.isfinite(vmax) and vmax > 0 else 1.0
        im = ax.contourf(lat_centers, depth, r.T, levels=np.linspace(0, vmax, 21),
                         cmap='YlOrRd', extend='max')
        ax.set_ylim(_dmax, 0)
        plt.colorbar(im, ax=ax, label=f'{V} cost per datum (ΔV/σ)²')
        ax.set_xlabel('Latitude (°N)')
        ax.set_ylabel('Depth (m)')
        ax.set_title(f'{source}: Analysis 7 — Zonal mean {V} cost per datum (unsplit, true)')
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
    fig.suptitle(f'{source}: Analysis 8 — Histogram of cost per datum')
    _savefig(fig, fig_dir / f'08_histogram_{source}.png')

    # ---- Analysis 10: depth-range summary ----
    fig, axes = plt.subplots(1, nv, figsize=(5*nv + 1, 5), squeeze=False)
    axes = axes.flatten()
    x = np.arange(N_DEPTH_RANGES)
    width = 0.3
    for ax, V in zip(axes, active_vars):
        for i, sp in enumerate(SPLIT_MODES):
            vals = res[('a10', V, sp, 'true')]
            ax.bar(x + i * width, vals, width, label=sp, alpha=0.8)
        ax.set_xticks(x + width / 2)
        ax.set_xticklabels(DEPTH_RANGE_LABELS, rotation=30, ha='right', fontsize=8)
        ax.set_ylabel(f'{V} mean cost per datum (ΔV/σ)²')
        ax.set_title(V)
        ax.legend(fontsize=7)
    fig.suptitle(f'{source}: Analysis 10 — Mean cost per datum by depth range')
    _savefig(fig, fig_dir / f'10_depth_ranges_{source}.png')


def _make_spatial_map(res, source, fig_dir, active_vars):
    """Analyses 5 and 6: spatial maps using geodesic bins."""
    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
    except ImportError:
        print('  [warn] cartopy not available; skipping spatial maps')
        return

    bl = common.BinLookup()
    raster, land, outside, extent = bl.projected_raster('global', npix=600)
    grid_a = raster['a']   # shape (ny, nx), values 0..10242

    for V in active_vars:
        for an, field, sp, wv, title in [
            ('05', 'total',     'unsplit', 'final', 'total cost'),
            ('06', 'per_datum', 'unsplit', 'true',  'cost per datum'),
        ]:
            k = (V, sp, wv)
            vals = res[('a05', *k)][field]  # (N_BINS,) for bins 1..N_BINS
            # map: per_bin[bin-1] -> raster pixel; bin_id 0 → 0 (unused)
            per_bin_full = np.concatenate([[0.0], vals])   # index 0 unused
            img = per_bin_full[grid_a]                      # (ny, nx)
            img = np.where(land | outside | (grid_a == 0), np.nan, img)

            with np.errstate(invalid='ignore'):
                vmax = np.nanpercentile(img[np.isfinite(img)], 98)
                if not (np.isfinite(vmax) and vmax > 0):
                    vmax = 1.0

            fig = plt.figure(figsize=(12, 5))
            ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
            ax.add_feature(cfeature.LAND, color='lightgray', zorder=1)
            ax.set_global()
            ax.imshow(img, origin='lower', extent=extent,
                      transform=ccrs.Robinson(), cmap='YlOrRd',
                      vmin=0, vmax=vmax, zorder=0)
            plt.colorbar(
                plt.cm.ScalarMappable(
                    norm=mcolors.Normalize(0, vmax), cmap='YlOrRd'),
                ax=ax, orientation='horizontal', pad=0.04, fraction=0.04,
                label=f'{V} {title} ({sp}, {wv})')
            ax.set_title(f'{source}: Analysis {an} — {V} {title} per geodesic bin')
            _savefig(fig, fig_dir / f'{an}_spatial_{V}_{source}.png')

# ---------------------------------------------------------------------------
# Per-source analysis
# ---------------------------------------------------------------------------

def _analyze_source(source, files, run_dir, out_dir, estim_dir, estim_dir_b, cost_types):
    """Run both passes for one obs class; write NetCDF and figures."""
    t0 = _time.time()
    n_files = len(files)
    print(f'  {source}: {n_files} files')

    # determine depth from the first readable file
    depth = None
    for f in files:
        fd = _load_file(f)
        if fd is not None and fd.get('depth') is not None:
            depth = fd['depth'].astype(np.float32)
            break
    if depth is None:
        print(f'  {source}: no readable files, skipped')
        return
    n_d = len(depth)

    est_index = EstimIndex(estim_dir, run_dir) if estim_dir else None
    est_index_b = EstimIndex(estim_dir_b, run_dir) if estim_dir_b else None
    if est_index is not None:
        found = sum(1 for f in files if (Path(estim_dir) / (Path(f).stem + '_model.nc')).exists() or Path(f).name in est_index.map)
        print(f'    model equivalents ({estim_dir}): for {found} of {n_files} files')
        if found == 0:
            print(f'    [warn] no model equivalents match these files; vs_model is skipped')

    cost_types_used = [ct for ct in cost_types
                       if ct == 'vs_clim'
                       or (ct == 'vs_model' and estim_dir)
                       or (ct == 'vs_run_b' and estim_dir_b)]

    # Accumulators: one per cost_type
    accs = {ct: _new_acc(n_d) for ct in cost_types_used}

    # ---- Pass 1: unsplit costs + bias accumulation ----
    print(f'    pass 1 (unsplit + bias)…', end='', flush=True)
    for i, f in enumerate(files):
        fd = _load_file(f, estim_index=est_index)
        if fd is not None and est_index_b is not None:
            fb = est_index_b.get(f, fd['n'], fd['depth'].size)
            for V in VARS:
                e = None if fb is None else fb[V]
                fd[f'estimb_{V}'] = None if e is None else np.where(e == -9999, np.nan, e)
        if fd is None:
            continue
        for ct in cost_types_used:
            _accumulate_file(accs[ct], fd, None, None, ct, 'unsplit')
        if (i + 1) % 50 == 0:
            print(f' {i+1}', end='', flush=True)
    print()

    # Compute biases from pass 1
    biases = {}
    for ct in cost_types_used:
        acc = accs[ct]
        biases[ct] = {}
        for V in VARS:
            s = acc[f'bias_sum_{V}']
            c = acc[f'bias_count_{V}']
            with np.errstate(invalid='ignore'):
                biases[ct][V] = np.where(c > 0, s / c, 0.0)   # shape (N_BINS+1, n_d)

    # ---- Pass 2: split costs ----
    print(f'    pass 2 (split)…', end='', flush=True)
    for i, f in enumerate(files):
        fd = _load_file(f, estim_index=est_index)
        if fd is not None and est_index_b is not None:
            fb = est_index_b.get(f, fd['n'], fd['depth'].size)
            for V in VARS:
                e = None if fb is None else fb[V]
                fd[f'estimb_{V}'] = None if e is None else np.where(e == -9999, np.nan, e)
        if fd is None:
            continue
        for ct in cost_types_used:
            _accumulate_file(accs[ct], fd,
                             biases[ct]['T'], biases[ct]['S'],
                             ct, 'split')
        if (i + 1) % 50 == 0:
            print(f' {i+1}', end='', flush=True)
    print()

    # Finalize and write
    for ct in cost_types_used:
        res = _finalize(accs[ct], depth)
        run_ts_name = common.run_timestamp(run_dir)
        fname = f'analysis_{source}_{run_ts_name}_{ct}.nc'
        out_nc = out_dir / source / fname
        _write_nc(res, source, run_ts_name, out_nc, [ct])
        fig_dir = out_dir / source / 'figures' / ct
        _make_figures(res, source, fig_dir)
        print(f'    {ct}: done  ({_time.time() - t0:.0f}s)')

# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def run(run_dir, out_dir, sources=None, estim_dir=None, estim_dir_b=None, n_workers=None):
    """Compute all standard analyses for every obs class (source) in ``run_dir``.

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
        Parallel workers (default 1; parallelism is over sources).

    Returns
    -------
    Path
        The output directory.
    """
    run_dir = common.yearly_dir(run_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cost_types = ['vs_clim']
    if estim_dir:
        cost_types.append('vs_model')
    if estim_dir_b:
        cost_types.append('vs_run_b')

    # group files by source
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

    for src, files in sorted(groups.items()):
        _analyze_source(src, files, run_dir, out_dir, estim_dir, estim_dir_b, cost_types)

    meta = {
        'run': str(run_dir), 'out': str(out_dir),
        'sources': sorted(groups.keys()),
        'cost_types': cost_types,
        'estim_dir': str(estim_dir) if estim_dir else None,
        'estim_dir_b': str(estim_dir_b) if estim_dir_b else None,
        'created': datetime.datetime.now().isoformat(timespec='seconds'),
    }
    (out_dir / 'meta.json').write_text(json.dumps(meta, indent=2))
    print(f'All done in {_time.time() - t0:.0f}s → {out_dir}')
    return out_dir


def cli(argv=None):
    """CLI entry point (``ecco-insitu profile-analyze``)."""
    p = argparse.ArgumentParser(
        description='Compute profile analyses (T/S profiles, cost distributions) '
                    'for each obs class in a processed run.')
    p.add_argument('--run', required=True,
                   help='Processed run directory (or its yearly/ sub-folder)')
    p.add_argument('--out', required=True, help='Output directory')
    p.add_argument('--sources', nargs='*', default=None,
                   help='Obs classes to process (default: all)')
    p.add_argument('--model-equivalents', default=None, dest='estim_dir',
                   help='Dir with <name>_model.nc files for vs-model cost')
    p.add_argument('--run-b', default=None, dest='estim_dir_b',
                   help='Second model-equivalent dir for vs-run-B cost')
    p.add_argument('-n', type=int, default=None,
                   help='Parallel workers (default: 1)')
    a = p.parse_args(argv)
    run(a.run, a.out, a.sources, a.estim_dir, a.estim_dir_b, a.n)
