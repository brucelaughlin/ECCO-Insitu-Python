"""
The assistant's tools: fixed, checked computations on one profile-analyze folder, each returning a few summary
numbers and links to the site's views that show them. The model only chooses tools and their settings; every
number in an answer comes from here (docs/planning/llm_assistant.md).

Conventions (as the report): costs are per datum with the true weights (1/sigma^2, no gamma or scaling factor);
'split' is the split with each class's own bin means (``CORE_SPLIT = 'own'``); misfits are reference - obs.
"""
import json
from functools import lru_cache
from pathlib import Path

import netCDF4
import numpy as np

from ecco_insitu.diagnostics.profile_analyses import compute as pa

PARTS = ('unsplit', 'anomaly', 'mean')     # unsplit; the split's two terms
VARS = pa.VARS
PEAK_CLUSTER_KM = 600                       # as the explorer's peak tool
HIGH_Q = 0.9


class ToolError(ValueError):
    """A tool was called with settings the analysis doesn't have; the message goes back to the model."""


class Analysis:
    """One profile-analyze folder: its classes, cost types and saved NetCDF files (read on demand)."""

    def __init__(self, analysis_dir, bin_csv=None):
        self.dir = Path(analysis_dir)
        self.meta = json.loads((self.dir / 'meta.json').read_text())
        self.sources = list(self.meta['sources'])
        self.cost_types = list(self.meta['cost_types'])
        self.bin_ll = np.loadtxt(bin_csv, delimiter=',')[:, :2] if bin_csv and Path(bin_csv).exists() else None
        self._open = {}

    def nc(self, src, ct):
        """The open analysis NetCDF of one obs class and cost type (ToolError for unknown names)."""
        if src not in self.sources:
            raise ToolError(f'unknown obs class {src!r}; the classes are {", ".join(self.sources)}')
        if ct not in self.cost_types:
            raise ToolError(f'unknown cost type {ct!r}; the cost types are {", ".join(self.cost_types)}')
        if (src, ct) not in self._open:
            f = sorted((self.dir / src).glob(f'analysis_{src}_*_{ct}.nc'))
            if not f:
                raise ToolError(f'no analysis file for {src} {ct}')
            ds = netCDF4.Dataset(f[-1])
            ds.set_auto_mask(False)
            self._open[(src, ct)] = ds
        return self._open[(src, ct)]

    @property
    def depth(self):
        """The depth levels (m)."""
        ds = self.nc(self.sources[0], self.cost_types[0])
        return np.asarray(ds['depth_m'][:], float)

    def has(self, src, ct, name):
        """Whether the class's file has variable ``name``."""
        return name in self.nc(src, ct).variables

    def get(self, src, ct, name):
        """Variable ``name`` of the class's file as an array (ToolError if it is missing)."""
        ds = self.nc(src, ct)
        if name not in ds.variables:
            raise ToolError(f'{src} {ct} has no {name} (e.g. no S for XBT)')
        return np.asarray(ds[name][:])


def _var(V):
    if V not in VARS:
        raise ToolError(f'variable must be T or S, not {V!r}')
    return V


def _part(part):
    if part not in PARTS:
        raise ToolError(f'part must be one of {", ".join(PARTS)}, not {part!r}')
    return part


def _pfx(V, part):
    return f'a__{V}__{"unsplit" if part == "unsplit" else "split_own"}__true'


def _levels(depth, z0, z1):
    lo, hi = (z0, z1) if z0 <= z1 else (z1, z0)
    k = np.nonzero((depth >= lo) & (depth <= hi))[0]
    if not len(k):
        raise ToolError(f'no depth levels between {lo} and {hi} m (levels run {depth[0]:g}-{depth[-1]:g} m)')
    return k


def _r(x, n=3):
    """Round for the model: n significant digits, NaN -> None."""
    if x is None or not np.isfinite(x):
        return None
    return float(f'{float(x):.{n}g}')


def _cost_depth(A, src, ct, V, part):
    """Cost per datum at each level (as analysis 02); the mean term is the mean over bins."""
    if part == 'mean':
        return (A.get(src, ct, f'meanterm__{V}__true__mean_depth').astype(float),
                A.get(src, ct, f'meanterm__{V}__true__nbins_depth').astype(float))
    p = _pfx(V, part)
    return (A.get(src, ct, f'a02__{p}__mean_depth').astype(float),
            A.get(src, ct, f'a02__{p}__count_depth').astype(float))


# ---------------------------------------------------------------------------------------------------------------
# the tools. Each returns {'result': ..., 'views': [{'page', 'params', 'caption'}]} (views become links)
# ---------------------------------------------------------------------------------------------------------------

def list_contents(A):
    """What the analysis has: classes, cost types, variables, depth levels, years."""
    out = {'obs_classes': A.sources, 'cost_types': A.cost_types, 'variables': list(VARS),
           'depth_levels_m': [_r(z, 4) for z in A.depth], 'parts': list(PARTS)}
    ds = A.nc(A.sources[0], A.cost_types[0])
    out['years'] = f'{ds.viewer_year0}-{ds.viewer_year0 + ds.viewer_n_blocks - 1}'
    return {'result': out, 'views': []}


def cost_vs_depth(A, src, cost_type, variable='T', part='unsplit', depth_from=None, depth_to=None):
    """The cost per datum against depth (analysis 02): the levels, their values, the largest value and where."""
    V, part = _var(variable), _part(part)
    z = A.depth
    v, n = _cost_depth(A, src, cost_type, V, part)
    k = _levels(z, depth_from if depth_from is not None else 0, depth_to if depth_to is not None else z[-1])
    vals = v[k]
    if not np.isfinite(vals).any():
        raise ToolError(f'no {part} values for {src} {cost_type} {V} in that depth range')
    kmax = k[int(np.nanargmax(vals))]
    step = max(1, len(k) // 25)                       # at most ~25 levels back to the model
    rows = [[_r(z[i], 4), _r(v[i]), int(n[i])] for i in k[::step] if np.isfinite(v[i])]
    res = {'obs_class': src, 'cost_type': cost_type, 'variable': V, 'part': part,
           'columns': ['depth_m', 'cost_per_datum' if part != 'mean' else 'mean_term_mean_over_bins',
                       'data' if part != 'mean' else 'bins'],
           'levels': rows, 'levels_shown': f'{len(rows)} of {len(k)} (every {step})' if step > 1 else 'all',
           'largest': {'depth_m': _r(z[kmax], 4), 'value': _r(v[kmax])},
           'mean_over_range': _r(np.nansum(v[k] * n[k]) / np.nansum(n[k])) if part != 'mean' else _r(np.nanmean(v[k]))}
    return {'result': res, 'views': [
        {'page': 'report', 'params': {'src': src, 'ct': cost_type, 'an': '02'},
         'caption': f'Analysis 02: {src} {cost_type} cost against depth'},
        {'page': 'explorer', 'params': {'src': src, 'ct': cost_type, 'V': V, 'part': part, 'depth': f'k{kmax}', 'years': 'all'},
         'caption': f'Cost map explorer: {part} at {z[kmax]:g} m (the largest level)'}]}


def cost_vs_time(A, src, cost_type, variable='T', part='unsplit'):
    """Cost per datum per month (analysis 03): yearly means, the highest and lowest years, the trend."""
    V, part = _var(variable), _part(part)
    if part == 'mean':
        raise ToolError('the mean term has no time dependence (it is defined over the whole record); use unsplit or anomaly')
    p = _pfx(V, part)
    m = A.get(src, cost_type, f'a03__{p}__mean_per_datum_time').astype(float)
    c = A.get(src, cost_type, f'a03__{p}__count_time').astype(float)
    ym = A.get(src, cost_type, 'yyyymm')
    ok = np.isfinite(m) & (c > 0)
    yr = ym[ok] // 100
    years = np.unique(yr)
    per = [(int(y), np.sum(m[ok][yr == y] * c[ok][yr == y]) / np.sum(c[ok][yr == y]), int(np.sum(c[ok][yr == y]))) for y in years]
    ys = np.array([p_[0] for p_ in per]); vs = np.array([p_[1] for p_ in per]); ns = np.array([p_[2] for p_ in per])
    trend = np.polyfit(ys, vs, 1, w=np.sqrt(ns))[0] if len(ys) > 2 else np.nan
    srt = np.argsort(vs)
    res = {'obs_class': src, 'cost_type': cost_type, 'variable': V, 'part': part,
           'columns': ['year', 'cost_per_datum', 'data'],
           'years': [[int(y), _r(v), int(n)] for y, v, n in per],
           'highest_years': [int(ys[i]) for i in srt[::-1][:3]], 'lowest_years': [int(ys[i]) for i in srt[:3]],
           'trend_per_decade': _r(10 * trend), 'note': 'trend: least squares on the yearly means, weighted by sqrt(data)'}
    return {'result': res, 'views': [
        {'page': 'report', 'params': {'src': src, 'ct': cost_type, 'an': '03'},
         'caption': f'Analysis 03: {src} {cost_type} cost per datum against time'}]}


def compare_classes(A, cost_type, variable='T', part='unsplit', depth_from=0, depth_to=2000):
    """The cost per datum of every obs class over one depth range, side by side."""
    V, part = _var(variable), _part(part)
    z = A.depth
    k = _levels(z, depth_from, depth_to)
    rows = []
    for src in A.sources:
        try:
            v, n = _cost_depth(A, src, cost_type, V, part)
        except ToolError:
            continue
        if part == 'mean':
            val, cnt = np.nanmean(v[k]), int(np.nansum(n[k]))
        else:
            cnt = np.nansum(n[k]); val = np.nansum(v[k] * n[k]) / cnt if cnt else np.nan
        rows.append([src, _r(val), int(cnt)])
    rows.sort(key=lambda r: -(r[1] or 0))
    return {'result': {'cost_type': cost_type, 'variable': V, 'part': part, 'depth_range_m': [depth_from, depth_to],
                       'columns': ['obs_class', 'cost_per_datum' if part != 'mean' else 'mean_term', 'data' if part != 'mean' else 'bin_levels'],
                       'classes': rows},
            'views': [{'page': 'report', 'params': {'src': rows[0][0] if rows else A.sources[0], 'ct': cost_type, 'an': '02'},
                       'caption': 'Analysis 02 (switch class with the Source tabs)'}]}


def compare_references(A, src, variable='T', part='unsplit', depth_from=0, depth_to=2000):
    """Climatology vs model: the cost per datum of one class over one depth range against each reference."""
    V, part = _var(variable), _part(part)
    k = _levels(A.depth, depth_from, depth_to)
    rows = []
    for ct in A.cost_types:
        v, n = _cost_depth(A, src, ct, V, part)
        val = np.nanmean(v[k]) if part == 'mean' else np.nansum(v[k] * n[k]) / np.nansum(n[k])
        rows.append([ct, _r(val)])
    return {'result': {'obs_class': src, 'variable': V, 'part': part, 'depth_range_m': [depth_from, depth_to],
                       'columns': ['cost_type', 'value'], 'references': rows},
            'views': [{'page': 'report', 'params': {'src': src, 'ct': ct, 'an': '02'},
                       'caption': f'Analysis 02, {src} {ct}'} for ct in A.cost_types]}


def _bin_cost(A, src, ct, V, part, k):
    """Per bin over the levels k: (value, weight). Costs: summed cost / data over all years (bins weighted by
    data, as analysis 02); mean term: the mean over the levels with a value."""
    if part == 'mean':
        mt = A.get(src, ct, f'meanterm__{V}__true__by_bin_depth').astype(float)[:, k]
        n = np.sum(np.isfinite(mt), axis=1).astype(float)
        with np.errstate(invalid='ignore', divide='ignore'):
            return np.where(n > 0, np.nansum(mt, axis=1) / n, np.nan), n
    ds = A.nc(src, ct)
    core = 'unsplit' if part == 'unsplit' else 'split_own'
    cell = np.asarray(ds[f'viewer__{V}__{core}__cell'][:], np.int64)
    s = np.asarray(ds[f'viewer__{V}__{core}__sum'][:], float)
    n = np.asarray(ds[f'viewer__{V}__{core}__n'][:], float)
    nd = len(A.depth)
    lev, b = cell % nd, (cell // nd) % (pa.N_BINS + 1)
    keep = np.isin(lev, k)
    tot = np.bincount(b[keep], s[keep], pa.N_BINS + 1)[1:]
    cnt = np.bincount(b[keep], n[keep], pa.N_BINS + 1)[1:]
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(cnt > 0, tot / cnt, np.nan), cnt


def _where(A, b):
    if A.bin_ll is None:
        return None
    lo, la = A.bin_ll[b - 1]
    lo = lo - 360 if lo > 180 else lo
    return f'{abs(la):.0f}°{"N" if la >= 0 else "S"} {abs(lo):.0f}°{"E" if lo >= 0 else "W"}'


def _hav(la1, lo1, la2, lo2):
    r = np.pi / 180
    a = np.sin((la2 - la1) * r / 2) ** 2 + np.cos(la1 * r) * np.cos(la2 * r) * np.sin((lo2 - lo1) * r / 2) ** 2
    return 12742 * np.arcsin(np.sqrt(a))


def where_high(A, src, cost_type, variable='T', part='unsplit', depth_from=0, depth_to=6000, top=10):
    """Where the cost is high over a depth range: the highest bins, and the top 10 % of bins grouped into areas
    (bins within 600 km chained together, as the explorer's peak tool), with each area's share of the summed cost."""
    V, part = _var(variable), _part(part)
    z = A.depth
    k = _levels(z, depth_from, depth_to)
    val, w = _bin_cost(A, src, cost_type, V, part, k)
    ok = np.nonzero(np.isfinite(val) & (w > 0))[0]
    if not len(ok):
        raise ToolError('no bins with data in that depth range')
    contrib = val[ok] * (1 if part == 'mean' else w[ok])           # additive part of the level's value
    order = ok[np.argsort(-val[ok])]
    thr = np.quantile(val[ok], HIGH_Q)
    high = [b for b in order if val[b] >= thr]
    clusters, seen = [], set()
    if A.bin_ll is not None:
        ll = A.bin_ll
        for b in high:
            if b in seen:
                continue
            c, i = [b], 0
            seen.add(b)
            while i < len(c):
                la, lo = ll[c[i], 1], ll[c[i], 0]
                for o in high:
                    if o not in seen and _hav(la, lo, ll[o, 1], ll[o, 0]) <= PEAK_CLUSTER_KM:
                        seen.add(o); c.append(o)
                i += 1
            clusters.append(c)
    total = contrib.sum()
    cmap = dict(zip(ok, contrib))
    cl = sorted(({'bins': len(c), 'center_of_highest': _where(A, c[0] + 1), 'share_of_bins_pct': _r(100 * len(c) / len(ok)),
                  'share_of_summed_cost_pct': _r(100 * sum(cmap[b] for b in c) / total), 'bin_ids': [int(b + 1) for b in c[:40]]}
                 for c in clusters), key=lambda d: -d['share_of_summed_cost_pct'])[:6]
    k_mid = int(k[len(k) // 2]) if len(k) > 1 else int(k[0])
    depth_param = _group_param(A, k)
    res = {'obs_class': src, 'cost_type': cost_type, 'variable': V, 'part': part,
           'depth_range_m': [_r(z[k[0]], 4), _r(z[k[-1]], 4)], 'bins_with_data': int(len(ok)),
           'highest_bins': [{'bin': int(b + 1), 'where': _where(A, b + 1), 'value': _r(val[b]), 'data': int(w[b])}
                            for b in order[:top]],
           'areas_of_the_top_10pct': cl,
           'note': 'bins are geodesic, ~220 km across; costs over all years'}
    region = '.'.join(str(b) for b in cl[0]['bin_ids']) if cl else ''
    views = [{'page': 'explorer', 'params': {'src': src, 'ct': cost_type, 'V': V, 'part': part, 'depth': depth_param,
                                             'years': 'all', **({'region': region} if region else {})},
              'caption': f'Cost map explorer: {part} over {z[k[0]]:g}–{z[k[-1]]:g} m' + (' (largest area outlined)' if region else '')}]
    if A.bin_ll is not None:
        b = order[0]
        views.append({'page': 'data', 'params': {'src': src, 'ref': 'clim' if cost_type == 'vs_clim' else cost_type[3:], 'V': V,
                                                 'bin': int(b + 1), 'lon': _r(A.bin_ll[b, 0], 5), 'lat': _r(A.bin_ll[b, 1], 5), 'k': k_mid},
                      'caption': f'In situ data viewer: the profiles of bin {b + 1} (the highest)'})
    return {'result': res, 'views': views}


def _group_param(A, k):
    """The explorer's depth setting for levels ``k``: one level 'k<n>', all levels 'all', else the range 'r<a>-<b>'."""
    k = sorted(int(x) for x in k)
    if len(k) == 1:
        return f'k{k[0]}'
    if len(k) == len(A.depth):
        return 'all'
    return f'r{k[0]}-{k[-1]}'


def offset_vs_variability(A, src, cost_type, variable='T', depth_from=0, depth_to=2000):
    """Is the misfit a steady offset or variability in time? Compares the split's mean term with its anomaly term,
    and gives the mean misfit (reference - obs) over the range."""
    V = _var(variable)
    k = _levels(A.depth, depth_from, depth_to)
    u, nu = _cost_depth(A, src, cost_type, V, 'unsplit')
    a, na = _cost_depth(A, src, cost_type, V, 'anomaly')
    m, nb = _cost_depth(A, src, cost_type, V, 'mean')
    res = {'obs_class': src, 'cost_type': cost_type, 'variable': V, 'depth_range_m': [depth_from, depth_to],
           'unsplit_per_datum': _r(np.nansum(u[k] * nu[k]) / np.nansum(nu[k])),
           'anomaly_term_per_datum': _r(np.nansum(a[k] * na[k]) / np.nansum(na[k])),
           'mean_term_mean_over_bins': _r(np.nanmean(m[k])),
           'note': 'the anomaly term measures variability in time (offsets removed); the mean term the steady offset '
                   'per bin, weighted by the bin\'s combined error, so it is not on the per-datum scale: compare its '
                   'size across depths and classes, not with the anomaly term directly'}
    name = f'bias_own_{V}'
    if A.has(src, cost_type, name):
        bias = -A.get(src, cost_type, name).astype(float)[:, k]          # saved as obs - ref
        cnt = A.get(src, cost_type, f'bias_own_count_{V}').astype(float)[:, k]
        ok = np.isfinite(bias) & (cnt > 0)
        res['mean_misfit_ref_minus_obs'] = _r(np.sum(np.where(ok, bias * cnt, 0)) / np.sum(np.where(ok, cnt, 0)))
        res['mean_misfit_units'] = 'degC' if V == 'T' else 'psu'
    return {'result': res, 'views': [
        {'page': 'report', 'params': {'src': src, 'ct': cost_type, 'an': '02'}, 'caption': 'Analysis 02: unsplit and the two split terms'},
        {'page': 'explorer', 'params': {'src': src, 'ct': cost_type, 'V': V, 'part': 'bias', 'depth': f'k{k[len(k) // 2]}', 'years': 'all'},
         'caption': 'Cost map explorer: the mean misfit (ref − obs) per bin, red–blue'}]}


def peak_attribution(A, src, cost_type, variable='T', part='mean', depth_from=2000, depth_to=2600):
    """Where a bump in the cost against depth comes from: the peak level in the range, how much of its excess over
    the neighboring levels comes from a different set of places being averaged at that depth, and the areas that carry it."""
    V, part = _var(variable), _part(part)
    z = A.depth
    v, n = _cost_depth(A, src, cost_type, V, part)
    k = [i for i in _levels(z, depth_from, depth_to) if np.isfinite(v[i])]
    if not k:
        raise ToolError('no values in that depth range')
    kp = max(k, key=lambda i: v[i])
    have = [i for i in range(len(z)) if np.isfinite(v[i])]
    ia = have.index(kp)
    nbrs = [have[j] for j in (ia - 1, ia + 1) if 0 <= j < len(have)]
    vals = {i: _bin_cost(A, src, cost_type, V, part, [i]) for i in [kp, *nbrs]}

    def level_value(i, keep=None):
        val, w = vals[i]
        ok = np.isfinite(val) & (w > 0)
        if keep is not None:
            ok &= keep
        if not ok.any():
            return np.nan
        return float(np.mean(val[ok])) if part == 'mean' else float(np.sum(val[ok] * w[ok]) / np.sum(w[ok]))

    base = np.nanmean([level_value(i) for i in nbrs]) if nbrs else np.nan
    keep = np.ones(pa.N_BINS, bool)
    for i in [kp, *nbrs]:
        keep &= np.isfinite(vals[i][0]) & (vals[i][1] > 0)
    sp, sb = level_value(kp, keep), np.nanmean([level_value(i, keep) for i in nbrs]) if nbrs else np.nan
    ex, sex = level_value(kp) - base, sp - sb
    sampling = float(np.clip(1 - sex / ex, 0, 1)) if np.isfinite(ex) and ex > 0 and np.isfinite(sex) else np.nan
    where = where_high(A, src, cost_type, V, part, z[kp], z[kp], top=5)
    res = {'obs_class': src, 'cost_type': cost_type, 'variable': V, 'part': part,
           'peak': {'depth_m': _r(z[kp], 4), 'value': _r(level_value(kp)), 'bins': int(np.sum(np.isfinite(vals[kp][0])))},
           'neighbor_levels': [{'depth_m': _r(z[i], 4), 'value': _r(level_value(i)), 'bins': int(np.sum(np.isfinite(vals[i][0])))} for i in nbrs],
           'same_places_only': {'places': int(keep.sum()), 'peak': _r(sp), 'depths_around_it_average': _r(sb)},
           'share_of_jump_from_a_different_set_of_places_pct': _r(100 * sampling) if np.isfinite(sampling) else None,
           'share_of_jump_from_a_worse_fit_at_that_depth_pct': _r(100 * (1 - sampling)) if np.isfinite(sampling) else None,
           'areas_at_peak': where['result']['areas_of_the_top_10pct'][:4],
           'method': 'Each depth\'s value is an average over the places (bins) with data at that depth, and that list of places '
                     'changes with depth (most Argo floats stop near 2000 m). The jump is the peak minus the average of the '
                     'two depths around it. It is computed twice: over all places, and over only the places with data at '
                     'all three depths. share from a different set of places = 1 - (jump over the same places) / (jump over '
                     'all places). A high share: the jump comes from averaging a different set of places at the peak depth '
                     'than at the depths around it, not from a worse fit at that depth in the same places'}
    views = [{'page': 'explorer', 'params': {'src': src, 'ct': cost_type, 'V': V, 'part': part,
                                             'peak': f'{z[min(k)]:g}-{z[max(k)]:g}'},
              'caption': f'Cost map explorer, peak tool on {z[min(k)]:g}–{z[max(k)]:g} m (outlines the largest area)'}]
    return {'result': res, 'views': views + where['views'][1:]}


def observation_counts(A, src, variable='T', depth_from=0, depth_to=6000):
    """How many observations there are: per year and over a depth range (counts of values with weight > 0)."""
    V = _var(variable)
    ct = A.cost_types[0]
    k = _levels(A.depth, depth_from, depth_to)
    if A.has(src, ct, f'tdc_{V}_count'):
        ym = A.get(src, ct, f'tdc_{V}_yyyymm'); c = A.get(src, ct, f'tdc_{V}_count').astype(float)[:, k].sum(1)
    else:
        ym = A.get(src, ct, 'yyyymm'); c = A.get(src, ct, f'a03__a__{V}__unsplit__true__count_time').astype(float)
    yr = ym // 100
    per = [[int(y), int(c[yr == y].sum())] for y in np.unique(yr) if c[yr == y].sum() > 0]
    tot = A.get(src, ct, f'a01_{V}_count_depth').astype(float)
    return {'result': {'obs_class': src, 'variable': V, 'depth_range_m': [depth_from, depth_to],
                       'total_values_in_range': int(tot[k].sum()), 'columns': ['year', 'values'], 'per_year': per,
                       'busiest_year': max(per, key=lambda r: r[1])[0] if per else None},
            'views': [{'page': 'report', 'params': {'src': src, 'ct': ct, 'an': '01'},
                       'caption': f'Analysis 01: {src} mean T and S and counts against depth'}]}


# ---------------------------------------------------------------------------------------------------------------
# regions: a named sea or ocean area (Natural Earth's marine areas), a common ocean region (the boxes below), or a
# lon/lat box; a region is the set of geodesic bins whose centers fall in it
# ---------------------------------------------------------------------------------------------------------------
# common ocean regions that are not one Natural Earth area: (name, aliases, lon0, lon1, lat0, lat1); lon0 > lon1
# crosses the dateline. Coarse boxes: descriptive only.
OCEAN_REGIONS = [
    ('eastern tropical Pacific', ('east pacific', 'eastern pacific', 'east tropical pacific'), -140, -70, -20, 20),
    ('western tropical Pacific', ('west pacific', 'western pacific', 'warm pool'), 120, -170, -20, 20),
    ('tropical Pacific', ('equatorial pacific',), 120, -70, -20, 20),
    ('North Pacific', ('north pacific',), 120, -100, 20, 66),
    ('South Pacific', ('south pacific',), 145, -70, -60, -20),
    ('subpolar North Atlantic', ('subpolar atlantic', 'north atlantic subpolar gyre'), -65, -5, 45, 66),
    ('subtropical North Atlantic', ('subtropical atlantic', 'sargasso'), -80, -15, 15, 45),
    ('North Atlantic', ('north atlantic',), -100, 20, 0, 70),
    ('tropical Atlantic', ('equatorial atlantic',), -60, 15, -15, 15),
    ('South Atlantic', ('south atlantic',), -70, 20, -60, 0),
    ('Indian Ocean', ('indian',), 20, 120, -60, 30),
    ('Southern Ocean', ('southern', 'antarctic', 'acc', 'antarctic circumpolar current'), -180, 180, -90, -50),
    ('Arctic', ('arctic ocean',), -180, 180, 66, 90),
    ('Gulf Stream region', ('gulf stream',), -80, -45, 32, 45),
    ('Kuroshio region', ('kuroshio', 'kuroshio extension'), 130, 175, 25, 42),
    ('Agulhas region', ('agulhas',), 10, 40, -45, -30),
]


def _norm(t):
    import re
    return re.sub(r'[^a-z0-9 ]', ' ', str(t).lower().replace('the ', ' ')).split()


def _in_box(lo, la, lon0, lon1, lat0, lat1):
    lo = np.where(lo > 180, lo - 360, lo)
    inlon = (lo >= lon0) & (lo <= lon1) if lon0 <= lon1 else (lo >= lon0) | (lo <= lon1)
    return inlon & (la >= lat0) & (la <= lat1)


@lru_cache(maxsize=1)
def _bin_area_names():
    """Each bin's area name (index = bin id - 1; '' if none), as the cost map explorer has them."""
    from ecco_insitu.diagnostics.profile_analyses.viewer import _bin_names
    return _bin_names() or []


def find_region(A, name=None, lon_min=None, lon_max=None, lat_min=None, lat_max=None):
    """(label, bin ids 1-based) of a region: a box if lon/lat are given, else the best match of ``name`` among the
    common ocean regions and Natural Earth's named marine areas."""
    if A.bin_ll is None:
        raise ToolError('no bin locations on this machine (paths.geodesic_bin_csv)')
    lo, la = A.bin_ll[:, 0], A.bin_ll[:, 1]
    if None not in (lon_min, lon_max, lat_min, lat_max):
        m = _in_box(lo, la, float(lon_min), float(lon_max), float(lat_min), float(lat_max))
        return f'the box {lat_min}–{lat_max}° lat, {lon_min}–{lon_max}° lon', np.nonzero(m)[0] + 1
    if not name:
        raise ToolError('give a region name or a lon/lat box')
    want = set(_norm(name)) - {'sea', 'ocean', 'region', 'area', 'basin', 'of'}
    full = ' '.join(_norm(name))
    for nm, aliases, a, b, c, d in OCEAN_REGIONS:               # common regions first (exact name or alias)
        if full in [' '.join(_norm(nm))] + [' '.join(_norm(x)) for x in aliases]:
            return f'{nm} (a box: {c}–{d}° lat, {a}–{b}° lon)', np.nonzero(_in_box(lo, la, a, b, c, d))[0] + 1
    # a named area: the bins the cost map explorer names after it (each bin carries the name of the smallest Natural
    # Earth marine area containing its center: profile_analyses.viewer._bin_names), so the two always agree
    names = _bin_area_names()
    uniq = sorted(set(n for n in names if n))
    best, score = None, 0
    for nm in uniq:
        words = set(_norm(nm)) - {'sea', 'ocean', 'of', 'gulf', 'bay'}
        sc = 2 if ' '.join(_norm(nm)) == full else len(want & words) / max(1, len(want | words))
        if sc > score:
            best, score = nm, sc
    if best is None or score < 0.5:
        for nm, aliases, a, b, c, d in OCEAN_REGIONS:           # a looser match on the common regions
            if want & set(_norm(nm) + [w for x in aliases for w in _norm(x)]):
                return f'{nm} (a box: {c}–{d}° lat, {a}–{b}° lon)', np.nonzero(_in_box(lo, la, a, b, c, d))[0] + 1
        raise ToolError(f'no region called {name!r}; try a sea or ocean name (e.g. "Labrador Sea", "eastern tropical '
                        f'Pacific", "Southern Ocean") or give a lon/lat box')
    return f'{best} (Natural Earth marine area)', np.asarray([i + 1 for i, n in enumerate(names) if n == best])


def region_cost(A, src, cost_type, region=None, variable='T', part='unsplit', depth_from=0, depth_to=2000,
                lon_min=None, lon_max=None, lat_min=None, lat_max=None):
    """The cost in a region (a named sea or ocean area, or a lon/lat box) against the rest of the world, over a depth
    range, all years: the region's cost per datum (or mean term), its share of the bins, the data and the summed
    cost, and its highest bins."""
    V, part = _var(variable), _part(part)
    label, ids = find_region(A, region, lon_min, lon_max, lat_min, lat_max)
    k = _levels(A.depth, depth_from, depth_to)
    val, w = _bin_cost(A, src, cost_type, V, part, k)
    ok = np.isfinite(val) & (w > 0)
    inr = np.zeros(pa.N_BINS, bool); inr[ids - 1] = True
    def summ(m):
        m = m & ok
        if not m.any():
            return None
        if part == 'mean':
            # pooled over every bin and level with a value, as the cost map explorer's region panel
            return {'bins_with_data': int(m.sum()), 'mean_term_per_bin_and_level': _r(np.sum(val[m] * w[m]) / np.sum(w[m]))}
        return {'bins_with_data': int(m.sum()), 'data': int(w[m].sum()), 'cost_per_datum': _r(np.sum(val[m] * w[m]) / np.sum(w[m]))}
    reg, rest = summ(inr), summ(~inr)
    if reg is None:
        raise ToolError(f'{label}: no {src} {V} data at {depth_from}–{depth_to} m')
    contrib = np.where(ok, val * w, 0)            # the summed cost (mean term: summed over bins and levels)
    res = {'region': label, 'region_bins': int(len(ids)), 'obs_class': src, 'cost_type': cost_type, 'variable': V,
           'part': part, 'depth_range_m': [_r(A.depth[k[0]], 4), _r(A.depth[k[-1]], 4)], 'years': 'all',
           'region_result': reg, 'rest_of_world': rest,
           'region_share_of_bins_with_data_pct': _r(100 * reg['bins_with_data'] / ok.sum()),
           'region_share_of_summed_cost_pct': _r(100 * contrib[inr].sum() / contrib.sum()) if contrib.sum() else None,
           'highest_bins_in_region': [{'bin': int(b + 1), 'where': _where(A, b + 1), 'value': _r(val[b])}
                                      for b in sorted(np.nonzero(inr & ok)[0], key=lambda b: -val[b])[:5]],
           'note': 'a share of the summed cost much bigger than the share of bins means the region weighs more than '
                   'its size; regions are sets of ~220 km geodesic bins, so small seas are coarse'}
    depth_param = _group_param(A, k)
    code = '.'.join(str(int(b)) for b in sorted(ids))
    return {'result': res, 'views': [{'page': 'explorer', 'params': {
        'src': src, 'ct': cost_type, 'V': V, 'part': part, 'depth': depth_param, 'years': 'all', 'region': code},
        'caption': f'{label.split(" (")[0]}: {src}, {cost_type.replace("_", " ")}, {V}, {"mean term" if part == "mean" else part} '
                   f'(cost map explorer with the region outlined)'}]}


def _text(html):
    """Plain text from the report's HTML snippets (tags dropped, entities decoded, spaces collapsed)."""
    import html as _h
    import re
    t = re.sub(r'<br\s*/?>|</p>|</tr>|</li>', '\n', html)
    t = re.sub(r'<[^>]+>', ' ', t)
    t = _h.unescape(t)
    return re.sub(r'[ \t]+', ' ', re.sub(r'\n\s*', '\n', t)).strip()


def explain(A, topic):
    """The report's own text on an analysis ('01'-'12') or a term: what it shows and how it is computed. Use it for
    any question about what something means, before or instead of computing numbers."""
    from ecco_insitu.diagnostics.profile_analyses import guide, report
    t = str(topic).strip().lower().replace('analysis', '').strip()
    if t.isdigit() and f'{int(t):02d}' in report._ANALYSIS_TITLES:
        an = f'{int(t):02d}'
        return {'result': {'analysis': an, 'title': report._ANALYSIS_TITLES[an],
                           'description': _text(report._ANALYSIS_DESC.get(an, ''))},
                'views': [{'page': 'report', 'params': {'src': 'PFL', 'an': an},
                           'caption': f'Analysis {an} in the report'}]}
    if t in ('analyses', 'list', 'all'):
        return {'result': {'analyses': report._ANALYSIS_TITLES}, 'views': []}
    if t in ('terms', 'terminology', 'glossary', 'weights', 'sigma', 'σ', 'cost', 'split', 'weight versions'):
        terms = next(body for key, _, body in guide.SECTIONS if key == 'terms')
        return {'result': {'report_terminology': _text(report._GLOSSARY_HTML)[:6000], 'guide_terms': _text(terms)[:4000]},
                'views': [{'page': 'guide', 'params': {}, 'caption': 'The guide (terms)'}]}
    for key, title, body in guide.SECTIONS:
        if t in (key, title.lower()) or t in title.lower():
            return {'result': {'guide_section': title, 'text': _text(body)[:6000]},
                    'views': [{'page': 'guide', 'params': {}, 'caption': f'The guide: {title}'}]}
    raise ToolError(f'no text on {topic!r}; topics: an analysis number 01-12, "analyses", "terms", or a guide section: '
                    + ', '.join(k for k, *_ in guide.SECTIONS))


def open_view(A, page, **params):
    """Open any view of the site directly (when no computation is needed)."""
    if page not in ('report', 'explorer', 'data', 'questions', 'guide'):
        raise ToolError('page must be report, explorer, data, questions or guide')
    return {'result': {'opened': page}, 'views': [{'page': page, 'params': {k: v for k, v in params.items() if v is not None},
                                                   'caption': params.pop('caption', None) or f'The {page} view'}]}


TOOLS = {f.__name__: f for f in (explain, region_cost, list_contents, cost_vs_depth, cost_vs_time, compare_classes, compare_references,
                                 where_high, offset_vs_variability, peak_attribution, observation_counts, open_view)}

# JSON schemas for the model (OpenAI tool format)
_SRC = {'type': 'string', 'description': 'obs class, e.g. PFL (Argo floats), CTD_WOD, XBT_WOD, GLD_WOD (gliders), MRB_WOD (moorings), MEOP (seals), ITP'}
_CT = {'type': 'string', 'description': 'reference: vs_clim (the WOA23 climatology) or vs_<model> for a model, e.g. vs_V4r6, '
       'vs_V5r1 (the analysis lists them: list_contents)'}
_V = {'type': 'string', 'enum': ['T', 'S']}
_PART = {'type': 'string', 'enum': list(PARTS), 'description': 'unsplit cost; or the split\'s anomaly term (variability in time) or mean term (steady offset per bin)'}
_Z = {'type': 'number', 'description': 'meters'}


def _schema(name, desc, props, required):
    return {'type': 'function', 'function': {'name': name, 'description': desc,
                                             'parameters': {'type': 'object', 'properties': props, 'required': required}}}


SCHEMAS = [
    _schema('region_cost', 'The cost in a region against the rest of the world, over a depth range (all years): a named sea '
            'or ocean area ("Labrador Sea", "Gulf of Mexico", "eastern tropical Pacific", "Southern Ocean", "subpolar '
            'North Atlantic", "Kuroshio region", ...) or a lon/lat box (lon -180..180; lon_min > lon_max crosses the dateline).',
            {'src': _SRC, 'cost_type': _CT, 'region': {'type': 'string'}, 'variable': _V, 'part': _PART,
             'depth_from': _Z, 'depth_to': _Z, 'lon_min': {'type': 'number'}, 'lon_max': {'type': 'number'},
             'lat_min': {'type': 'number'}, 'lat_max': {'type': 'number'}}, ['src', 'cost_type']),
    _schema('explain', 'The report\'s own text: what an analysis (01-12) shows and how it is computed, the list of analyses '
            '("analyses"), the definitions of terms ("terms": cost, weights, sigma, the split), or a guide section '
            '(e.g. "attribution", "peak", "explorer", "viewer"). Use it for questions about meaning or method.',
            {'topic': {'type': 'string'}}, ['topic']),
    _schema('list_contents', 'The obs classes, cost types, variables, depth levels and years the analysis has.', {}, []),
    _schema('cost_vs_depth', 'Cost per datum against depth (analysis 02) for one class and reference; the largest level.',
            {'src': _SRC, 'cost_type': _CT, 'variable': _V, 'part': _PART, 'depth_from': _Z, 'depth_to': _Z}, ['src', 'cost_type']),
    _schema('cost_vs_time', 'Cost per datum per year (analysis 03): highest and lowest years, trend per decade.',
            {'src': _SRC, 'cost_type': _CT, 'variable': _V, 'part': {**_PART, 'enum': ['unsplit', 'anomaly']}}, ['src', 'cost_type']),
    _schema('compare_classes', 'Cost per datum of every obs class over one depth range.',
            {'cost_type': _CT, 'variable': _V, 'part': _PART, 'depth_from': _Z, 'depth_to': _Z}, ['cost_type']),
    _schema('compare_references', 'Cost per datum of one class against the climatology and against the model.',
            {'src': _SRC, 'variable': _V, 'part': _PART, 'depth_from': _Z, 'depth_to': _Z}, ['src']),
    _schema('where_high', 'Where the cost is high over a depth range: the highest bins and the areas holding the top 10 % of bins.',
            {'src': _SRC, 'cost_type': _CT, 'variable': _V, 'part': _PART, 'depth_from': _Z, 'depth_to': _Z}, ['src', 'cost_type']),
    _schema('offset_vs_variability', 'Is the misfit a steady offset or variability in time: the split terms and the mean misfit (reference - obs).',
            {'src': _SRC, 'cost_type': _CT, 'variable': _V, 'depth_from': _Z, 'depth_to': _Z}, ['src', 'cost_type']),
    _schema('peak_attribution', 'Where a bump in cost against depth comes from: peak level, share due to which bins have data, areas.',
            {'src': _SRC, 'cost_type': _CT, 'variable': _V, 'part': _PART, 'depth_from': _Z, 'depth_to': _Z}, ['src', 'cost_type', 'depth_from', 'depth_to']),
    _schema('observation_counts', 'How many observations a class has, per year, over a depth range.',
            {'src': _SRC, 'variable': _V, 'depth_from': _Z, 'depth_to': _Z}, ['src']),
    _schema('open_view', 'Open a page of the site with settings, when no numbers are needed. report: src, ct, an (01-12); '
            'explorer: src, ct, V, part (unsplit|anomaly|mean|count|bias|nmis|sig), depth (all|k<level>), years (all|Y0-Y1); '
            'data: src, ref (clim|model), V, k (level index); questions; guide.',
            {'page': {'type': 'string', 'enum': ['report', 'explorer', 'data', 'questions', 'guide']},
             'src': {'type': 'string'}, 'ct': {'type': 'string'}, 'an': {'type': 'string'}, 'V': {'type': 'string'},
             'part': {'type': 'string'}, 'depth': {'type': 'string'}, 'years': {'type': 'string'},
             'ref': {'type': 'string'}, 'k': {'type': 'integer'}, 'caption': {'type': 'string'}}, ['page']),
]


def call(A, name, args):
    """Run tool ``name`` with ``args`` (dict); errors come back as a result the model can read."""
    if name not in TOOLS:
        return {'error': f'unknown tool {name}'}, []
    try:
        out = TOOLS[name](A, **(args or {}))
        return out['result'], out['views']
    except (ToolError, TypeError) as e:
        return {'error': str(e)}, []
