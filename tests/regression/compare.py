"""
Compare two chain output folders (``yearly/``, ``ecco_arrays/``, ``documentation/``) for the
regression test.

NetCDF files must have the same set of files, dimensions, variables, dtypes, attributes and
values (exactly; NaN equals NaN). Attributes and JSON/text content that legitimately change
from run to run (creation time, absolute paths) are listed in VOLATILE and ignored.
"""
import json
import os
import re

import netCDF4
import numpy as np

VOLATILE_ATTRS = {'created'}                       # global attributes that hold the run time
VOLATILE_JSON_KEYS = {'profiles_dir', 'arrays_dir', 'grid_dir'}   # absolute paths in the manifest


def _files(root):
    return sorted(os.path.relpath(os.path.join(d, f), root)
                  for d, _, fs in os.walk(root) for f in fs if not f.startswith('.'))


def _compare_nc(a, b, rel, problems):
    with netCDF4.Dataset(a) as x, netCDF4.Dataset(b) as y:
        for ds in (x, y):
            ds.set_auto_maskandscale(False)
        gx = {k: x.getncattr(k) for k in x.ncattrs() if k not in VOLATILE_ATTRS}
        gy = {k: y.getncattr(k) for k in y.ncattrs() if k not in VOLATILE_ATTRS}
        if gx.keys() != gy.keys() or any(str(gx[k]) != str(gy[k]) for k in gx):
            problems.append(f'{rel}: global attributes differ')
        dx = {k: len(v) for k, v in x.dimensions.items()}
        dy = {k: len(v) for k, v in y.dimensions.items()}
        if dx != dy:
            problems.append(f'{rel}: dimensions differ: {dx} vs {dy}')
            return
        if set(x.variables) != set(y.variables):
            problems.append(f'{rel}: variables differ: only in new {sorted(set(x.variables) - set(y.variables))}, '
                            f'only in reference {sorted(set(y.variables) - set(x.variables))}')
        for v in sorted(set(x.variables) & set(y.variables)):
            vx, vy = x[v], y[v]
            if vx.dtype != vy.dtype or vx.dimensions != vy.dimensions:
                problems.append(f'{rel}:{v}: dtype/dims differ ({vx.dtype}{vx.dimensions} vs {vy.dtype}{vy.dimensions})')
                continue
            ax = {k: str(vx.getncattr(k)) for k in vx.ncattrs()}
            ay = {k: str(vy.getncattr(k)) for k in vy.ncattrs()}
            if ax != ay:
                problems.append(f'{rel}:{v}: attributes differ')
            p, q = np.asarray(vx[:]), np.asarray(vy[:])
            same = np.array_equal(p, q, equal_nan=True) if p.dtype.kind == 'f' else np.array_equal(p, q)
            if not same:
                if p.dtype.kind == 'f':
                    d = np.abs(p - q); d = d[np.isfinite(d)]
                    n = int(np.sum(~((p == q) | (np.isnan(p) & np.isnan(q)))))
                    problems.append(f'{rel}:{v}: {n} values differ (max abs diff {d.max() if d.size else "nan"})')
                else:
                    problems.append(f'{rel}:{v}: values differ')


def _strip_json(o):
    if isinstance(o, dict):
        return {k: _strip_json(v) for k, v in o.items() if k not in VOLATILE_JSON_KEYS}
    if isinstance(o, list):
        return [_strip_json(v) for v in o]
    return o


def compare_runs(new, ref):
    """Return a list of differences between output folders ``new`` and ``ref`` (empty if identical)."""
    problems = []
    fn, fr = _files(new), _files(ref)
    if fn != fr:
        problems.append(f'file lists differ: only in new {sorted(set(fn) - set(fr))}, '
                        f'only in reference {sorted(set(fr) - set(fn))}')
    for rel in sorted(set(fn) & set(fr)):
        a, b = os.path.join(new, rel), os.path.join(ref, rel)
        if rel.endswith('.nc'):
            _compare_nc(a, b, rel, problems)
        elif rel.endswith('.json'):
            if _strip_json(json.load(open(a))) != _strip_json(json.load(open(b))):
                problems.append(f'{rel}: content differs')
        else:
            ta, tb = (re.sub(r'/\S+', '<path>', open(p).read()) for p in (a, b))
            if ta != tb:
                problems.append(f'{rel}: content differs')
    return problems
