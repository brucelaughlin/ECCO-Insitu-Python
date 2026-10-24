"""
Cut the small regression-test input files in tests/fixtures/inputs/ from full NCEI input files.

Each fixture keeps the source file's format, variables and attributes, and a subset of its
profiles chosen to exercise the chain:
  - a contiguous block (keeps near-duplicate profiles together, for step 09);
  - profiles north of 60°N and south of 60°S (high-latitude sigma, step 06 tests);
  - profiles with out-of-range HHMMSS (date/time normalization).
Selection is deterministic (fixed seed), so rerunning gives the same files.

Run once; the fixtures are committed. Usage:
    python tests/make_fixtures.py --src /path/to/Interp_Profiles
"""
import argparse
import os

import netCDF4
import numpy as np

# (path under --src, output subfolder, contiguous block start, block length, n north, n south)
PLAN = [
    ('CTD_WOD/WOD_WO_1992_CTD_OSD.nc', 'CTD_WOD', 20000, 150, 40, 40),
    ('PFL/ARGO_WO_2015_PFL_D.nc',      'PFL',     60000, 150, 30, 30),
    ('XBT_WOD/WOD_WO_1992_XBT.nc',     'XBT_WOD',  5000, 150, 20, 20),
    ('GLD_WOD/WOD_WO_2013_GLD.nc',     'GLD_WOD', 90000, 150, 20, 20),
    ('MRB_WOD/WOD_WO_2025_MRB.nc',     'MRB_WOD',  3000, 150,  0,  0),
    ('MEOP/MEOP_WO_2024_CTD.nc',       'MEOP',      500, 100,  0, 60),
    ('ITP/L2/ITP2_WO_2023_CTD.nc',     'ITP/L2',    100, 120,  0,  0),
    ('ITP/L3/ITP_WO_2004_CTD.nc',      'ITP/L3',     50, 100,  0,  0),
]
N_BADTIME = 10


def pick(ds, start, length, n_north, n_south, rng):
    lat = np.asarray(ds['prof_lat'][:], float)
    t = np.asarray(ds['prof_HHMMSS'][:], float)
    n = len(lat)
    keep = set(range(min(start, max(0, n - length)), min(n, start + length)))
    for mask, k in ((lat > 60, n_north), (lat < -60, n_south),
                    ((t % 100 >= 60) | ((t // 100) % 100 >= 60) | (t >= 240000), N_BADTIME)):
        idx = np.flatnonzero(mask)
        if len(idx) and k:
            keep.update(rng.choice(idx, min(k, len(idx)), replace=False).tolist())
    return np.array(sorted(keep))


def copy_subset(src, dst, idx):
    with netCDF4.Dataset(src) as s, netCDF4.Dataset(dst, 'w', format=s.data_model) as d:
        s.set_auto_maskandscale(False)
        d.setncatts({a: s.getncattr(a) for a in s.ncattrs()})
        for name, dim in s.dimensions.items():
            d.createDimension(name, len(idx) if name == 'iPROF' else (None if dim.isunlimited() else len(dim)))
        for name, v in s.variables.items():
            fill = v.getncattr('_FillValue') if '_FillValue' in v.ncattrs() else None
            o = d.createVariable(name, v.dtype, v.dimensions, fill_value=fill)
            o.setncatts({a: v.getncattr(a) for a in v.ncattrs() if a != '_FillValue'})
            o.set_auto_maskandscale(False)
            o[:] = v[:][idx] if v.dimensions and v.dimensions[0] == 'iPROF' else v[:]


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--src', required=True, help='folder of full input files (Interp_Profiles)')
    p.add_argument('--out', default=os.path.join(os.path.dirname(__file__), 'fixtures', 'inputs'))
    a = p.parse_args()
    rng = np.random.default_rng(20261023)
    for rel, sub, start, length, nn, ns in PLAN:
        src = os.path.join(a.src, rel)
        with netCDF4.Dataset(src) as ds:
            idx = pick(ds, start, length, nn, ns, rng)
        os.makedirs(os.path.join(a.out, sub), exist_ok=True)
        dst = os.path.join(a.out, sub, os.path.basename(rel))
        copy_subset(src, dst, idx)
        print(f'{dst}: {len(idx)} profiles, {os.path.getsize(dst) / 1e6:.1f} MB')


if __name__ == '__main__':
    main()
