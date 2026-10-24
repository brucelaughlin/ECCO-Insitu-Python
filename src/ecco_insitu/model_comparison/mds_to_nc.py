"""
**Work in progress** (see :mod:`ecco_insitu.model_comparison`). Verified: bit-identical to
MITprof_gcm2nc.m on all 21 ECCO v4r6 iter47 datasets.

Convert MITgcm profiles-package output (<name>.<tile>.001.equi.data, one file per tile) into
NetCDF: a copy of the input profile file with prof_<V>estim filled with the model equivalents.
Python replacement for MITprof_gcm2nc.m (MITgcm/MITprof, MITprof_IO/), which the legacy
reconstruct_profile_nc_from_model_data_binary.m called.

Binary layout (from MITprof_gcm2nc.m): big-endian float64, per profile a block of
(nr+1) x 2 x nvar values (Fortran order), nr = number of depth levels::

  [0:nr, 0, v]  model value of variable v at each level
  [0:nr, 1, v]  mask (0 = no model value)
  [nr,   0, v]  1-based profile index into the input file (<= 0 for unused records)

Variables are those of T, S, U, V, ptr, ssh, OXY that the input file has, in that order.

Differences from the MATLAB:

- every variable of the input file is kept (MITprof_write kept only a subset, which drops
  e.g. prof_bin_id_a and the interpolation fields);
- checks that each profile is found at most once and reports how many were never found;
- datasets are processed in parallel.

Usage::

    ecco-insitu mds-to-nc --mds <dir with *.equi.data> --input <dir with input .nc>
                          --out <output dir> [--only NAME ...] [-n 8]

Output files are ``<input name>_model.nc`` in ``--out``.
"""
import argparse
import glob
import os
import shutil
import sys
from concurrent.futures import ProcessPoolExecutor

import netCDF4
import numpy as np

VAR_ORDER = ('T', 'S', 'U', 'V', 'ptr', 'ssh', 'OXY')
FILL = -9999.0


def read_tiles(tile_files, nr, nvar, nprof):
    """Model values (nvar, nprof, nr) with FILL where masked or missing, and a count of
    how many times each profile appeared."""
    est = np.full((nvar, nprof, nr), FILL)
    seen = np.zeros((nvar, nprof), np.int32)
    block = (nr + 1) * 2 * nvar
    for path in tile_files:
        a = np.fromfile(path, '>f8')
        if a.size % block:
            raise ValueError(f'{path}: {a.size} values is not a multiple of the record size {block}')
        a = a.reshape(-1, nvar, 2, nr + 1)          # profile, var, (value|mask), level
        for v in range(nvar):
            idx = a[:, v, 0, nr]
            keep = idx > 0
            ip = idx[keep].astype(np.int64) - 1
            if ip.size and ip.max() >= nprof:
                raise ValueError(f'{path}: profile index {ip.max() + 1} > {nprof} profiles in the input file')
            val, msk = a[keep, v, 0, :nr], a[keep, v, 1, :nr]
            est[v, ip] = np.where(msk == 0, FILL, val)
            np.add.at(seen[v], ip, 1)
    return est, seen


def convert(name, mds_dir, input_dir, out_dir):
    """Write ``<out_dir>/<name>_model.nc``: a copy of ``<input_dir>/<name>.nc`` with ``prof_<V>estim``
    filled from the tile files ``<mds_dir>/<name>.<tile>.001.equi.data``.

    Returns
    -------
    str
        A one-line summary: profiles, tiles, variables, and how many profiles had no model record.

    Raises
    ------
    ValueError
        If a tile file's size doesn't fit the record layout, a profile index exceeds the input file,
        or a profile appears in more than one tile.
    """
    src = os.path.join(input_dir, name + '.nc')
    tiles = sorted(glob.glob(os.path.join(mds_dir, name + '.[0-9][0-9][0-9].[0-9][0-9][0-9].equi.data')))
    if not tiles:
        return f'{name}: no model files found, skipped'
    with netCDF4.Dataset(src) as ds:
        nr = len(ds.dimensions['iDEPTH'])
        nprof = len(ds.dimensions['iPROF'])
        vars_ = [v for v in VAR_ORDER if f'prof_{v}' in ds.variables]
    est, seen = read_tiles(tiles, nr, len(vars_), nprof)
    if (seen > 1).any():
        raise ValueError(f'{name}: {int((seen > 1).sum())} profile records appear in more than one tile')

    dst = os.path.join(out_dir, name + '_model.nc')
    tmp = dst + '.part'
    shutil.copyfile(src, tmp)
    with netCDF4.Dataset(tmp, 'r+') as ds:
        ds.set_auto_mask(False)
        for v, vv in enumerate(vars_):
            key = f'prof_{vv}estim'
            if key not in ds.variables:
                ds.createVariable(key, 'f8', ds[f'prof_{vv}'].dimensions, fill_value=FILL)
            ds[key][:] = est[v]
        ds.setncattr('model_equivalents_from', os.path.abspath(mds_dir))
    os.replace(tmp, dst)
    missing = ', '.join(f'{vv} {int((seen[v] == 0).sum())}' for v, vv in enumerate(vars_))
    return f'{name}: {nprof} profiles, {len(tiles)} tiles, vars {"/".join(vars_)}; profiles with no model record: {missing}'


def cli(argv=None):
    """Command-line entry point (``ecco-insitu`` subcommand)."""
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--mds', required=True, help='folder with the <name>.<tile>.001.equi.data files')
    p.add_argument('--input', required=True, help='folder with the input profile files <name>.nc given to the model')
    p.add_argument('--out', required=True, help='output folder (created if needed)')
    p.add_argument('--only', nargs='+', help='convert only these dataset names')
    p.add_argument('-n', type=int, default=4, help='parallel workers (default 4)')
    a = p.parse_args(argv)

    names = a.only or sorted({os.path.basename(f).rsplit('.', 4)[0]
                              for f in glob.glob(os.path.join(a.mds, '*.equi.data'))})
    names = [n for n in names if os.path.exists(os.path.join(a.input, n + '.nc'))
             or print(f'{n}: no input file {n}.nc, skipped', file=sys.stderr)]
    os.makedirs(a.out, exist_ok=True)
    with ProcessPoolExecutor(a.n) as ex:
        for msg in ex.map(convert, names, [a.mds] * len(names), [a.input] * len(names), [a.out] * len(names)):
            print(msg, flush=True)


if __name__ == '__main__':
    sys.exit(cli())
