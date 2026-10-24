"""
Set up the climatology the NCEI chain compares profiles with (step 03): check the one in the
settings, or make one.

The chain needs a monthly potential-temperature and salinity climatology on a regular lon/lat
grid (the contract is in docs/input_contract.md). There are two ways to provide it:

1. **Use your own file.** Any NetCDF file that meets the contract works: set
   ``paths.climatology_file`` in config/local.yaml (or a --config file) and run
   ``ecco-insitu climatology --check``.
2. **Build the standard one** (what this project uses): WOA23 1991-2020, 1-degree, monthly
   fields to 1500 m spliced onto seasonal fields to 5500 m, in-situ T converted to potential T
   with the chain's own UNESCO routine. ``ecco-insitu climatology --build`` downloads about
   1.4 GB from NOAA and writes it to the folder of ``paths.climatology_file``.

Either way, pre-baked copies (the climatology already interpolated onto each input depth grid,
which makes step 03 much faster) can then be made with ``--prebake``. The chain only uses
pre-baked files made from the climatology currently set, so they can't get out of step.

With no option, this checks the current setup and says what to do.

Usage:
    ecco-insitu climatology [--check | --build | --prebake] [--config FILE] [--set key=value]
"""
import argparse
import sys
from pathlib import Path

import ecco_insitu.config as config_loader
import ecco_insitu.validate as validate_inputs


def _status(cfg):
    clim = cfg['paths']['climatology_file']
    pre_dir = cfg['paths']['prebaked_clim_dir']
    print(f'paths.climatology_file  = {clim}')
    print(f'paths.prebaked_clim_dir = {pre_dir}')
    if not clim:
        print('\nNo climatology is set. Either set paths.climatology_file in config/local.yaml to a file of '
              'your own, or to where the standard one should go, then run:  ecco-insitu climatology --build')
        return False
    if not Path(clim).exists():
        print(f'\nThe climatology file does not exist. To build the standard WOA23 one there, run:\n'
              f'    ecco-insitu climatology --build\n'
              f'Or point paths.climatology_file at a climatology of your own (see docs/input_contract.md).')
        return False
    rep = validate_inputs.check_climatology_file(clim)
    print(f'\nclimatology: {"OK" if rep.ok else "NOT USABLE"}')
    for line in rep.lines():
        print(f'  {line}')
    if not rep.ok:
        return False
    if pre_dir and Path(pre_dir).is_dir():
        import netCDF4
        n = 0
        for pb in sorted(Path(pre_dir).glob('*_prebaked.nc')):
            with netCDF4.Dataset(pb) as d:
                src = d.getncattr('source_file') if 'source_file' in d.ncattrs() else None
            match = src is not None and Path(src).resolve() == Path(clim).resolve()
            n += match
            print(f'  pre-baked {"used   " if match else "IGNORED"}  {pb.name}'
                  + ('' if match else f'  (made from {src})'))
        if n == 0:
            print('  no pre-baked files for this climatology (optional; faster). Make them with:  '
                  'ecco-insitu climatology --prebake')
    return True


def cli(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    g = ap.add_mutually_exclusive_group()
    g.add_argument('--check', action='store_true', help='check the climatology in the settings (default)')
    g.add_argument('--build', action='store_true', help='download WOA23 and build the standard climatology')
    g.add_argument('--prebake', action='store_true', help='make pre-baked copies for the input depth grids')
    config_loader.add_arguments(ap)
    a = ap.parse_args(argv)
    try:
        cfg = config_loader.from_args(a)
    except config_loader.ConfigError as e:
        sys.exit(f'config error: {e}')
    clim = cfg['paths']['climatology_file']
    if a.build:
        if not clim:
            sys.exit('set paths.climatology_file (where the file should go) in config/local.yaml first')
        if Path(clim).exists():
            sys.exit(f'{clim} already exists; remove it or set another paths.climatology_file to build a new one')
        import ecco_insitu.climatology.build_woa23 as build_woa23_climatology
        out = build_woa23_climatology.build(Path(clim).parent)
        if Path(out).resolve() != Path(clim).resolve():
            print(f'\nNOTE: the file was written as {out}; set paths.climatology_file to that path.')
        return 0 if _status(cfg) else 1
    if a.prebake:
        if not (clim and Path(clim).exists()):
            sys.exit('no climatology yet; run --build or set paths.climatology_file first')
        if not cfg['run']['input_dir']:
            sys.exit('set run.input_dir (the profile files whose depth grids to pre-bake for)')
        import ecco_insitu.climatology.prebake as prebake_woa23_climatology
        prebake_woa23_climatology.prebake(clim, cfg['run']['input_dir'],
                                          cfg['paths']['prebaked_clim_dir'] or Path(clim).parent)
        return 0 if _status(cfg) else 1
    return 0 if _status(cfg) else 1


if __name__ == '__main__':
    sys.exit(cli())
