"""
Check input files against what the NCEI chain needs, before processing them.

Two contracts are defined here, and documented in ``docs/input_contract.md``:

**Profile files** (the chain's inputs, MITprof format). :func:`check_profile_file` reads only
the metadata and the small 1-D variables, so a whole input folder is checked in seconds. Each
problem is reported as

- ``ERROR``: the file can't be processed correctly (e.g. no ``prof_lat``, or a 1-D variable
  on the wrong dimension); the file is skipped, with the reason in the run summary.
- ``WARNING``: the file is processed, but something is unusual and worth a look (e.g. dates
  before 1950, or no label saying whether ``prof_T`` is potential or in-situ temperature).

**Climatology files** (``paths.climatology_file``, and the pre-baked ones). The chain needs
monthly potential temperature and salinity on a regular lon/lat grid; any file with the
expected variable names, dimensions and value ranges can be used, whether made by
``ecco_insitu/climatology/build_woa23.py`` or by the user. :func:`check_climatology_file` checks it.

Command line::

    ecco-insitu validate <input_dir>            # check every *.nc profile file under it
    ecco-insitu validate --climatology <file>   # check a climatology file
"""
import argparse
import datetime
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import netCDF4

# ---------------------------------------------------------------------------------------------
# Profile-file contract
# ---------------------------------------------------------------------------------------------

#: Variables every profile file must have: name -> required dimensions.
REQUIRED_PROFILE_VARS = {
    'prof_lon': ('iPROF',),
    'prof_lat': ('iPROF',),
    'prof_YYYYMMDD': ('iPROF',),
    'prof_HHMMSS': ('iPROF',),
    'prof_depth': ('iDEPTH',),
}
#: Observed variables: at least one must be present, on (iPROF, iDEPTH), each with its weight.
OBS_VARS = ('prof_T', 'prof_S')
#: Fill value the chain's input files use (also accepted: NaN).
FILL = -9999.0

# plausible ranges for the warnings (values outside are flagged, not changed)
LAT_RANGE = (-90.0, 90.0)
LON_RANGE = (-180.0, 360.0)
LON_MISSING = 99999.0           # MITprof "no position" marker; step 06 zeroes those weights
DEPTH_RANGE = (0.0, 11000.0)
T_RANGE = (-3.0, 40.0)          # deg C
S_RANGE = (0.0, 42.0)           # psu
FIRST_YEAR = 1950               # step 06 rejects earlier dates


@dataclass
class FileReport:
    """Result of checking one file.

    Attributes
    ----------
    path : str
        The file checked.
    errors : list of str
        Problems that stop the file being processed.
    warnings : list of str
        Unusual things that don't stop it.
    n_profiles : int
        Number of profiles (``iPROF`` length), when it could be read.
    """
    path: str
    errors: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    n_profiles: int = 0

    @property
    def ok(self):
        """True if the file can be processed (warnings allowed)."""
        return not self.errors

    def lines(self):
        """Human-readable lines, errors first: ``ERROR   ...`` / ``WARNING ...``.

        Returns
        -------
        list of str
        """
        return [f'ERROR   {e}' for e in self.errors] + [f'WARNING {w}' for w in self.warnings]


def _valid(x):
    x = np.asarray(x, float)
    return x[np.isfinite(x) & (x != FILL)]


def _t_label(var):
    text = ' '.join(str(var.getncattr(a)) for a in ('long_name', 'standard_name') if a in var.ncattrs()).lower()
    if 'potential' in text:
        return 'potential'
    if any(k in text for k in ('in situ', 'in-situ', 'insitu', 'sea_water_temperature')):
        return 'in-situ'
    return None


def check_profile_file(path, variables=OBS_VARS):
    """Check one input profile file against the chain's contract.

    Parameters
    ----------
    path : str or Path
        A MITprof-format NetCDF file.
    variables : sequence of str
        The observed variables the chain processes (``chain.variables`` in the settings).

    Returns
    -------
    FileReport
        ``report.ok`` is False if the file must be skipped; ``report.errors`` and
        ``report.warnings`` say why.
    """
    rep = FileReport(str(path))
    try:
        ds = netCDF4.Dataset(path)
    except Exception as e:                      # not NetCDF, truncated, permissions, ...
        rep.errors.append(f'cannot be opened as NetCDF: {e}')
        return rep
    with ds:
        ds.set_auto_mask(False)
        dims = {k: len(v) for k, v in ds.dimensions.items()}
        for d in ('iPROF', 'iDEPTH'):
            if d not in dims:
                rep.errors.append(f'no "{d}" dimension (dimensions are: {", ".join(dims) or "none"})')
        if rep.errors:
            return rep
        rep.n_profiles = dims['iPROF']
        if dims['iPROF'] == 0:
            rep.errors.append('no profiles (iPROF has length 0)')
            return rep

        for name, want in REQUIRED_PROFILE_VARS.items():
            if name not in ds.variables:
                rep.errors.append(f'required variable "{name}" is missing')
            elif ds[name].dimensions != want:
                rep.errors.append(f'"{name}" has dimensions {ds[name].dimensions}, expected {want}')
        present = [v for v in variables if v in ds.variables]
        if not present:
            rep.errors.append(f'none of the observed variables {", ".join(variables)} is present')
        for v in present:
            if ds[v].dimensions != ('iPROF', 'iDEPTH'):
                rep.errors.append(f'"{v}" has dimensions {ds[v].dimensions}, expected (iPROF, iDEPTH)')
            if f'{v}weight' not in ds.variables:
                rep.errors.append(f'"{v}" has no "{v}weight" (the chain sets weights from sigma, '
                                  f'but needs the variable to exist; add it filled with 1)')
            elif ds[f'{v}weight'].dimensions != ds[v].dimensions:
                rep.errors.append(f'"{v}weight" dimensions {ds[f"{v}weight"].dimensions} differ from "{v}"')
        if rep.errors:
            return rep

        # --- values (warnings: the chain handles these, but they're worth knowing about) ---
        depth = np.asarray(ds['prof_depth'][:], float)
        good_d = depth[np.isfinite(depth) & (depth != FILL)]
        if good_d.size < depth.size:
            rep.errors.append(f'prof_depth has {depth.size - good_d.size} missing levels')
        elif np.any(np.diff(depth) <= 0):
            rep.errors.append('prof_depth is not strictly increasing')
        elif depth.min() < DEPTH_RANGE[0] or depth.max() > DEPTH_RANGE[1]:
            rep.errors.append(f'prof_depth outside {DEPTH_RANGE} m: {depth.min()}..{depth.max()} '
                              '(is it pressure, or negative-down?)')

        lat = np.asarray(ds['prof_lat'][:], float)
        lon = np.asarray(ds['prof_lon'][:], float)
        bad_lat = int(np.sum(~np.isfinite(lat) | (lat < LAT_RANGE[0]) | (lat > LAT_RANGE[1])))
        lon_missing = int(np.sum(lon == LON_MISSING))
        bad_lon = int(np.sum(~np.isfinite(lon) | (((lon < LON_RANGE[0]) | (lon > LON_RANGE[1])) & (lon != LON_MISSING))))
        if bad_lat:
            rep.warnings.append(f'{bad_lat} profiles with latitude outside -90..90 (their weights will be zeroed)')
        if bad_lon:
            rep.warnings.append(f'{bad_lon} profiles with longitude outside -180..360 (their weights will be zeroed)')
        if lon_missing:
            rep.warnings.append(f'{lon_missing} profiles with no position (lon = 99999; their weights will be zeroed)')

        ymd = np.asarray(ds['prof_YYYYMMDD'][:], float)
        hms = np.asarray(ds['prof_HHMMSS'][:], float)
        y, m, d = ymd // 10000, (ymd // 100) % 100, ymd % 100
        no_day = int(np.sum((m >= 1) & (m <= 12) & (d == 0)))
        bad_date = int(np.sum((m < 1) | (m > 12) | (d > 31) | ~np.isfinite(ymd) | ((d == 0) & ((m < 1) | (m > 12)))))
        if no_day:
            rep.warnings.append(f'{no_day} profiles with the day missing (YYYYMM00: year and month known, day 00); '
                                'step 06 zeroes their weights')
        if bad_date:
            rep.warnings.append(f'{bad_date} profiles with an impossible date (their weights will be zeroed)')
        old = int(np.sum(y < FIRST_YEAR))
        future = int(np.sum(y > datetime.date.today().year))
        if old or future:
            rep.warnings.append(f'{old + future} profiles dated before {FIRST_YEAR} or in the future '
                                '(their weights will be zeroed)')
        roll = int(np.sum((hms % 100 >= 60) | ((hms // 100) % 100 >= 60) | (hms >= 240000)))
        if roll:
            rep.warnings.append(f'{roll} profiles with times like SS=60 or 24:00:00 (rolled forward to valid times)')

        # --- observed values and labels ---
        for v in present:
            vals = _valid(ds[v][:])
            if vals.size == 0:
                rep.warnings.append(f'"{v}" has no valid values')
                continue
            lo, hi = T_RANGE if v == 'prof_T' else S_RANGE if v == 'prof_S' else (-np.inf, np.inf)
            n_out = int(np.sum((vals < lo) | (vals > hi)))
            if n_out:
                rep.warnings.append(f'{n_out} "{v}" values outside {lo}..{hi} (min {vals.min():.3g}, max {vals.max():.3g}); '
                                    'check units and fill values')
            fv = [ds[v].getncattr(a) for a in ('missing_value', '_FillValue') if a in ds[v].ncattrs()]
            if fv and all(float(f) != FILL for f in fv):
                rep.warnings.append(f'"{v}" uses fill value {fv[0]}, not {FILL:g}; values equal to it are still read as '
                                    'missing')
        if 'prof_T' in present:
            label = _t_label(ds['prof_T'])
            if label is None:
                rep.warnings.append('prof_T has no label saying potential or in-situ temperature: step 05 will decide '
                                    'with a climatology check (see the log)')
    return rep


def check_input_dir(input_dir, variables=OBS_VARS):
    """Check every ``*.nc`` file under ``input_dir`` (recursively) with :func:`check_profile_file`.

    Parameters
    ----------
    input_dir : str or Path
    variables : sequence of str
        The observed variables the chain processes.

    Returns
    -------
    list of FileReport
        In sorted path order.
    """
    return [check_profile_file(f, variables) for f in sorted(Path(input_dir).rglob('*.nc'))]


# ---------------------------------------------------------------------------------------------
# Climatology-file contract
# ---------------------------------------------------------------------------------------------

#: Variables a climatology file must have. ``depth`` (full-depth file) or ``obs_depth``
#: (pre-baked file) is the vertical coordinate.
CLIM_VARS = ('potential_T_monthly', 'S_monthly')
CLIM_DIMS = ('month', None, 'lat', 'lon')        # None: depth or obs_depth


def check_climatology_file(path, prebaked=False):
    """Check a climatology file against what step 03 needs.

    Required:

    - ``potential_T_monthly`` and ``S_monthly``, dimensions (month, depth, lat, lon) with 12
      months (pre-baked files: (month, obs_depth, lat, lon));
    - 1-D coordinates ``lon``, ``lat``, ``depth`` (or ``obs_depth``), each evenly spaced for lon
      and lat (step 03 interpolates bilinearly on a regular grid), increasing, lat within
      -90..90, depth in meters (positive down);
    - temperature that is **potential** temperature (checked from the attributes when given,
      and by its value range), salinity in psu.

    Parameters
    ----------
    path : str or Path
        The climatology NetCDF file.
    prebaked : bool
        True for a pre-baked file (vertical coordinate ``obs_depth``).

    Returns
    -------
    FileReport
    """
    rep = FileReport(str(path))
    zname = 'obs_depth' if prebaked else 'depth'
    try:
        ds = netCDF4.Dataset(path)
    except Exception as e:
        rep.errors.append(f'cannot be opened as NetCDF: {e}')
        return rep
    with ds:
        ds.set_auto_mask(False)
        for c in ('lon', 'lat', zname, 'month'):
            if c not in ds.variables:
                rep.errors.append(f'coordinate "{c}" is missing')
        for v in CLIM_VARS:
            if v not in ds.variables:
                rep.errors.append(f'variable "{v}" is missing (names must be potential_T_monthly and S_monthly)')
            elif ds[v].dimensions != ('month', zname, 'lat', 'lon'):
                rep.errors.append(f'"{v}" has dimensions {ds[v].dimensions}, expected (month, {zname}, lat, lon)')
        if rep.errors:
            return rep
        if len(ds.dimensions['month']) != 12:
            rep.errors.append(f'{len(ds.dimensions["month"])} months; expected 12 (a monthly climatology)')
        lon, lat, z = (np.asarray(ds[c][:], float) for c in ('lon', 'lat', zname))
        for name, c in (('lon', lon), ('lat', lat)):
            dc = np.diff(c)
            if c.size < 2 or np.any(dc <= 0):
                rep.errors.append(f'"{name}" must be increasing')
            elif np.ptp(dc) > 1e-3 * abs(dc.mean()):
                rep.errors.append(f'"{name}" is not evenly spaced (step 03 needs a regular grid)')
        if lat.min() < -90 or lat.max() > 90:
            rep.errors.append(f'lat outside -90..90: {lat.min()}..{lat.max()}')
        if np.any(np.diff(z) <= 0):
            rep.errors.append(f'"{zname}" must be increasing (meters, positive down)')
        if z.min() < 0 or z.max() > 11000:
            rep.errors.append(f'"{zname}" outside 0..11000 m: {z.min()}..{z.max()} (meters, positive down?)')
        if not prebaked and z.max() < 5000:
            rep.warnings.append(f'deepest level is {z.max():.0f} m; deeper profile data will use the deepest '
                                'climatology value (nearest-depth fallback)')
        if rep.errors:
            return rep
        # value ranges, from a few months to keep it quick. Read as step 03 does (xarray): values
        # equal to a declared _FillValue / missing_value are missing; any other value is used as
        # data, so an undeclared fill value (-9999, 1e20) shows up here as out of range.
        ds.set_auto_mask(True)
        for v, (lo, hi), what in (('potential_T_monthly', T_RANGE, 'temperature'), ('S_monthly', S_RANGE, 'salinity')):
            x = ds[v][::4]
            x = np.asarray(x.compressed() if np.ma.isMaskedArray(x) else x, float).ravel()
            x = x[np.isfinite(x)]
            if x.size == 0:
                rep.errors.append(f'"{v}" has no valid values')
                continue
            if x.min() < lo - 1 or x.max() > hi + 1:
                rep.errors.append(f'"{v}" values {x.min():.3g}..{x.max():.3g} are outside the plausible {what} range '
                                  f'{lo}..{hi}; check units, and that missing values are NaN or declared as '
                                  f'_FillValue / missing_value')
        text = ' '.join(str(ds.getncattr(a)) for a in ds.ncattrs()).lower() + ' ' + \
            ' '.join(str(ds['potential_T_monthly'].getncattr(a)) for a in ds['potential_T_monthly'].ncattrs()).lower()
        if 'in situ' in text or 'in-situ' in text:
            if 'potential' not in text:
                rep.errors.append('the attributes say the temperature is in-situ; the chain needs potential temperature')
        elif 'potential' not in text:
            rep.warnings.append('nothing in the attributes says the temperature is potential temperature; '
                                'make sure it is (WOA t_an is in-situ and must be converted)')
    return rep


# ---------------------------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------------------------

def summarize(reports, root=None):
    """Lines summarizing a list of reports: the counts, then each file with problems and its lines.

    Parameters
    ----------
    reports : list of FileReport
    root : str or Path, optional
        If given, file names are shown relative to it.

    Returns
    -------
    list of str
    """
    n_err = sum(not r.ok for r in reports)
    n_warn = sum(r.ok and bool(r.warnings) for r in reports)
    out = [f'{len(reports)} files checked: {len(reports) - n_err} can be processed '
           f'({n_warn} with warnings), {n_err} will be skipped']
    for r in reports:
        if r.errors or r.warnings:
            name = str(Path(r.path).relative_to(root)) if root else r.path
            out.append(f'  {name}{"  -> SKIPPED" if not r.ok else ""}')
            out.extend(f'      {line}' for line in r.lines())
    return out


def cli(argv=None):
    """Command-line entry point (``ecco-insitu validate``): check a folder of input files, or one
    climatology file with ``--climatology FILE`` (``--prebaked`` for a pre-baked one).

    Parameters
    ----------
    argv : list of str, optional
        Arguments (default: ``sys.argv[1:]``).

    Returns
    -------
    int
        0 if every file can be used, 1 otherwise.
    """
    ap = argparse.ArgumentParser(description='Check NCEI chain input files (see docs/input_contract.md).')
    ap.add_argument('input_dir', nargs='?', help='folder of input profile files (searched recursively)')
    ap.add_argument('--climatology', metavar='FILE', help='check a climatology file instead')
    ap.add_argument('--prebaked', action='store_true', help='the climatology file is a pre-baked one')
    a = ap.parse_args(argv)
    if a.climatology:
        rep = check_climatology_file(a.climatology, prebaked=a.prebaked)
        print('\n'.join([f'{a.climatology}: {"OK" if rep.ok else "NOT USABLE"}'] + [f'  {x}' for x in rep.lines()]))
        return 0 if rep.ok else 1
    if not a.input_dir:
        ap.error('give an input folder or --climatology FILE')
    reps = check_input_dir(a.input_dir)
    print('\n'.join(summarize(reps, a.input_dir)))
    return 0 if all(r.ok for r in reps) else 1


if __name__ == '__main__':
    sys.exit(cli())
