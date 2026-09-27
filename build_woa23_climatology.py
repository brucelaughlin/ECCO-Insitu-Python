"""
Build a WOA23 monthly T/S climatology for the NCEI chain.

Boss's spec (Sept 2026):
  - Source: World Ocean Atlas 2023, 1991-2020 normal ('decav91C0'), 1-degree,
    monthly, temperature + salinity.
  - WOA23 't_an' is IN-SITU temperature; the chain compares against POTENTIAL
    temperature, so step one is converting T to potential temperature.
  - Output: a single dataset with two variables, monthly potential temperature
    and salinity, on the WOA grid.

Output variable names (potential_T_monthly, S_monthly, lon, lat, depth) match
what step03.py already reads, so the built file drops into the existing chain.

Potential-temperature conversion reuses step06.calculate_potential_T verbatim
(by reshaping each WOA (lat,lon) cell into a pseudo-profile), so the climatology
is converted with the EXACT same UNESCO algorithm as the observations.

Downloads ~24 files (~1.4 GB) to woa23_climatology/raw/ (kept for provenance;
delete manually if you want the space back).

Run:
    python build_woa23_climatology.py
"""

import numpy as np
import xarray as xr
from pathlib import Path
import urllib.request
import urllib.error
import time

import step06  # reuse calculate_potential_T verbatim

# ==============================================================================
# Config
# ==============================================================================

DECADE     = 'decav91C0'      # 1991-2020 climate normal
RES        = '1.00'           # 1-degree grid
GRID_CODE  = '01'             # filename suffix for 1-degree
OUT_DIR    = Path('/Users/brucel/ecco/yip/woa23_climatology')
RAW_DIR    = OUT_DIR / 'raw'
OUT_FILE   = OUT_DIR / f'woa23_{DECADE}_TS_clim_potential_T_1deg_fulldepth.nc'

BASE = 'https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA'

# WOA "time period" codes: 01-12 = months; 13-16 = seasons (JFM, AMJ, JAS, OND).
# Monthly fields reach 1500 m (57 levels); seasonal reach 5500 m (102 levels).
SEASON_CODES = [13, 14, 15, 16]
# Season center day-of-year (mid-season), for interpolating seasons -> months.
# JFM~mid-Feb, AMJ~mid-May, JAS~mid-Aug, OND~mid-Nov. Wraps OND<->JFM.
_SEASON_CENTER_DOY = np.array([46.0, 135.0, 227.0, 319.0])
_MID_MONTH_DOY = np.array([15, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349],
                          dtype=float)
_YEAR_LEN = 365.0

def _url(var_letter, code):
    """code: 1-12 monthly, 13-16 seasonal. Filename pattern is identical."""
    folder = 'temperature' if var_letter == 't' else 'salinity'
    fname = f'woa23_{DECADE}_{var_letter}{code:02d}_{GRID_CODE}.nc'
    return f'{BASE}/{folder}/netcdf/{DECADE}/{RES}/{fname}', fname

# ==============================================================================
# Download
# ==============================================================================

def _expected_size(url):
    """Content-Length from a HEAD-like request, or None if unavailable."""
    try:
        req = urllib.request.Request(url, method='HEAD')
        with urllib.request.urlopen(req, timeout=60) as r:
            cl = r.headers.get('Content-Length')
            return int(cl) if cl else None
    except Exception:
        return None


def _download_one(url, dest, tries=4):
    """Download with retries; verify against Content-Length; keep only a
    complete file. Resumes across script runs by skipping already-complete files."""
    expected = _expected_size(url)
    if dest.exists() and expected and dest.stat().st_size == expected:
        print(f"  have {dest.name} ({expected/1e6:.0f} MB)")
        return
    for attempt in range(1, tries + 1):
        try:
            print(f"  downloading {dest.name} (try {attempt}/{tries}) ...", flush=True)
            tmp = dest.with_suffix(dest.suffix + '.part')
            urllib.request.urlretrieve(url, tmp)
            size = tmp.stat().st_size
            if expected and size != expected:
                raise IOError(f"size mismatch: got {size} expected {expected}")
            # final integrity check: does it open as NetCDF with t_an/s_an?
            ds = xr.open_dataset(tmp, decode_times=False); ds.close()
            tmp.replace(dest)
            print(f"    ok ({size/1e6:.0f} MB)")
            return
        except Exception as e:
            print(f"    attempt {attempt} failed: {e}")
            try: tmp.unlink()
            except FileNotFoundError: pass
            if attempt < tries:
                time.sleep(5 * attempt)   # backoff
    raise RuntimeError(f"failed to download {dest.name} after {tries} tries")


def download_all():
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    for var in ('t', 's'):
        for code in list(range(1, 13)) + SEASON_CODES:   # months + seasons
            url, fname = _url(var, code)
            _download_one(url, RAW_DIR / fname)
    print("  all WOA23 files present")

# ==============================================================================
# Load + stack fields
# ==============================================================================

def _load_field(var_letter, code, an_name):
    """Load one WOA field: return (arr (depth,lat,lon), depth, lat, lon)."""
    _, fname = _url(var_letter, code)
    ds = xr.open_dataset(RAW_DIR / fname, decode_times=False)
    arr = ds[an_name].isel(time=0).values
    depth = ds['depth'].values.astype(float)
    lat = ds['lat'].values.astype(float)
    lon = ds['lon'].values.astype(float)
    ds.close()
    return arr, depth, lat, lon


def load_monthly_stack(var_letter, an_name):
    """Return (12, depth, lat, lon) monthly stack + depth/lat/lon (from month 1)."""
    months, depth, lat, lon = [], None, None, None
    for month in range(1, 13):
        arr, d, la, lo = _load_field(var_letter, month, an_name)
        months.append(arr)
        if depth is None:
            depth, lat, lon = d, la, lo
    return np.stack(months, axis=0), depth, lat, lon


def load_seasonal_stack(var_letter, an_name):
    """Return (4, depth_full, lat, lon) seasonal stack (JFM,AMJ,JAS,OND) +
    the FULL-depth axis (102 levels to 5500 m)."""
    seasons, depth = [], None
    for code in SEASON_CODES:
        arr, d, la, lo = _load_field(var_letter, code, an_name)
        seasons.append(arr)
        if depth is None:
            depth = d
    return np.stack(seasons, axis=0), depth


def seasonal_to_monthly(seasonal_stack):
    """Interpolate a (4, depth, lat, lon) seasonal stack to (12, depth, lat, lon)
    monthly, linearly in day-of-year between season centers, wrapping OND<->JFM.
    Deep seasonality is weak, so this mostly just smooths season boundaries."""
    a = _SEASON_CENTER_DOY
    nseas, nd, nlat, nlon = seasonal_stack.shape
    out = np.full((12, nd, nlat, nlon), np.nan)
    for mi in range(12):
        doy = _MID_MONTH_DOY[mi]
        if doy < a[0]:                      # before mid-Feb -> OND..JFM wrap
            s0, s1 = 3, 0
            span = (_YEAR_LEN - a[3]) + a[0]
            w1 = (doy + (_YEAR_LEN - a[3])) / span
        elif doy >= a[3]:                   # after mid-Nov -> OND..JFM wrap
            s0, s1 = 3, 0
            span = (_YEAR_LEN - a[3]) + a[0]
            w1 = (doy - a[3]) / span
        else:
            s1 = int(np.searchsorted(a, doy, side='right'))
            s0 = s1 - 1
            w1 = (doy - a[s0]) / (a[s1] - a[s0])
        f0, f1 = seasonal_stack[s0], seasonal_stack[s1]
        both = np.isfinite(f0) & np.isfinite(f1)
        blended = np.where(np.isfinite(f0), f0, f1).astype(float)
        blended[both] = (1 - w1) * f0[both] + w1 * f1[both]
        out[mi] = blended
    return out

# ==============================================================================
# Potential-temperature conversion (reuse step06 verbatim)
# ==============================================================================

def to_potential_T(T_insitu, S, depth, lat, lon):
    """Convert in-situ T (12, nd, nlat, nlon) to potential T using the SAME
    UNESCO routine the chain applies to observations (step06.calculate_potential_T).

    Each (lat, lon) cell is treated as a pseudo-profile: iPROF = nlat*nlon,
    iDEPTH = nd. calculate_potential_T needs prof_depth (1-D over iDEPTH),
    prof_lat (1-D over iPROF), prof_T (iPROF, iDEPTH) and an explicit S arg.
    """
    nmon, nd, nlat, nlon = T_insitu.shape
    lat2d = np.broadcast_to(lat[:, None], (nlat, nlon))
    prof_lat = lat2d.reshape(-1)                       # (nlat*nlon,)

    pot = np.full_like(T_insitu, np.nan)
    for m in range(nmon):
        # (nd, nlat, nlon) -> (iPROF=nlat*nlon, iDEPTH=nd)
        T_prof = np.moveaxis(T_insitu[m], 0, -1).reshape(-1, nd)
        S_prof = np.moveaxis(S[m],        0, -1).reshape(-1, nd)
        ds = xr.Dataset(
            {
                'prof_T': (['iPROF', 'iDEPTH'], T_prof),
                'prof_lat': (['iPROF'], prof_lat),
                'prof_depth': (['iDEPTH'], depth),
            }
        )
        S_da = xr.DataArray(S_prof, dims=['iPROF', 'iDEPTH'])
        pot_prof = step06.calculate_potential_T(ds, S_da)          # (iPROF, iDEPTH)
        pot[m] = np.moveaxis(np.asarray(pot_prof).reshape(nlat, nlon, nd), -1, 0)
    return pot

# ==============================================================================
# Main
# ==============================================================================

def _splice_full_depth(monthly, seas_monthly, depth_m, depth_full):
    """Combine monthly (12, 57, ...) [0-1500m] with seasonal-derived-monthly
    (12, 102, ...) [0-5500m] into a full-depth (12, 102, ...) field: use monthly
    where available (top 1500 m), seasonal-derived below. The two share identical
    depth values for the first 57 levels (verified), so this is a clean splice at
    the 1500 m seam with no re-interpolation."""
    n_m = depth_m.size                                   # 57
    assert np.array_equal(depth_full[:n_m], depth_m), "seasonal top levels must match monthly depths"
    full = seas_monthly.copy()                           # (12, 102, ...) deep source everywhere
    full[:, :n_m, :, :] = monthly                        # overwrite top 1500 m with monthly
    return full


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("1/5 downloading WOA23 decav91C0 1-deg monthly + seasonal T and S ...")
    download_all()

    print("2/5 stacking monthly (0-1500m) fields ...")
    T_mon, depth_m, lat, lon = load_monthly_stack('t', 't_an')
    S_mon, depth_ms, lat_s, lon_s = load_monthly_stack('s', 's_an')
    assert np.array_equal(depth_m, depth_ms) and np.array_equal(lat, lat_s) and np.array_equal(lon, lon_s), \
        "monthly T and S grids differ"

    print("    stacking seasonal (0-5500m) fields and interpolating seasons -> months ...")
    T_seas, depth_full = load_seasonal_stack('t', 't_an')
    S_seas, depth_full_s = load_seasonal_stack('s', 's_an')
    assert np.array_equal(depth_full, depth_full_s), "seasonal T and S depth axes differ"
    T_seas_mon = seasonal_to_monthly(T_seas)             # (12, 102, lat, lon)
    S_seas_mon = seasonal_to_monthly(S_seas)

    print("3/5 splicing monthly (top 1500m) onto seasonal-derived deep (to 5500m) ...")
    T_insitu = _splice_full_depth(T_mon, T_seas_mon, depth_m, depth_full)
    S        = _splice_full_depth(S_mon, S_seas_mon, depth_m, depth_full)
    print(f"    full-depth T,S shape = {T_insitu.shape}  (month, depth={depth_full.size}, lat, lon), "
          f"max depth {depth_full.max():.0f} m")

    print("4/5 converting in-situ T -> potential T (UNESCO, via step06) ...")
    pot_T = to_potential_T(T_insitu, S, depth_full, lat, lon)
    delta = (T_insitu - pot_T)
    finite = np.isfinite(delta)
    print(f"    in-situ minus potential T: mean={np.nanmean(delta[finite]):.4f} "
          f"max={np.nanmax(delta[finite]):.4f} degC (expect small, >=0 at depth)")

    print("5/5 writing full-depth climatology NetCDF ...")
    out = xr.Dataset(
        {
            'potential_T_monthly': (['month', 'depth', 'lat', 'lon'], pot_T.astype(np.float32)),
            'S_monthly':           (['month', 'depth', 'lat', 'lon'], S.astype(np.float32)),
        },
        coords={
            'month': np.arange(1, 13),
            'depth': depth_full.astype(np.float32),
            'lat':   lat.astype(np.float32),
            'lon':   lon.astype(np.float32),
        },
        attrs={
            'title': 'WOA23 1991-2020 (decav91C0) monthly climatology, 1-degree, full depth',
            'source': 'World Ocean Atlas 2023, objectively-analyzed mean (t_an, s_an)',
            'depth_note': 'monthly fields 0-1500 m; below 1500 m, seasonal fields '
                          '(JFM/AMJ/JAS/OND) interpolated to monthly by day-of-year',
            'potential_T_note': 'in-situ t_an converted to potential temperature via '
                                'UNESCO 1983 (step06.calculate_potential_T)',
            'salinity_note': 's_an, PSS-78',
        },
    )
    enc = {v: {'zlib': True, 'complevel': 4} for v in ('potential_T_monthly', 'S_monthly')}
    out.to_netcdf(OUT_FILE, encoding=enc)
    print(f"    wrote {OUT_FILE}  ({OUT_FILE.stat().st_size/1e6:.1f} MB)")
    print("done.")


if __name__ == '__main__':
    main()
