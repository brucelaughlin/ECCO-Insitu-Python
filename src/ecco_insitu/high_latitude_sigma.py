"""High-latitude T / S uncertainties (sigma) for step 04.

Poleward of ``high_latitude_sigma.min_lat`` (60 degrees), step 04 uses these instead of the CTD
sigma. Step 06 then applies its climatology-cost tests to those profiles, which it otherwise skips
poleward of 60 degrees. The fields (e.g. ``Theta_rmse_se`` / ``S_rmse_se``, "expected rmse of in situ
T / S from the larger of the WOA23-based and AOTD-based estimates") are on an EASE2 100 km polar
grid, ``(depth, y, x)`` with 2-D lat / lon, and are defined only in sampled cells.

For each profile with ``|lat| > min_lat``:

- find the nearest grid cell that has data (great-circle distance);
- if it is within ``max_dist_km``, interpolate that cell's sigma profile linearly to the profile's
  depths (constant beyond the cell's shallowest / deepest value, as step 04 does for the CTD sigma
  files);
- otherwise leave the profile on the CTD sigma (not "covered").

The files used today cover 59.4-89.4 N only, so southern profiles are never covered here; step 04
gives them a constant sigma instead (``high_latitude_sigma.south_constant`` in the settings). A file
with southern values would cover southern profiles near its cells the same way as in the north.
"""

from pathlib import Path

import numpy as np
import xarray as xr
from scipy.spatial import cKDTree

_R_EARTH_KM = 6371.0
_CACHE = {}


def _xyz(lat, lon):
    la, lo = np.radians(np.asarray(lat, float)), np.radians(np.asarray(lon, float))
    return np.column_stack((np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)))


class HighLatSigma:
    """One sigma field (one variable of one file), loaded once.

    Parameters
    ----------
    path : str or Path
        NetCDF file with the field on ``(depth, y, x)`` and 2-D ``lat`` / ``lon``.
    var : str
        Variable name, e.g. ``Theta_rmse_se``.

    Attributes
    ----------
    depth : numpy.ndarray
        The field's depths (m).
    profiles : numpy.ndarray, (n_cells, n_depth)
        Sigma profiles of the cells that have any value.
    cell_lat : numpy.ndarray
        Latitude of those cells.
    name : str
        ``<file name>:<var>``, for messages.
    """

    def __init__(self, path, var):
        with xr.open_dataset(path) as ds:
            field = ds[var].values.astype(np.float64)          # (depth, y, x)
            self.depth = ds['depth'].values.astype(np.float64)
            lat, lon = ds['lat'].values, ds['lon'].values      # (y, x)
        has = np.isfinite(field).any(axis=0)
        self.profiles = field[:, has].T                       # (n_cells, depth)
        self.cell_lat = lat[has]
        self._tree = cKDTree(_xyz(lat[has], lon[has]))
        self.name = f'{Path(path).name}:{var}'

    def _interp(self, cell, prof_depth):
        p = self.profiles[cell]
        ok = np.isfinite(p)
        return np.interp(prof_depth, self.depth[ok], p[ok])  # constant beyond the ends

    def sigma(self, lat, lon, prof_depth, min_lat=60.0, max_dist_km=300.0):
        """Sigma at each profile, from the nearest cell with data.

        Parameters
        ----------
        lat, lon : array_like, (n_prof,)
            Profile positions (degrees).
        prof_depth : array_like, (n_depth,)
            Profile depths (m; the sign is ignored).
        min_lat : float
            Only profiles with ``|lat| > min_lat`` are considered.
        max_dist_km : float
            The nearest cell with data must be within this distance.

        Returns
        -------
        sigma : numpy.ndarray, (n_prof, n_depth)
            Sigma (same units as the file), NaN for profiles not covered.
        covered : numpy.ndarray of bool, (n_prof,)
            Which profiles got a value.
        """
        lat, lon = np.asarray(lat, float), np.asarray(lon, float)
        prof_depth = np.abs(np.asarray(prof_depth, float))
        out = np.full((lat.size, prof_depth.size), np.nan)
        covered = np.zeros(lat.size, bool)
        want = np.flatnonzero(np.abs(lat) > min_lat)
        if want.size == 0 or len(self.profiles) == 0:
            return out, covered
        chord, cell = self._tree.query(_xyz(lat[want], lon[want]))
        dist_km = 2 * _R_EARTH_KM * np.arcsin(np.clip(chord / 2, 0, 1))
        ok = dist_km <= max_dist_km
        for c in np.unique(cell[ok]):
            rows = want[ok & (cell == c)]
            out[rows] = self._interp(c, prof_depth)
        covered[want[ok]] = True
        return out, covered


def get(path, var):
    """The :class:`HighLatSigma` for ``(path, var)``, loaded once per process and cached (each chain
    worker process loads a field once).

    Parameters
    ----------
    path : str or Path
    var : str

    Returns
    -------
    HighLatSigma
    """
    key = (str(path), var)
    if key not in _CACHE:
        _CACHE[key] = HighLatSigma(path, var)
    return _CACHE[key]
