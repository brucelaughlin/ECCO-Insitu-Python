"""
Optimal climatology interpolation for the NCEI chain (step03).

Interpolates a monthly gridded climatology (WOA23: potential_T_monthly /
S_monthly on a regular lon/lat/depth grid) onto observed profile
locations / depths / dates.

Scheme (each axis independently switchable via the config knobs below):
  - TIME  : linear between adjacent monthly fields by day-of-year, each month
            centered mid-month, wrapping Dec<->Jan. ('nearest' = pick month.)
  - SPACE : bilinear in lon/lat on the regular grid. ('nearest' = pick cell.)
  - DEPTH : linear onto each profile's own depth levels. ('nearest' = snap.)

Edge-case FALLBACK (config):
  Where the primary interpolation yields NaN — a profile depth beyond the
  climatology's deepest level, or a location whose surrounding cells are all
  land/NaN — fall back to the NEAREST valid climatology value. Set
  FALLBACK = 'none' to leave such points NaN (which zeros their weight later).

Public entry point:
  interpolate_climatology(field_months, lon_grid, lat_grid, clim_depths,
                          plon, plat, pdoy, prof_depths) -> (nprof, ndepth)

All behavior is parameterized so the method can change later without editing
logic: flip TIME_METHOD / SPACE_METHOD / DEPTH_METHOD / FALLBACK.
"""

import numpy as np

# ==============================================================================
# Config knobs — change these, not the code below.
# ==============================================================================

TIME_METHOD  = 'linear'     # 'linear' | 'nearest'
SPACE_METHOD = 'bilinear'   # 'bilinear' | 'nearest'
DEPTH_METHOD = 'linear'     # 'linear' | 'nearest'
FALLBACK     = 'nearest'    # 'nearest' | 'none'

# Mid-month day-of-year anchors (non-leap; fine for a monthly climatology).
_MID_MONTH_DOY = np.array([15, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349],
                          dtype=float)
_YEAR_LEN = 365.0


# ==============================================================================
# TIME: blend the 12 monthly fields into one field per profile
# ==============================================================================

def _time_bracket(doy):
    """For one day-of-year, return (m0, m1, w1): bracketing month indices (0..11)
    and the weight on m1 so value = (1-w1)*field[m0] + w1*field[m1]. Wraps."""
    a = _MID_MONTH_DOY
    if doy < a[0]:
        span = (_YEAR_LEN - a[11]) + a[0]
        return 11, 0, (doy + (_YEAR_LEN - a[11])) / span
    if doy >= a[11]:
        span = (_YEAR_LEN - a[11]) + a[0]
        return 11, 0, (doy - a[11]) / span
    b = int(np.searchsorted(a, doy, side='right'))
    return b - 1, b, (doy - a[b - 1]) / (a[b] - a[b - 1])


def _blended_field(field_months, doy):
    """Return a single (depth, nlat, nlon) field for a given day-of-year, blending
    the two bracketing monthly fields (NaN-aware). 'nearest' picks one month."""
    m0, m1, w1 = _time_bracket(doy)
    if TIME_METHOD == 'nearest':
        return field_months[m1 if w1 >= 0.5 else m0]
    f0, f1 = field_months[m0], field_months[m1]
    # NaN-aware blend: where one month is NaN, use the other; where both, blend.
    both = np.isfinite(f0) & np.isfinite(f1)
    out = np.where(np.isfinite(f0), f0, f1).astype(float)      # start from whichever exists
    out[both] = (1 - w1) * f0[both] + w1 * f1[both]
    return out


# ==============================================================================
# SPACE: interpolate a (depth, nlat, nlon) field onto profile lon/lat
# ==============================================================================

def _space_interp(field3d, lon_grid, lat_grid, plon, plat):
    """Return (nprof, ndepth_clim): field sampled at each (plon, plat).
    Bilinear (NaN-aware, renormalized over valid corners) or nearest."""
    nlon, nlat = lon_grid.size, lat_grid.size
    dlon = lon_grid[1] - lon_grid[0]
    dlat = lat_grid[1] - lat_grid[0]
    lon0 = lon_grid[0]
    plon_w = ((plon - lon0) % 360.0) + lon0     # into grid's lon frame
    nd = field3d.shape[0]
    nprof = plon.size

    if SPACE_METHOD == 'nearest':
        ix = np.clip(np.round((plon_w - lon_grid[0]) / dlon).astype(int) % nlon, 0, nlon - 1)
        iy = np.clip(np.round((plat - lat_grid[0]) / dlat).astype(int), 0, nlat - 1)
        return np.moveaxis(field3d[:, iy, ix], 0, -1)

    fx = (plon_w - lon_grid[0]) / dlon
    fy = (plat - lat_grid[0]) / dlat
    ix0 = np.floor(fx).astype(int)
    iy0 = np.floor(fy).astype(int)
    tx, ty = fx - ix0, fy - iy0
    ix0m, ix1 = ix0 % nlon, (ix0 + 1) % nlon       # wrap lon seam
    iy0c = np.clip(iy0, 0, nlat - 1)
    iy1c = np.clip(iy0 + 1, 0, nlat - 1)
    w = ((1 - tx) * (1 - ty), tx * (1 - ty), (1 - tx) * ty, tx * ty)

    out = np.full((nprof, nd), np.nan)
    for d in range(nd):
        f = field3d[d]
        corners = (f[iy0c, ix0m], f[iy0c, ix1], f[iy1c, ix0m], f[iy1c, ix1])
        num = np.zeros(nprof); den = np.zeros(nprof)
        for c, wc in zip(corners, w):
            v = np.isfinite(c)
            num[v] += c[v] * wc[v]
            den[v] += wc[v]
        good = den > 0
        out[good, d] = num[good] / den[good]
    return out


# ==============================================================================
# DEPTH: (nprof, ndepth_clim) -> (nprof, ndepth_prof) along a shared depth axis
# ==============================================================================

def _depth_interp(clim_cols, clim_depths, prof_depths):
    nprof = clim_cols.shape[0]
    out = np.full((nprof, prof_depths.size), np.nan)
    for i in range(nprof):
        col = clim_cols[i]
        v = np.isfinite(col)
        if not v.any():
            continue
        cd, cv = clim_depths[v], col[v]
        if DEPTH_METHOD == 'nearest':
            in_range = (prof_depths >= cd.min()) & (prof_depths <= cd.max())
            idx = np.abs(prof_depths[:, None] - cd[None, :]).argmin(axis=1)
            out[i, in_range] = cv[idx][in_range]
        else:
            out[i] = np.interp(prof_depths, cd, cv, left=np.nan, right=np.nan)
    return out


# ==============================================================================
# FALLBACK: fill remaining NaN with nearest valid climatology
# ==============================================================================

def _fallback_fill(prof_clim, blended_cols_all_depths, clim_depths, prof_depths):
    """For any NaN left in prof_clim, use nearest valid depth from that profile's
    own (already space-interpolated) clim column. This covers depths beyond the
    clim's deepest valid level (extends the deepest/shallowest valid value) and
    profiles where linear depth interp left interior gaps. Space-level gaps
    (all-land neighborhood) are handled by _fallback_nearest_cell upstream."""
    if FALLBACK == 'none':
        return prof_clim
    nprof = prof_clim.shape[0]
    for i in range(nprof):
        need = np.isnan(prof_clim[i])
        if not need.any():
            continue
        col = blended_cols_all_depths[i]
        v = np.isfinite(col)
        if not v.any():
            continue
        cd, cv = clim_depths[v], col[v]
        idx = np.abs(prof_depths[need][:, None] - cd[None, :]).argmin(axis=1)
        prof_clim[i, need] = cv[idx]
    return prof_clim


def _fallback_nearest_cell(blended_cols, field3d, lon_grid, lat_grid, plon, plat):
    """For profiles whose space interpolation returned an all-NaN column (their
    grid neighborhood is entirely land/NaN), replace with the nearest grid cell
    that has any valid data, searched on the unit sphere. Modifies/returns
    blended_cols (nprof, ndepth_clim)."""
    if FALLBACK == 'none':
        return blended_cols
    bad = np.where(~np.isfinite(blended_cols).any(axis=1))[0]
    if bad.size == 0:
        return blended_cols

    # valid surface cells (any depth finite) -> unit-sphere coords for KD-tree
    valid_surf = np.isfinite(field3d).any(axis=0)          # (nlat, nlon)
    la = np.deg2rad(lat_grid); lo = np.deg2rad(lon_grid)
    LO, LA = np.meshgrid(lo, la)
    xs = (np.cos(LA) * np.cos(LO))[valid_surf]
    ys = (np.cos(LA) * np.sin(LO))[valid_surf]
    zs = (np.sin(LA))[valid_surf]
    flat_idx = np.argwhere(valid_surf)                     # (k, 2) -> (iy, ix)
    from scipy.spatial import cKDTree
    tree = cKDTree(np.column_stack([xs, ys, zs]))

    pla, plo = np.deg2rad(plat[bad]), np.deg2rad(plon[bad])
    q = np.column_stack([np.cos(pla) * np.cos(plo),
                         np.cos(pla) * np.sin(plo),
                         np.sin(pla)])
    _, nn = tree.query(q, k=1)
    for j, i in enumerate(bad):
        iy, ix = flat_idx[nn[j]]
        blended_cols[i] = field3d[:, iy, ix]
    return blended_cols


# ==============================================================================
# Public entry point
# ==============================================================================

def interpolate_climatology(field_months, lon_grid, lat_grid, clim_depths,
                            plon, plat, pdoy, prof_depths):
    """
    field_months : (12, ndepth_clim, nlat, nlon) monthly climatology
    lon_grid, lat_grid, clim_depths : 1-D grid coordinate vectors
    plon, plat   : (nprof,) profile longitudes/latitudes (deg)
    pdoy         : (nprof,) profile day-of-year (1..365/366)
    prof_depths  : (ndepth_prof,) shared profile depth axis
    returns      : (nprof, ndepth_prof) interpolated climatology
    """
    plon = np.asarray(plon, float); plat = np.asarray(plat, float)
    pdoy = np.asarray(pdoy, float); prof_depths = np.asarray(prof_depths, float)
    nprof = plon.size

    # Group profiles by their (m0, m1, rounded-w1) time bracket so we build each
    # blended field once instead of per profile. Round w1 to 0.01 for grouping.
    keys = np.empty(nprof, dtype=object)
    brackets = {}
    for i in range(nprof):
        m0, m1, w1 = _time_bracket(pdoy[i])
        if TIME_METHOD == 'nearest':
            w1 = 1.0 if w1 >= 0.5 else 0.0
        key = (m0, m1, round(float(w1), 2))
        keys[i] = key
        brackets.setdefault(key, []).append(i)

    clim_cols = np.full((nprof, clim_depths.size), np.nan)
    for key, idxs in brackets.items():
        m0, m1, w1 = key
        if TIME_METHOD == 'nearest':
            field = field_months[m1 if w1 >= 0.5 else m0]
        else:
            f0, f1 = field_months[m0], field_months[m1]
            both = np.isfinite(f0) & np.isfinite(f1)
            field = np.where(np.isfinite(f0), f0, f1).astype(float)
            field[both] = (1 - w1) * f0[both] + w1 * f1[both]
        idxs = np.array(idxs)
        clim_cols[idxs] = _space_interp(field, lon_grid, lat_grid, plon[idxs], plat[idxs])
        # space-level fallback for all-land neighborhoods within this group
        clim_cols[idxs] = _fallback_nearest_cell(
            clim_cols[idxs], field, lon_grid, lat_grid, plon[idxs], plat[idxs])

    prof_clim = _depth_interp(clim_cols, clim_depths, prof_depths)
    prof_clim = _fallback_fill(prof_clim, clim_cols, clim_depths, prof_depths)
    return prof_clim
