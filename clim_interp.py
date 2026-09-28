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
                          plon, plat, pdoy, prof_depths, *, prebaked=False)
  -> (nprof, ndepth)

  PREBAKED MODE: when field_months has already been depth-interpolated onto the
  observation grid (via prebake_woa23_climatology.py), pass prebaked=True.
  Time-blend + space interpolation run as normal; _depth_interp and _fallback_fill
  are skipped.  clim_depths is unused in this mode.

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


def _time_brackets_vec(pdoy):
    """Vectorized version of _time_bracket over an array of day-of-year values.
    Returns (m0, m1, w1) arrays: bracketing month indices and weight on m1."""
    a = _MID_MONTH_DOY
    pdoy = np.asarray(pdoy, float)
    m0 = np.empty(pdoy.shape, dtype=int)
    m1 = np.empty(pdoy.shape, dtype=int)
    w1 = np.empty(pdoy.shape, dtype=float)

    wrap_span = (_YEAR_LEN - a[11]) + a[0]
    # wrap region: before mid-Jan or on/after mid-Dec -> bracket (Dec, Jan)
    before = pdoy < a[0]
    after = pdoy >= a[11]
    wrap = before | after
    m0[wrap] = 11
    m1[wrap] = 0
    w1[before] = (pdoy[before] + (_YEAR_LEN - a[11])) / wrap_span
    w1[after] = (pdoy[after] - a[11]) / wrap_span

    # interior: searchsorted gives the upper bracket month
    mid = ~wrap
    b = np.searchsorted(a, pdoy[mid], side='right')     # 1..11
    aidx = b - 1
    m0[mid] = aidx
    m1[mid] = b
    w1[mid] = (pdoy[mid] - a[aidx]) / (a[b] - a[aidx])
    return m0, m1, w1


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
    """Interpolate (nprof, ndepth_clim) -> (nprof, ndepth_prof) along depth.

    Vectorized: prof_depths and clim_depths are shared 1-D axes, so the bracket
    indices and weights are the same for every profile and are computed once,
    then applied to all profiles via fancy indexing. Interior NaNs in a clim
    column propagate to NaN here and are handled by the nearest-valid fallback
    downstream (the deepest/shallowest valid value), matching the previous
    per-profile behavior after fallback.
    """
    clim_depths = np.asarray(clim_depths, float)
    prof_depths = np.asarray(prof_depths, float)
    dmin, dmax = clim_depths.min(), clim_depths.max()

    if DEPTH_METHOD == 'nearest':
        idx = np.abs(prof_depths[:, None] - clim_depths[None, :]).argmin(axis=1)  # (ndepth_prof,)
        out = clim_cols[:, idx]                                                    # (nprof, ndepth_prof)
    else:
        hi = np.clip(np.searchsorted(clim_depths, prof_depths), 1, clim_depths.size - 1)
        lo = hi - 1
        w = (prof_depths - clim_depths[lo]) / (clim_depths[hi] - clim_depths[lo])   # (ndepth_prof,)
        lo_v = clim_cols[:, lo]                                                     # (nprof, ndepth_prof)
        hi_v = clim_cols[:, hi]
        out = (1.0 - w)[None, :] * lo_v + w[None, :] * hi_v

    # out-of-range prof depths -> NaN (no extrapolation; fallback handles them)
    oor = (prof_depths < dmin) | (prof_depths > dmax)
    out = out.copy()
    out[:, oor] = np.nan
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
    clim_depths = np.asarray(clim_depths, float)
    prof_depths = np.asarray(prof_depths, float)
    cols = blended_cols_all_depths                      # (nprof, ndepth_clim)

    # Distance from every prof depth to every clim depth (shared axes, computed
    # once): (ndepth_prof, ndepth_clim).
    dist = np.abs(prof_depths[:, None] - clim_depths[None, :])

    valid = np.isfinite(cols)                           # (nprof, ndepth_clim)
    # For each (profile, prof_depth), pick the nearest clim depth that is VALID
    # for that profile: mask invalid clim depths to +inf distance, argmin.
    # dist broadcast to (nprof, ndepth_prof, ndepth_clim) would be large, so loop
    # over the (few) prof depths that actually need filling instead of profiles.
    need_any = np.isnan(prof_clim)                      # (nprof, ndepth_prof)
    cols_safe = np.where(valid, cols, np.nan)
    for k in np.where(need_any.any(axis=0))[0]:         # loop over prof-depth cols (<=102), not profiles
        rows = np.isnan(prof_clim[:, k])                # profiles needing fill at this depth
        if not rows.any():
            continue
        # masked distance for these profiles: invalid clim depths -> inf
        md = np.where(valid[rows], dist[k][None, :], np.inf)   # (nrows, ndepth_clim)
        nn = md.argmin(axis=1)                                 # nearest valid clim depth idx
        picked = cols_safe[rows, nn]
        # rows whose entire column is invalid stay NaN (argmin over all-inf -> 0,
        # but cols_safe there is nan, so picked is nan -> correct)
        prof_clim[rows, k] = picked
    return prof_clim


def _build_valid_cell_tree(field3d, lon_grid, lat_grid):
    """Build (once) a KD-tree over valid surface cells of a field, on the unit
    sphere, plus the (iy, ix) index of each. Returns (tree, flat_idx)."""
    from scipy.spatial import cKDTree
    valid_surf = np.isfinite(field3d).any(axis=0)          # (nlat, nlon)
    la = np.deg2rad(lat_grid); lo = np.deg2rad(lon_grid)
    LO, LA = np.meshgrid(lo, la)
    xs = (np.cos(LA) * np.cos(LO))[valid_surf]
    ys = (np.cos(LA) * np.sin(LO))[valid_surf]
    zs = (np.sin(LA))[valid_surf]
    flat_idx = np.argwhere(valid_surf)                     # (k, 2) -> (iy, ix)
    return cKDTree(np.column_stack([xs, ys, zs])), flat_idx


def _fallback_nearest_cell(blended_cols, field3d, lon_grid, lat_grid, plon, plat,
                           cache=None, key=None):
    """For profiles whose space interpolation returned an all-NaN column (their
    grid neighborhood is entirely land/NaN), replace with the nearest grid cell
    that has any valid data, searched on the unit sphere. Modifies/returns
    blended_cols. The KD-tree is expensive to build, so it is cached per field
    (keyed by month index) when `cache`/`key` are supplied."""
    if FALLBACK == 'none':
        return blended_cols
    bad = np.where(~np.isfinite(blended_cols).any(axis=1))[0]
    if bad.size == 0:
        return blended_cols

    if cache is not None and key in cache:
        tree, flat_idx = cache[key]
    else:
        tree, flat_idx = _build_valid_cell_tree(field3d, lon_grid, lat_grid)
        if cache is not None:
            cache[key] = (tree, flat_idx)

    pla, plo = np.deg2rad(plat[bad]), np.deg2rad(plon[bad])
    q = np.column_stack([np.cos(pla) * np.cos(plo),
                         np.cos(pla) * np.sin(plo),
                         np.sin(pla)])
    _, nn = tree.query(q, k=1)
    iy = flat_idx[nn, 0]; ix = flat_idx[nn, 1]
    blended_cols[bad] = np.moveaxis(field3d[:, iy, ix], 0, -1)   # vectorized assign
    return blended_cols


# ==============================================================================
# Public entry point
# ==============================================================================

def interpolate_climatology(field_months, lon_grid, lat_grid, clim_depths,
                            plon, plat, pdoy, prof_depths, *, prebaked=False):
    """
    field_months : (12, ndepth_clim, nlat, nlon) monthly climatology.
                   In prebaked mode: (12, ndepth_obs, nlat, nlon) — depth interp
                   already applied; ndepth_clim == ndepth_obs.
    lon_grid, lat_grid, clim_depths : 1-D grid coordinate vectors.
                   clim_depths is unused when prebaked=True.
    plon, plat   : (nprof,) profile longitudes/latitudes (deg)
    pdoy         : (nprof,) profile day-of-year (1..365/366)
    prof_depths  : (ndepth_prof,) shared profile depth axis (unused when prebaked=True)
    prebaked     : if True, skip _depth_interp and _fallback_fill (depth axis
                   of field_months already matches the observation grid)
    returns      : (nprof, ndepth_prof) interpolated climatology
    """
    plon = np.asarray(plon, float); plat = np.asarray(plat, float)
    pdoy = np.asarray(pdoy, float); prof_depths = np.asarray(prof_depths, float)
    nprof = plon.size

    # Per-profile time brackets, vectorized (no Python loop over profiles).
    m0, m1, w1 = _time_brackets_vec(pdoy)
    if TIME_METHOD == 'nearest':
        w1 = (w1 >= 0.5).astype(float)

    # Group by the (m0, m1) MONTH PAIR only — at most 12 distinct pairs, not one
    # per rounded weight. For each pair, space-interpolate the profiles against
    # BOTH monthly fields once, then blend the two interpolated columns using
    # each profile's own w1 (exact time interp, no rounding). Space-level
    # nearest-cell fallback is applied per field (KD-tree built once per field,
    # cached), not per group.
    # In prebaked mode field_months.shape[1] == ndepth_obs (not ndepth_clim),
    # so derive the working-array width from the array, not from clim_depths.
    ndepth_work = field_months.shape[1]
    clim_cols = np.full((nprof, ndepth_work), np.nan)
    pairs = {}
    for i in range(nprof):
        pairs.setdefault((int(m0[i]), int(m1[i])), []).append(i)

    field_cache = {}   # month index -> (space-interp fn is cheap; cache fallback tree)
    for (a, b), idxs in pairs.items():
        idxs = np.array(idxs)
        pl, pa = plon[idxs], plat[idxs]
        # interpolate against each of the two bracketing monthly fields
        c0 = _space_interp(field_months[a], lon_grid, lat_grid, pl, pa)
        c0 = _fallback_nearest_cell(c0, field_months[a], lon_grid, lat_grid, pl, pa, cache=field_cache, key=a)
        if b == a:
            blended = c0
        else:
            c1 = _space_interp(field_months[b], lon_grid, lat_grid, pl, pa)
            c1 = _fallback_nearest_cell(c1, field_months[b], lon_grid, lat_grid, pl, pa, cache=field_cache, key=b)
            wcol = w1[idxs][:, None]                      # (ngroup, 1), broadcasts over depth
            both = np.isfinite(c0) & np.isfinite(c1)
            # blend where both valid; else take whichever exists
            blended = np.where(np.isfinite(c0), c0, c1)
            blended = np.where(both, (1 - wcol) * c0 + wcol * c1, blended)
        clim_cols[idxs] = blended

    if prebaked:
        # Depth interp and fallback were applied offline; clim_cols is already
        # on the observation depth grid.  Return directly.
        return clim_cols
    prof_clim = _depth_interp(clim_cols, clim_depths, prof_depths)
    prof_clim = _fallback_fill(prof_clim, clim_cols, clim_depths, prof_depths)
    return prof_clim
