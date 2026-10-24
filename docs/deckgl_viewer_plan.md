# Plan: deck.gl profile viewer (1-1 web copy of `geodesic_binning/` plots)

Status: planned 2026-10-23, nothing built yet.

## Goal
Make a web viewer that reproduces the offline tool in `geodesic_binning/plotting` (and its
binning step), extended from one yearly file to any time window, up to the whole dataset. It has no
hover and no time series. It sits alongside `profile_analysis/app_profile_map.py` and does not
replace it, which keeps the time series and source filtering.

## What to reproduce from the original
- **Bins:** a T or S anomaly map at a chosen depth (value − clim).
  - fill = bin mean;
  - edge color = bin std;
  - edge width grows with profile count.
- **Ocean clipping:** bins are clipped to ocean with the 10 m land mask; bins under 2% ocean are dropped; for
  MultiPolygons, the largest piece is kept.
- **Past a zoom threshold:**
  - each bin splits into Voronoi cells, one per profile;
  - the values are sorted ascending and assigned south to north by cell centroid;
  - every profile is a red dot;
  - each bin shows its mean as a label;
  - thin outlines are drawn around the parent bins.
- **Controls** (originally keys):
  - depth down/up (1/2), variable back/forward (3/4), reset zoom (9), reset depth and zoom (0);
  - all of these become UI buttons.

## Why Plotly/Dash couldn't do it
- Plotly geo choropleths don't use WebGL, so they slow down after a few thousand shapes.
- One cell per profile means tens to hundreds of thousands of cells in view.
- A WebGL map library (deck.gl) draws that many smoothly.

## Architecture
**1. Preprocessing**, run once per NCEI run (roughly 10–20 min). It writes to its own output folder:
- **Bin polygons:**
  - clipped to the ocean;
  - projected into World, Arctic and Southern Ocean coordinates in Python, because deck.gl has no
    polar stereographic view, so the browser draws flat x/y in every projection;
  - coastlines projected the same way;
  - saved as small binary files.
- **Per-profile zarr cache**, about 2.3 GB, for all 5.9 M profiles:
  - lon, lat, time, bin ID and source;
  - T and S anomalies at all 97 depths (float16, chunked by depth);
  - profiles ordered by bin, so any bin + time window is one slice lookup.

**2. Python server** (Flask). Endpoints:
- **files and runs:** the available runs, sources and the year range.
- **positions:** profile positions for a time window, as a binary array.
  - 5.9 M dots are about 50 MB, sent once per window.
- **bin stats:** mean, std and count per bin for a variable, depth and window, computed with
  `np.bincount` in well under 1 s.
- **cells:** Voronoi cells and their values for a list of bins in view.
  - reuses `profile_analysis/subpolygons.py` (gnomonic projection, bounded loops, south-to-north
    ordering);
  - cached per (window, variable, depth, bin).

**3. Browser page** (one HTML file + deck.gl, about 400–600 lines of JS):
- **Layers:**
  - bin fill;
  - std edges (width by count);
  - sub-cells and red dots past the zoom threshold;
  - parent outlines;
  - coastlines;
  - bin-mean labels.
- **Controls:** run, T/S, a depth slider with ◀/▶, year range, projection, reset.
- **Colorbars:** on both sides (mean and std).

## Time windows and the cell budget
- Use a budget for the whole view, about 150,000 cells, instead of a fixed cap per bin (the
  original capped at 1000).
  - **Within budget** (zoomed in enough): every bin gets one cell per profile.
  - **Over budget:** show the bin outlines and dots only, or quantile cells as the current app does,
    until you zoom in further.
- The busiest bin has ~45,000 profiles over the whole dataset. Its cells take seconds the first
  time and are cached after that.

## Sizes, for reference
- 241 yearly files, 5,938,524 profiles (run 20261021_100417).
- Per file: median ~17,000 profiles, largest ~152,000 (Argo 2016).
- The original tool's pickle was 5 GB for one CTD year. It stored matplotlib objects for every
  variable and depth; the new design computes only what's in view.

## Effort
- The single-file version is a few working sessions; most of it is the JS page, colorbars and line
  widths.
- Time windows up to the whole dataset add about 30–50%:
  - the preprocessing script;
  - a year-range control;
  - the budget rule.

## Open decisions
1. Cell budget value (~150k suggested).
2. When over budget: outlines + dots only, or quantile cells.
3. Whether to filter by source too, or only by time window.

## Hosting (to-do, discussed 2026-10-23, nothing set up)
Hosting is optional; the local version needs no cloud at all.
- **Static hosting:** the data sits in object storage, and the browser reads chunks and computes the
  cells itself (Voronoi ported to JS with d3-delaunay), so no server is needed.
- **Never host from the team's shared buckets, including `my-bucket`.** Opening even one folder means
  editing that bucket's policy. Use a separate bucket that holds only the viewer files, copied in by
  hand, with no link to the other buckets.
- **Isolation, best first:**
  1. a separate AWS account; a member account under the team's AWS Organization would be billed to the
     team and overseen by its admins;
  2. Cloudflare R2: free up to 10 GB, no transfer fees, needs a custom domain (about $10 a year)
     and Cloudflare Access for a login; it sits outside the team's oversight;
  3. a separate bucket in the team account.
- **AWS cost:** about 5 cents a month for storage; transfer is likely within CloudFront's free tier
  for a small team.
- **Security:** CloudFront with Origin Access Control, block public access on the bucket, a
  WAF rate limit, CORS limited to the viewer's page, a billing alert, and optionally signed
  cookies or a Cognito login.
- **Ask the team first:**
  - Is a new bucket or account allowed?
  - May the processed profiles be shared publicly or behind a login?
  - Who pays for and looks after the hosting?
  - Is there a preferred route, such as an institutional server or data portal?
