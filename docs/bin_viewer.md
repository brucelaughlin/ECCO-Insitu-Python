# Geodesic bin viewer (web)

A web copy of the original `geodesic_binning` plotting app (`geodesic_binning/plotting/`), with time windows of
any length. Code: `src/ecco_insitu/diagnostics/bin_viewer/` (`export.py` writes the data and the page,
`viewer.html` is the page template, `vendor/` holds deck.gl 9.0.38 and d3-delaunay 6.0.4).
Hosting: Cloudflare R2 (a static bucket; the data are open).

## What matches the original

Everything that the original draws is drawn here in the same way:
- the bin polygons (the same construction: 0.1° grid, convex hull, the >100° dateline drop, clipped to 10 m land,
  dropped below 2 % ocean, largest piece);
- the colors (PRGn cut so white is 0, cividis for the std);
- the color limits (5/95 % of all values; 25/95 % of the bins' std; both over all depth levels);
- the face opacity sqrt(n/max) clipped to 0.4–1;
- the edge width (0.01 pt per profile, capped by the smallest bin);
- the zoom factor sqrt(global area / view area) and the threshold of 10;
- the cells: one Voronoi cell per profile, y-stratified seeds, 2 Lloyd steps, values sorted south to north, at
  most 1000 per bin;
- at zoom: the bin borders, the red profile dots, and the bin-mean labels;
- the legend, the title lines, the caption, and the keys 1/2/3/4/8/9/0.

Check (2026-11-02): GLD_WOD 2017, T, level 1. The page's bin statistics agree with a direct Python computation
from the yearly file: the same 143 bins, 0 count differences, mean and std differences below 0.0007 (float16
storage).

## Deliberate differences

1. **Only values the chain kept (weight > 0).** The original used every non-missing value. Restricting the
   original's logic to weight > 0 reproduces the page exactly.
2. **Color limits per class over the whole record**, so colors stay fixed when the time window changes. The
   original computed them from the one file it showed.
3. **Time windows.** "All years", or any range of whole years (slider, the ◀ ▶ buttons, `[` `]` to shift).
4. **Cell budget.** If more than 200,000 cells would be in view, bins stay whole until you zoom further (the
   title says so). Also a "max cells/bin" selector (the original's 1000 is the default).
5. Small: land is drawn from Natural Earth 10 m simplified to 0.02°; gridline labels are simpler than cartopy's.

## Build and update

```
# data + page (about 5 min; ~1.5 GB, of which PFL ~1 GB; each level file is 1-7 MB, loaded on demand)
python -m ecco_insitu.diagnostics.bin_viewer.export \
    --run <run>/yearly --out ~/ecco/yip/bin_viewer_data/<run name>

# after editing viewer.html: only rewrite index.html
python -m ecco_insitu.diagnostics.bin_viewer.export --out <folder> --page-only

# look at it locally
cd <folder> && python -m http.server 8765      # then open http://localhost:8765/
```

## Cloudflare R2: one-time setup (by hand, in the Cloudflare dashboard)

1. Create a Cloudflare account (free). Under **R2 Object Storage**, enable R2. This needs a card on file; the free
   tier is 10 GB of storage with no transfer fees, so this site costs nothing.
2. **Create bucket**, e.g. `ecco-bin-viewer`.
3. Public access, pick one:
   - Quick: bucket → Settings → **Public Development URL** → Allow. This gives a `https://pub-….r2.dev` address.
     It is rate-limited, which is fine for a few users.
   - Proper: bucket → Settings → **Custom Domains** → connect a domain you have on Cloudflare (about $10 a year).
4. **API token**: R2 → **Manage API tokens** → Create → "Object Read & Write", limited to this bucket. Note the
   Access Key ID, the Secret Access Key, and the endpoint `https://<account id>.r2.cloudflarestorage.com`.
5. Store the token as its own aws CLI profile. Keep it apart from the team's AWS credentials:
   ```
   aws configure --profile r2        # key id, secret; region: auto; output: json
   ```
CORS is not needed, because the page and its data are served from the same bucket.

## Upload

```
python -m ecco_insitu.diagnostics.bin_viewer.export --out <folder> \
    --upload ecco-bin-viewer --prefix <run name> \
    --endpoint https://<account id>.r2.cloudflarestorage.com [--dry-run]
```
The page is then at `<public URL>/<run name>/index.html`.

The data files are uploaded as plain gzip bytes, which the page decompresses itself. The page goes up last, with
`no-cache`, so a new version shows at once. The data are cached for a day.
