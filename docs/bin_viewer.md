# Geodesic bin viewer (web)

A web copy of the original `geodesic_binning` plotting app (`geodesic_binning/plotting/`), with time windows of
any length. Code: `src/ecco_insitu/diagnostics/bin_viewer/` (`export.py` writes the data and the page,
`viewer.html` is the page template, `vendor/` holds deck.gl 9.0.38 and d3-delaunay 6.0.4).
Hosting (since 2026-11-20): the page and its main data on GitHub Pages, everything bigger on S3 + CloudFront:
- **`brucelaughlin/ecco-bin-viewer`** (GitHub Pages): the page, bins, maps, climatology anomalies and platforms.
  **Live:** https://brucelaughlin.github.io/ecco-bin-viewer/run_20261027_142247/
- **S3 bucket `ecco-insitu-site`** (the user's own AWS account, us-east-2, private) behind the CloudFront
  distribution **`d1r1rp2b701s7u.cloudfront.net`** (Origin Access Control: only CloudFront reads the bucket):
  the profile details and the model anomalies (V4r6, V5r1), under `bin_viewer/`. The page reads them through
  `detail_urls` and `ref_urls` in `index.json`. New `*.cloudfront.net` addresses get through the facility's
  secure network (tested 2026-11-20); `*.r2.dev` and new domains don't.

CORS: the distribution uses its own response-headers policy (`ecco-insitu-cors-any-header`: any origin, any
request header, GET/HEAD/OPTIONS). AWS's managed SimpleCORS policy is not enough: it leaves out
`Access-Control-Allow-Origin` whenever the request carries a header beyond the simple ones, and Chrome always
sends `priority`.

One run at a time: each publish replaces what is there. Older runs stay in their local folders.

History: 2026-11-02 first on Cloudflare R2 (blocked on the secure network), then GitHub Pages only, split over
seven repositories to stay under the 1 GB per-site limit; 2026-11-20 the details and the model levels moved to
S3 + CloudFront, so only `ecco-bin-viewer` remains on GitHub Pages.

## Where the data are (since 2026-11-20)

| place | holds | about |
|---|---|---|
| GitHub Pages `ecco-bin-viewer` | the page, bins, maps, climatology anomalies, platforms | 830 MB |
| S3/CloudFront `bin_viewer/run_<name>/detail/` | profile details, all classes | 1.85 GB |
| S3/CloudFront `bin_viewer/model/run_<name>/` | anomalies vs V4r6 | 774 MB |
| S3/CloudFront `bin_viewer/V5r1/run_<name>/` | anomalies vs V5r1 | 776 MB |

The former sites `ecco-bin-viewer-model`, `-v5r1` and `-detail-a/b/c` are no longer read.

Profile details (`export --details`): every profile's observation, climatology, model, σ and QC code (step 06's
`prof_Vweight_code`), date, time and description, to 0.001 (σ 0.0001), in chunks of 2000 profiles, loaded only when a
profile is opened. Checked against the yearly files for 400 random profiles per class (178,000 values; script
`tests/tools/check_bin_viewer_detail.py`): all within rounding, QC codes, dates and descriptions exact.

**Publish (since 2026-11-20):** the details and model levels to S3, then the main site (the S3 step writes their
addresses into `index.json`, which the main site carries). Needs the AWS CLI signed in to the bucket's account
(`aws login --profile ecco-site`; that profile needs `region = us-east-2` in `~/.aws/config`).
```
python -m ecco_insitu.diagnostics.bin_viewer.export --out <folder> \
    --s3 s3://ecco-insitu-site/bin_viewer https://d1r1rp2b701s7u.cloudfront.net/bin_viewer --aws-profile ecco-site \
    --github git@github.com:brucelaughlin/ecco-bin-viewer.git
```
Add `--dry-run` to see what would be uploaded. A new model: `--add-model NAME DIR <folder>_NAME_site` writes its
level files there, and the next `--s3` uploads them and adds it to `ref_urls`.

Before 2026-11-20 the details went to three GitHub Pages sites (no longer used):
```
python -m ecco_insitu.diagnostics.bin_viewer.export --out <folder> \
    --github-detail git@github.com:brucelaughlin/ecco-bin-viewer-detail-a.git=PFL:0-653 \
                    git@github.com:brucelaughlin/ecco-bin-viewer-detail-b.git=PFL:654- \
                    git@github.com:brucelaughlin/ecco-bin-viewer-detail-c.git=CTD_WOD,GLD_WOD,ITP,MEOP,MRB_WOD,XBT_WOD \
    --github git@github.com:brucelaughlin/ecco-bin-viewer.git
```

**The profile panel:** zoomed in, click a red dot (or a cyan platform ring). Tabs: *this profile* (observation,
climatology, model, gray neighbors from the same bin within ±2 years and other platforms; obs − reference with ±1σ
and ±2σ bands; rejected levels ×; the values at the map's depth), *this platform over time* (all its profiles as a
depth-time picture of obs − reference) and *largest misfits here* (the bin's profiles ranked by |Δ|/σ at the
map's depth; click a row to open one).

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
3. **Time windows.** "All years", or any range of whole years (the slider and the ◀ ▶ buttons).
4. **Cell budget.** If more than 200,000 cells would be in view, bins stay whole until you zoom further (the
   title says so). Also a "max cells/bin" selector (the original's 1000 is the default).
5. **Edge widths** (2026-11-02): 0.1 px × √(profile count), growing with √(zoom factor), at most 25 % of the
   bin's own size. The original's 0.01 pt per profile, capped by the smallest bin on the map, made every bin
   above ~60 profiles look the same on the world map; now busy regions (e.g. XBT shipping lanes) stand out.
6. **Polar maps** below the world map (Arctic, Southern Ocean; polar stereographic, the 40° circle at the
   edges), with their own bins built the original's way in polar coordinates, so bins at the poles and across
   the dateline are whole.
7. **Depth ranges** (standard bands, 100 m bands, 15 or 30 groups), as in the cost map explorer: each profile
   counts once, with the mean of its values at the range's levels. Checked against Python (ITP, 250–500 m band:
   0 count differences, mean differences < 10⁻⁷).
8. **Anomaly vs the model** as well as the climatology, from the model equivalents the cost analyses use; model
   values of 0 or missing are left out (as in the costs).
9. **Isobaths** (Natural Earth, 200–5000 m) with depth labels; a "← report" button; a home button; the view
   cannot be panned off the map. No keyboard shortcuts.
10. Small: land is drawn from Natural Earth 10 m simplified to 0.02°; gridline labels are simpler than cartopy's.

The zoom threshold for cells, dots and labels is the original's: zoom factor √(world area / view area) above 10
(on the world map, a view of about 36° × 18°); the polar maps use the same factor, so the switch happens at the
same on-screen bin size.

## Build and update

```
# data + page (about 9 min with model values; ~1.55 GB, split over two sites below; level files 1-4 MB, loaded on demand)
python -m ecco_insitu.diagnostics.bin_viewer.export \
    --run <run>/yearly --out ~/ecco/yip/bin_viewer_data/<run name>

# a folder written before 2026-11-02 (format 1): rewrite as format 2 (lossless, about half the size)
python -m ecco_insitu.diagnostics.bin_viewer.export --out <folder> --convert

# after editing viewer.html: only rewrite index.html
python -m ecco_insitu.diagnostics.bin_viewer.export --out <folder> --page-only

# look at it locally
cd <folder> && python -m http.server 8765      # then open http://localhost:8765/
```

## Publish to GitHub Pages

```
B=~/ecco/yip/bin_viewer_data; R=<run name>
# the data, with model anomalies; the report it belongs to (its "← report" button)
python -m ecco_insitu.diagnostics.bin_viewer.export --run <run>/yearly --out $B/$R \
    --model-equivalents <model equivalents>/netCDF
python -m ecco_insitu.diagnostics.bin_viewer.export --out $B/$R \
    --report https://brucelaughlin.github.io/ECCO-Insitu-Python/<report name>/
# move the model files to their own site's folder ($B/${R}_model_site/), and point the page at that site
python -m ecco_insitu.diagnostics.bin_viewer.export --out $B/$R \
    --split-model model https://brucelaughlin.github.io/ecco-bin-viewer-model/$R
# publish both sites (add --dry-run to check sizes first)
python -m ecco_insitu.diagnostics.bin_viewer.export --out $B/$R \
    --github git@github.com:brucelaughlin/ecco-bin-viewer.git \
    --github-model git@github.com:brucelaughlin/ecco-bin-viewer-model.git
```
Each publish replaces the repository's `main` branch with one commit (no history, so nothing piles up). GitHub
Pages is already switched on for both repositories (Settings → Pages → Deploy from a branch → `main`, `/ (root)`);
a new repository would need that once. Then rebuild the report with
`--bin-viewer https://brucelaughlin.github.io/ecco-bin-viewer/<run name>/` so its button points at the run.

GitHub's limits (docs.github.com, "GitHub Pages limits", checked 2026-11-02): published site at most 1 GB; soft
bandwidth limit 100 GB a month per site; 10 builds an hour; no stated limit on the number of project sites; not
for commercial use. Split one viewer over a few sites at most.

