# Input contract: what the NCEI chain needs from its input files

The chain checks every input file against this contract before processing any of them
(`validate_inputs.py`; run it alone with `python validate_inputs.py <input_dir>`).

- An **error** means the file can't be processed correctly. The file is **skipped**, the rest of
  the run continues, and the reason is listed in `<run>/documentation/run_summary.txt`.
- A **warning** means the file is processed but something is unusual. Warnings are listed in the
  same summary and at the top of the log.

The climatology file has its own contract (section 2).

## 1. Profile files

Format: NetCDF (classic or NetCDF-4) in the MITprof layout, as written by the MITprof MATLAB
toolbox and `csv_to_nc.py`. Every file in `Interp_Profiles` meets this.

### 1.1 Required (an error if missing or wrong)

| Item | Requirement |
|---|---|
| Dimensions | `iPROF` (profiles, length > 0) and `iDEPTH` (depth levels) |
| `prof_lon`, `prof_lat` | on `(iPROF,)`; degrees east / north |
| `prof_YYYYMMDD`, `prof_HHMMSS` | on `(iPROF,)`; date and time as numbers, e.g. 19920615 and 134500 |
| `prof_depth` | on `(iDEPTH,)`; meters, positive down, strictly increasing, within 0–11,000 m, no missing levels. One depth grid per file, shared by all its profiles. |
| Observed variables | at least one of `prof_T`, `prof_S` (the `chain.variables` setting), on `(iPROF, iDEPTH)`. A file with only T, such as XBT, is fine. |
| Weights | each observed variable `prof_V` needs `prof_Vweight` with the same dimensions. Its values are replaced by step 04 (1/σ²), unless `chain.respect_existing_zero_weights` is set, which keeps the input's zeros; but the variable must exist. A file without it can be fixed by adding `prof_Vweight` filled with 1. |

### 1.2 Conventions the chain assumes

| Item | Convention |
|---|---|
| Missing values | `-9999` (as `missing_value` / `_FillValue`) or NaN. A different fill value gives a warning. |
| `prof_T` | **potential temperature**, labeled in `long_name` (or `standard_name`) as "potential temperature". If it is labeled in-situ ("in situ temperature", `sea_water_temperature`), step 05 converts it. With **no label**, step 05 decides with a climatology check and says so loudly in the log. The decision is recorded in each output file as `prof_T:T_conversion`. To force it, use `chain.convert_insitu_T_to_potential: true / false`. |
| `prof_S` | practical salinity (psu) |
| Longitude | −180–180 or 0–360; `99999` means no position known |

### 1.3 Warnings (the file is processed; the affected profiles get zero weight or are corrected)

| Warning | What the chain does |
|---|---|
| Latitude outside −90–90, longitude outside −180–360, or no position (`lon = 99999`) | step 06 zeroes those profiles' weights |
| Date before 1950, in the future, impossible (month 13), or day missing (`YYYYMM00`: year and month known, day 00) | step 06 zeroes those profiles' weights. The day-missing ones are 18,374 CTD profiles, 1.6% of all CTD, all 1992–2002; assigning mid-month is a possible future option, not implemented. |
| Time like `SS = 60`, `MM = 60`, `24:00:00` | rolled forward to a valid date and time (e.g. 23:59:60 → 00:00:00 the next day) |
| T outside −3–40 °C or S outside 0–42 psu | processed normally; often a unit or fill-value problem, so check |
| A variable with no valid values (e.g. a CTD file with an empty `prof_S`) | processed normally; nothing from that variable survives |
| `prof_T` has no potential / in-situ label | step 05 decides by the climatology check (see above) |

On the full `Interp_Profiles` set (246 files) the check takes about 12 seconds and every file passes:

- 120 files have warnings, the most common being empty `prof_S` (77 files), times rolled forward
  (28) and dates with the day missing (12 files, 18,377 profiles);
- no file is skipped.

## 2. Climatology files

Step 03 compares every profile with a monthly climatology (`paths.climatology_file`). **Any file
that meets this contract can be used**, whether it's the project's standard one or a user's own.
The chain checks it before processing anything, and stops with the reason if it isn't usable.

| Item | Requirement |
|---|---|
| Variables | `potential_T_monthly` and `S_monthly`, dimensions `(month, depth, lat, lon)` |
| `month` | 12 (a monthly climatology) |
| `lon`, `lat` | 1-D, increasing, **evenly spaced** (step 03 interpolates bilinearly on a regular grid); lat within −90–90; lon in any 360° range (it wraps) |
| `depth` | 1-D, meters, positive down, increasing. Profile depths below the deepest level use the deepest value (a warning if it is shallower than 5000 m). |
| Temperature | **potential** temperature, °C. If the attributes say in-situ (and not potential), the file is rejected; if they say nothing, there's a warning. WOA's `t_an` is in-situ and must be converted. |
| Salinity | psu |
| Missing values | NaN, or a value declared as the variable's `_FillValue` / `missing_value` attribute. An **undeclared** fill value (e.g. −9999 or 10²⁰ with no attribute) would be used as a real temperature, so the check rejects values outside the plausible range. Gaps (land, below the bottom) are filled from the nearest valid value. |

Check a file with `python validate_inputs.py --climatology <file>`, or `python setup_climatology.py`
(which also checks the pre-baked files).

### 2.1 Getting a climatology

- **Your own file:** set `paths.climatology_file` to it, then run `python setup_climatology.py`.
- **The standard one** (what this project uses): `python setup_climatology.py --build`. It
  downloads WOA23 (1991–2020, 1°) from NOAA, about 1.4 GB, and builds it like this:
  - monthly fields to 1500 m, spliced onto seasonal fields (interpolated to months) down to
    5500 m;
  - in-situ T converted to potential T with the chain's own UNESCO routine (step 05).

  The result is written where `paths.climatology_file` points. The method is in
  `build_woa23_climatology.py`; the decade and resolution are constants at its top.
- If the climatology is missing when a run starts, the run stops before processing and names
  both options.

### 2.2 Pre-baked climatologies (optional, faster)

`python setup_climatology.py --prebake` writes copies of the climatology already interpolated onto
each depth grid found in the inputs (`*_prebaked.nc` in `paths.prebaked_clim_dir`). Step 03 then
skips the vertical interpolation.

Each pre-baked file records the climatology it was made from (`source_file`). **The chain only uses
pre-baked files made from the current `paths.climatology_file`**; others are ignored with a message.
So switching to a different climatology can't silently mix old pre-baked values with the new file.

## 3. Other data files

The grid, the geodesic bins and the σ files are checked for existence before a run (section
"Settings" of the README). Their formats are fixed by the model setup (LLC90) and are not meant to
be replaced per run.
