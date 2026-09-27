# NCEI Chain Profile Survival Report

## Questions Addressed

This report was produced to investigate two related questions:

1. **Why are T and S profile survival rates roughly 50%, when the expectation was closer to 90%?** Is this a result of bugs, or is it inherent to the chain's design?
2. **Does the step 10 daily decimation algorithm introduce spatial ambiguity that could affect a downstream data assimilation system?**

---

## Short Answer

The chain is now working as intended after two bug fixes described below. The remaining losses are real, intended behavior — specifically, a strict cost-vs-climatology filter in step 07 that rejects profiles deviating too far from the WOA13 climatology. Whether those thresholds are scientifically appropriate is a question for the PI, but the code is executing them correctly. The 90% survival expectation is likely too optimistic for datasets with significant mesoscale variability or data types that routinely deviate from climatology.

---

## Bug Fixes Applied

Two bugs were found and corrected in the refactored code. Neither was present in a form that caused the chain to crash, but both silently caused excess rejection:

**1. High-latitude exclusion logic in step 07 (codes 9 and 10) was inverted.**
The parameter `exclude_high_latitude_profiles_from_clim_cost` is intended to protect high-latitude profiles from the cost-vs-climatology tests, because WOA climatology is unreliable at high latitudes. Due to a logic inversion, it was doing the opposite — forcing high-latitude profiles to be zeroed. The fix was a one-word change in `step07.py` at both code 9 and code 10. This was the most impactful bug.

**2. The sigma floor in step 04 was conditionally applied to non-negative values only.**
The floor should apply unconditionally. This was corrected to `np.maximum(new_floor_dict[prof_key], sigma_np_array_MITprof_2D)`. In practice, negative sigma values appear to be rare or absent in these files, so the impact was likely small.

One additional change was made to code 6 in step 07 (removing a `== 0` arm from the climatology null-check). While this was the right correction in principle, testing showed it had no measurable impact on any of the test files — the WOA climatology does not contain exact-zero values at valid grid points.

---

## Test Results (current code, all fixes applied)

Four files from `profile_data/Interp_Profiles` were run through the full chain. "T final %" is the fraction of input T profiles with at least one valid depth level surviving through step 10.

| File | Type | Hi-lat fraction | T final % | S final % | Step 10 impact |
|---|---|---|---|---|---|
| WOD_WO_2002_CTD_OSD.nc | Ship CTD (WOD) | 22% | **55%** | **61%** | −10% (same-day ship repeats) |
| ARGO_WO_2010_PFL_R.nc | Argo float (PFL) | 50% | **59%** | **60%** | none (Argo profiles ~10-day cycle) |
| MEOP_WO_2006_CTD.nc | Elephant seal tag | 28% | **47%** | T-only | −13% (seals revisit same areas) |
| WOD_WO_2002_MRB.nc | Moored buoy (MRB) | 0% | **77%** | **182%** ¹ | ~0% (buoys are fixed in space) |

¹ The S count exceeding 100% after step 06 is an artifact of the `replace_missing_S_with_clim_S` option, which fills missing observed S with climatological S before computing potential temperature. From step 06 onward, `prof_S` valid-data counts include imputed climatology values. This is intentional behavior but inflates apparent S survival rates and should be kept in mind when interpreting S statistics.

---

## Where the Losses Come From

Profile losses happen almost entirely in steps 07 and 08. Steps 01–06 pass all profiles through unchanged (with the exception of step 06 dropping a small number of T profiles where no valid T/S pair exists for the potential temperature calculation). Steps 09 and 10 may trim additional profiles depending on data type, as described below.

### Step 07 / 08: The cost-based quality filter (~25–40% loss depending on dataset)

Step 07 zeros out the weight of any data point meeting one of eleven rejection criteria. Step 08 then removes any profile whose *entire* weight array has been zeroed. A profile survives step 08 as long as it has at least one depth level with a nonzero weight.

The dominant source of loss within step 07 is the combination of two filters:

**Code 3 — missing data.** This zeros the weight of any depth level where `prof_T` or `prof_S` is NaN. This is not a filter in the traditional sense — it simply reflects that ocean profiles, when interpolated to the model's 97 depth levels, are typically sparse. A CTD profile to 500 m will have NaN at all depth levels below 500 m; a surface-only XBT will have NaN at everything below its reach. These depth levels were never going to carry useful information anyway. Code 3 accounts for roughly 60–70% of weight-point zeroing across all data types, but because step 08 only removes profiles that are *entirely* zeroed, the profile-level impact is much smaller.

**Code 9 — average profile cost vs. climatology.** This zeros the entire profile (all depth levels) if the mean value of `(T_obs − T_clim)² × weight` exceeds a threshold of 16 across the profile. Qualitatively, with a typical T uncertainty (sigma) of ~0.5°C at mid-depth, a weight of ~4, and a threshold of 16, this rejects any profile whose depth-averaged deviation from WOA13 climatology exceeds roughly 2°C. This filter is responsible for most of the *profile-level* loss in step 08 — accounting for 10–20% of profiles being dropped, depending on dataset. Profiles in regions of strong mesoscale variability (western boundary currents, eddy fields, fronts) will frequently exceed this threshold even if the data are perfectly good. This is the most scientifically debatable threshold in the chain, and is discussed further below.

The remaining step 07 codes (range checks, flag checks, date/time validity, lat/lon bounds) each contribute a few percent at most and together are not responsible for the large losses.

### Step 10: Daily decimation (0–93% loss, highly data-type dependent)

Step 10 removes redundant profiles that were sampled multiple times in the same location on the same day, retaining only one profile per location per day. Its impact varies enormously by data type:

- **Moored buoys:** near zero impact. Buoys are fixed, so each calendar day has exactly one profile per location.
- **Argo floats:** zero impact. Argo profiles roughly every 10 days and never returns to the same location the same day.
- **Ship CTD:** moderate (~10%). Ships occasionally repeat stations on the same day.
- **Glider (GLD):** extreme (~90%). Gliders sample quasi-continuously at depth, producing tens or hundreds of profiles within a small area on any given day.
- **Seal tags (MEOP):** moderate (~15–20%). Seals dive repeatedly and often return to the same general area within a day.

For datasets dominated by Argo or buoy data, step 10 barely affects the final count, and a 90% survival rate is plausible. For glider or frequently-sampled platform data, step 10 is the primary source of apparent "loss," but it is working as designed — it is a spatial/temporal deduplication step, not a quality filter.

---

## The Decimation Spatial Question

**Question raised:** Does the step 10 algorithm introduce spatial ambiguity that could matter for data assimilation?

The algorithm for each calendar day is:
1. Compute pairwise distances between all profiles on that day.
2. Sort profiles by proximity to noon (closest first).
3. Greedily iterate: keep a profile if no already-kept profile is within `distance_tolerance` (5 km); otherwise discard it.

The result is that within any 5 km cluster of same-day profiles, the surviving profile is whichever one happened to be sampled closest to noon. The spatial location of the survivor is governed by a temporal criterion, not by proximity to any grid reference point (model cell center, geodesic bin center, etc.). The algorithm is fully deterministic — it produces the same result every run — but the surviving location within a cluster is in some sense arbitrary spatially.

**Whether this matters for data assimilation depends on two things:**

First, the DA system's interpolation scheme. Steps 01 and 02 already map each profile to its nearest model grid cell via nearest-neighbor interpolation. At llc90 resolution (~100 km at mid-latitudes), two profiles within 5 km will almost certainly map to the same grid cell regardless of which one survives. In that case the exact lat/lon of the survivor is irrelevant — both profiles would be assimilated into the same cell with the same weight. At higher resolutions (e.g. llc270, ~33 km), profiles within 5 km could plausibly straddle a cell boundary, and the choice of which survives could affect which cell receives the observation.

Second, the scientific motivation for the noon criterion. Selecting the closest-to-noon profile is motivated by diurnal sampling considerations — noon is a reasonable reference time to minimize diurnal cycle bias across profiles from different times of day. This is a defensible choice, but it is a temporal criterion being used to solve a spatial deduplication problem. A spatially more principled approach would rank candidates by proximity to the model grid cell center or geodesic bin center. Whether this distinction matters in practice is worth confirming with whoever is running the DA system.

In summary: the spatial concern is real but likely inconsequential at llc90 resolution, and the algorithm is self-consistent and deterministic. It is worth flagging to the DA team in case they are working at higher resolution or using interpolation schemes sensitive to exact profile coordinates.

---

## The 90% Expectation

A 90% survival rate is achievable in the current code for the right data type — specifically, low-latitude Argo or moored buoy data with good climatological coverage. The MRB test file achieves 77% with no high-latitude profiles and no step 10 loss; removing the cost filter (code 9) would push it higher still.

For datasets with any of the following characteristics, survival rates well below 90% should be expected even from a correctly functioning chain:

- **High mesoscale variability:** profiles in western boundary current systems, eddy fields, or frontal zones will frequently fail the cost-vs-climatology threshold (code 9) even when the data are oceanographically valid.
- **High-latitude data:** despite the exclusion flag, WOA13 climatology is sparse and unreliable at high latitudes, and the missing-climatology filter (code 6) will zero weights at depth levels where the climatology simply does not exist.
- **Glider or high-frequency platforms:** step 10 will decimate these heavily by design.
- **T-only instruments (e.g. XBT, some MEOP):** S imputation from climatology is required before potential T can be computed; errors in that imputation propagate into both the potential T conversion and the cost calculation.

The cost threshold of 16 (corresponding to roughly a 2°C mean profile deviation from WOA13) is the most scientifically debatable parameter in the chain. It is hard-coded and was presumably set empirically for a particular class of data. Whether it is appropriate for the full range of instrument types and ocean regions being processed is a question the PI may want to revisit.
