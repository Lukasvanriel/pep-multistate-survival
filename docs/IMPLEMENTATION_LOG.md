# Implementation log (final pipeline)

One section per work package, appended by the agent that did it (newest at the bottom). Spec: `docs/PIPELINE_FIX_PLAN_2026-11-03.md`
(read-only). Status: `docs/IMPLEMENTATION_CHECKLIST.md`.

Template:

## WP N -- items ... (YYYY-MM-DD)
- Commit: <hash>
- Files / functions changed:
- Checks run (command + key output):
- Checks the user must run on the Mac (exact commands):
- Deviations from the plan (each one asked and approved by the user? say when):
- Open issues for the next agent:

---

## WP0 -- baseline (2026-11-04)
- Commit: 2f4382a (baseline of the code and docs before any change of the fix plan).
- Not committed (left untracked on purpose): data-output/simulation/, data-output/main/ (non-ignored parts), Manuscript/,
  Present/, "Claude outputs/", docs/*.html and the hidden docs/.*backup* files.

## WP1 -- A1 (2026-11-04) -- CODE DONE, RUN CHECKS PENDING
- Commit: "WP1: A1 (code; run checks pending)" (hash: see `git log`; the final entry below fills in the run results and the hash).
- Files / functions changed:
  - `R/03c_build_climate_covariates.R`: step 4 is plot-aware (`need` = every plot-year the visit windows require, `miss` = `need` minus the
    checkpoint; the loop runs over years with missing plots and extracts ONLY those plots: `pep_yr`, `pts_ll` rebuilt per year, `id_pe = pep_yr$id_pe`).
    Before step 4: one-time copies `climate_by_year_raw.pre_A1.rds` and `climate_covariates.pre_A1.rds` (never overwritten); checks that PEP ids are
    unique and every visit's plot has a location. After step 4b: stop unless every needed plot-year is present and has no NA in gdd_3, cmi_gs,
    start_gs, end_gs (user decision: also stop on NA in a used plot-year). End: the FLAG on `n_years_used == 0` is replaced by `stop()` unless every
    visit has `n_years_used == win_end - win_start + 1` and `>= 1`; it runs BEFORE `saveRDS`, so a bad file is never written.
  - `R/03z_combine_covariates.R`: before `saveRDS(combined)`: prints NA counts and `stopifnot(!anyNA(...))` on gdd_3, cmi_gs, soil_pc1 (all visits,
    including plots 04 later drops; if this fails on such a plot, report, do not loosen).
  - `R/04_prepare_inla.R`: after `iv` is built: the same `stopifnot` on `iv`, and A1.5 on the first visit of every plot in `dt2` (done now, not in WP4;
    no new `plot_info` columns -- WP4 adds them). Prints the NA counts.
  - New `R/tests/wp1_A1_check.R` (read-only): new `climate_covariates.rds` vs `climate_covariates.pre_A1.rds`: NA counts, visits that differ
    (must be 0 outside the 5 plots), the 5 plots' values with same-year percentiles, the pair 7209002301 / 7209002302.
  - Not touched: 04x, the 04 header ("expected" counts left as they are, user decision), download_zip().
- Checks run: no R in the workspace and no access to NRCan, so no R code was run; bracket/quote balance of the three scripts equals the committed
  versions. Cell check (Python, 1961 CMI raster cmi60_06.tif, EPSG:4269, 0.01667 deg cells): none of the 5 plots shares a cell with a checkpoint plot
  (7109604001 shares its cell only with 7109604002, which is not in the checkpoint either), so the values cannot be taken from existing plots.
- Checks the user must run on the Mac, in this order, from the project root (zips already placed by hand, see below):
  `Rscript R/03c_build_climate_covariates.R` -> `Rscript R/tests/wp1_A1_check.R` -> `Rscript R/03z_combine_covariates.R` ->
  `Rscript R/04_prepare_inla.R` -> `Rscript R/04x_exclusion_flow.R`. Compare against the A1 Check items (checklist, WP1).
- Deviations from the plan (all asked and approved by the user on 2026-11-04):
  - A1.5 NA check implemented in WP1 on the first rows of `dt2`, not deferred to WP4.
  - Extra 03c guard (NA in gdd_3 / cmi_gs / start_gs / end_gs of used plot-years; missing plot-years) beyond the plan text.
  - `pre_A1` backup copies made by 03c itself; extra file `R/tests/wp1_A1_check.R`.
  - Plan A1 point 6 is WRONG in its explanation: the 5 plots are NOT removed by the 03z join (they have rows with NA climate and pass the inner
    join). The 11,281 / 36,242 / 8,766 counts come from an old cloud build. 04x is to be compared with the actual 04 output (11,286 plots /
    36,251 intervals / 8,770 changes expected unchanged by A1).
- Open issues for the next agent / later:
  - `download_zip()` (03c) is unchanged: if a download is interrupted, a partial `data.zip` stays on disk, a rerun sees the file, skips the
    download and then tries to unzip it. Proposed fix (not implemented, user decision): download to `data.zip.part` with `curl -fL -C - --retry`, check
    the size against Content-Length and `unzip(list = TRUE)`, rename only if all pass, delete the part file after repeated failure.
  - The zips for A1 (sg and cmi, 1961-1996) were fetched MANUALLY by the user with `~/Desktop/get_nrcan.sh` (curl with resume and retries,
    `unzip -t` check, moved to `data/raw/climate_raw/<sg|cmi>/<year>/data.zip`); 03c's `download_zip()` finds them and skips the download.
  - The extracted folder `data/raw/climate_raw/cmi/1961/1961/` is incomplete (no `cmi60_sum.tif`; the zip itself is complete). `download_zip()`
    unzips again on each call, so it is completed when 03c runs for that year.
  - A direct-read (`/vsizip//vsicurl/`) test was proposed and DROPPED by the user (zip entries are deflated, so random access saves little);
    nothing of it is in the project.
