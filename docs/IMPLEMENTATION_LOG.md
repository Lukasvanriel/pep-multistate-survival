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

## WP2 -- A2-A11, B7 (2026-11-04) -- CODE DONE, FINAL CHECK PENDING WP1 RUN
- Commit: "WP2: A2-A11, B7 (final check pending WP1 run)" (hash: see `git log`; a follow-up commit records it in the checklist).
- Files / functions changed:
  - `R/config.R`: new switches `WEIBULL_VARIANT <- 0L`, `COEF_PRIOR_SD <- 5`, `TSD_FORM <- "log1p"`, `MIN_EVENTS <- 15L`, `XREGEN_EXIT_RULE <- TRUE`
    (WEIBULL_VARIANT, COEF_PRIOR_SD, XREGEN_EXIT_RULE are defined only; used from WP3 / WP5). `LOW_EVENTS_MAX` NOT added (WP7, user decision).
  - `R/04_prepare_inla.R`: local `MIN_EVENTS` removed (config). TSD kept in YEARS: `tsd_fire_y`, `tsd_pest_y`, `tsd_harv_y` (0 where unflagged); old
    decade columns `tsd_fire/pest/harv` no longer exist in `visit_intervals`. A9: `D_harv_heavy`, `D_fire`, `D_heavy` (= earlier) in `visit_intervals`;
    prints the X->regen count and the C(i) edge-case count. A2: `harv_term` = "3level" iff kept and `n_harv1 + n_harv2 >= MIN_EVENTS`, else "none".
    A11: new column `n_undist` in `transition_set` (events with none of the transition's disturbance terms, computed AFTER the new term set; end-visit
    flags; pest/fire/harvest events in a transition without that term count as undisturbed). 04 now reads `data/build/disturbance_events.rds` (C(i) count only).
  - `R/covariate_spec.R`: `add_model_covariates(d, scaling, tsd_form = TSD_FORM)` builds `tsd_fire`, `tsd_pest`, `tsd_harv_light`, `tsd_harv_heavy` from the
    `_y` columns (log1p, or linear = decades); `harv_any`, `gdd_x_harv_any`, `tsd_harv` removed. `all_terms()` = 15 terms. `term_mask()` (needs `n_harv1`,
    `n_harv2`, `n_undist` in the transition set): A6 per-level GDD x / TSD only with >= 1 event of the level; A11 `z_gdd_3` FALSE where `n_undist == 0` and
    the transition has a disturbance term. The config must be sourced before these functions are used (TSD_FORM default).
  - `R/04b_build_rows.R` `recount_events()`: default `min_events = MIN_EVENTS`; harvest condition `harv_term != "3level" | n_harv1 + n_harv2 >= min_events`.
  - `R/tests/04c_check_rows.R`: `harv_any` count dropped; new CHECK 6 (stop() on failure): B7 switches, A2, A3, A6, A7, A11, A4, A5 (log1p and linear), A9.
  - Not touched (as agreed): 05a (calls `all_terms()`/`term_mask()`; `prepare_fit_inputs()` verified to run on the new tables), 09, 12, tests 14/15/16/22/23.
- Checks run (R 4.x in the cloud workspace; scratch copy of `combined_covariates.rds` WITHOUT plots 7109604001, 7209002301, 7209002302, 7209605702,
  7501600102; baseline = 04 at 7f5a1b6 on the same scratch table; nothing of the scratch is in this folder):
  - Scratch totals: 11,281 plots / 36,242 intervals / 8,766 changes (plan: 11,286 / 36,251 / 8,770 = 5 plots, 9 intervals, 4 changes more). 90 kept transitions.
  - New vs baseline: `kept`, `n_events`, `n_fire/pest/harv0/1/2` (end and start), `fire_term`, `pest_term` identical; visit_intervals differ only in the TSD columns
    (old decades x 10 == new years) and the new `D_*` columns. Harvest terms 69 (17 "3level" kept, 21 "any" and 31 "none" -> "3level", 21 "none").
  - 04c: ALL CHECKS PASSED (checks 1-5 unchanged, CHECK 6 new).
  - A9 on scratch: X->regen intervals with D_heavy 1,535 (harvest 1,221, fire 315, both 1); plan 1,537 / 1,223 / 315. The 2 missing are the X->regen heavy-harvest events
    of the dropped plots 7209002301 / 7209002302 (verified on the full table). C(i) edge case (user definition 2026-11-04: dated heavy record with year in
    (year_start, year_end], or undated heavy record first seen at the end visit): 9 intervals, 0 of them X->regen.
- Checks the user must run on the Mac, after the WP1 run (on the new combined table): `Rscript R/04_prepare_inla.R` then `Rscript R/tests/04c_check_rows.R`.
  Expect: 69 harvest terms, fire 14, pest 28, transitions/events/fire/pest identical to the pre-WP2 04 output, 1,537 X->regen intervals with D_heavy
  (1,223 / 315), `regen->sugar_maple` the only transition without `z_gdd_3`, `other_temperate->yellow_birch` without `gdd_x_harv_heavy` / `tsd_harv_heavy`,
  04c ALL CHECKS PASSED. Then tick the two "User ran on the Mac" / final-check items.
- Deviations from the plan (all asked and approved by the user on 2026-11-04):
  - `LOW_EVENTS_MAX <- MIN_EVENTS + 4L` not done here (09 is WP7); noted next to B7 in the checklist.
  - WP2 checks added as a permanent section of `tests/04c`; grep check limited to the WP2 files plus 05a and R/sensitivity (clean; the only hit is 04c's own check naming the removed terms).
  - C(i) definition as stated above.
- Open issues for the next agent / later:
  - `09_hazard_ratios.R` (`harv_any`, `tsd = "tsd_harv"`) and `12_cross_transition_synthesis.R` still use the OLD term names until WP7. Older diagnostics
    (`tests/14`, `15`, `16`, `22`, `23`, `04t`) belong to the old build and were left untouched; `term_mask()` fails on transition sets without `n_undist`.
  - GDD main term resting on 1 undisturbed event (plan rule: only 0 is a boundary, kept; user decision): `white_birch->regen`, `jack_pine->regen`,
    `sugar_maple->regen` (scratch table).
  - `05a_fit_funs.R` still uses its earlier Weibull variant / priors until WP3; `WEIBULL_VARIANT`, `COEF_PRIOR_SD`, `XREGEN_EXIT_RULE` are not read anywhere yet.
  - A5: P1/P2 get the new TSD columns through `add_model_covariates()` as planned; nothing else reads the decade columns of `iv`.
