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
  - Not touched: 04x, download_zip(). (04's printed "expected" line was updated in the review fixes below.)
- Checks run: no R in the workspace and no access to NRCan, so no R code was run; bracket/quote balance of the three scripts equals the committed
  versions. Cell check (Python, 1961 CMI raster cmi60_06.tif, EPSG:4269, 0.01667 deg cells): none of the 5 plots shares a cell with a checkpoint plot
  (7109604001 shares its cell only with 7109604002, which is not in the checkpoint either), so the values cannot be taken from existing plots.
- Checks the user must run on the Mac, in this order, from the project root (zips already placed by hand, see below):
  `Rscript R/03c_build_climate_covariates.R` -> `Rscript R/tests/wp1_A1_check.R` -> `Rscript R/03z_combine_covariates.R` ->
  `Rscript R/04_prepare_inla.R` -> `Rscript R/04x_exclusion_flow.R`. Compare against the A1 Check items (checklist, WP1).
- Deviations from the plan (all asked and approved by the user on 2026-11-04):
  - A1.5 NA check implemented in WP1 on the first rows of `dt2`, not deferred to WP4.
  - Extra 03c guard (NA in gdd_3 / cmi_gs / start_gs / end_gs of used plot-years; missing plot-years) beyond the plan text (widened to all averaged columns in the review fixes).
  - Stricter window guard than plan A1.2/A1.4: besides `n_years_used == win_end - win_start + 1` the 03c guard also requires `n_years_used >= 1`
    (an empty window can never pass; the user confirmed none exist, earliest visit 1970).
  - `pre_A1` backup copies made by 03c itself; extra file `R/tests/wp1_A1_check.R`.
  - Plan A1 point 6 is WRONG in its explanation: the 5 plots are NOT removed by the 03z join (they have rows with NA climate and pass the inner
    join). The 11,281 / 36,242 / 8,766 counts come from an old cloud build. 04x is to be compared with the actual 04 output (11,286 plots /
    36,251 intervals / 8,770 changes expected unchanged by A1).
- Open issues for the next agent / later:
  - A1.5 was done in WP1 on the first rows of `dt2` (checks gdd_3, cmi_gs, soil_pc1 at each plot's first visit). WP4 must re-add the check on the
    first-visit columns written to `plot_info` (C1.1), because that is what P1 reads.
  - `download_zip()` (03c) is unchanged: if a download is interrupted, a partial `data.zip` stays on disk, a rerun sees the file, skips the
    download and then tries to unzip it. Proposed fix (not implemented, user decision): download to `data.zip.part` with `curl -fL -C - --retry`, check
    the size against Content-Length and `unzip(list = TRUE)`, rename only if all pass, delete the part file after repeated failure.
  - The zips for A1 (sg and cmi, 1961-1996) were fetched MANUALLY by the user with `~/Desktop/get_nrcan.sh` (curl with resume and retries,
    `unzip -t` check, moved to `data/raw/climate_raw/<sg|cmi>/<year>/data.zip`); 03c's `download_zip()` finds them and skips the download.
  - The extracted folder `data/raw/climate_raw/cmi/1961/1961/` is incomplete (no `cmi60_sum.tif`; the zip itself is complete). `download_zip()`
    unzips again on each call, so it is completed when 03c runs for that year.
  - A direct-read (`/vsizip//vsicurl/`) test was proposed and DROPPED by the user (zip entries are deflated, so random access saves little);
    nothing of it is in the project.

### WP1 review fixes (2026-11-04, commit "WP1: A1 review fixes"; review of 7f5a1b6, all approved by the user)
- `R/config.R`: `RUN_TAG <- "final"` (was "main"). 04 now writes `data/tables/final/`; `data/tables/main` and `data-output/main` stay untouched.
  Consequences: (1) the WP6 checklist item "RUN_TAG set" is done early in config but left unticked until WP6; (2) 05, 07, 08, 09, 12 and
  `sensitivity/10`, `11` read `data-output/final/`, which does not exist until WP6; (3) `R/tests/22_cap_test_new_rows.R`, `23_p1_onset_bench.R`,
  `08f_failed_vs_valid.R` and `08g_bs_regen_drop_plot.R` read the MAIN run's stored draws (`extras_<RUN_TAG>_m01.rds`, `mi_raw_<RUN_TAG>`, ...) through
  `RUN_TAG` and now stop with file-not-found. WP3 (B6 check: reproduce `cap_test_new_rows_m01_nocut`) and WP5 (comparison against it) must read the main run
  explicitly (e.g. a separate source tag for the extras, tables of main for the old rows). Nothing was changed in those tests now.
  The simulation scripts name `main` explicitly and are unaffected. `data/build/` has no tag (shared; differs from the main run only in the 5 plots).
- `R/03z_combine_covariates.R`: one-time copy `combined_covariates.pre_A1.rds` before `saveRDS` (never overwritten).
- `R/03c_build_climate_covariates.R`: the NA guard covers all averaged columns (16 sg, cmi_1-12, cmi_sum, cmi_gs; per-column NA counts printed) and `out` is
  checked for NA in all averaged columns before `saveRDS`. Read-only check before the edit: the checkpoint (717,835 plot-years) has no NA in its 29 raw columns.
- `R/04_prepare_inla.R`: the printed "expected" line now says: pre-A1 table, unchanged by A1: 11,286 plots, 36,251 intervals, 8,770 state changes.
- `R/tests/wp1_A1_check.R` rewritten. Hard stops: 0 NA in the averaged columns, 0 differing visits outside the 5 plots (exact, all columns), every window full;
  part B (only after 03z, needs `combined_covariates.pre_A1.rds`): same rows/columns, identical outside the 5 plots, 0 NA in gdd_3 / cmi_gs / soil_pc1,
  only climate columns change in the 5 plots. FLAGS only (printed, never stop): each of the 5 plots vs its 3 nearest plots of the pre-A1 checkpoint (PEP
  coordinates, same window) and the pair 7209002301 / 7209002302, bound = 99th percentile of the same statistic over the pool plots, same window
  (neighbours: |value - median of the 3 nearest|; pair: |difference| of existing nearest-neighbour pairs <= 0.6 km). The check holds a COPY of `cmi_gs_matched()`
  (the checkpoint has no cmi_gs); it validates the copy on 200 random old visits and skips the flags if the copy does not reproduce them.
  Reference values (Python, windows 1961-70 / 1987-96): neighbour 99th percentile gdd_3 144 / 148, cmi_gs 5.4 / 4.5; pair gdd_3 86 / 86, cmi_gs 5.2 / 3.1.
- Not run (no R in the workspace); bracket/quote balance only. The user runs: 03c -> check -> 03z -> check (adds part B) -> 04 -> 04x.

## WP2 -- A2-A11, B7 (2026-11-04) -- CODE DONE, FINAL CHECK PENDING WP1 RUN
- Commit: "WP2: A2-A11, B7 (final check pending WP1 run)" (hash: fa1d432).
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

## WP2 -- review fixes (2026-11-04)
- Commit: "WP2: review fixes" (b2c314d) and "WP2: log correction" (this entry).
- Changes: WP7 checklist gets the item `09: LOW_EVENTS_MAX <- MIN_EVENTS + 4L (moved from B7)`; `tests/04c` CHECK 6 compares the new `transition_set`
  (TABLE_DIR) with `data/tables/main/transition_set.rds` (hard check: kept set, n_events, n_fire, n_pest, n_harv0/1/2, fire_term, pest_term; prints the
  changed harvest terms and the transitions without `z_gdd_3`; skipped with a printed note if the main table is missing or is the table under test);
  the `covariate_spec.R` header comment on "main effects both present" is corrected for A11 transitions.
- Correction: the review chat said TABLE_DIR equals `data/tables/main`, so 04 would overwrite the baseline and a copy was needed. That was wrong. Since
  `e4e1f41` `RUN_TAG` is "final", so TABLE_DIR is `data/tables/final` and `data/tables/main` stays untouched (config.R line 9). No copy is needed and
  no code change. The "table under test IS the baseline" skip in `tests/04c` only fires if RUN_TAG is set back to "main".
- Not run: the new 04c block has only been parse-checked (no data.table in the cloud workspace). On a scratch table without the 5 WP1 plots it is
  expected to differ by those plots; on the full table it must pass.

## WP3 -- B1-B6 + B7 item "recount_events() in run_stage1()" (2026-11-04) -- CODE DONE, MAC CHECKS PENDING
- Commit: "WP3: B1-B6" (hash: see `git log`; recorded in the next commit).
- Files / functions changed:
  - `R/05a_fit_funs.R`: `fit_inla_models_pc` (B1/B2): `variant = WEIBULL_VARIANT`, `control.fixed = list(mean = 0, prec = 1/COEF_PRIOR_SD^2, mean.intercept = 0, prec.intercept = 0)`.
    `tidy_fit` (B3): intercept `est` = INLA mean (decade scale), no `- log(10)`. New `v0_to_p1p2()` (B4; stops unless variant 0): intercept = b0/alpha - log 10, beta = beta/alpha.
    `extract_theta` (posterior means) and `draw_theta_m` (per draw, own alpha) return the converted values plus the raw `intercept_dec`, `beta_v0`; `ok` also needs alpha > 0; a draw with
    alpha <= 0 / NA falls back to the seed. `extras$theta_m` now has `alpha, intercept_dec, <raw beta columns>` (was `intercept_years` + converted beta). New: `CAP`, `fit_map_v0`, `H_at` (moved from
    tests/22, SIGMA_B -> COEF_PRIOR_SD), `CAP_COLS`, `cap_diagnostics`, `add_cap_columns` (called in `run_one_imputation`; columns `maxH_inla, n_above_inla, maxH_map, n_above_map, z_map_max, z_map_alpha`,
    one value per fit repeated over its rows; NA for a failed fit). `prepare_fit_inputs` returns `subsampled`. `run_stage1`: `recount_events(b$rows, inp$tr)` after `build_rows()`; skipped with a printed note when subsampling.
  - `R/05_fit_inla.R`, `R/07_redraw.R`: `SIM_DIR` and `source(model_fitting.R)` removed; `run_meta` records `WEIBULL_VARIANT, COEF_PRIOR_SD, TSD_FORM, MIN_EVENTS, XREGEN_EXIT_RULE` (and keeps `inla_args` with control.fixed / control.family). 07 prints the full screen-reference table (all kept transitions, 90) before any candidate is fitted, then the n_ref < 10 list. Screen logic unchanged.
  - `R/08_final_pool.R`: new `intercept_years_<tag>.csv` (per used fit: b0 - alpha log 10 from the `mi_raw` means, variance var(b0) + log10^2 var(alpha) - 2 log10 cov(b0, alpha) from the stored covariance; Rubin-pooled) and `cap_check_<tag>.csv`
    (per transition: n_used, n_map, share maxH_map > 22.85, share maxH_inla > 22.85, medians of maxH_map, z_map_max, z_map_alpha). Pooling and screen stay on the decade scale.
  - `R/config.R`: comments only. New tests: `R/tests/wp3_B4_check.R`, `R/tests/wp3_B6_reference_check.R` (scratch-tree harness, header explains the layout).
- Checks run (cloud; R 4.x, data.table; no INLA):
  - B4 mock: `Rscript R/tests/wp3_B4_check.R`: 4 transitions with different term sets, 3 random rows x t = 5/20/50: max |S0 - S1| 2.8e-16 (extract_theta), 4.4e-16 (draw_theta_m, with a fallback transition) -> ALL CHECKS PASSED.
  - B6 reproduction: rows of tests/22 (`CAP_CUT=FALSE`, imputation 1 of the main run) built with the BASELINE 2f4382a files (covariate_spec, 04b, config RUN_TAG "main", 05a, blob hashes verified) in a scratch tree outside the project, then the NEW `fit_map_v0` / `H_at`
    on those rows against `cap_test_new_rows_m01_nocut/results.csv` and the stored `fits/*.rds` (INLA estimates; INLA not rerun): 180 fits (arms A and B x 90): `n_above_map`, `n_above_inla`, `n_above_ev_inla` identical, no NA mismatch,
    max relative difference `maxH_map` 6.3e-7 (arm B; arm A 4.1e-8), `alpha_map` 1.3e-8, `maxH_inla` 2.8e-13; the full MAP coefficient / SE tables of the stored fits agree to max 5.0e-6 / 9.3e-7 relative (coefficients near 0 amplify; not part of the 1e-6 criterion).
    The tolerance was NOT loosened; the largest `maxH_map` difference is within a factor 1.6 of it.
  - `add_cap_columns` on mock `tidy_fit` output from 3 stored INLA fits + 1 failed fit (black_spruce->regen, regen->fir, sugar_maple->regen): all six columns equal the values computed from the stored cap-test tabs; failed fit -> NA.
  - 08 on synthetic `mi_raw` / extras (2 transitions, 31 imputations, one screened fit): `intercept_years` pooled estimate and SE equal an independent calculation; `cap_check` share 0.4 as constructed.
  - `codetools` free-variable check of 05a + 04b + covariate_spec: nothing from model_fitting.R is needed. `run_stage1` wiring (build_rows -> recount_events -> fit_rows; skip when subsampled) tested with stubs.
  - B7 grep (`MIN_EVENTS *<-|variant = |prec = ` in R/*.R): only config.R and `fit_inla_models_pc`.
- Checks the user must run on the Mac (project root):
  - `Rscript R/tests/wp3_B4_check.R inla` -> expect "ALL CHECKS PASSED", `INLA call: variant 0`, `control.fixed ... prec 0.04`. Closes the B4 unit test (one real fit, real `inla.posterior.sample` names) and shows B1/B2 on a real fit.
  - B1/B2/B5 run-time checks come with the WP5 runs (`run_meta_final.rds`; 07 reference printout).
- Deviations from the plan (asked and approved by the user 2026-11-04):
  - `z_map_alpha` added to the B6 columns (alpha collapses in cap fits); `z_map_max` covers the term coefficients only. `cap_check` also has the `maxH_inla` share and medians. `run_meta` records all config switches. 07 prints the full reference table.
  - B6 reference build from the baseline code in a scratch tree (WP2 changed the main-table columns, so tests/22 cannot run at HEAD); sensitivity/11 left for WP7.
  - `recount_events()` is skipped when `SUBSAMPLE_N` is set.
- Open issues for the next agent / later:
  - NOT ticked until the Mac run: B1, B2, B4 unit test (one real fit), B5.
  - The WP5 smoke test (SUBSAMPLE_N = 500) skips `recount_events()`, so its first real run is Stage 1 of the full imputation.
  - `08` existing behaviour (not changed): `dcast` of the invalid-fits table stops when there are no invalid fits at all.
  - Older diagnostics (tests 06a, 08f, 17, 23, test_fit_mock) read `extras$theta_m$intercept_years` of the main run; they were left untouched and still match the old files.
  - `v0_to_p1p2` stops for `WEIBULL_VARIANT != 0` (a variant-1 sensitivity would need its own conversion).
  - MAP cost per fit (about 1 s) is added to every fit in each imputation; `n_above_ev` is computed by `H_at` but not stored in `mi_raw`.
