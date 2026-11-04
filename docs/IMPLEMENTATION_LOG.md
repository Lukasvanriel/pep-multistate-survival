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
- Commit: "WP3: B1–B6" (hash 4237012). Committed with the repository's existing author identity (the VM has no git identity); trailers name Claude.
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

## WP3 -- review fixes (2026-11-04)
- Commit: "WP3: review fixes" (hash: see `git log`).
- Correction: the first WP3 commit (4237012) did NOT contain the `control.fixed` assertions in `R/tests/wp3_B4_check.R` (the write-back of the edited test raced with the copy and the older version was committed), so the B2 checklist note and the "Checks the user must run" text overstated what was asserted. Now `wp3_B4_check.R inla` checks, on the real fit, the settings INLA was given (`.args`): `control.family$variant` 0, `control.fixed$prec` == 1/COEF_PRIOR_SD^2 (0.04), `mean` 0, `mean.intercept` 0, `prec.intercept` 0. It prints the table of the settings first.
- Changes: (1) the test's failure policy: a failed check prints what it has (settings table; estimate, SE, bound and true value; or the S0/S1 difference) and stops with "bounds were NOT changed, report to the user"; bounds unchanged (alpha 0.15, intercept 0.3, identity 1e-10). Dry-run of the INLA block with a stub model (no INLA in the cloud): passes with a model near the truth; a shifted alpha and a wrong `prec.intercept` each stop with the table. The real-INLA run itself is still the user's.
  (2) 08: the `stopifnot` on the B6 cap columns now sits right after `mi_raw` is read, before any computation or output (tested: a `mi_raw` without a column stops with the message and nothing else printed). (3) New B6 column `z_map_intercept` (|INLA - MAP| / MAP SE of the intercept, like `z_map_alpha`) in `cap_diagnostics` / `mi_raw` and `median_z_map_intercept` in `cap_check`; `add_cap_columns` rechecked against the stored cap-test tabs (black_spruce->regen 15.17, regen->fir 0.020, sugar_maple->regen 0.246) and 08 on synthetic data.
  (4) Checklist: the WP5 item "B1/B2/B5 run-time checks from the first 05/07 runs" with the commands; B1/B2 notes point to it.
- B6 reproduction: the largest relative difference (maxH_map 6.3e-7, tolerance 1e-6) is accepted by the user; its cause is UNEXPLAINED (a platform difference Linux vs Mac was not tested).
- Mac command: `Rscript R/tests/wp3_B4_check.R inla` (project root). If it stops, report the printed table; do not change bounds.

## WP3 -- fix variant check (2026-11-04)
- Commit: "WP3: fix variant check" (hash: see `git log`).
- Diagnosis (by the user on the Mac, `Rscript R/tests/wp3_B4_check.R inla` failed): INLA was given the right settings. `.args$control.family` is a LIST OF ONE per-likelihood list, so the variant is at `.args$control.family[[1]]$variant` (= 0L); `control.fixed` mean 0 / prec 0.04 / mean.intercept 0 / prec.intercept 0 confirmed in `.args` and in `all.hyper$fixed`. The test read `.args$control.family$variant`, which is NULL; `c()` dropped it and the settings table shifted by one row, so the failure was the test's, not the fit's.
- The variant assertion in 4237012 (`mods[[1]]$.args$control.family$variant == 0L` in `stopifnot`) was VACUOUS: NULL == 0L is logical(0), which passes `stopifnot`. The review-fixes version (7538746; 0b58fb6 was the pre-amend commit, now orphaned) built the table with the same wrong path. Variant 0 is confirmed by the user's Mac fit via the corrected path (`control.family[[1]]$variant` = 0L); it was not confirmed by any earlier test run.
- Changes: `R/tests/wp3_B4_check.R` reads `control.family[[1]]$variant`; the settings are collected in a named list and every value must have length 1, be numeric and not NA before the table is built, otherwise the check stops with a message naming the missing setting (no silent drop). Dry-run with a stub model shaped like INLA's `.args` (list of one): passes; a missing variant stops naming it; variant 1 stops with the table. Checklist lines 45 and 70 (B1 note, WP5 run_meta command) use `control.family[[1]]$variant`.
- Done (user, Mac, 2026-11-04): `Rscript R/tests/wp3_B4_check.R inla` rerun after this fix: ALL CHECKS PASSED, variant 0 confirmed. B4 unit test ticked.

## WP3 -- close (2026-11-04)
- Commit: "WP3: close" (docs only, no code).
- Plan error (plan is read-only, not edited): the B1 check string `run_meta$inla_args$control.family$variant == 0` is NULL in practice; `.args$control.family` is a list of one, so the correct path is `control.family[[1]]$variant`. The checklist (B1 note, WP5 run_meta command) and `wp3_B4_check.R` already use it.
- B4 unit test ticked (Mac run 2026-11-04, ALL CHECKS PASSED, variant 0 confirmed). B1 and B2 stay unticked until the WP5 run_meta check; B5 until the first 07 run.
- WP3 review ticked.

## WP4 -- C1, C5 (2026-11-04) -- CODE DONE, MAC CHECKS PENDING
- Commit: "WP4: C1, C5" (hash c2b2b46). Author Lukas Van Riel as in the earlier commits.
- User decisions before the work (2026-11-04): (1) yearly grid from the first onset, last piece cut at W_j with the midpoint of the shortened piece; (2) EXTRA breakpoint at the exact
  onset of a later type (piece splits in two, each with TSD at its own midpoint) -- DEVIATION from "yearly pieces" in the plan (C1.5); (3) time since entry stays the clock when lambda
  changes; (4) the "same candidates" check = new code vs the WP3 `p1_candidates` fed the same base-row lambda; (5) the user's own toy for C5; (6) the real-data shift as a test script run in WP5
  (checklist item added under WP5); (7) fv_* columns and defaults as proposed. (The message numbered the clock answer "2. A" and the toy under "3"; read as above.)
- Files / functions changed:
  - `R/04_prepare_inla.R`: `plot_info` gets `fv_gdd_3, fv_cmi_gs, fv_soil_pc1, fv_fire, fv_harv, fv_pest, fv_tsd_fire_y, fv_tsd_pest_y, fv_tsd_harv_y` (raw, first row of `dt2`; TSD per TSD_SOURCE, 0 where the flag is 0).
    New `stopifnot` on these columns (A1.5 on what P1 reads: no NA; flags 0/1 and 0/1/2; TSD >= 0 and 0 where unflagged); prints NA counts and flagged / unflagged plot counts. The WP1 check on `dt2` stays.
  - `R/covariate_spec.R`: new `tsd_fun(tsd_form)` (the TSD transform of `add_model_covariates`, now shared with P1); no change of results.
  - `R/05a_fit_funs.R`: `prepare_fit_inputs` also returns `X1` (first-visit model row), `X0` (same without disturbance), `on` (onset = W_j - first-visit TSD per type, harvest level), `iv_first`, `iv_n`; stops with a clear message if
    `plot_info.rds` has no fv_* columns (rerun 04) or the intervals of a plot are not contiguous from W_j. `inp$X` (end-visit values of interval 1) kept as LEGACY, no longer read by P1.
    New: `p1_pieces` (pieces + lambda per piece; observed intervals after W_j up to the end of the weighting bracket), `p1_next_exit` (exact piece-by-piece inversion; earliest competing exit), `p1_candidates_pw`,
    `p1_H_between`, `log1mexp`, `p1_weights_pw` (C5: [S(lo - e) - S(hi - e)] / S(W_j - e) in log space = exp(-H(W,lo)) (1 - exp(-H(lo,hi))); no change: exp(-H(W,hi))), `p1_fallback_by_state`.
    `p1_entries` now uses them (Stage 2 mean mode and the imputation loop). The WP3 `p1_candidates` / `p1_weights` are kept unchanged as the reference of the identity test. `run_one_imputation` returns `p1_state`.
  - `R/05_fit_inla.R`, `R/07_redraw.R`: P1 fallbacks by first state printed (Stage 2 table; one line per imputation) and appended to `p1_fallback_by_state_<tag>.csv` / `p1_fallback_by_state_redraw_<tag>.csv` (first_state, n_plots, n_fallback, share, m).
  - New tests: `R/tests/wp4_C1_check.R` (synthetic tables through the real `prepare_fit_inputs`; `real` argument for the real tables), `R/tests/wp4_C5_toy.R`, `R/tests/wp4_C5_shift.R` (WP5).
- Checks run (cloud, R 4.3.3, data.table; no INLA, no real data):
  - `Rscript R/tests/wp4_C1_check.R` -> ALL CHECKS PASSED: (1) X1 = scaled first-visit values, X0 without disturbance, onset = W_j - TSD, 40 plots; (2) pieces: plot with onsets 25.4 and 27.5 gets breakpoints 25.4, 26.4, ..., 39.4, 40 plus 27.5, last piece [39.4, 40] with midpoint 39.7; onset -5 -> grid from 0 and TSD 5.5 at the first midpoint; lambda per piece equals an independent build;
    (3) identity: 26 plots, max difference 4.3e-14; (4) KS: see checklist (p = 0.013 is the lowest of four; 1e6-draw rerun p = 0.17 / 0.63 / 0.06); (5) weights equal an independent numerical integral (max 1.4e-16); (6) `p1_entries` reproducible with the seed; no-ok-transition case; (7) `c5_shift` runs.
  - `Rscript R/tests/wp4_C5_toy.R` -> exact 30.52 | current 36.59 | corrected 30.56 | production functions: current 36.52, corrected 30.43 (-1.3 SE). ALL CHECKS PASSED.
  - `Rscript R/tests/wp3_B4_check.R` (mock) still passes (max |S0 - S1| 4.4e-16).
  - The new 04 block was run on a synthetic `dt2` (flags set / unset, nonzero TSD where unflagged -> 0; an NA first-visit gdd_3 stops with "A1.5").
  - Two diagnosed bugs of my own on the way (fixed before the run above): rounding of W_j to 9 decimals could round it up past W_j and drop the last breakpoint (non-integer W_j); a stray name on the breakpoint vector.
- Checks the user must run on the Mac (project root, in this order):
  1. `Rscript R/04_prepare_inla.R`  (writes the fv_* columns; expect NA counts 0 in the new "plot_info first-visit columns" line, 11,286 plots; compare "none flagged" with the plan's 6,330)
  2. `Rscript R/tests/wp4_C1_check.R real`  (synthetic part again, then: real onsets = W_j - TSD, X1 vs scaled first-visit values for 300 plots, pieces contiguous from 0 to the end of the weighting bracket for 300 flagged plots; prints the piece counts and the number of onsets <= 0)
  3. optional, `Rscript R/tests/04c_check_rows.R` (04 changes only plot_info; 04c should still pass)
  4. in WP5, after the first full imputation: `Rscript R/tests/wp4_C5_shift.R 1 1000` (C5 real-data shift)
- Deviations from the plan (each asked and approved by the user 2026-11-04):
  - C1.5: extra breakpoint at the exact onset of every type that starts after the first onset (the plan has only yearly pieces after the first onset).
  - C1 identity check redefined (the plan's wording cannot hold because C1.2 moves the climate/soil to the first visit).
  - C5 toy: the plan's numbers (27.3 / 37.3 / 27.4) came from a version that cannot be found; the user's toy gives 30.52 / 36.59 / 30.56. Same pattern (current later by about 6 y, corrected = exact), different numbers.
  - C5 real-data shift moved to a script run in WP5 (the main-run theta has the old terms and cannot be evaluated with the new X). It compares old and corrected weights on the SAME piecewise hazards and candidates, so it isolates the C5 effect.
- Open issues / things to watch:
  - NOT verified: the real 04 output, the real `prepare_fit_inputs`, the P1 fallbacks by first state, P1 run time. Cloud timing on the synthetic plots with a random theta: 0.011 s per unflagged and 0.021 s per flagged plot (M = 200, one core, W_j median 64): not comparable with the plan's 1.19 min per imputation on the Mac; the real figure comes with the WP5 full imputation.
  - Pieces after W_j use each observed interval's end-visit row, so the hazard jumps at W_j from the first-visit climate to the interval-1 end-visit climate (as decided in the plan).
  - A flagged plot with TSD 0 (onset = W_j) has its flag off for the whole of P1's history (C1.3); its P1 row is X0.
  - `tests/23_p1_onset_bench.R` (benchmark prototype) and older diagnostics that read `inp$X` were not touched.
  - A plot starting in regen with no change observed keeps entry 0 even when no transition is ok (same as the WP3 code).
- Git note: the WP4 commits left `.git/HEAD.lock` / `.git/index.lock` and a few `tmp_obj_*` files that the sandbox may not delete; the stale locks were moved to `.git/stale-*-wp4*` (safe to delete), the `tmp_obj_*` files are unreferenced leftovers (safe to delete).
