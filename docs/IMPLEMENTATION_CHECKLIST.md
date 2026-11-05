# Implementation checklist (final pipeline)

**Spec:** `docs/PIPELINE_FIX_PLAN_2026-11-03.md` -- READ-ONLY. Never edit it. If it is unclear, wrong or needs a methodological
choice, stop and ask the user with an example.
**Details of each work package:** `docs/IMPLEMENTATION_LOG.md` (append one section per work package).
**This file:** tick the boxes (`[ ]` -> `[x]`) and fill in the commit hash. Do not change the wording of the items.

Rules for every agent:
- Read this file, the log, the plan sections of your items, "FINAL PIPELINE LIST 2026-11-03" in `docs/PIPELINE_NOTES.md`, and the
  code you touch, before changing anything.
- Do ONLY your work package. Do not start the next one.
- Diagnose before fixing. No methodological or data-coding decision on your own.
- Commit at the end ("WP N: items ..."). Tick an item only when its plan **Check** passed; checks that need INLA or the full data
  on the Mac are ticked by the agent after the user reports the output.
- A review agent (no edits) checks each work package's diff against the plan before the next package starts.

Order is strict: WP1 -> WP2 -> WP3 -> WP4 -> WP5 -> WP6 -> WP7. WP3-WP5 all edit `05a_fit_funs.R`.

---

## WP1 -- A1 (03c, 03z, 04 guards)   Commit: 7f5a1b6, e4e1f41, (ticks) see log   Review: [ ]
- [x] A1.1 03c step 4 plot-aware: missing plot-years extracted (re-download of the needed years)
- [x] A1.2/A1.4 03c `stop()` unless every visit has `n_years_used == win_end - win_start + 1`
- [x] A1.3 `stopifnot(!anyNA(...))` on gdd_3, cmi_gs, soil_pc1 at the end of 03z and in 04 (`iv`)
- [x] A1.5 04: same NA check on the first-visit covariates in `plot_info` (needs C1 point 1 columns; add the check in WP4 if the columns do not exist yet -- note it in the log)
- [x] User ran on the Mac: 03c -> 03z -> 04 -> 04x
- [x] Check: 0 NA climate; all other visits identical to the previous `climate_covariates.rds`; 04 plots/intervals/events unchanged; 5 plots plausible (7209002301 ~ 7209002302); 04x matches 04

## WP2 -- A2-A11, B7 (04, covariate_spec.R, config.R, 04b, tests/04c)   Commit: fa1d432   Review: [ ]
- [x] A2 harvest: light + heavy only, never "any" (04, covariate_spec, recount_events, tests/04c)  -- *69 of 69 harvest terms on scratch (17 + 21 + 31 -> "3level"); no `harv_any`/`tsd_harv` left in the WP2 files or 05a*
- [x] A3 fire/pest rule unchanged (MIN_EVENTS from config)  -- *fire 14, pest 28 on scratch; transition_set counts and terms identical to the 7f5a1b6 baseline on the same table*
- [x] A4 tsd_harv_light / tsd_harv_heavy  -- *04c: at most one harvest TSD non-zero per interval; 0 where the level flag is 0*
- [x] A5 TSD in years in 04; model columns via TSD_FORM (log1p main, linear = S24)  -- *max TSD 174 y fire / 157 pest / 165 harvest; max log1p column 5.16 (plan ~5.2); linear = old decades x 10 exactly (checked vs baseline)*
- [x] A6 per-level main + GDD x + TSD terms; 0-event level = main term only  -- *0-event level: only other_temperate->yellow_birch (heavy); GDD x / TSD heavy terms in 68 transitions, light and main in 69*
- [x] A7 other_temperate->yellow_birch has harv_term "3level"  -- *3level on scratch; terms: z_gdd_3 + z_cmi_gs + z_soil_pc1 + harv_light + harv_heavy + gdd_x_harv_light + tsd_harv_light*
- [x] A8 no TSD cap (nothing to change; confirm)  -- *confirmed, no TSD cap in 03a/03z/04*
- [x] A9 D_harv_heavy, D_fire, D_heavy (= earlier) in visit_intervals; C(i) edge-case count printed  -- *1,535 X->regen intervals with D_heavy on scratch (harvest 1,221, fire 315) vs plan 1,537/1,223/315: the 2 missing are the X->regen events of dropped plots 7209002301/7209002302; dates all in (year_start, year_end]; C(i): 9 intervals, 0 of them X->regen*
- [x] A10 counts on end-visit flags (unchanged; confirm)  -- *confirmed unchanged; recount_events passes on the final rows*
- [x] A11 n_undist in transition_set; no z_gdd_3 where n_undist == 0 and a disturbance term exists  -- *only regen->sugar_maple (n_undist 0) loses z_gdd_3 on scratch*
- [x] B7 switches in config.R (WEIBULL_VARIANT, COEF_PRIOR_SD, TSD_FORM, MIN_EVENTS, XREGEN_EXIT_RULE); LOW_EVENTS_MAX <- MIN_EVENTS + 4L  -- *switches in config.R; `LOW_EVENTS_MAX <- MIN_EVENTS + 4L` is done in WP7 (09), not here (user decision 2026-11-04)*
- [x] User ran on the Mac: 04 and tests/04c  -- *2026-11-04: 04 and 04c ALL CHECKS PASSED on the full table (11,286 plots / 36,251 intervals / 8,770 changes)*
- [x] Check (implementation step 2): only harvest terms (A2, A4, A6) and regen->sugar_maple's GDD main term (A11) change; transitions, events, fire and pest terms identical; 69 harvest terms; 1,537 X->regen intervals with D_heavy  -- *04c baseline check vs data/tables/main: 0 differences in kept set, n_events, n_fire, n_pest, n_harv0/1/2, fire_term, pest_term (90 kept); harvest term changed in 52 (21 any + 31 none -> 3level); 69 harvest terms, fire 14, pest 28; 1,537 (1,223 / 315); only regen->sugar_maple lacks z_gdd_3*

## WP3 -- B1-B6 (05a fit + theta, 05, 07, 08)   Commit: 4237012   Review: [x]
- [x] B1 variant from config in every INLA call; SIM_DIR / model_fitting.R source removed from 05 and 07; run_meta records it  -- *code done (05a, 05, 07; sensitivity/11 is WP7); check `run_meta$inla_args$control.family[[1]]$variant == 0` needs a 05 run on the Mac (see the WP5 item "B1/B2/B5 run-time checks"); `wp3_B4_check.R inla` asserts the variant given to INLA on one real fit*  *Mac 2026-11-04 (step 12, run_meta of the first smoke-test attempt in `fits_subsample`): user reported it passed (variant 0); the output was not pasted into the chat, no values recorded here*
- [x] B2 control.fixed: N(0, COEF_PRIOR_SD^2), flat intercept  -- *code done; check (run_meta shows the values) pending the Mac run (WP5 item "B1/B2/B5 run-time checks"); `wp3_B4_check.R inla` asserts prec == 1/COEF_PRIOR_SD^2 (0.04), mean == 0, mean.intercept == 0, prec.intercept == 0 in the `control.fixed` given to INLA on one real fit (assertions added in the review fixes; the first WP3 commit did not have them)*  *Mac 2026-11-04 (step 12): user reported it passed (control.fixed as expected); output not pasted, no values recorded here*
- [x] B3 intercept kept on INLA's (decade) scale in tidy_fit / mi_raw / screen / pooling; years-scale intercept computed in 08 (reporting only)  -- *no `- log(10)` left on the INLA estimates (only `v0_to_p1p2` for P1/P2); 08 `intercept_years_<tag>.csv` matches an independent calculation on synthetic fits*
- [x] B4 variant-0 -> P1/P2 conversion in extract_theta and draw_theta_m; extras$theta_m stores the variant-0 draws  -- *`v0_to_p1p2()`; theta_m columns: alpha, intercept_dec, raw beta*
- [x] B4 unit test: S(t) from variant-0 parameters == S(t) from converted lambda_1, alpha (1e-10), 3 rows x t = 5, 20, 50 y  -- *mock fits passed (max difference 4e-16, extract_theta, draw_theta_m, fallback); the one-real-fit part needs INLA: `Rscript R/tests/wp3_B4_check.R inla` on the Mac* *Mac run 2026-11-04: `wp3_B4_check.R inla` ALL CHECKS PASSED, variant 0 confirmed.*
- [ ] B5 screen unchanged (both rules; refs from the new m1-30)  -- *logic untouched; 07 now prints the full reference table; check is the printout in the first 07 run on the Mac*
- [x] B6 fit_map_v0 + H_at moved to 05a; per-fit cap columns in mi_raw; 08 writes cap_check_<tag>.csv  -- *columns maxH_inla, n_above_inla, maxH_map, n_above_map, z_map_max, z_map_alpha, z_map_intercept (z_map_alpha and z_map_intercept added by user decision; also medians in cap_check); `add_cap_columns` and 08 tested on mock / synthetic data; columns filled in a real mi_raw is a WP5 check*
- [x] B6 check: the moved functions reproduce `cap_test_new_rows_m01_nocut/results.csv` on the test's rows  -- *180 fits (arms A, B): n_above_* identical; max relative difference maxH_map 6.3e-7, alpha_map 1.3e-8, maxH_inla 2.8e-13; tolerance 1e-6 met, cause of the 6.3e-7 unexplained (platform difference untested), accepted by the user*
- [x] B7 recount_events() called in run_stage1() after build_rows()  -- *plan grep finds only config.R and the fit function; wiring tested with stubs; skipped (printed) when SUBSAMPLE_N is set, so the first real run is the full Stage 1*

## WP4 -- C1, C5 (P1; 04 plot_info first-visit values)   Commit: c2b2b46   Review: [ ]
- [ ] C1.1 04 writes first-visit covariates to plot_info (+ A1.5 NA check if not done in WP1)  -- *code done (`fv_*` columns, NA/coding stopifnot on plot_info; block tested on a synthetic dt2: structural zero TSD, NA stops); the check on the real 04 output is the Mac run below*
- [ ] C1.2-C1.5 P1 base row = first visit; flags switch on at onset = W_j - first-visit TSD; yearly pieces after the first onset, TSD at mid-year; exact piecewise draw  -- *code done (`p1_pieces`, `p1_next_exit`, `p1_candidates_pw`; extra breakpoint at the exact onset of a later type = user decision, see log); synthetic-table checks 1-2 passed; the sample-of-real-plots check is the Mac run below*
- [ ] C1.6 + C5 weights: same piecewise hazards to the first visit, observed-interval covariates after it, divided by S(W_j - e)  -- *code done (`p1_weights_pw`, log space); equals an independent numerical integral (max difference 1.4e-16, 5 plots); the real-data shift is the WP5 item below*
- [x] C1 unit test: 100,000 piecewise draws match the exact piecewise survival (KS p > 0.01)  -- *`R/tests/wp4_C1_check.R` check 4, rewritten in the WP4 review fixes (exact side = closed-form piecewise Weibull survival written independently of the production code, pieces of every case rebuilt from the raw first-visit values; a separate seed per case): 4 cases of 1e5 draws, KS p / exit-type chi-square p -- case 1 plot 1, entry 0 (seed 20261105): 0.228 / 0.352; case 2 plot 1, entry 5.3 (20261106): 0.605 / 0.945; case 3 plot 2, entry 0 (20261107): 0.080 / 0.077; case 4 plot 3, entry 11.7 (20261108): 0.267 / 0.455. Bound 0.01 unchanged; all four pass. (Superseded first version: shared seed, spline CDF, p = 0.553 / 0.976 / 0.487 / 0.013; that 0.013 was chance, 1e6-draw reruns gave p = 0.17 / 0.63 / 0.06.)*
- [x] C1 check: unflagged plots give the same candidates as the current P1 for the same seed  -- *user decision 2026-11-04 (A): reference = the WP3 `p1_candidates` fed the same base-row lambda (the climate/soil of the first visit, C1.2, so HEAD's interval-1 covariates cannot give identical candidates); 26 synthetic plots without an onset before W_j (unflagged, and flagged with onset = W_j): same seed, max abs difference 4.3e-14*
- [x] C5 check: toy example reproduces the exact posterior mean  -- *`R/tests/wp4_C5_toy.R` (user's toy, 2026-11-04): exact 30.52 | current 36.59 | corrected 30.56 (0.9 SE); production functions on the same toy: current 36.52 (99 SE away), corrected 30.43 (-1.3 SE). The plan's 27.3 / 37.3 / 27.4 came from a version that cannot be found*
- [ ] P1 fallbacks printed per imputation and by first state  -- *code done (`p1_fallback_by_state`; printed in 05 Stage 2 and in the 05/07 loops, `p1_fallback_by_state[_redraw]_<tag>.csv`); first real output with the WP5 full imputation*
- [ ] User ran on the Mac: `Rscript R/04_prepare_inla.R`, then `Rscript R/tests/wp4_C1_check.R real` (first Mac run failed at plot 7609800202 with W_j = 0: test fixed, see the log "WP4 -- follow-up")  -- *(added in WP4) expect the 04 NA counts 0 and 11,286 plots; the real section: onsets = W_j - TSD, X1 = scaled first-visit values (300 plots), pieces contiguous 0 -> end of the weighting bracket (300 flagged plots)*

## WP5 -- C2, C3 (P2; run_stage1) + first test runs   Commit: 44bde15   Review: [ ]
- [ ] C2 exit rule (D_heavy) in Stage 1, Stage 2 and every imputation; exit_source recorded in extras$p2  -- *code done (`p2_rule_T`, `make_p2_stage1`, `make_p2_moment`, `make_p2_draw`; `stage1_rows`); cloud (`R/tests/wp5_C23_check.R`, synthetic): rule date = t_start + D_heavy - year_start in all three, for ok / failed / non-kept transitions, X->regen only, clamp > 1e-9 stops, XREGEN_EXIT_RULE = FALSE off; plan Check (1,537 rule exits, every rule exit = D_heavy, floored count) needs the real data: `R/tests/wp5_C23_real_check.R` (Mac list in the log)*
- [ ] C3 make_p2_draw (closed-form truncated draw; uniform where no usable fit); Stage 2 keeps the moment-matched point  -- *code done; Stage 2 identical to the pre-WP5 function where there is no D_heavy (460 synthetic events, difference 0); real-data checks (spread over 5 imputations > 0, draw mean vs moment point) are in `wp5_C23_real_check.R` (Mac)*
- [x] C3 unit test: 100,000 draws match the truncated Weibull CDF (KS p > 0.01)  -- *`R/tests/wp5_C23_check.R` check 1 (cloud, R 4.3.3): production `make_p2_draw` vs the CDF written independently (conditional CDF with expm1), own seed per case (20261110 + case): KS p = 0.181 / 0.198 / 0.451 / 0.465 / 0.918 for cases 1-5 (a, b, lam, alpha in the log; case 5 has H(a) = 800, S(a) underflows). Bound 0.01 unchanged. Negative control: the same draws against a CDF with alpha + 2% give p = 3e-4*
- [ ] User ran on the Mac: `Rscript R/tests/wp5_recount_stage1.R` (added in WP5; full data, no fitting: `prepare_fit_inputs` without subsampling, then `stage1_rows(inp)` = exactly the build_rows call of `run_stage1`, then `recount_events`; run BEFORE the smoke test)  -- *tested in the cloud on a synthetic table (400 plots); pass = "recount_events(): OK" and ALL CHECKS PASSED*
- [x] User ran on the Mac: `Rscript R/tests/wp5_p1_guard_check.R` (added in WP5 "P1 guards"; needs `data-output/final/fits_subsample/stage2.rds` from the first smoke-test attempt)  -- *pass = p1_entries on the Stage 2 theta runs (< 60 s, diagnosed 0.5 s) and p1_entries on the drawn theta of m = 1 stops within 60 s with the guard-A message ("GUARD OK: ..."); cloud: tested with a fake sampler (one drawn beta of 54, alpha 0.17)*  *Mac 2026-11-04 (step 11): user reported it passed (stage-2 theta and the drawn-theta guard stop); the output was not pasted into the chat when this was ticked, so no numbers are recorded here*
- [ ] User ran on the Mac: smoke test 05 with SUBSAMPLE_N = 500 (code runs)  -- *`mkdir -p data-output/final && PEP_SUBSAMPLE_N=500 PEP_M_MI=6 Rscript R/05_fit_inla.R 2>&1 | tee data-output/final/smoke_console.log` (environment overrides added in WP5, defaults unchanged); writes to `fits_subsample`* -- *DROPPED (user decision 2026-11-04): with 500 plots the smoke test cannot finish an imputation by design (sparse fits give drawn betas up to 54, infinite hazards; first attempt hung in the MI loop, see the log "WP5 -- P1 guards"). Replaced by the guard test above; run_meta of that attempt is checked in the B1/B2 item*
- [ ] B1/B2/B5 run-time checks from the first 05/07 runs
  - run_meta (smoke test writes to `fits_subsample`; the full run to `fits`): `Rscript -e 'm <- readRDS("data-output/final/fits_subsample/run_meta_final.rds"); print(m$inla_args$control.family[[1]]$variant); str(m$inla_args$control.fixed); print(unlist(m[c("WEIBULL_VARIANT","COEF_PRIOR_SD","TSD_FORM","MIN_EVENTS","XREGEN_EXIT_RULE")]))'`  -- expect variant 0; control.fixed mean 0, prec 0.04, mean.intercept 0, prec.intercept 0; switches 0 / 5 / log1p / 15 / TRUE (B1, B2)
  - 07 reference table (after the full 05 run): `Rscript R/07_redraw.R 2>&1 | tee data-output/final/redraw_console.log`, then `grep -n -A95 "Screen reference table" data-output/final/redraw_console.log`  -- expect all 90 kept transitions (k, icpt_ref, icpt_sd_ref, n_ref, label) printed before the first "Batch m=" line (B5)
- [ ] User ran on the Mac: `Rscript R/tests/wp4_C5_shift.R 1 1000` after the first full imputation (C5 real-data shift, from WP4)  -- pass = corrected entries earlier on average (mean old - corrected > 0); size vs the plan's 2.9-3.8 y for information only
- [ ] User ran on the Mac: `Rscript R/tests/wp5_C23_real_check.R 5` after the first batch (m1-6) of `PEP_M_MI=6 Rscript R/05_fit_inla.R` (added in WP5)  -- *pass = ALL CHECKS PASSED: 1,537 rule exits per imputation, each = D_heavy, 0 non-rule events with zero spread over m1-5; prints exit_source per label / transition, uniform_nonfinite per transition, draw mean vs the Stage 2 moment point (information only)*
- [ ] User ran on the Mac: one full imputation; checks: 1,537 rule exits; B6 columns filled; P1 fallbacks by first state vs main run's 124-196 (watch the 681 unharvested sugar-maple-first plots)

## WP6 -- full run (user)   Run tag: ______
- [ ] RUN_TAG set to the new tag (e.g. "final"); `main` untouched
- [ ] 05 -> 07_redraw -> 08 done
- [ ] Watch: cap_check (caveat list comes from here), failure counts of regen->sugar_maple and other_temperate->yellow_birch, P1 fallbacks

## WP7 -- D1-D8 (09, 12, 08, sensitivity/10, 11, tests/08f)   Commit: ______   Review: [ ]
- [ ] D1 HR = exp(beta); no alpha in the delta method; alpha_cov columns out of the HR table
- [ ] D2 tau grid 1/10/30/60 y for every transition; range columns; _tsd_1_vs_10y, _tsd_10_vs_30y, _tsd_coef
- [ ] D3 alpha_report_<tag>.csv
- [ ] D4 one-sided bounds incl. the main-term-only contrast; regen->sugar_maple at z = 0 only
- [ ] D5 12: no harv_any; block D on _tsd_coef; block B at each tau; robust without alpha_cov; A11 transitions out of A-GDD, B, C
- [ ] D6 tests/08f: fit_map_v0 / H_at, new mask, no NA-row dropping, B6 columns
- [ ] D7 sensitivity/10 (n_dist per level, MIN_EVENTS_S1 from config, no alpha_cov_low_share) and 11 (variant 0, control.fixed, intercept scale, exp(beta), A11)
- [ ] D8 09 for transitions without z_gdd_3
- [ ] 09: `LOW_EVENTS_MAX <- MIN_EVENTS + 4L` (moved from B7)
- [ ] User ran on the Mac: 09 -> 12 -> 10 -> 11; every 12 block has k > 0

## Outside the code
- [ ] C4 Methods text (exit rule incl. undated heavy disturbances and "fire"; exits drawn; limitations list in the plan)
- [ ] C5 simulation rerun with corrected weights (copy of the simulation code inside this project)
- [ ] Refit script for S21, S22, S24 (before those sensitivities)
