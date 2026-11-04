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

## WP1 -- A1 (03c, 03z, 04 guards)   Commit: ______   Review: [ ]
- [ ] A1.1 03c step 4 plot-aware: missing plot-years extracted (re-download of the needed years)
- [ ] A1.2/A1.4 03c `stop()` unless every visit has `n_years_used == win_end - win_start + 1`
- [ ] A1.3 `stopifnot(!anyNA(...))` on gdd_3, cmi_gs, soil_pc1 at the end of 03z and in 04 (`iv`)
- [ ] A1.5 04: same NA check on the first-visit covariates in `plot_info` (needs C1 point 1 columns; add the check in WP4 if the columns do not exist yet -- note it in the log)
- [ ] User ran on the Mac: 03c -> 03z -> 04 -> 04x
- [ ] Check: 0 NA climate; all other visits identical to the previous `climate_covariates.rds`; 04 plots/intervals/events unchanged; 5 plots plausible (7209002301 ~ 7209002302); 04x matches 04

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
- [ ] User ran on the Mac: 04 and tests/04c
- [ ] Check (implementation step 2): only harvest terms (A2, A4, A6) and regen->sugar_maple's GDD main term (A11) change; transitions, events, fire and pest terms identical; 69 harvest terms; 1,537 X->regen intervals with D_heavy

## WP3 -- B1-B6 (05a fit + theta, 05, 07, 08)   Commit: 4237012   Review: [ ]
- [ ] B1 variant from config in every INLA call; SIM_DIR / model_fitting.R source removed from 05 and 07; run_meta records it  -- *code done (05a, 05, 07; sensitivity/11 is WP7); check `run_meta$inla_args$control.family$variant == 0` needs a 05 run on the Mac (wp3_B4_check.R inla checks the same on one real fit)*
- [ ] B2 control.fixed: N(0, COEF_PRIOR_SD^2), flat intercept  -- *code done; check (run_meta shows the values) pending the Mac run; wp3_B4_check.R inla asserts them on one real fit*
- [x] B3 intercept kept on INLA's (decade) scale in tidy_fit / mi_raw / screen / pooling; years-scale intercept computed in 08 (reporting only)  -- *no `- log(10)` left on the INLA estimates (only `v0_to_p1p2` for P1/P2); 08 `intercept_years_<tag>.csv` matches an independent calculation on synthetic fits*
- [x] B4 variant-0 -> P1/P2 conversion in extract_theta and draw_theta_m; extras$theta_m stores the variant-0 draws  -- *`v0_to_p1p2()`; theta_m columns: alpha, intercept_dec, raw beta*
- [ ] B4 unit test: S(t) from variant-0 parameters == S(t) from converted lambda_1, alpha (1e-10), 3 rows x t = 5, 20, 50 y  -- *mock fits passed (max difference 4e-16, extract_theta, draw_theta_m, fallback); the one-real-fit part needs INLA: `Rscript R/tests/wp3_B4_check.R inla` on the Mac*
- [ ] B5 screen unchanged (both rules; refs from the new m1-30)  -- *logic untouched; 07 now prints the full reference table; check is the printout in the first 07 run on the Mac*
- [x] B6 fit_map_v0 + H_at moved to 05a; per-fit cap columns in mi_raw; 08 writes cap_check_<tag>.csv  -- *columns maxH_inla, n_above_inla, maxH_map, n_above_map, z_map_max, z_map_alpha (z_map_alpha added by user decision); `add_cap_columns` and 08 tested on mock / synthetic data; columns filled in a real mi_raw is a WP5 check*
- [x] B6 check: the moved functions reproduce `cap_test_new_rows_m01_nocut/results.csv` on the test's rows  -- *180 fits (arms A, B): n_above_* identical; max relative difference maxH_map 6.3e-7, alpha_map 1.3e-8, maxH_inla 2.8e-13*
- [x] B7 recount_events() called in run_stage1() after build_rows()  -- *plan grep finds only config.R and the fit function; wiring tested with stubs; skipped (printed) when SUBSAMPLE_N is set, so the first real run is the full Stage 1*

## WP4 -- C1, C5 (P1; 04 plot_info first-visit values)   Commit: ______   Review: [ ]
- [ ] C1.1 04 writes first-visit covariates to plot_info (+ A1.5 NA check if not done in WP1)
- [ ] C1.2-C1.5 P1 base row = first visit; flags switch on at onset = W_j - first-visit TSD; yearly pieces after the first onset, TSD at mid-year; exact piecewise draw
- [ ] C1.6 + C5 weights: same piecewise hazards to the first visit, observed-interval covariates after it, divided by S(W_j - e)
- [ ] C1 unit test: 100,000 piecewise draws match the exact piecewise survival (KS p > 0.01)
- [ ] C1 check: unflagged plots give the same candidates as the current P1 for the same seed
- [ ] C5 check: toy example reproduces the exact posterior mean
- [ ] P1 fallbacks printed per imputation and by first state

## WP5 -- C2, C3 (P2; run_stage1) + first test runs   Commit: ______   Review: [ ]
- [ ] C2 exit rule (D_heavy) in Stage 1, Stage 2 and every imputation; exit_source recorded in extras$p2
- [ ] C3 make_p2_draw (closed-form truncated draw; uniform where no usable fit); Stage 2 keeps the moment-matched point
- [ ] C3 unit test: 100,000 draws match the truncated Weibull CDF (KS p > 0.01)
- [ ] User ran on the Mac: smoke test 05 with SUBSAMPLE_N = 500 (code runs)
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
