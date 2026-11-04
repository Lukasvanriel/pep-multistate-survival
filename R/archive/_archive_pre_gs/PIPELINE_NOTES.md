# PEP multistate pipeline: decisions, order, status
Last updated 2026-10-20. Simulation folder `p1p2_v2_test` is FINISHED and READ-ONLY: never modify it.
All new work lives in this project (`pep-multistate-survival`).

## Goal
Adapt the validated simulation pipeline (P1 entry imputation, P2 exit imputation, per-transition INLA
Weibull fits, Rubin pooling; orchestrated by `fit_p1p2_combined.R`) to the real PEP data.
No truth exists: no bias/coverage evaluation, no bias correction, no oracle columns (`true_*`), no hidden hops.

## Order (strict)
A. Data preparation (no imputation) - script `R/04_prepare_inla.R`
   1. W_j per plot  2. visit-interval table  3. transition set + covariate terms
B. Imputation + fitting (simulation pipeline, adapted)
   4. Builder: one sequential pass per plot over intervals; P1 gives first entry, P2 gives exits; emits split rows.
      Stage 1 (naive: entry W/2, midpoint exits) -> Stage 2 (P1 weighted-mean entry + moment-matched P2) -> MI loop.
   5. INLA fit per transition, Rubin pooling (`05_fit_inla.R`), diagnostics (`06_inspect_plot.R`).

## Decisions made
- Regen state = "state 0" at time 0. Hidden hops ignored. Number of states irrelevant (transitions decomposed). Interval gaps irrelevant.
- Clock reset: time since entry. Absolute time is not a concern.
- Data: `data/processed/combined_covariates.rds` (48,012 visits, 11,393 plots), keys id_pe + no_mes.
- W_j (stand age at FIRST visit, fixed per plot, used only by P1): from `cl_age`.
  Source: `pee_ori_sond.cl_age`, gaps filled from `station_pe.cl_age` (both in data/raw/pep/PEP.gpkg).
  Class -> number: VIN/VIR=100, JIN/JIR=50, else leading 120/90/70/50/30/10.
  Rules: 1 first visit; 2 later visit back-calculated (minus elapsed years; if <= 0: 0 when the first state is regen, otherwise the later class value itself, decided 2026-10-14); 3 min(time since fire, harvest) at first
  visit from raw origine/perturb; 4 first state regen AND later exits regen -> 0; else EXCLUDE.
  RESULT of 04_prepare_inla.R (dataset visits): rule1 9758, rule2 1584, rule3 0, rule4 2, excluded 49 of 11,393 (44 regen, 3 black_spruce, 1 jack_pine, 1 yellow_birch);
  11,344 plots included. (Raw-visit estimate was 47: two plots only had cl_age at raw visits that are not in the combined dataset.)
  Layer comparison: pee more accurate where checkable (80% vs 66% exact class vs age from origin year), station far more complete at visit 1.
- Transition filter: a transition needs >= 15 actual events (observed state changes, not N at risk), else dropped as a whole.
  Covariate terms per kept transition: continuous: none extra; fire needs >= 15 fire events; pest >= 15 pest events;
  harvest (0 none / 1 light / 2 heavy) needs >= 15 in EACH level; if a level fails collapse to any-harvest and need >= 15; else drop harvest.
  Counted at the event interval.
  RESULT (04_prepare_inla.R, 11,344 plots, 36,505 intervals, 8,879 state changes): 146 observed pairs of 156, 91 kept (>= 15) holding 8,525 events;
  354 events (4.0%, 312 plots) fall on dropped transitions; all 13 states reachable from regen via kept pairs.
  Terms with end-visit flags (91 kept): fire 14, pest 19, harvest 3level 17 / any 32 / none 42.
  ERRATUM: an earlier count of 6,382 events / 77 kept was wrong (pandas string-dtype bug in the run grouping); 8,881 events matches data-output/data/transition_counts.csv.
- 2026-10-19: the covariate terms per kept transition now ALSO include GDD x disturbance interactions and TSD terms; see section "Status 2026-10-19" (earlier results in data/processed/fit/ are main-effects only).
- Covariates change over time (flags AND climate), so rows are SPLIT per visit interval (no collapsing into runs).
  `rle` is not needed; a running `entry_used` carries the clock origin.
- Row rules (time since entry, in decades = years/10), per sojourn with entry E:
  truncation of first segment of first sojourn = max(W - E1, 0)/10; first segment of a later sojourn = 0 (starts at imputed exit);
  any later segment = previous visit - E (end of previous segment). Tstop = segment end (visit) or imputed exit T for the event segment.
  Every segment has status-0 rows for all destinations; the event segment has one status-1 row (realized destination) + status-0 for the others.
  Stub: the segment from imputed entry to the first visit in the new state is its own segment.
- Dropped transitions: hazard 0 in P1; their events still get an imputed exit (needed for competing censoring and the next entry), by midpoint.
- Plots must never be dropped because P1 finds no candidate path: needs a fallback (proposal: naive Stage-1 entry W/2, counted and reported). TO IMPLEMENT.
- P1 uses the covariates of the plot's first interval. P2 uses the event interval's covariates.
- Final recount of events per transition/level right before INLA (last check).

## Verified (do not redo)
- Collapsed vs split rows give identical INLA results (static covariates): `R/test_split_vs_collapsed.R`
  (events 1281 = 1281, one per id; log-lik diff 5e-12; INLA diffs <= 3e-7).
  The test's "true value" labels used the wrong Weibull parametrization (pipeline fits variant 1: S = exp(-(lambda t)^alpha));
  estimates equal true/alpha as expected. Rules also checked in Python on 2000 random two-sojourn plots (max diff 4e-16).
- Step D thresholds re-read from the simulation header: numbers reproduce, but event-count thresholds cannot be re-derived from the saved rds.
  Red flag: imputed continuous-term coverage degrades at large N (gdd ~0.47-0.50, cmi 0.33-0.77 at N=2000), alpha coverage 0 at N=2000.

## Findings about the data that matter
- Disturbance flags `fire/harv/pest` in combined_covariates are CARRIED FORWARD within an epoch (03_build_covariates.R step 7);
  an origine-level event starts a new epoch and resets. ~99% of flagged starts stay flagged at the next visit.
  New onsets between consecutive visits: fire 287, pest 1000, harv 2848 of 36,619 pairs. tsf/tsh/sfp are time since disturbance.
  Harvest coding: 2 = origine-level (heavy), 1 = perturb-level (light).
- Start vs end flag on the 8,880 state-change intervals: they differ in 199 (fire), 482 (pest), 1682 (harv). Example black_spruce->regen (n=809): heavy harvest start 44 vs end 423; fire 70 vs 141.
- Sim details: sojourn covariates were static per individual (cov1); D drawn per sojourn; multi-covariate Step D/E generation code is not in the folder.

- W_j data quality: 96 plots have W_j = 0 (94 via rule 2 flooring at 0), 35 of them with a NON-regen first state; non-regen first state with W_j <= 5: 75, <= 10: 275.
  P1 (forward simulation from regen to W_j) cannot reach those first states -> they need the no-candidate fallback.

## Open / to decide
- Resolved 2026-10-19: continuous set is final (gdd_3, cmi_sum, soil_pc1); no domain random effect for now.
- Interval covariates: flags at the END visit (default, carried in table) vs START visit (also carried as *_start).
- Which continuous covariates enter the model; scaling.
- Fallback for plots with no P1 candidate path (see above).
- Adapt P1 (covariate row = first interval) and P2 (event interval covariates) for per-transition formulas (beta matrix with zeros for dropped terms).
- No R/INLA in the cloud environment: tests are run by the user locally.

## Status 2026-10-14 (B1, B2, B8, B7 implemented and run in the cloud with R on the real intervals)
- `R/covariate_spec.R` (B1): COV_CONT is a PLACEHOLDER (gdd_3, cmi_sum, soil_pc1), z-scored; fire/pest/harvest dummies from END-visit flags;
  `term_mask()` / `rhs_formulas()` give per-kept-transition formulas (zeros in the beta matrix for dropped terms).
- `R/04b_build_rows.R` (B2): `build_rows()` split-row builder, `p2_midpoint()` hook, `recount_events()` (B7).
- `R/04c_check_rows.R` (B8): Stage-1 build (entry W_j/2, midpoint exits) -> 544,608 rows, 45,384 segments; all checks pass:
  status-1 rows 8,879 = state changes; kept events 8,525; time at risk per plot = last - first visit; segments contiguous;
  row log-lik == per-sojourn log-lik (diff 1e-11); recount matches transition_set.
- BUG FOUND AND FIXED while building: the stub segment was first only emitted when another interval followed the event, which dropped the
  time from the imputed exit to the last visit when the event interval was a plot's last one. Stub is now emitted after EVERY event.
- Caveat: check 4 compares rows with a collapsed likelihood built from the same segments; the independent INLA equivalence test is `R/test_split_vs_collapsed.R`.

## Status 2026-10-15: fitting code written (05a_fit_funs.R, 05_fit_inla.R); NOT yet run with INLA
- Decisions: M = simulation values (SIM_M 200 paths, M_MI 30, 6 workers). P1 with no candidate -> entry W_j/2, marked and counted, NO retry (assess afterwards).
  Stage-1 fit failed / NA for a kept transition -> hazard 0 in P1 and midpoint in P2 for that transition only, counted (assess afterwards).
- Stages: 1 naive (W_j/2, midpoint) -> theta_stage1; 2 P1 weighted mean + moment-matched P2 with theta_stage1 -> seed models;
  MI loop (draw theta_m from seed posterior, P1 draw, chained P2, fit); Rubin pooling. Checkpoints in data/processed/fit/ (delete to redo).
- Copied/adapted from the simulation (folder untouched): sim_entry_candidates_2b, compute_weights_2b, moment_match_weibull_naive, draw_theta_m, pool_rubin.
  fit_inla_models is SOURCED read-only from the simulation folder.
- P1 bracket (real): first change interval (last visit in first state, first visit in new state); no change -> right-censored at last visit. No hidden-hop flag.
- build_rows: exit imputed exactly at the last visit would give Tstop == truncation; Tstop is floored at truncation + EPSILON and counted (n_floored).
- Mock test (R/test_fit_mock.R, fake INLA): flow works on all 11,344 plots; Stage 2 about 136 s with 2 cores (P1 dominates).
  Mock theta is meaningless, so the mock P1 fallback counts (30-60 plots, mostly jack_pine first states) say nothing about the real rate.
- To check on the first real run: posterior-sample row names printed by 05_fit_inla.R, number of failed Stage-1 fits, P1 fallback plots (p1_fallback_plots_stage2.csv).

## Status 2026-10-15: first real run (INLA, all 11,344 plots) and diagnostics rerun
First run (outputs mi_raw.csv, mi_diag.csv, pooled.csv in data/processed/fit/; kept untouched):
- Stage 1 about 1 min, Stage 2 about 2 min, MI loop 30 imputations about 25 min. Posterior-sample name matching works (n_theta_fallback 0 in all 30).
- Stage 2: 0/91 fits failed, alpha range 0.82-2.69, events per fit min 15 / median 48. P1: 124 of 11,344 plots without a usable candidate (124-196 per MI run),
  concentrated in young non-regen plots (W_j quality, upstream); entry W_j/2 for those, marked. P2: 8,879 exits, 354 by midpoint (dropped transitions).
- Problems seen in the MI results (cause NOT understood yet):
  * 5 failed fits (black_spruce->regen m=5,11,15; other_temperate->white_birch m=24,28); the error text was not stored (tidy_fit discarded it).
  * 9 degenerate fits (alpha posterior SD exactly 0) for black_spruce->regen: pooled alpha 1.846 vs 1.888 without them, fmi 0.82.
  * red_maple->yellow_birch m=1 outlier (alpha 3.34, SD 0.028) doubles the pooled SE; 7 more fits with unusually small alpha SD (mostly alpha ~1.05).
  * Rubin's B is not robust to outliers. Do not trim because values are extreme; drop only fits shown to be computational failures (non-convergence, SD 0),
    keep spread that is imputation-driven. Cannot be decided without the error/diagnostic output -> rerun below.
  * alpha > 2 for 6 transitions; Stage 1 -> Stage 2 alpha shift. No truth, so bias cannot be judged.
- Seed (stage 2) models: all 91 ok, mode.status 0, no warnings, so INLA flags alone will not separate the bad fits. Several seed fits have alpha mode 1.000, SD ~0.05
  (prior or sparse data? not checked).
- COV_CONT is still a placeholder (supervisor discussion pending). Bioclimatic domain random effect (earlier spec) is NOT in the model.
- Pooling still lacks Barnard-Rubin df and log-scale pooling of alpha (offered, not implemented).
- Reproducibility: the first run's theta draws used INLA's default (clock) sampler seed, so imputations cannot be reproduced.

Diagnostics rerun (code in R/05a_fit_funs.R, R/05_fit_inla.R; tested with the mock only):
- RUN_TAG "diag": outputs mi_raw_diag.csv, mi_diag_diag.csv, mi_errors_diag.log, pooled_diag.csv, so the first-run files stay.
- Sampler seed now taken from R's RNG per call (seed > 0, num.threads "1:1") -> reproducible after set.seed(SEED*1000+m) (not yet verified with real INLA).
- Per fit in mi_raw: error_msg, mode_status, mlik, log_warn (warning-like INLA log lines).
- extras_diag_mXX.rds per imputation: theta_m (drawn parameters), p1 (entry, n_cand, fallback per plot), p2 (imputed exit per event),
  fits[[label]] = summary.fixed / summary.hyperpar (quantiles, mode, kld), internal hyperpar, marginals (fixed, alpha), mlik, mode theta/status, log posterior mode,
  optimiser trace, nfunc, warnings, cpu, n_rows, last 40 log lines, and the joint posterior mean + covariance of (Intercept, betas, alpha) from N_COV_SAMPLES (1000)
  posterior samples. Intercept on INLA's decade scale (subtract log 10 for years); covariance unaffected. Failed fit -> list(error = text).
- models_diag_m01.rds: full INLA objects for imputation 1 (SAVE_MODELS_FOR), to measure size (first-run stage2.rds: 217 MB for 91 fits, ~9 MB per fit in memory).
  Unneeded parts: per-row summaries, model.matrix, .args$data (about 2.3 + 2.6 MB per fit); misc$configs (3.9 MB) only for later sampling.
- run_meta_diag.rds (once): INLA version, sessionInfo, INLA control settings (priors) of a seed fit, source of fit_inla_models, formulas, transition table.
- Unknown until the run: time for the 1000 posterior samples per fit (lower N_COV_SAMPLES, or 0, if too slow) and the size of extras per imputation.
- Next: read error_msg / log_warn / opt_trace of the failed and degenerate fits; compare their imputed inputs (p1, p2, theta_m) with the normal fits; then decide drop vs keep.

## Status 2026-10-15 (later): result of the "diag" run and further changes
Run result (data/processed/fit/mi_*_diag.csv, pooled_diag.csv, models_diag_m01.rds 211 MB):
- The device copies of 05_fit_inla.R / 05a_fit_funs.R were the earlier diagnostics version, so extras_diag_mXX.rds and run_meta_diag.rds were NOT written
  (cause unknown; possibly an older editor buffer saved over the committed file). Extras version re-committed; check the files after reloading the editor.
- Saved diagnostics: error_msg is only INLA's generic "inla-program exited with an error" (uninformative); mode_status is 0 for every fit that ran (INLA's flags do not catch
  the bad fits); log_warn matched a harmless header line (fixed now).
- Problem fits: black_spruce->regen 5 failed (m=1,4,6,17,20) + 4 with alpha SD < 1e-3 (m=7,14,15,23); other_temperate->white_birch 3 failed (m=8,9,11).
  Different imputations than in the first run -> imputation-dependent, not a fixed-seed effect. red_maple->yellow_birch outlier recurs (alpha 3.37 at m=7; first run 3.34 at m=1); its SD is normal so no SD screen catches it.
  other_temperate->white_birch (21 events): alpha near 1.04 in most imputations, 1.4 in at least 5 (large SD) -> weak data, not a failure (not confirmed).
- Stability between the two runs (pooled alpha): SE ratio new/old 0.49-1.56; median change in pooled estimate 0.02 SE, max 0.43 SE (trembling_aspen->sugar_maple 1.163 -> 1.080, SE halved).
  P1 fallbacks 133-191 per imputation (stage 2: 124).
Changes made after this run (tested with the mock and on the real mi_raw_diag.csv; not yet run with INLA):
- fit data of every failed or near-degenerate fit (alpha SD < SD_MIN = 1e-3) saved per imputation as baddata_<tag>_mXX.rds (data, rhs, label) -> rerun with verbose = TRUE.
- log_warn pattern fixed (no longer matches INLA's "Relative/Absolute error" header lines).
- Pooling: pool_all() keeps the regular Rubin columns unchanged (checked identical to pooled_diag.csv) and adds r, fmi, Barnard-Rubin df (complete-data df = events - parameters), t-based 95% interval,
  and for alpha pooling on the log scale (log_est, log_se, est_from_log, ci_*_log). Written as pooled_<tag>.csv (all fits that ran) and pooled_<tag>_excl.csv (flagged fits removed:
  INLA error or alpha SD < SD_MIN, removed for all parameters of that fit); flagged_fits_<tag>.csv lists them.

## Status 2026-10-19: failed-fit redraw done; interactions and TSD added to the model (implemented, NOT yet run with INLA)
Failed-fit handling (07_redraw_failed.R, run on the main-effects model, tag diag):
- black_spruce->regen 8 of 30 fits collapsed (m=2,8,14,15,19,23,25,29), other_temperate->white_birch 3 INLA errors (m=20,24,29). Per-transition redraw from candidate imputations m31-42:
  replacements black_spruce->regen m31,33,34,35,37,39,40,41; other_temperate->white_birch m31,32,33. Screen: ratio = alpha_SD/(0.78*alpha/sqrt(events)); a fit fails on INLA error, NA estimate or ratio < 0.1.
  Healthy fits in the original 30 had ratio >= 0.265. Root cause of the collapse unknown (INLA's inner latent solve fails intermittently).
- Pooled alpha black_spruce->regen 1.884 (SE 0.0669); excluded 1.886; as-is 1.855. 14 of 554 series changed; m_eff 30 everywhere.
  Unexplained: candidates scatter from the originals more than random subsets do (|z| > 2 share 9.8% vs 5.3%, p = 0.013).

Omission found: the real-data fits so far were MAIN-EFFECTS ONLY. This was not deliberate; the model should have had interactions from the start.
Model specification (decided by the user 2026-10-19):
- INLA weibullsurv variant 1: S(t) = exp(-(lambda t)^alpha), lambda = exp(eta), t in decades; hazard h(t|x) = alpha t^(alpha-1) exp(alpha eta). Intercept is in eta.
  Hazard ratio = exp(alpha * beta), NOT exp(beta). The ecology summaries written earlier used exp(beta) and understate effects (signs and rankings unaffected); recompute.
- GDD x disturbance interactions, GDD only (not CMI): ":" products with both main effects present: gdd_x_fire, gdd_x_pest, gdd_x_harv_light, gdd_x_harv_heavy, gdd_x_harv_any = z_gdd_3 * flag.
  Main effect of a disturbance = effect at mean GDD (z_gdd_3 is centred over all intervals of the final plot set).
- TSD: one term per disturbance TYPE (tsd_fire, tsd_pest, tsd_harv = any harvest). RAW time since the disturbance in DECADES; not centred, not averaged, not capped. 0 = at the disturbance.
  Structural zero where the flag is 0 (nested "x" structure; the flag main term carries the on/off effect).
- Term mask: fire_term -> gdd_x_fire + tsd_fire; pest_term -> gdd_x_pest + tsd_pest; harvest 3level -> gdd_x_harv_light + gdd_x_harv_heavy + tsd_harv; harvest any -> gdd_x_harv_any + tsd_harv.
- Thresholds: unchanged (>= 15 events per kept transition / per level). NO additional threshold for interactions or TSD: allowed wherever the disturbance main term exists (about 53 transitions).
  Basis: Step D (stepD_oracle_vs_imputed.R) needs no stricter threshold for interactions under the oracle; under real imputation interactions are systematically attenuated toward zero
  (not fixable with more imputations). The older rule "min(n_dist, n_undist) >= 20" came from an earlier power-curve analysis (derivation script not found) and was superseded.
  Planned reporting: flag interactions with 15-19 disturbed events as lower confidence; >= 20 subset as a sensitivity check; state the attenuation limitation.
- COV_CONT final (gdd_3, cmi_sum, soil_pc1). No domain random effect for now. Full 30-imputation pipeline to be rerun under RUN_TAG "int", outputs in data/processed/fit_int/.

Missing TSD (decisions):
- TSD source (decision A): `*_imputed` (minimum-plausible TSD, a LOWER bound; an undated disturbance first seen at a visit gets that visit's year). `*_imputed_prev` (upper bound) not used.
  Consequence: before exclusions about 24% (fire), 45% (pest), 45% (harvest) of flagged intervals have TSD exactly 0 (22 / 42 / 36% of event intervals). Sensitivity run with `imputed_prev` not done.
- 204 intervals in 72 plots (65 event intervals) have a flagged disturbance without TSD; all because the disturbance is already active at the plot's first visit with no year (NA carries through the episode).
- Origine-level disturbances (10 plots: 1 BR, 9 CT): TSD = W_j + (visit year - first-visit year) = t_abs ("W_j trick"), user-approved.
- Perturb-level (partial) disturbances (62 plots: 44 CP, 1 CE, 17 EL): cannot be dated -> plots EXCLUDED (tsd_excluded_plots.rds). Cost: 76 of 8,879 events (0.9%); three terms lost, each by one event
  (fir->trembling_aspen pest term 15->14; regen->other_pioneer harvest 3level->any, light level 15->14; white_birch->trembling_aspen harvest any->none). 91 kept transitions unchanged.
- Rescue by assigning ages to some perturb plots was examined and rejected: W_j (stand age class) bounds a partial cut only from above; rescuing all three terms needs plots with W_j 30-50.
  Cored-tree ages (dendro_arbres_etudes): 7609604002 and 7501401202 show one surviving cohort and no information on the CP date. 7501702402 shows a stand-replacing event 1975-77 and station_pe
  records origine CT there, i.e. a clear-cut (heavy), so it would not rescue the light level either.
- Open observation, NOT checked: 03_build_covariates.R reads origine/perturb from pee_ori_sond only (line 52). station_pe also has origine/perturb codes (e.g. CT at visits 2-3 of plot 7501702402
  where pee_ori_sond has none). Whether undated or station-only clear-cuts are mis-coded as light harvest or missed is unknown.
- Open interpretation risk, NOT tested: flags and TSD are taken at the END visit of the interval, so negative tsd coefficients may partly reflect the timing convention. A start-visit TSD variant exists as an option.

Files changed 2026-10-19 (originals copied outside the project):
- covariate_spec.R: 8 new terms in add_model_covariates / all_terms / term_mask, checks on the TSD columns (no NA, >= 0, zero where the flag is 0).
- 04_prepare_inla.R: TSD_SOURCE <- "imputed", W_j trick for origine plots, exclusion of plots with undated flagged intervals (writes tsd_excluded_plots.rds), tsd_fire/pest/harv in decades in the interval table.
- 05_fit_inla.R and 07_redraw_failed.R: RUN_TAG "int", OUT_DIR data/processed/fit_int (the old untagged stage1.rds, stage2.rds, scaling.rds in fit/ must not be resumed or overwritten).
- Rerunning 04 overwrites visit_intervals.rds, transition_set.rds, plot_info.rds; results in fit/ stay but no longer match those tables.
Checked on staged data (R in the cloud): files parse; 62 plots / 183 flagged intervals excluded; 91 kept transitions, 3-14 terms each; only the three terms above change. NOT checked: the .gpkg read in 04, any INLA fit with the final spec.
Run order: 04_prepare_inla.R -> test_interaction_inla_v2.R (one transition, real spec) -> 05_fit_inla.R -> 07_redraw_failed.R if needed.
After the run: recompute hazard ratios as exp(alpha*beta); consider the end-visit TSD timing test; Methods_draft.docx not yet updated.

## Status 2026-10-20: "int" run complete (30 valid fits for all 91 transitions); final pooling done
Decisions (user, final; do not reopen):
- INLA weibullsurv variant 1 for ALL transitions. No MLE or Poisson-offset engine as the analysis fit; no transition is fitted differently from another.
- Failed/invalid fits are replaced by fresh imputations (per-transition replacement). Rubin pooling is per parameter, so transitions may be pooled over different imputation indices.
- 30 valid fits per transition are pooled (the same number for every transition).
- Model: GDD x disturbance interactions (":" products) and one TSD term per disturbance type (see "Status 2026-10-19").

What was run (tag int, data/processed/fit_int/):
1. 05_fit_inla.R: imputations m1-30, all 91 transitions (2730 fits).
2. 07_redraw_failed.R (MAX_CAND 80): candidate imputations m31-114 (84), all 91 transitions, to replace failed fits (ratio screen only at that stage).
3. Screen revision. The ratio rule alone let partial collapses through (black_spruce->regen m21, m37, m73: intercept stuck or shifted). FROZEN SCREEN (below) applied to every fit.
4. 07b_bsregen_sweep.R: black_spruce->regen only (the one transition without 30 valid fits), candidates m115-156 (42), same data build (identity check against a saved dataset), same screen. 14 valid, 7 INLA errors.
   The sweep was only a PREFILTER: it fitted one transition per imputation and stored no posterior samples/extras.
5. 07c_bsregen_cov.R: regular run_one_imputation() (all 91 transitions, all standard outputs: extras, covariance, diagnostics) for the 14 sweep-valid imputations (m116,117,120,122,123,127,133,140,145,148,150,151,153,156).
   Rule: the REGULAR full-pipeline fit decides validity. m145 and m148 (valid in the sweep) gave INLA errors in the regular run, so 12 of 14 count. Their fits for the other 90 transitions are stored but not used.
6. 08_final_pool.R: pulls all regular-pipeline fits, applies the frozen screen, takes the FIRST 30 valid fits per transition in increasing m (selection by index only), Rubin-pools estimates and the joint covariance.

FROZEN SCREEN (identical for every transition, fixed before applying it to new fits). A fit is valid if ALL hold:
- INLA returned no error and est/sd are not NA;
- alpha-SD ratio sd_alpha / (0.78 * alpha / sqrt(events)) >= 0.1 (0.78*alpha/sqrt(n) is the asymptotic Weibull-shape MLE SE; the rule catches collapsed alpha SDs of about 2e-4);
- |intercept - ref| / ref_sd <= 4, where ref and ref_sd are the median intercept and median intercept SD of the transition's ratio-passing fits in m1-114
  (black_spruce->regen: -4.985135842 and 0.06137571849; 08_final_pool.R recomputes them and checks this). Across ~8,800 ratio-passing fits only m21, m37, m73 exceed |z| = 4 (z = 11 to 44; all other fits < 1.6).
The screen detects collapse signatures. It does NOT prove a fit is correct (no MLE comparison exists for healthy fits; an MLE comparison was done only for the saved collapsed datasets of the main-effects run: 06c).

Failure statistics (all 91 transitions, frozen screen, fit_validity_int.csv):
- m1-30: 36 of 2730 fits invalid (1.3%); m31-114: 107 of 7644 (1.4%); m115-156 (14 pre-selected imputations): 7 of 1274 (0.55%).
- Reasons: alpha-SD ratio 115, INLA error/NA 32, intercept 3.
- Invalid fits concentrate in 8 transitions: black_spruce->regen 96, regen->fir 18, fir->sugar_maple 18, other_temperate->white_birch 9, regen->red_maple 6, jack_pine->regen 1, trembling_aspen->regen 1, white_spruce->red_maple 1.
- Outside black_spruce->regen: 16 of 2700 invalid in m1-30.
- black_spruce->regen: valid 10 of 30 (33%) in m1-30; 10 of 84 (12%) in m31-114; sweep m115-156 14 of 42 (33%) before the regular-run check, 12 of those 14 valid in the regular run. Total valid 32; the first 30 by m are used (m1-114 plus m116 ... m151; m153 and m156 unused).
  Why the rate fell from 33% to 12% and recovered to 33% in the single-transition sweep: UNEXPLAINED (not tested; ideas: 91 sequential fits per worker vs one, resource contention).
- Reproducibility: refits of the main-effects collapsed datasets under default settings reproduced the collapse in 7 of 8 black_spruce->regen cases (06c). But m145 and m148 fitted in the sweep and errored in the regular run, so INLA is not perfectly deterministic across run contexts.
- Selection concern: replacing/pooling only valid fits selects on INLA convergence. Argument for ignorability: failure does not depend on the imputed entry summaries (tested for black_spruce->regen only, weak power) plus the failure being numerical. Not proven.
- Low alpha-SD ratio fits remain in the pool (they pass the frozen screen): black_spruce->regen m30 (ratio 0.36), m99 (0.18), m109 (0.14) [m145 (0.15) not used]; regen->fir m36 (alpha 1.945, ratio 0.21); many fits at alpha 1.0-1.1 in sparse transitions with ratio 0.3-0.5 (not interpreted).
  Sensitivity for black_spruce->regen without m30/m99/m109 (29 fits of the earlier 32-fit pool): alpha 1.830 -> 1.838, every coefficient moved < 0.1 SE. 10_sensitivity_posthoc.R repeats this for all transitions.

Final pooled results (data/processed/fit_int/): pooled_int_final.csv (727 rows; estimates, SE, CI, fmi, m_eff = 30 everywhere), pooled_int_final_cov.rds (named list by transition: Rubin-pooled joint covariance of intercept, betas, alpha), fit_validity_int.csv (every fit, screen result, used yes/no).
- black_spruce->regen (30 fits): alpha 1.828 (SE 0.070, fmi 0.53); heavy harvest 1.548 (0.078); fmi for alpha is high (half of its uncertainty is imputation).
- The other 90 transitions are identical to pooled_int_redraw.csv. pooled_int_redraw.csv's black_spruce->regen row (23 fits including m21/m37/m73) and the 07b/07c black_spruce->regen pooled files are SUPERSEDED.

Known limitation (documented, not corrected): joint covariance and alpha.
- The joint covariance comes from 1000 draws of inla.posterior.sample(). Per the INLA documentation, the hyperparameter (alpha) is drawn only from the configurations of the numerical integration (int.strategy "grid", dz 0.75, diff.logdens 6 = INLA defaults; model_fitting.R of the simulation uses the same call). The intercept/betas are drawn from their Gaussian conditional on alpha.
- When alpha is poorly determined (sparse transitions, alpha near 1, long right tail) the sampled alpha is narrower than INLA's own alpha marginal. Sampled alpha SD < 0.9 x INLA's in 91 of 2730 used fits (cedar->yellow_birch, fir->sugar_maple, red_maple->sugar_maple, red_maple->yellow_birch, trembling_aspen->red_maple, black_spruce->white_birch, ...; lowest pooled ratio 0.78 in black_spruce->white_birch). All beta/intercept diagonals agree with the pooled SEs (0.99-1.015).
- Pooled estimates and SEs (from INLA's own summaries) are NOT affected. The alpha variance and alpha-beta covariances used for hazard-ratio CIs are somewhat too small in those transitions. We stay with INLA's documented procedure; flag these transitions in the HR table (column alpha_cov_ratio_med in hazard_ratios_int.csv).
- A post-hoc correction (law of total variance, INLA alpha mean/variance substituted) was tried and rejected: it made 124 of 19,080 beta SDs deviate by > 10% from INLA's (2 before) or produced non-positive-definite matrices (87 of 2730). The proper fix would be refitting with a finer grid (not done). INLA's own fixed-effect marginals use the same integration grid, so they are not an independent check of the tails.

Hazard ratios (09_hazard_ratios.R) and sensitivity analyses (10_sensitivity_posthoc.R) were written 2026-10-20; NOT yet run on the user's machine.
- HR = exp(alpha * c'beta) for a contrast vector c over the model terms. Per imputation: log HR = alpha * c'beta (INLA posterior means from mi_raw), variance by the delta method with the stored joint covariance (gradient: alpha*c for beta, c'beta for alpha). Rubin-pool the log HR (Barnard-Rubin df), then exponentiate. Contrasts: per-SD GDD/CMI/soil; disturbance vs none at GDD z in {-1,0,1} and TSD in {0,1,3} decades; GDD effect within a disturbance; GDD x disturbance ratio of HRs; HR per +1 decade since disturbance.
- Sensitivity (post-hoc, no refits): (S1) GDD x disturbance interactions restricted to >= 20 disturbed events (15-19 flagged as lower confidence); (S2) pooling without fits with alpha-SD ratio < 0.5; (S3) pooling all valid fits (about 128 imputations per transition) instead of the first 30.
- Sensitivities that need new INLA runs and are NOT scripted: TSD upper bound (`imputed_prev`; 04 has TSD_SOURCE but rerunning 04 overwrites visit_intervals.rds, transition_set.rds, plot_info.rds, so a separate output tag is required first), and start-visit flags/TSD instead of end-visit (timing convention may partly drive negative TSD coefficients). Each is a full 30-imputation run plus redraw plus the frozen screen.

Still open: station_pe vs pee_ori_sond origine/perturb codes (undated or station-only clear-cuts possibly mis-coded as light harvest); the 62 excluded plots (76 events) to be reported; Methods_draft.docx not updated; MLE comparison for healthy fits does not exist.
