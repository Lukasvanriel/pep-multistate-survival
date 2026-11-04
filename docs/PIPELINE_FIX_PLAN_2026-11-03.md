# Pipeline fix plan (2026-11-03): specific changes per item, for review before implementation

Source list: `docs/PIPELINE_NOTES.md`, section "FINAL PIPELINE LIST 2026-11-03" (items A1–D7).
This document gives one entry per item, in four parts:
- **Goal:** what the item must achieve in the whole pipeline.
- **Now:** what the code does today.
- **Change:** the specific edit, by file and function.
- **Check:** how to verify it after implementation.

Items with an open choice say so under **Choice**, with a suggestion. Nothing in this document has been implemented.

**Review of 2026-11-04 (user + Claude):** every item was reviewed and every choice is now decided (marked **Decided 2026-11-04**). New items: **A11**, **C5** (bug fix in the P1 weights), **D8**. Additions to A1, A9, D2, D4, D5, D6, D7. **C1 rewritten** (disturbance flags switch on at their onset in P1's history, TSD yearly). Summary of decisions and the new limitations at the end. Still nothing implemented.

Conventions used throughout:
- Survival time in the rows stays in **decades** (`TIME_SCALE = 10`, `04b_build_rows.R`).
- TSD (time since disturbance) is kept in **years**.
- Weibull variant 0 means S(t) = exp(−e^η t^α). β is then a log hazard ratio.
- The P1/P2 code works with the variant-1 form S(t) = exp(−(λt)^α) in years.

Pipeline order: 03c → 03z → 04 → 05 → 07_redraw → 08 → 09 → 12 → sensitivity/10, 11.

---

## A. Data preparation

### A1 Missing GDD/CMI for 5 later-added plots

**Goal:** every interval that reaches INLA has real climate values. INLA silently treats an NA covariate as 0, i.e. the mean climate.

**Now (diagnosed 2026-11-03):**
- 14 visits in 5 plots have `gdd_3` and `cmi_gs` NA in `data/build/climate_covariates.rds`, all with `n_years_used = 0`. The plots are 7109604001, 7209002301, 7209002302, 7209605702 and 7501600102.
- These visits produce 9 intervals in `visit_intervals.rds`.
- Cause: none of the 5 plots is in the extraction checkpoint `data/build/climate_by_year_raw.rds`. That file was saved 2026-09-10 and holds 11,469 plots and 63 years.
- 03c step 4 skips every year already in the checkpoint (`if (yr %in% done_years) next`). The checkpoint is per year, not per plot-year, so plots added to the analysis after 10 Sep were never extracted.
- Step 5 then returns NA (`nrow(sub) == 0`). The only signal is a printed FLAG.

**Change:**
1. `03c` step 4, made plot-aware:
   - compute the plot-years needed (every visit's trailing window) that are missing from `long_dt`;
   - for every year with missing plots, download and extract that year **only for the missing plots** (`pts_ll` restricted to them) and append the rows to `long_dt`;
   - keep the checkpoint save and the raster deletion as now.
   - For these 5 plots the windows span 1961–1996, so 36 years of sg and cmi zips must be downloaded again (the raw rasters were deleted after the first extraction).
2. `03c` end: replace the printed FLAG with `stop()` when any visit inside the 1960–2022 data range has `n_years_used == 0`.
3. `03z` end and `04` (after building `iv`): `stopifnot(!anyNA(iv[, .(gdd_3, cmi_gs, soil_pc1)]))`.
4. (2026-11-04) The 03c guard in point 2 checks the **full** window: `stop()` unless every visit has `n_years_used == win_end - win_start + 1`. The windows are already clipped to 1960-2022, so visits in 2023-2025 with shorter windows pass; a visit extracted only in part would pass a `== 0` check.
5. (2026-11-04) 04 also checks the first-visit covariates written to `plot_info.rds` for C1: no NA in `gdd_3`, `cmi_gs`, `soil_pc1`. A plot's first visit is never an end visit, so the `iv` check does not cover it.
6. (2026-11-04) Rerun `04x_exclusion_flow.R` after 04. These 5 plots / 9 intervals / 4 changes are exactly the gap between the current table (11,286 plots / 36,251 intervals / 8,770 changes) and the 04x counts in the notes (11,281 / 36,242 / 8,766, attributed there to the 03z join).

**Check:**
- The new `climate_covariates.rds` has 0 NA `gdd_3`/`cmi_gs`.
- All other visits are **identical** to the current file (the checkpoint rows for the old plots are unchanged).
- 04 reports the same plots, intervals and events as now.
- The 5 plots' values are plausible next to their neighbours. 7209002301 and 7209002302 are a plot pair and should be near-identical.
- 04x matches 04 again, and its 03z-join step no longer removes these 5 plots.

### A2 Harvest term: light and heavy together, never "any"

**Goal:** harvest is always modelled as light and heavy (each vs none), so the two levels can always be compared. It is included whenever harvest matters for the transition; the reference level (no harvest) has no minimum.

**Now:**
- `04_prepare_inla.R` (lines 140–144): `harv_term` is "3level" (every level ≥ 15), else "any" (none ≥ 15 and light + heavy ≥ 15), else "none".
- `covariate_spec.R` builds `harv_any`, `gdd_x_harv_any` and `tsd_harv`.
- `recount_events()` (04b) checks the old rule.

**Change:**
1. 04, `harv_term := fcase(!kept, "none", n_harv1 + n_harv2 >= MIN_EVENTS, "3level", default = "none")`. The label "3level" is kept to limit edits; it now means "light + heavy vs none". Update the comment block (lines 115–121).
2. `covariate_spec.R`:
   - delete `harv_any`, `gdd_x_harv_any` and `tsd_harv` from `add_model_covariates()`, `all_terms()` and `term_mask()`;
   - `term_mask()` keeps `harv_light`, `harv_heavy`, `gdd_x_harv_light` and `gdd_x_harv_heavy` for "3level";
   - the per-level TSD terms are covered under A4 and A6.
3. `04b` `recount_events()`: the harvest condition becomes `harv_term != "3level" | (n_harv1 + n_harv2 >= min_events)`.
4. `tests/04c_check_rows.R` line 95: drop the `harv_any` count.

**Check:**
- 04 prints 69 transitions with a harvest term on the main build (17 "3level" + 21 "any" today, plus the ones that had "none" only because of the reference).
- `grep -n "harv_any\|tsd_harv\b"` over `R/*.R`, `R/sensitivity` and `R/tests/04c` returns nothing.

### A3 Fire and pest unchanged

**Goal:** fire or pest term at ≥ 15 events under that disturbance; transitions need ≥ 15 events.

**Now:** already the rule (04 lines 137–139).

**Change:** none. `MIN_EVENTS` moves to `config.R` (B7).

**Check:** 04 prints fire 14 and pest 28 terms (main build), unless A1 changes a count. It should not, since A1 only fills climate.

### A4 Separate harvest TSD per level

**Goal:** time since light harvest and time since heavy harvest are separate terms, so the two levels' recovery curves are not mixed.

**Now:**
- one `tsd_harv` term for any harvest;
- 04 takes `tsh_imputed_mid` (per `TSD_SOURCE`) and stores it in decades.
- Checked in `03a visit_values()`: per visit, the harvest level and the harvest TSD come from the **same event**, so splitting by level is consistent.

**Change:**
1. 04 stores the raw harvest TSD in years as `tsd_harv_y` (0 where `harv == 0`), alongside `tsd_fire_y` and `tsd_pest_y` (see A5).
2. `covariate_spec.R add_model_covariates()`:
   - `tsd_harv_light = harv_light * f(tsd_harv_y)`;
   - `tsd_harv_heavy = harv_heavy * f(tsd_harv_y)`;
   - `f()` is the TSD form (A5).
3. `all_terms()` and `term_mask()` gain `tsd_harv_light` and `tsd_harv_heavy`, present when `harv_term == "3level"` (subject to A6).

**Check:**
- In every interval at most one of `tsd_harv_light` and `tsd_harv_heavy` is non-zero.
- Each is 0 wherever its level flag is 0.

### A5 TSD form log(TSD in years + 1)

**Goal:** all TSD terms follow the a-priori form chosen on 2026-11-01. The effect at TSD = 0 is the main term; the change is fast early and levels off later. Linear TSD stays available as S24.

**Now:** 04 converts TSD to decades (lines 94–97). `covariate_spec.R` checks and uses the decade values as is.

**Change:**
1. 04: keep TSD in **years**. Rename the columns `tsd_fire_y`, `tsd_pest_y` and `tsd_harv_y`, with a structural 0 where unflagged. Remove the `/10`.
2. `covariate_spec.R add_model_covariates()` builds the model columns from a config switch `TSD_FORM` (B7):
   - `"log1p"`: f(x) = log(x + 1);
   - `"linear"`: f(x) = x / 10 (decades, S24).
   - Model columns: `tsd_fire = fire_f * f(tsd_fire_y)`, `tsd_pest = pest_f * f(tsd_pest_y)`, plus A4.
   - Keep the existing checks: no NA, ≥ 0, 0 where unflagged.
3. Because P1 and P2 take their covariate matrices from `add_model_covariates()` (via `prepare_fit_inputs`), they automatically use the same form. No separate change is needed there.

**Check:**
- `max(tsd_*)` ≈ log(176) ≈ 5.2.
- With `TSD_FORM = "linear"`, the old decade columns are reproduced exactly (except the A4 split).

### A6 Per-level GDD × disturbance and TSD terms

**Goal:** each harvest level has its own main term, GDD interaction and TSD term, even when the level has few events (decided: the prior keeps them finite and the intervals honest).

**Now:** the GDD × level terms exist for "3level"; TSD is shared.

**Change:** `term_mask()`: for "3level", include `harv_light`, `gdd_x_harv_light` and `tsd_harv_light`, and the same three for heavy.

**Decided 2026-11-04: (a), main term only for a level with 0 events.** 0 vs >= 1 event is a real boundary: with 0 events the estimate has no finite maximum, so one side of its interval comes only from the prior. No other threshold. Options considered: On the main build only other_temperate→yellow_birch, heavy (0 events) is concerned.
- **(a) Main term only (suggested).** With 0 events the GDD and TSD terms carry no information. They would only widen the main term's interval through the prior, and the bound in D4 would then depend on the prior of the TSD slope.
- **(b) All three.**

Implementation of (a) in `term_mask()`: `gdd_x_harv_heavy` and `tsd_harv_heavy` require `n_harv2 > 0`; likewise for light and `n_harv1`.

**Check:** `term_mask()` printed for other_temperate→yellow_birch shows the expected set.

### A7 other_temperate→yellow_birch keeps the heavy-harvest term

**Goal:** the heavily harvested rows stay in the transition with their own term, and are not counted as unharvested.

**Now:** no harvest term, so 303 heavily harvested intervals (3,455 plot-years, 0 events) count as unharvested.

**Change:** none beyond A2 (light + heavy = 22 ≥ 15 includes the term) and A6. Reporting is in D4.

**Check:** after 04, `harv_term == "3level"` for this transition.

### A8 Long TSDs kept

**Goal:** no cap on TSD. The log form already compresses long TSDs.

**Now:** no cap anywhere.

**Change:** none. State the TSD range and the stand-age dating in the Methods.

**Check:** none needed.

### A9 Dates of heavy disturbances inside each interval (input for C2)

**Goal:** P2 needs, per change interval, the date of a heavy harvest or fire that falls inside the interval (recorded year, or imputed date for undated events).

**Now:** 04 carries end-visit flags, TSD and `*_dating`, but no date.

**Change:** 04 adds two columns to `visit_intervals`:
- `D_harv_heavy = year_end − tsd_harv_y` if `harv == 2` and `tsd_harv_y < year_end − year_start`, else NA;
- `D_fire = year_end − tsd_fire_y` if `fire == 1` and `tsd_fire_y < year_end − year_start`, else NA;
- (2026-11-04) `D_heavy = pmin(D_harv_heavy, D_fire, na.rm = TRUE)`: the earlier date (C2 decision). C2 uses `D_heavy`.

TSD follows `TSD_SOURCE`, so S11 (lower/upper bound) moves the date consistently. The same logic was used and checked in `tests/22_cap_test_new_rows.R` (1,537 X→regen exits on the main build).

Edge case to log, not solve: under rule C(i), a heavy harvest followed by a light harvest inside the same interval shows the light one at the end visit, so the heavy date is not seen. 04 counts and prints the intervals where `harv_start != 2`, the end visit shows light, and a heavy record exists in `data/build/disturbance_events.rds` within the interval. Fire isn't affected (fire is 0/1).

Methods wording (2026-11-04): "fire", not "stand-replacing fire". Partial fires are merged into the 0/1 fire flag, so they also trigger the exit rule (C2).

**Check:**
- Number of X→regen change intervals with a date = 1,537 (heavy harvest 1,223, fire 315) on the main build.
- Every date lies in (`year_start`, `year_end`].

### A10 Term counts and TSD value per row

**Goal:** the ≥ 15 rules are evaluated on a count that is the same in every imputation, so every imputation has the same formula.

**Now:** counts use the end-visit flag of the event interval. That stays valid because rows are not cut.

**Change:** none.

**Decided 2026-11-04: (a), end-visit value.** Options considered:
- **(a) End-visit value for every row (now; suggested).** It is consistent with the flag covering the whole interval, and it is what the cap test (`22_cap_test_new_rows.R`, `CAP_CUT=FALSE`) validated: 0 screen failures, 2 transitions above the cap.
- **(b) The value at the row's end:** the end visit for observed and stub rows, the exit time for event rows. This sets TSD = 0 on the X→regen event rows placed at the disturbance (C2), which raises their linear predictor. It must be cap-tested before adoption.

### A11 No GDD main term where a transition has 0 undisturbed events (new 2026-11-04)

**Goal:** no coefficient that only the prior determines. The GDD main term (`z_gdd_3`, beta_gdd) is the GDD effect in plots with
none of the transition's disturbances. With 0 events there, the data say nothing about it.

**Now:** every kept transition has `z_gdd_3`.

**Change:**
1. 04 counts, per transition, the events with none of its disturbance terms (fire = 0 if it has a fire term, pest = 0 if a pest
   term, harv = 0 if a harvest term) and stores it in `transition_set` as `n_undist`. Main build: 0 only for regen->sugar_maple
   (67 events: 0 unharvested, 18 light, 49 heavy; no fire or pest term).
2. `term_mask()`: `z_gdd_3` is FALSE where `n_undist == 0` and the transition has at least one disturbance term. All other terms
   stay: intercept, CMI, soil, harvest main terms, GDD x harvest, TSD.
3. Meaning for these transitions: `gdd_x_harv_light` / `gdd_x_harv_heavy` become the GDD slopes WITHIN light / heavy harvest.
   The unharvested GDD slope is assumed 0 (untestable: no events). The baseline stays "no harvest". Harvest vs none therefore holds
   at mean GDD only (D4 point 5).

**Check:** `term_mask()` for regen->sugar_maple has no `z_gdd_3`; no other transition loses it on the main build.

---

## B. Model and fitting

### B1 Weibull variant 0 everywhere

**Goal:** coefficients are log hazard ratios. sim02 showed that variant 1's alpha is biased low with few events per coefficient, and variant 0's is not. sim03 showed variant 0 is as usable on the real data.

**Now:** `05a fit_inla_models_pc()` line 48 has `variant = 1`; `sensitivity/11` line 106 has `variant = 1`.

**Change:**
1. Both calls use `variant = WEIBULL_VARIANT` (B7, value 0).
2. `05_fit_inla.R` and `07_redraw.R`: delete the unused `source(file.path(SIM_DIR, "model_fitting.R"))` and `SIM_DIR`, removing the last dependency on the simulation folder.
3. `run_meta` records `WEIBULL_VARIANT`.

**Check:** `run_meta_<tag>.rds$inla_args$control.family$variant == 0`.

### B2 Priors

**Goal:** N(0, 5²) on all coefficients, flat intercept, pc.alphaw(0.1) on alpha.

**Now:** coefficients use INLA's default (precision 0.001); alpha is already pc.alphaw(0.1).

**Change:** `fit_inla_models_pc()` and `sensitivity/11` add `control.fixed = list(mean = 0, prec = 1 / COEF_PRIOR_SD^2, mean.intercept = 0, prec.intercept = 0)`.

**Check:** `run_meta` shows the `control.fixed` values.

### B3 Intercept scale

**Goal:** pooling and the screen use the intercept exactly as INLA reports it, with a correct SD.

**Now:** `tidy_fit()` (line 301) and `extract_theta()` / `draw_theta_m()` (lines 107, 138) subtract `log(10)`. That converted variant 1 from decades to years. Under variant 0 the conversion is − α·log 10, so the subtraction is wrong, and the SD of a converted intercept would also need α's variance.

**Change:**
1. `tidy_fit()`: the intercept `est` is INLA's mean, without the subtraction. Add a comment "decade scale, variant 0".
2. `mi_raw`, the screen references (07, 08) and pooling all stay on this scale.
3. Reporting only, in **08** (decided 2026-11-04: 08 already loads each fit's joint covariance, needed for cov(b0, α)): intercept in years per fit = b0 − α·log 10, with variance var(b0) + (log 10)²·var(α) − 2·log 10·cov(b0, α) from the stored joint covariance; then Rubin-pooled.

**Check:** the 07 screen references are computed from the new m1–30, so the scale is consistent by construction.

### B4 Variant-0 → P1/P2 conversion

**Goal:** P1 (`lambda_matrix`, `p1_candidates`, `p1_weights`) and P2 (`moment_match_weibull_naive`, the new draw C3) keep working in the variant-1 form in years.

**Now:** `extract_theta()` and `draw_theta_m()` give the intercept in years and β unchanged (variant-1 scale).

**Change:** in both functions, for a variant-0 fit with intercept b0 (decades), coefficients β and shape α:
- `theta$intercept[k] = b0 / α − log(10)`
- `theta$beta[k, ] = β / α`
- `theta$alpha[k] = α`

Derivation: S = exp(−e^η (t/10)^α) = exp(−(λ₁t)^α) with log λ₁ = η/α − log 10.

- **`draw_theta_m()`:** per posterior draw, using that draw's own α.
- **`extract_theta()`:** posterior means (b0/α evaluated at the means). This is used only for Stage-2 seeds and is acceptable.
- **`extras$theta_m`:** store the **variant-0 draws** (α, b0, β), so the stored values match the fitted model and the conversion is visible.

**Check:** unit test with one fit. For 3 random covariate rows, compare S(t) from INLA's variant-0 parameters with S(t) from the converted λ₁, α at t = 5, 20 and 50 y. They must agree to 1e-10.

### B5 Screen and replacement

**Goal:** the same frozen screen as the main run, with references from this run.

**Now:** 07 and 08 apply (1) no error/NA, (2) alpha-SD ratio ≥ 0.1, (3) |intercept z| ≤ 4, with references from m1–30. Failures are replaced by new imputations.

**Change:** none in the logic (decided 2026-11-03: both rules kept). The references are recomputed automatically because `screen_ref_<tag>.csv` belongs to the new run tag.

**Check:** 07 prints the reference table before any candidate is fitted.

### B6 Cap diagnostic per fit

**Goal:** detect fits distorted by INLA's truncation cap, which the screen misses. In the cap test, black_spruce→regen passed the screen with heavy-harvest log HR 1.69 against 4.58 from the cap-free fit.

**Now:** nothing per fit (only `tests/08f`, after the fact, for transitions with failed fits).

**Change:**
1. Move the penalised-ML fit (`fit_map_v0()` in `tests/22_cap_test_new_rows.R`: the same variant-0 likelihood, N(0, `COEF_PRIOR_SD`²) penalty, no cap) and the helper `H_at()` into `05a_fit_funs.R`.
2. In `run_one_imputation()`, after `tidy_fit()`, for every fitted transition compute:
   - `maxH_inla`: max over rows of e^η·truncation^α at INLA's posterior means;
   - `n_above_inla`: rows above 22.85;
   - `maxH_map`;
   - `n_above_map`;
   - `z_map_max`: largest |INLA − MAP| / MAP SE over the coefficients.

   Add them as columns to `mi_raw`, one value per fit, repeated over its parameter rows like `mlik`.
3. 08 writes `cap_check_<tag>.csv` (per transition: share of used fits with `maxH_map > 22.85`, median `z_map_max`).

These columns are **diagnostics, not a validity gate** (user decision 2026-10-22: maximum likelihood is not the gate).

Cost: about 1 s per fit (cap test: 180 penalised fits in 0.9 min on 4 workers).

**Check:**
- (Revised 2026-11-04) The functions moved into 05a (`fit_map_v0()`, `H_at()`) reproduce `cap_test_new_rows_m01_nocut/results.csv` when given the test's rows (arm B, no cut). The values of the new imputation 1 will differ, because C1, C3 and C5 change the entries and exits.
- The caveat list of cap transitions comes from B6 after the final run, not from the cap test (see "What to expect").

### B7 Switches in `config.R`

**Goal:** sensitivity runs only change `config.R` plus the run tag.

**Change:** add to `config.R`:
- `WEIBULL_VARIANT <- 0L`
- `COEF_PRIOR_SD <- 5` (S21 = 2.5; S22 = `1/sqrt(0.001)`)
- `TSD_FORM <- "log1p"` (S24 = "linear")
- `MIN_EVENTS <- 15L`, moved from 04 (S17 = 20, S23 = 10)
- `XREGEN_EXIT_RULE <- TRUE` (C2; FALSE only as a possible check)

`recount_events()` default reads `MIN_EVENTS`; 09: `LOW_EVENTS_MAX <- MIN_EVENTS + 4L` (the 15–19 band; 20–24 under S17, 10–14 under S23).

Addition 2026-11-04: `recount_events()` is today called only in `tests/04c` (line 90), although 04b's header calls it "the very last check before INLA". Call it in `run_stage1()` (05a) right after `build_rows()`, before the first fit.

**Check:** `grep -n "MIN_EVENTS <-\|variant = \|prec = "` finds only `config.R` and the fit function.

---

## C. Imputation

### C1 P1 covariates: first visit, disturbance flags switching on at their onset (REWRITTEN 2026-11-04)

**Goal:** P1 simulates the history before the first visit with the covariates that applied at each time: climate and soil of the
first visit; each disturbance only from its onset, with TSD growing.

**Now:** `prepare_fit_inputs()` line 80: `f1 <- iv[interval == 1L]`, so `X` holds end-visit values of interval 1 (climate window,
flags, TSD), constant over the whole simulated history.

**Change (decided 2026-11-04):**
1. 04 writes the first-visit covariates to `plot_info.rds`: `gdd_3`, `cmi_gs`, `soil_pc1`, `fire`, `harv`, `pest` and `tsd_*_y`
   at visit 1 (from `dt2`, the first row per plot). NA check: A1 point 5.
2. Climate and soil: first-visit values, scaled with the interval scaling, constant over P1's history.
3. Each disturbance type flagged at the first visit switches on at its onset, on P1's clock: onset = W_j - TSD at the first visit
   (TSD per `TSD_SOURCE`). Before the onset: flag, GDD x flag and TSD = 0. From the onset: flag on, TSD = time since the onset.
   Onset <= 0 (disturbance at or before stand origin): on from time 0, TSD = t - onset.
4. Heavy disturbance dated more recently than the stand age (main build: 282 heavy-harvest and 141 fire plots more than 5 y after
   stand origin, all dated; median 22 y and 15 y): handled by point 3 like any other, option (i). No forced move to regen, no
   change to W_j.
5. Draw: before the first onset, the current one-formula Weibull draw. From the first onset, yearly pieces with the covariates
   constant within a piece and TSD at the piece MIDPOINT. TSD depends only on calendar time and the onset, so it is known before any
   draw (not circular). The exit is drawn exactly for the piecewise-constant hazard (cumulative-hazard inversion piece by piece);
   the earliest competing exit wins.
6. The P1 weights use the same piecewise hazards up to the first visit and, after it, the covariates of the observed intervals
   (as in the fitted rows). Weight formula: C5.
7. theta converted per draw from variant 0 (B4).
8. Plots without a flag at the first visit (6,330 on the main build) are unchanged.

Options considered and rejected: (a) first-visit flags/TSD over the whole history (applies e.g. a 1996 light harvest to
1900-1996); (b) flags off (regen exits would use the undisturbed reference, which has 0 events for regen->sugar_maple); (c) only
stand-origin flags kept; fixed TSD after the onset (about 7 min faster per full run on the Mac, cruder).

**Cost (measured 2026-11-04 with `R/tests/23_p1_onset_bench.R`; Mac, one core, 300 flagged + 300 unflagged plots, theta of m01):**
P1 per imputation 0.10 min now -> 1.19 min (yearly TSD; fixed TSD 0.85 min). INLA per imputation 4.1 min. Whole imputation about
x1.26; full run (about 120 imputations, 6 workers) about +22 min. The benchmark code is a prototype (plain R loops).

**Check:**
- The base row equals the first-visit values for a sample of plots; onsets equal W_j - first-visit TSD.
- Unit test: for one plot with a fixed theta, 100,000 drawn exits from a piecewise sojourn match the exact piecewise survival
  (Kolmogorov-Smirnov p > 0.01).
- Plots without a flag get the same candidates as the current P1 for the same seed (unchanged draw).
- **P1 fallbacks** (no usable candidate) printed per imputation **and by first state**, and compared with the main run's 124-196. Watch the 681 plots that start in sugar maple without a harvest flag at the first visit (median W_j 70 y, IQR 50-90): P1 has essentially no direct regen->sugar_maple route before a harvest onset (0 unharvested events, A11), only indirect routes. The benchmark sample gave
  4 (current) vs 4-8 (onset versions) of 300 flagged plots: too few to judge. Check on the first full imputation.

### C2 Exit at the heavy disturbance for X→regen

**Goal:** implement the Methods rule. A transition to regen in an interval that contains a heavy harvest or fire (dated or undated) happened at that disturbance. Pest is not a heavy disturbance.

**Now:** P2 (`make_p2_moment()`) places every exit by moment matching, ignoring disturbances.

**Change:**
1. The P2 function passed to `build_rows()` (`make_p2_moment()` and its successor in C3) checks first: if `state_end == "regen"` and `D_heavy` (A9) is not NA for interval `r`, return T = t_start[r] + (D_heavy − year_start[r]) − E, clamped to [a_dur, b_dur].
2. This applies for every theta, including transitions whose fit failed, and in Stage 1 (added 2026-11-04), Stage 2 and every imputation. Stage 1: `run_stage1()` passes a P2 function that returns the rule date when the interval has `D_heavy` and the midpoint otherwise (today it uses the default `p2_midpoint`). The rule needs no fitted parameters, and Stage 1 seeds Stage 2, so its X→regen fits must use the same exit times.
3. Record the exit source per event in `extras$p2` (`exit_source = "rule" | "draw" | "midpoint"`) and count it in the diagnostics.
4. Not changed: P1's weighting bracket for plots whose first change is such an exit. It still uses the probability of leaving within the interval, not the known date. That is a small loss of information, stated here so a reviewer can object.

**Decided 2026-11-04: (a), the earlier date** (affects 1 interval on the main build). Options considered:
- (a) The earlier date (suggested; the first stand-replacing event makes the stand regen). The cap test used this.
- (b) The later date.

**Check:**
- `extras$p2` shows 1,537 `exit_source == "rule"` on the main build.
- Every rule exit equals its interval's `D_heavy`.
- `build_rows()` bracket checks pass.
- Floored segments rise for dates equal to the end visit; print the count.

### C3 P2 as a proper imputation (drawn exit)

**Goal:** the imputed exits vary between imputations as much as the data allow, so Rubin's between-imputation variance is not understated. Today's moment-matched point is the same in every imputation given θ_m.

**Now:** `make_p2_moment()` uses `moment_match_weibull_naive()`, a deterministic point.

**Change:**
1. New `make_p2_draw(theta, inp, diag)` in 05a, used in `run_one_imputation()`. `run_stage2()` keeps `make_p2_moment()`, the point for the seed models. Both apply C2 first.
2. Draw for the realised transition k, durations a < b from entry, H(t) = (λt)^α:

   U ~ Uniform(0, 1), Δ = H(b) − H(a), T = (H(a) − log1p(−U·(−expm1(−Δ))))^(1/α) / λ

   This is the inverse of the survival function conditional on (a, b]. It is written in the cumulative-hazard form so it stays exact when S(a) is tiny, where the probability form would hit the same underflow as INLA. Clamp T to [a, b].
3. Transitions without a usable fit (dropped, or θ not ok): a uniform draw on (a, b) instead of the midpoint. This is suggested so that every exit is drawn.
4. Random numbers come from R's RNG after `set.seed(SEED * 1000 + m)`, as now, so imputations are reproducible.
5. Methods: "drawn from the distribution … conditional on the transition occurring within the interval". The text matches the rule; only "dated" is removed (C4).
6. Limitation: this departs from the simulation-validated point (SIM4 was not run).

**Decided 2026-11-04: (a), realised transition only; uniform draw for transitions without a usable fit (point 3).** Options considered:
- (a) The realised transition only, f ∝ h_k(t)·S_k(t), as above and as your Methods say.
- (b) The correct conditional density given that the exit was to k, f ∝ h_k(t)·S_total(t), with S_total over all kept exits from the state. This has no closed form; numeric inverse CDF on 200 grid points in (a, b).

Suggested: **(a)**. It matches the Methods text and the simulation design. With low competing hazards over one interval the difference is small. (b) can be checked on one imputation if a reviewer asks.

**Check:**
- Unit test: 100,000 draws for fixed (a, b, λ, α) reproduce the truncated Weibull CDF (Kolmogorov–Smirnov p > 0.01).
- On the real data, the spread of exits over 5 imputations is > 0 for every non-rule event.
- The mean of the draws is close to the moment-matched point (sanity check, not exact).

### C4 Methods sentence

**Goal:** the Methods text describes what the code does.

**Change** (in `Manuscript/Methods_*.docx`, not in the pipeline): "we assumed that the transition time coincided with any heavy disturbance (heavy harvest or fire, dated or with an imputed date) …". Also state that exits are **drawn**, with the C3 limitation.

**Check:** text read against C2 and C3.

### C5 P1 weights: divide by survival to the first visit (new 2026-11-04, BUG FIX)

**Goal:** each candidate entry e is weighted by the probability of what was observed AFTER the first visit, given that the stand is
in the observed state s_j at W_j.

**Now:** `p1_weights()` (05a) uses S(lo - e) - S(hi - e) when a change is observed in (lo, hi], and S(hi - e) when not. The
candidates are already conditioned on being in s_j at W_j (`p1_candidates()` keeps only paths in s_j at W_j), so survival to W_j is
counted twice. This favours late entries.
The simulation has the same form (read-only folder, checked 2026-11-04): `compute_weights_2b()` in
`Sim4stateExact/p1p2_v2_test/fit_p1_rung4_mi.R` (the parent folder's copy differs only in whitespace), with candidates from
`sim_entry_candidates_2b()` (`fit_p1_rungs_2a_2b.R`, also kept only if in s_j at W_j). So the simulation's P1 validation carried
the same double counting. Whether it explains part of the simulation's alpha / coverage problems is unproven.

**Evidence (2026-11-04):**
- Toy (regen -> A -> B, exact posterior by integration): mean entry exact 27.3; current weights 37.3; divided by S(W - e) 27.4.
- Real data (theta of m01, current P1, 1,466 plots with a non-regen first state, 1,000 paths each): the current weights put the
  entry later by mean 3.8 y (median 2.9, p90 8.1, max 21.4) when a change is observed later, and by 2.9 y (2.3, 6.8, 17.2) when
  not. Mean time in the first state at the first visit 18.5 -> 22.3 y and 24.0 -> 26.9 y (about 15% too short now). By first
  state: 1.9 y (jack_pine) to 5.5 y (other_temperate, yellow_birch). This sets the truncation times of every first sojourn, hence
  alpha and the intercepts of transitions out of first states, in all imputations of the main run.

**Change:** weight = [S(lo - e) - S(hi - e)] / S(W_j - e) with a change observed, S(hi - e) / S(W_j - e) without; in log space;
S from the same hazards as the candidate paths (C1 point 6). Implemented together with C1.

**Check:** the toy example reproduces the exact posterior mean; on the real data the shift above is reproduced.

**Decided 2026-11-04 (user):** rerun the simulation with the corrected weights, in a copy of the simulation code inside this project (the simulation folder stays read-only); the reporting follows from that rerun.

---

## D. Outputs

### D1 Hazard ratios = exp(β)

**Goal:** under variant 0, a contrast c over the coefficients has log HR = c'β, with no α.

**Now:** `09_hazard_ratios.R` lines 63–64 compute `q = α·c'β`, with the delta method over (β, α). `sensitivity/11` lines 213–215 use `exp(al * …)`.

**Change:**
1. 09:
   - `q = c'β`, `u = c'V_ββ c`, with V_ββ the β block of the stored joint covariance;
   - drop the α gradient;
   - drop `alpha_cov_ratio_med` and `alpha_cov_low_share` from the hazard-ratio output (they only mattered because α entered the HRs), or move them to a separate alpha diagnostics table.
2. 11: `exp(β)` and its interval directly.
3. Update the header comments of 09 and 11.

**Check:** 09's printed check becomes "pooled log HR per SD CMI == pooled β of `z_cmi_gs`". It should be identical up to pooling arithmetic.

### D2 TSD effects at fixed times, within the observed range

**Goal:** report the disturbance effect b + c·log(τ + 1) at τ = 1, 10, 30 and 60 y (at GDD z = −1, 0, 1) for every transition and disturbance level (decided 2026-11-04; was: only where the data have that TSD).

**Now:** `TAU_GRID <- c(0, 1, 3)` in decades; there is a `_tsd_per_decade` contrast.

**Change:**
1. 09:
   - `TAU_GRID_Y <- c(1, 10, 30, 60)`;
   - the contrast vector for τ is c = (1, z, f(τ)) on (main, `gdd_x`, tsd), with f from `TSD_FORM`;
   - `DIST` maps `harv_light → tsd_harv_light` and `harv_heavy → tsd_harv_heavy`; drop `harv_any`;
   - replace `_tsd_per_decade` with `_tsd_1_vs_10y` and `_tsd_10_vs_30y` (differences between fixed times; meaningful for log TSD);
   - (2026-11-04) add `_tsd_coef`: the TSD coefficient c itself (per unit of log(TSD + 1)), input to 12's block D (D5).
2. (Decided 2026-11-04) **No range filter**: every transition and level is reported at all four τ, for consistency. The at-risk range (min, 5th and 95th percentile, max of `tsd_*_y` over the at-risk flagged intervals of the origin state) and the event range are reported as extra columns. Main build: every fire, pest and heavy-harvest pair reaches 60 y; light harvest falls short in 45 transitions (max 47-59 y), where 60 y is an extrapolation along the log curve: caveat in the manuscript (PIPELINE_NOTES, E). τ = 1 y is below the 5th percentile in all 179 pairs.
3. **Decided 2026-11-04 (user): keep both** `_tsd_1_vs_10y` and `_tsd_10_vs_30y` (they show the size of the change early and later). Note: they are both fixed multiples of the TSD coefficient (log(11/2) and log(31/11) times it), so they share one test (same sign and p-value); state this so they are not read as independent evidence.

**Check:** every transition × level has rows for the four τ; the range columns are filled.

### D3 Alpha reporting

**Goal:** alpha is reported for every transition. Where the state starts at the disturbance, it is shown with the TSD effect, since part of the duration dependence moves into the TSD term there.

**Now:** 08 pools alpha, also on the log scale; 12 synthesises log alpha.

**Change:** 09 writes `alpha_report_<tag>.csv`, one row per transition, with:
- pooled alpha and its CI (log-scale pooling);
- events;
- class (regen_exit / to_regen / between);
- flag `starts_at_disturbance` = class == "regen_exit";
- the pooled TSD coefficients present in that transition;
- the cap flags from B6.

**Decided 2026-11-04:** the transitions that count as "state starts at the disturbance" are regen→X, where regen begins at the stand-replacing event; the cap test showed alpha shifts mainly there and in jack_pine→regen.

**Check:** 90 rows; alpha equals `pooled_<tag>_final.csv`.

### D4 One-sided bounds

**Goal:** a level with 0 events in a contrast gets a one-sided statement, because only one side of its interval comes from the data:
- regen→sugar_maple: 0 events without harvest, so harvest HR > X;
- other_temperate→yellow_birch: 0 heavy-harvest events, so heavy HR < X.

**Change:** 09 adds a column `bound`:
- if the reference has 0 events: report the lower one-sided 95% limit exp(log HR − t₀.₉₅,df·se);
- if the level has 0 events: report the upper limit exp(log HR + t₀.₉₅,df·se);
- the two-sided CI for these rows is set to NA in the summary tables.

Event counts come from `transition_set` (`n_harv0`, `n_harv1`, `n_harv2`).

Additions 2026-11-04:
4. `make_contrasts()` must also build the disturbance-vs-none contrast for a level that has ONLY its main term (A6: heavy harvest in other_temperate→yellow_birch); today it skips any level without main, GDD and TSD terms together, so that level would get no row to bound. The contrast is then the main term alone.
5. regen→sugar_maple has no GDD main term (A11): its harvest-vs-none bound is reported at mean GDD (z = 0) only.

**Check:** exactly these two transitions are flagged on the main build. 09 prints them.

### D5 Script 12 (cross-transition synthesis)

**Goal:** the synthesis uses the new terms.

**Now:**
- `DIST` includes `harv_any`;
- line 84 de-duplicates harvest TSD ("light/heavy/any share one TSD term");
- the robust subset uses `alpha_cov_low_share` (line 69).

**Change:**
1. Remove `harv_any`.
2. TSD block D1: one estimate per transition × type × level (no de-duplication). **Decided 2026-11-04:** the estimate is the TSD coefficient c (`_tsd_coef`, with its pooled SE). Every fixed-time contrast is c times a constant, so the meta-analysis gives the same tests whichever is used; the pooled results are reported scaled to windows in the text and tables: 1 vs 10 y = c × log(11/2) ≈ 1.70c, 10 vs 30 y = c × log(31/11) ≈ 1.04c, 30 vs 60 y = c × log(61/31) ≈ 0.68c (intervals scale the same way). Not both contrasts in one model: each transition would enter twice with the same information.
3. `robust := n_events >= MIN_EVENTS_ROBUST & !low_events` (no `alpha_cov`).
4. Header comments.
5. (Decided 2026-11-04) Disturbance block B (B1, B2): the old selection `"_vs_none z=0 tau=0$"` matches nothing under the new τ grid. Run block B once at each τ = 1, 10, 30 and 60 y (z = 0).
6. (2026-11-04) Transitions without a GDD main term (A11; main build: regen→sugar_maple) are left out of block A for GDD (CMI and soil stay), out of block C, and out of block B (harvest vs none is only a one-sided bound, D4).

The shade-tolerance values (v5 R4) are a separate, later item.

**Check:** 12 runs on the new 09 output, and every block has k > 0.

### D6 `tests/08f` (supplement diagnostic D2: are replacements neutral?)

**Goal:** compare failed and valid imputations under the final model.

**Now:** a variant-1 maximum-likelihood fit (`fit_mle`); old term set; rows rebuilt from `extras$p1` and `extras$p2`.

**Change:**
1. Use `fit_map_v0()` and `H_at()` from 05a (B6) instead of `fit_mle`.
2. Terms from the new mask.
3. The rebuild stays as is, because `extras$p2` holds the final exits (rule or draw).
4. Add the B6 per-fit columns to the comparison.
5. (2026-11-04) After A1, remove the dropping of rows with NA climate (`complete.cases`); there are none.

**Check:** the rebuild check against `baddata_*` still matches 100%.

### D7 `sensitivity/10` and `11`

**Goal:** the sensitivity scripts match the final model.

**Change:**
- **10:**
  - S1 uses `n_dist` per level (provided by 09);
  - the threshold is read from config (`MIN_EVENTS_S1 = 20`);
  - S2 and S3 unchanged.
- **11:**
  - the inla call gets variant 0 (B1) and `control.fixed` (B2);
  - intercept without `− LOG_TIME_SCALE` (B3, line 117);
  - HRs as `exp(β)` (D1, lines 213–215);
  - (2026-11-04) transitions without a GDD main term (A11): drop both `gdd_btw` and `gdd_wth` (the CMI split stays).
- **10** (2026-11-04): remove `alpha_cov_low_share` from the output columns (line 47); 09 no longer produces it (D1).

**Check:** 11's step 0 (reproducing `mi_raw` for `CHECK_M`) passes under the new run tag.

### D8 09 for transitions without a GDD main term (new 2026-11-04)

**Goal:** 09 reports only quantities the model contains.

**Change:** for a transition without `z_gdd_3` (A11): `<d>_gdd_per_sd` = the `gdd_x_<d>` coefficient alone (GDD slope within that
level); no `per_sd_gdd_3` and no `<d>_gdd_interaction` rows; `<d>_vs_none` at z = 0 only (D4 point 5).

**Check:** regen->sugar_maple has `harv_light_gdd_per_sd` and `harv_heavy_gdd_per_sd` equal to exp(pooled `gdd_x_harv_*`), and
none of the excluded rows.

---

## Implementation order and run checks

1. **A1** (03c, 03z): rerun 03c, 03z and 04, then 04x. Check the A1 items.
2. **A2–A11 and B7** (04, `covariate_spec.R`, `config.R`, 04b): rerun 04 and `tests/04c`. Compare transitions, terms and events with the current main build: only the harvest terms (A2, A4, A6) and the GDD main term of regen→sugar_maple (A11) may change; transitions, events and the fire and pest terms must be identical (wording revised 2026-11-04).
3. **B1–B6, C1–C3, C5** (05a, 05, 07): run unit tests (B4, C1, C3, C5), then a smoke test of 05 with `SUBSAMPLE_N = 500` (checks that the code runs, not the results). Then one full imputation, comparing against `cap_test_new_rows_m01_nocut`. Expect differences, since C1, C3 and C5 change P1 and P2. Check B6, the C2 counts and the P1 fallbacks (C1).
4. Full run under a **new run tag** (e.g. `final`), so `main` stays intact: 05 → 07_redraw → 08.
5. **D1–D8**: 09 → 12 → 10 → 11.

## What to expect in the first run (2026-11-04)

- **The cap list may grow.** C5 puts first entries earlier (time in the first state at the first visit about 15% longer on average),
  and H at truncation scales with that time to the power alpha. Transitions just below the cap (22.85) in the cap test (m01, no cut):

  | Transition | max H(tr), reference fit | max H(tr), INLA | alpha |
  |---|---|---|---|
  | sugar_maple->regen | 20.6 | 22.0 | 0.96 |
  | fir->regen | 18.2 | 18.2 | 1.60 |
  | jack_pine->regen | 16.1 | 16.0 | 1.42 |

  Scaled by 1.15^alpha: about 23.5, 22.7 and 19.6. Rough only: the mean shift does not apply directly to the largest-H row, and the
  refit re-estimates alpha and the intercept. The caveat list comes from B6 after the run.
- **Transitions with a 0-event level** (regen->sugar_maple, other_temperate->yellow_birch): check their failure counts in 07
  (variant 0 failed more often with a 0-event reference in sim01; figure not re-checked).
- **P1 fallbacks by first state**, in particular plots starting in sugar maple without harvest (C1 check).

## Outside this plan

- S21, S22 and S24 are planned as refits on the stored imputations, but no refit script exists (already on the v5 "missing scripts"
  list); the config switches only allow full new runs. Not needed for the main run; needed before these sensitivities.

## Decisions (all taken 2026-11-04)

| Item | Decision |
|---|---|
| A6 | A level with 0 events: main term only (0 vs >= 1 event is a real boundary; no other threshold) |
| A10 | TSD per row: end-visit value |
| A11 (new) | Drop the GDD main term where a transition has 0 undisturbed events (main build: regen->sugar_maple) |
| C1 | Flags switch on at their onset in P1's history; TSD yearly, at the midpoint of each year; heavy disturbances more recent than the stand age get no special rule |
| C2 | Heavy harvest and fire in one interval: the earlier date |
| C3 | Draw with the realised transition only; uniform draw where there is no usable fit |
| C5 (new) | P1 weights divided by S(W_j - e) (bug fix, also present in the simulation) |
| D2 | All four τ (1, 10, 30, 60 y) for every transition; no range filter; caveat for light harvest at 60 y |
| D3 | "Starts at the disturbance" = regen->X |
| D5 | Block B at τ = 1, 10, 30, 60 y; block D uses the TSD coefficient c, results scaled to windows in the text; A11 transitions out of block A (GDD), B and C |
| B3 | Intercept in years computed in 08 |
| B7 | `LOW_EVENTS_MAX <- MIN_EVENTS + 4L`; `recount_events()` called in `run_stage1()` |
| A9 | Also `D_heavy` = the earlier of the two dates |

D2 point 3: keep both TSD-difference contrasts (decided 2026-11-04). C5: rerun the simulation with the corrected weights (decided 2026-11-04).

## Limitations to state (additions 2026-11-04; also in PIPELINE_NOTES, E)

- P1 represents only the most recent disturbance of each type before the first visit.
- P1's hazards change within a sojourn (onset, yearly TSD), which departs from the simulation-validated P1 (static covariates).
- P1 uses real-time TSD, the fitted model the end-visit TSD of each interval (A10). For fading effects (negative TSD slope), P1's
  hazards just after a disturbance are higher than the fitted model implies.
- Fire in the X->regen exit rule includes partial fires.
- regen->sugar_maple: no GDD effect in unharvested plots (0 events); its harvest effects hold at mean GDD only.
- Light harvest at 60 y is extrapolated in 45 transitions.
- P1 weighting corrected relative to the simulation-derived code (C5).
