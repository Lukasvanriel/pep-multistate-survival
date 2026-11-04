# ==============================================================================
# R/config.R -- run configuration and project paths. Sourced, not run.
# T0 2026-10-24; folder layout reorganised 2026-10-24 (see docs/PIPELINE_NOTES.md "Folder layout").
# Rule: everything up to the INLA input lives in data/; everything produced by INLA and
#       used afterwards lives in data-output/<RUN_TAG>/.
# Sourced by: 02, 03, 03b, 03c, 03z, 04, 05, 07_redraw, 08, 09, 12, R/sensitivity/10, 11, R/tests/04c.
# Switches for the full-run sensitivities (TSD_SOURCE, prior strength, TSD form, MIN_EVENTS, ...) are in this file; see the end.
# ==============================================================================
RUN_TAG <- "final"  # run tags: noint (archived), inttest, main (previous main run: keep data/tables/main and data-output/main untouched), final (2026-11-04: final pipeline) -- see PIPELINE_NOTES "Run tags"

# ---- prep side (data/) --------------------------------------------------------
# data/raw/ = inputs, never written by the pipeline (except the 03c download scratch data/raw/climate_raw/)
BUILD_DIR     <- here::here("data", "build")                  # 02-03z outputs: state, disturbance, soil, climate, combined covariates
BUILD_SUM_DIR <- file.path(BUILD_DIR, "summaries")            # 03z summary CSVs + figure
TABLE_DIR     <- here::here("data", "tables", RUN_TAG)        # 04 output = INLA input (+ transition_set.csv)
DIAG_DATA_DIR <- here::here("data", "diagnostics")            # prep-side diagnostics (e.g. code export of both layers)

# ---- INLA side (data-output/) ---------------------------------------------------
RUN_OUT_DIR   <- here::here("data-output", RUN_TAG)
FIT_DIR       <- file.path(RUN_OUT_DIR, "fits")               # 05 + 07_redraw: stage1/2, scaling, extras_m*, mi_raw, mi_diag, screen_ref
RESULTS_DIR   <- file.path(RUN_OUT_DIR, "results")            # 08, 09, 12: fit_validity, pooled_final(+cov), cov_check, hazard ratios, synthesis
SENS_DIR      <- file.path(RUN_OUT_DIR, "sensitivity")        # R/sensitivity/10, 11 and later S-refits
FIG_DIR       <- file.path(RUN_OUT_DIR, "figures")
DIAG_OUT_DIR  <- here::here("data-output", "diagnostics")     # INLA-side diagnostics not tied to one run (H(tr) simulation)

# ---- earlier run used by the diagnostics R/tests/13-21 ---------------------------
INTTEST_TABLE_DIR <- here::here("data", "tables", "inttest")
INTTEST_FIT_DIR   <- here::here("data-output", "inttest", "fits")   # files inside keep their "_int" suffix

# ---- thresholds shared by 02 and 03 (2026-10-26) --------------------------------------------------------------
BA_REGEN       <- 5      # m2/ha live basal area: below = regeneration state (02); also the "drop below regen" BA event (03)
BA_LOSS_THRESH <- 0.25   # BA loss since the previous visit that counts as a BA event (03 steps 5, E0)

# ---- disturbance TSD settings (2026-10-26) -------------------------------------------------------------------
TSD_SOURCE      <- "mid"   # TSD of undated events in 04: "mid" (main), "lower" / "upper" (sensitivity S11)
FV_AGE_QUANTILE <- 0.5     # undated PARTIAL harvest/pest code at a plot's first visit: age = this quantile of the ages of
                           # DATED first-visit partial codes of that type (0.5 main; 0.1 / 0.9 = sensitivity S19)

# ---- station_pe disturbance codes (2026-10-26, user: full merge) ----------------------------------------------
USE_STATION_CODES <- TRUE  # TRUE (main): station_pe origine/perturb codes merged with pee_ori_sond in 03; FALSE = sensitivity S20
                           # (pee_ori_sond only). Merge rules: header of R/03a_disturbance_events.R (merge_station_records).

# ---- INLA fitting (2026-10-26, T5/T6) -------------------------------------------------------------------------
ALPHA_PRIOR <- 0.1   # PC prior strength on the Weibull shape alpha, pc.alphaw(ALPHA_PRIOR), in ALL fits (stage 1, stage 2, every
                     # imputation, redraw, script 11). 5 = INLA's default (checked 2026-10-26: identical fits) = sensitivity S4.
                     # Fits use one thread each (num.threads "1:1", T6).

# ---- final pipeline switches (WP2 / B7, 2026-11-04; spec: docs/PIPELINE_FIX_PLAN_2026-11-03.md) ----------------------------
# Sensitivity runs change only these lines plus RUN_TAG.
WEIBULL_VARIANT  <- 0L        # Weibull variant 0: S(t) = exp(-e^eta t^alpha); beta = log hazard ratio (used from WP3, B1)
COEF_PRIOR_SD    <- 5         # SD of the N(0, sd^2) prior on the coefficients (S21 = 2.5; S22 = 1/sqrt(0.001)); used from WP3 (B2)
TSD_FORM         <- "log1p"   # TSD model columns: "log1p" = log(TSD in years + 1) (main); "linear" = TSD in decades (S24)
MIN_EVENTS       <- 15L       # transition kept with >= this many state changes; also the per-term minimum for fire / pest / harvest
                              # (light + heavy together) (S17 = 20, S23 = 10). Read by 04, 04b recount_events(), tests/04c.
XREGEN_EXIT_RULE <- TRUE      # X->regen exit placed at the heavy disturbance (C2); FALSE only as a check; used from WP5
