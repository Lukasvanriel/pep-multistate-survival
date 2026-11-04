# ==============================================================================
# R/config.R -- run configuration and project paths. Sourced, not run.
# T0 2026-10-24; folder layout reorganised 2026-10-24 (see docs/PIPELINE_NOTES.md "Folder layout").
# Rule: everything up to the INLA input lives in data/; everything produced by INLA and
#       used afterwards lives in data-output/<RUN_TAG>/.
# Sourced by: 02, 03, 03b, 03c, 03z, 04, 05, 07_redraw, 08, 09, 12, R/sensitivity/10, 11, R/tests/04c.
# Later switches for the full-run sensitivities (TSD_SOURCE, prior strength, ...) will be added here.
# ==============================================================================
RUN_TAG <- "main"   # run tags: noint (archived), inttest, main -- see PIPELINE_NOTES "Run tags"

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
