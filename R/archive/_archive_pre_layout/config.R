# ==============================================================================
# R/config.R -- run configuration (T0, 2026-10-24). Sourced, not run.
# Sourced by: 04_prepare_inla.R, 05_fit_inla.R, 07_redraw.R, 08_final_pool.R,
#   09_hazard_ratios.R, 12_cross_transition_synthesis.R, R/sensitivity/10, 11,
#   R/tests/04c_check_rows.R.
# RUN_TAG names the run: 04 writes its tables to TABLE_DIR, 05-12 write fits and
# results to FIT_DIR. Run tags are documented in docs/PIPELINE_NOTES.md ("Run tags").
# Later switches for the full-run sensitivities (TSD_SOURCE, prior strength, ...)
# will be added here.
# ==============================================================================
RUN_TAG   <- "main"
TABLE_DIR <- here::here("data", "processed", paste0("tables_", RUN_TAG))   # 04 output: visit_intervals, plot_info, transition_set, wj_table, tsd_excluded_plots
FIT_DIR   <- here::here("data", "processed", paste0("fit_", RUN_TAG))      # 05-12 output
INT_TABLE_DIR <- here::here("data", "processed", "tables_int")           # backup of the int-run 04 tables (diagnostics 13-21)
