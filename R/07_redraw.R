# ==============================================================================
# 07_redraw.R   RUN ON THE MACHINE WITH INLA (after 05_fit_inla.R; resumable; batches of N_WORKERS)
# Replaces 07_redraw_failed.R + 07b_bsregen_sweep.R + 07c_bsregen_cov.R (those were specific to the "int" run and are kept for that archive).
#
# PURPOSE: give every transition TARGET (30) valid fits, with ONE rule for all transitions.
# SCREEN (user decision 2026-11-05; same rule as 08_final_pool.R and sensitivity/11; replaces the frozen intercept-reference rule, z <= 4):
#   valid if  (1) INLA returned no error and est/sd of alpha and intercept are not NA
#             (2) alpha-SD ratio  sd_alpha / (0.78 * alpha / sqrt(events)) >= RATIO_MIN (0.1)   [collapsed alpha SD ~ 2e-4]
#             (3') |z_map_alpha| <= Z_MAP_MAX (3) and |z_map_intercept| <= Z_MAP_MAX (3): INLA mode vs. direct maximisation of the same
#                  log-posterior (fit_map_v0, B6 columns of mi_raw / mi_raw_redraw)               [non-converged fits]
#   No reference table. Resume: the screen is applied to ALL stored fits (m1-30 and earlier candidates) at every start, so a restart after a change of
#   the screen picks up exactly the transitions that are short under the new rule (candidates are numbered on from the largest m tried).
# CANDIDATES: fresh imputations m = 31, 32, ... with the SAME code and seed rule as the MI loop (set.seed(SEED * 1000 + m), the same stage-2
#   models). Each candidate imputation is built in full (theta draw for all transitions, P1, P2: identical to what a full run would build), but
#   only the transitions that still lack valid fits are FITTED (only_k; transitions are fitted independently, so the fit of a transition is the
#   same INLA call either way). Full outputs (extras with joint covariance, diagnostics) are written for the fitted transitions.
# STOP: when every transition has >= TARGET valid fits, or MAX_CAND candidates have been tried (then 08 stops and reports).
# Outputs (FIT_DIR = data-output/<tag>/fits/): mi_raw_redraw_<tag>.csv, mi_diag_redraw_<tag>.csv,
#   extras_<tag>_mXX.rds / baddata_<tag>_mXX.rds for candidates, redraw_progress_<tag>.csv (valid counts after every batch).
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(INLA) })
source(here::here("R", "config.R"))   # RUN_TAG, TABLE_DIR, FIT_DIR, ALPHA_PRIOR
INLA::inla.setOption(num.threads = "1:1")   # T6: one thread per INLA call (also set in every fit call)
SEED <- 2026L; N_WORKERS <- 6L; M_MI <- 30L; SIM_M <- 200L; SD_MIN <- 1e-3; N_COV_SAMPLES <- 1000L
TARGET <- 30L; MAX_CAND <- 500L; RATIO_MIN <- 0.1; Z_MAP_MAX <- 3
OUT_DIR <- FIT_DIR; options(width = 200)
source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
# WP5 review: 07 only runs on fits made by the same code as 05 used (md5 of the files 05 sources, stored in run_id.rds by 05); a difference stops
check_resume_guard(OUT_DIR, run_id_now(RUN_TAG, NA_integer_, NULL, here::here()), show_dir = sub(paste0("^", here::here(), "/"), "", OUT_DIR), only = c("code_md5", "RUN_TAG"), write = FALSE)
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED); fit_fun <- fit_inla_models_pc   # T5 (2026-10-26): same fit function and prior as 05
cat(sprintf("INLA fits: Weibull variant %d, N(0, %g^2) on coefficients, alpha prior pc.alphaw(%s), one thread per fit\n", WEIBULL_VARIANT, COEF_PRIOR_SD, format(ALPHA_PRIOR)))
inp <- prepare_fit_inputs(TABLE_DIR, subsample_n = NULL, seed = SEED)
s2 <- readRDS(file.path(OUT_DIR, "stage2.rds"))

screen_fits <- function(raw) {                         # one row per (k, m): validity under the screen
  stopifnot("mi_raw lacks z_map_alpha / z_map_intercept (B6 columns)" = all(c("z_map_alpha", "z_map_intercept") %in% names(raw)))
  a <- raw[parameter == "alpha", .(k, label, m, alpha = est, a_sd = sd, n_events, fit_error, z_map_alpha, z_map_intercept)]
  i <- raw[parameter == "(Intercept)", .(k, m, icpt = est, i_sd = sd)]
  x <- merge(a, i, by = c("k", "m"), all.x = TRUE); x[, ratio := a_sd / (0.78 * alpha / sqrt(n_events))]
  x[, reason := fifelse(fit_error %in% TRUE | is.na(alpha) | is.na(a_sd) | is.na(icpt), "error/NA",
                fifelse(ratio < RATIO_MIN, "alpha_sd_ratio",
                fifelse(is.na(z_map_alpha) | is.na(z_map_intercept) | abs(z_map_alpha) > Z_MAP_MAX | abs(z_map_intercept) > Z_MAP_MAX, "z_map", "ok")))]
  x[, valid := reason == "ok"][order(k, m)] }

raw0 <- fread(file.path(OUT_DIR, sprintf("mi_raw_%s.csv", RUN_TAG)))
stopifnot(setequal(unique(raw0$m), seq_len(M_MI)))
raw_path  <- file.path(OUT_DIR, sprintf("mi_raw_redraw_%s.csv", RUN_TAG)); diag_path <- file.path(OUT_DIR, sprintf("mi_diag_redraw_%s.csv", RUN_TAG))
prog_path <- file.path(OUT_DIR, sprintf("redraw_progress_%s.csv", RUN_TAG))
state_path <- file.path(OUT_DIR, sprintf("p1_fallback_by_state_redraw_%s.csv", RUN_TAG))
p2src_path <- file.path(OUT_DIR, sprintf("p2_source_redraw_%s.csv", RUN_TAG))
all_raw <- function() if (file.exists(raw_path)) rbind(raw0, fread(raw_path), fill = TRUE) else raw0
err_path  <- file.path(OUT_DIR, sprintf("redraw_runerrors_%s.csv", RUN_TAG))   # candidates whose whole imputation failed to run (counted as tried)
tried   <- function() c(if (file.exists(diag_path)) fread(diag_path)$m else integer(0), if (file.exists(err_path)) fread(err_path)$m else integer(0))

t0 <- proc.time()["elapsed"]
repeat {
  v <- screen_fits(all_raw())
  cnt <- merge(inp$td[, .(k, label)], v[valid == TRUE, .(n_valid = .N), by = k], by = "k", all.x = TRUE)[is.na(n_valid), n_valid := 0L]
  need <- cnt[n_valid < TARGET]
  fwrite(data.table(time = format(Sys.time()), candidates = length(tried()), need[, .(label, n_valid)]), prog_path, append = file.exists(prog_path))
  if (!nrow(need)) { cat("All transitions have >=", TARGET, "valid fits.\n"); break }
  if (length(tried()) >= MAX_CAND) { cat("!! MAX_CAND reached; still short:\n"); print(need); break }
  b0 <- max(c(M_MI, tried())) + 1L; batch <- b0:(b0 + N_WORKERS - 1L)
  cat(sprintf("Batch m=%d-%d, fitting %d transition(s): %s\n", min(batch), max(batch), nrow(need), paste(sprintf("%s (%d)", need$label, need$n_valid), collapse = ", ")))
  out <- mclapply(batch, function(m) tryCatch({ set.seed(SEED * 1000L + m)
      run_one_imputation(m, inp, s2$models, s2$theta, fit_fun, M = SIM_M, keep_models = FALSE, only_k = need$k) },
    error = function(e) list(error = conditionMessage(e), m = m)), mc.cores = N_WORKERS)
  for (r in out) {
    if (is.null(r$res)) { cat(sprintf("  m=%s FAILED to run: %s\n", r$m, r$error))
                          fwrite(data.table(m = r$m, error = substr(as.character(r$error), 1, 300)), err_path, append = file.exists(err_path)); next }
    if (!is.null(r$bad_data)) saveRDS(r$bad_data, file.path(OUT_DIR, sprintf("baddata_%s_m%02d.rds", RUN_TAG, r$diag$m)))
    saveRDS(r$extras, file.path(OUT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, r$diag$m)))
    fwrite(r$res, raw_path, append = file.exists(raw_path)); fwrite(r$diag, diag_path, append = file.exists(diag_path))
    fwrite(r$p1_state, state_path, append = file.exists(state_path))     # WP4: P1 fallbacks by first state, one block per imputation
    fwrite(r$p2_src, p2src_path, append = file.exists(p2src_path))      # WP5: exit_source counts per transition
    cat(sprintf("  m=%d: P2 exits by source: %s\n", r$diag$m, p2_source_line(r$p2_src)))
    cat(sprintf("  m=%d: P1 max jumps of any path: %d (plot %s); limit %d\n", r$diag$m, r$diag$p1_max_jumps, r$diag$p1_max_jumps_plot, P1_MAX_JUMPS))
    cat(sprintf("  m=%d: P1 fallbacks %d of %d plots; by first state (n_fallback / n_plots):", r$diag$m, sum(r$p1_state$n_fallback), sum(r$p1_state$n_plots)),
        paste(sprintf("%s %d/%d", r$p1_state$first_state, r$p1_state$n_fallback, r$p1_state$n_plots), collapse = "; "), "\n") }
  cat(sprintf("  done (%.0f s elapsed)\n", proc.time()["elapsed"] - t0))
}
v <- screen_fits(all_raw())
cat("\nInvalid fits per transition and reason (m1-30 and candidates):\n"); inv <- v[valid == FALSE, .N, by = .(label, reason)]
if (nrow(inv)) print(dcast(inv, label ~ reason, value.var = "N", fill = 0)) else cat("  none\n")
cat("\nFailure rate in m1-", M_MI, " (all 91 transitions fitted there): ", round(v[m <= M_MI, mean(!valid)], 4), "\n", sep = "")
