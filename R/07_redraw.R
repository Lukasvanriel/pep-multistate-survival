# ==============================================================================
# 07_redraw.R   RUN ON THE MACHINE WITH INLA (after 05_fit_inla.R; resumable; batches of N_WORKERS)
# Replaces 07_redraw_failed.R + 07b_bsregen_sweep.R + 07c_bsregen_cov.R (those were specific to the "int" run and are kept for that archive).
#
# PURPOSE: give every transition TARGET (30) valid fits, with ONE rule for all transitions.
# FROZEN SCREEN (same rules as before; the reference values are computed ONCE from this run's m1-30 and saved before any candidate is fitted):
#   valid if  (1) INLA returned no error and est/sd are not NA
#             (2) alpha-SD ratio  sd_alpha / (0.78 * alpha / sqrt(events)) >= RATIO_MIN (0.1)   [collapsed alpha SD ~ 2e-4]
#             (3) |intercept - ref| / ref_sd <= Z_MAX (4); ref, ref_sd = median intercept and median intercept SD of the transition's
#                 ratio-passing fits among m1-30                                                  [intercept stuck / shifted]
# CANDIDATES: fresh imputations m = 31, 32, ... with the SAME code and seed rule as the MI loop (set.seed(SEED * 1000 + m), the same stage-2
#   models). Each candidate imputation is built in full (theta draw for all transitions, P1, P2: identical to what a full run would build), but
#   only the transitions that still lack valid fits are FITTED (only_k; transitions are fitted independently, so the fit of a transition is the
#   same INLA call either way). Full outputs (extras with joint covariance, diagnostics) are written for the fitted transitions.
# STOP: when every transition has >= TARGET valid fits, or MAX_CAND candidates have been tried (then 08 stops and reports).
# Outputs (FIT_DIR = data-output/<tag>/fits/): screen_ref_<tag>.csv (frozen references), mi_raw_redraw_<tag>.csv, mi_diag_redraw_<tag>.csv,
#   extras_<tag>_mXX.rds / baddata_<tag>_mXX.rds for candidates, redraw_progress_<tag>.csv (valid counts after every batch).
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(INLA) })
source(here::here("R", "config.R"))   # RUN_TAG, TABLE_DIR, FIT_DIR, ALPHA_PRIOR
INLA::inla.setOption(num.threads = "1:1")   # T6: one thread per INLA call (also set in every fit call)
SEED <- 2026L; N_WORKERS <- 6L; M_MI <- 30L; SIM_M <- 200L; SD_MIN <- 1e-3; N_COV_SAMPLES <- 1000L
TARGET <- 30L; MAX_CAND <- 500L; RATIO_MIN <- 0.1; Z_MAX <- 4
OUT_DIR <- FIT_DIR; options(width = 200)
source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED); fit_fun <- fit_inla_models_pc   # T5 (2026-10-26): same fit function and prior as 05
cat(sprintf("INLA fits: Weibull variant %d, N(0, %g^2) on coefficients, alpha prior pc.alphaw(%s), one thread per fit\n", WEIBULL_VARIANT, COEF_PRIOR_SD, format(ALPHA_PRIOR)))
inp <- prepare_fit_inputs(TABLE_DIR, subsample_n = NULL, seed = SEED)
s2 <- readRDS(file.path(OUT_DIR, "stage2.rds"))

screen_fits <- function(raw, ref) {                    # one row per (k, m): validity under the frozen screen
  a <- raw[parameter == "alpha", .(k, label, m, alpha = est, a_sd = sd, n_events, fit_error)]
  i <- raw[parameter == "(Intercept)", .(k, m, icpt = est, i_sd = sd)]
  x <- merge(a, i, by = c("k", "m"), all.x = TRUE); x[, ratio := a_sd / (0.78 * alpha / sqrt(n_events))]
  x <- merge(x, ref, by = "k", all.x = TRUE); x[, z_icpt := (icpt - icpt_ref) / icpt_sd_ref]
  x[, reason := fifelse(fit_error %in% TRUE | is.na(alpha) | is.na(a_sd) | is.na(icpt), "error/NA",
                fifelse(ratio < RATIO_MIN, "alpha_sd_ratio", fifelse(is.na(z_icpt) | abs(z_icpt) > Z_MAX, "intercept", "ok")))]
  x[, valid := reason == "ok"][order(k, m)] }

raw0 <- fread(file.path(OUT_DIR, sprintf("mi_raw_%s.csv", RUN_TAG)))
stopifnot(setequal(unique(raw0$m), seq_len(M_MI)))
ref_path <- file.path(OUT_DIR, sprintf("screen_ref_%s.csv", RUN_TAG))
if (file.exists(ref_path)) { ref <- fread(ref_path); cat("Screen references: loaded (frozen)\n") } else {
  a0 <- raw0[parameter == "alpha", .(k, m, ratio = sd / (0.78 * est / sqrt(n_events)), err = fit_error %in% TRUE | is.na(est) | is.na(sd))]
  ok0 <- a0[err == FALSE & !is.na(ratio) & ratio >= RATIO_MIN, .(k, m)]
  ref <- raw0[parameter == "(Intercept)"][ok0, on = .(k, m)][, .(icpt_ref = median(est), icpt_sd_ref = median(sd), n_ref = .N), by = k]
  if (nrow(ref) < nrow(inp$td)) stop("No ratio-passing fit among m1-", M_MI, " for: ", paste(setdiff(inp$td$label, inp$td$label[ref$k]), collapse = ", "), " -- decide before continuing.")
  fwrite(ref, ref_path); cat("Screen references computed from m1-", M_MI, " and FROZEN in ", basename(ref_path), "\n", sep = "")
}
ref_tab <- merge(ref, inp$td[, .(k, label)], by = "k")[order(k)]
cat("Screen reference table (intercept on INLA's decade scale, variant 0), printed before any candidate is fitted:\n"); print(ref_tab, nrows = nrow(ref_tab))
cat("Transitions whose reference rests on < 10 fits:\n"); print(ref_tab[n_ref < 10])

raw_path  <- file.path(OUT_DIR, sprintf("mi_raw_redraw_%s.csv", RUN_TAG)); diag_path <- file.path(OUT_DIR, sprintf("mi_diag_redraw_%s.csv", RUN_TAG))
prog_path <- file.path(OUT_DIR, sprintf("redraw_progress_%s.csv", RUN_TAG))
state_path <- file.path(OUT_DIR, sprintf("p1_fallback_by_state_redraw_%s.csv", RUN_TAG))
all_raw <- function() if (file.exists(raw_path)) rbind(raw0, fread(raw_path), fill = TRUE) else raw0
err_path  <- file.path(OUT_DIR, sprintf("redraw_runerrors_%s.csv", RUN_TAG))   # candidates whose whole imputation failed to run (counted as tried)
tried   <- function() c(if (file.exists(diag_path)) fread(diag_path)$m else integer(0), if (file.exists(err_path)) fread(err_path)$m else integer(0))

t0 <- proc.time()["elapsed"]
repeat {
  v <- screen_fits(all_raw(), ref)
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
    cat(sprintf("  m=%d: P1 fallbacks %d of %d plots; by first state (n_fallback / n_plots):", r$diag$m, sum(r$p1_state$n_fallback), sum(r$p1_state$n_plots)),
        paste(sprintf("%s %d/%d", r$p1_state$first_state, r$p1_state$n_fallback, r$p1_state$n_plots), collapse = "; "), "\n") }
  cat(sprintf("  done (%.0f s elapsed)\n", proc.time()["elapsed"] - t0))
}
v <- screen_fits(all_raw(), ref)
cat("\nInvalid fits per transition and reason (m1-30 and candidates):\n")
print(dcast(v[valid == FALSE, .N, by = .(label, reason)], label ~ reason, value.var = "N", fill = 0))
cat("\nFailure rate in m1-", M_MI, " (all 91 transitions fitted there): ", round(v[m <= M_MI, mean(!valid)], 4), "\n", sep = "")
