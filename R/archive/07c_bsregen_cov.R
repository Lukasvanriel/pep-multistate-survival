# ARCHIVE (2026-10-22): used for the 'int' run only. For the 'gs' run and later use 07_redraw.R.
# ==============================================================================
# 07c_bsregen_cov.R   RUN ON THE MACHINE WITH INLA (14 full imputations: ~3 batches of the 28-min full run; resumable)
# Purpose: the black_spruce->regen fits from the sweep (07b, m >= 115) were fitted for that one transition only, so they lack everything the regular pipeline
#   stores per imputation. For the 14 VALID sweep imputations this script runs the REGULAR run_one_imputation() (same function and seed rule as 05/07): all 91 transitions,
#   mi_raw rows, diagnostics, extras (theta draws, P1 entries, P2 exits, per-fit summaries/marginals/traces, joint covariance from N_COV_SAMPLES samples), baddata of failed fits.
#   No new imputations are drawn; no pooled table other than the black_spruce->regen covariance is changed (the other 90 transitions keep their existing valid fits;
#   the new fits of those 90 transitions are stored, not used).
# Steps: (1) full run for each valid sweep imputation; outputs written like 07: extras_int_mXXX.rds, baddata_int_mXXX.rds, and (own files, no overwriting)
#            mi_raw_sweepfull_int.csv, mi_diag_sweepfull_int.csv
#        (2) compare the black_spruce->regen fit of each full run with the stored sweep fit (all 16 parameters, difference in SD units) and re-apply the frozen screen;
#        (3) validity = the frozen screen applied to the regular full-run fit (no other rule). Rubin-pool the joint covariance and the estimates over the 20 valid fits of m1-114
#            plus the valid full runs. The comparison with the sweep fit in (2) is informational only.
# Outputs also: bsregen_refit_compare_int.csv, pooled_int_bsregen_cov.rds/.csv, pooled_int_bsregen_regular.csv, pooled_int_bsregen_regular_nolowratio.csv
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(INLA) })
SIM_DIR <- "/Users/lukas/Desktop/RProjects/Quebec-forest-succession/R-scripts/Sim4stateExact/p1p2_v2_test"
SEED <- 2026L; N_WORKERS <- 6L; M_MI <- 30L; SIM_M <- 200L; RUN_TAG <- "int"; SD_MIN <- 1e-3; N_COV_SAMPLES <- 1000L; LAST_M_EXISTING <- 114L
RATIO_MIN <- 0.1; ICPT_REF <- -4.985135842; ICPT_SD_REF <- 0.06137571849; Z_MAX <- 4; REPRO_TOL <- 0.01     # frozen screen (same as 07b); reproduction tolerance in SD units
LABEL <- "black_spruce->regen"; OUT_DIR <- here("data", "processed", "fit_int"); options(width = 200)
source(file.path(SIM_DIR, "model_fitting.R")); source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED); fit_fun <- fit_inla_models
inp <- prepare_fit_inputs(here("data", "processed"), subsample_n = NULL, seed = SEED)
s2 <- readRDS(file.path(OUT_DIR, "stage2.rds"))
K <- match(LABEL, inp$td$label); stopifnot(!is.na(K))
tms_k <- colnames(inp$mask)[inp$mask[K, ]]
inp_k <- list(td = inp$td[K], mask = inp$mask[K, , drop = FALSE])
cols <- c("from", "to", "status", "Tstop", "truncation", inp$terms)

build_data_k <- function(m) {                                                                     # identical to 07b
  set.seed(SEED * 1000L + m)
  th <- draw_theta_m(s2$models, inp, s2$theta)
  p1 <- p1_entries(th, inp, SIM_M, "draw", 1L)
  dg <- new.env()
  b <- build_rows(inp$iv, setNames(p1$entry, p1$id_pe), inp$states, inp$tr, p2_fun = make_p2_moment(th, inp, dg), keep_cols = inp$keep_cols)
  d <- b$rows[kept == TRUE]; idx <- match(paste(d$from, d$to, sep = "->"), inp$td$label); stopifnot(!anyNA(idx))
  as.data.frame(d[idx == K, ..cols])
}
screen <- function(raw) {                                                                         # identical rule to 07b
  a <- raw[parameter == "alpha", .(m, alpha = est, a_sd = sd, n_events, fit_error)]; a[, ratio := a_sd / (0.78 * alpha / sqrt(n_events))]
  i <- raw[parameter == "(Intercept)", .(m, icpt = est)]; x <- merge(a, i, by = "m", all.x = TRUE); x[, z_icpt := (icpt - ICPT_REF) / ICPT_SD_REF]
  x[, reason := fifelse(fit_error %in% TRUE | is.na(alpha) | is.na(a_sd) | is.na(icpt), "error/NA",
                fifelse(ratio < RATIO_MIN, "alpha_sd_ratio", fifelse(abs(z_icpt) > Z_MAX, "intercept", "ok")))]
  x[, valid := reason == "ok"]; x[order(m)] }
pool_cov <- function(covs) {                                                                      # covs: list of list(mean, cov); Rubin for a vector parameter
  M <- length(covs); nm <- colnames(covs[[1]]$cov); stopifnot(all(vapply(covs, function(x) identical(colnames(x$cov), nm), logical(1))))
  Q <- do.call(rbind, lapply(covs, function(x) x$mean[nm])); Wbar <- Reduce(`+`, lapply(covs, function(x) x$cov)) / M
  B <- stats::cov(Q); Tm <- Wbar + (1 + 1 / M) * B; dimnames(Tm) <- list(nm, nm)
  list(M = M, mean = colMeans(Q), W = Wbar, B = B, T = Tm) }

# ---- sweep fits to redo: the valid ones among m >= 115 (frozen screen applied to the stored sweep fits) -------------------------------------------------
allr <- rbind(fread(file.path(OUT_DIR, sprintf("mi_raw_%s.csv", RUN_TAG))), fread(file.path(OUT_DIR, sprintf("mi_raw_redraw_%s.csv", RUN_TAG))),
              fread(file.path(OUT_DIR, sprintf("mi_raw_sweep_%s.csv", RUN_TAG))), fill = TRUE)[label == LABEL]
sc <- screen(allr); old_valid <- sc[valid == TRUE & m <= LAST_M_EXISTING, m]; todo <- sc[valid == TRUE & m > LAST_M_EXISTING, m]
cat(sprintf("Valid fits: %d (m1-114: %d with stored covariance; sweep: %d to refit: %s)\n", sum(sc$valid), length(old_valid), length(todo), paste(todo, collapse = ",")))

raw_path <- file.path(OUT_DIR, sprintf("mi_raw_sweepfull_%s.csv", RUN_TAG)); diag_path <- file.path(OUT_DIR, sprintf("mi_diag_sweepfull_%s.csv", RUN_TAG))
ex_path <- function(m) file.path(OUT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, m)); ex_path_old <- ex_path
todo_run <- todo[!file.exists(vapply(todo, ex_path, ""))]                                           # resumable
cat(sprintf("Full runs needed: %d (%s)\n", length(todo_run), paste(todo_run, collapse = ",")))
t0 <- proc.time()["elapsed"]
for (i in if (length(todo_run)) seq(1L, length(todo_run), by = N_WORKERS) else integer(0)) {
  batch <- todo_run[i:min(i + N_WORKERS - 1L, length(todo_run))]
  res <- mclapply(batch, function(m) tryCatch({ set.seed(SEED * 1000L + m)
      run_one_imputation(m, inp, s2$models, s2$theta, fit_fun, M = SIM_M, keep_models = FALSE) }, error = function(e) list(error = conditionMessage(e), m = m)), mc.cores = N_WORKERS)
  for (r in res) {
    if (is.null(r$res)) { cat(sprintf("m=%s FAILED to run: %s\n", r$m, r$error)); next }
    if (!is.null(r$bad_data)) saveRDS(r$bad_data, file.path(OUT_DIR, sprintf("baddata_%s_m%02d.rds", RUN_TAG, r$diag$m)))
    saveRDS(r$extras, ex_path(r$diag$m)); fwrite(r$res, raw_path, append = file.exists(raw_path)); fwrite(r$diag, diag_path, append = file.exists(diag_path)) }
  cat(sprintf("  batch done: m = %s (%.0f s)\n", paste(batch, collapse = ","), proc.time()["elapsed"] - t0)) }

# ---- compare the black_spruce->regen fit of each full run with the stored sweep fit, re-screen ------------------------------------------------------------
full <- fread(raw_path)[label == LABEL]
cmp <- rbindlist(lapply(todo, function(mi) {
  nw <- full[m == mi]; if (!nrow(nw)) return(data.table(m = mi, max_diff_sd = NA_real_, refit_valid = FALSE, reason = "no full-run fit"))
  old <- allr[m == mi, .(parameter, est_old = est, sd_old = sd)]                                  # stored sweep fit (mi, not m: data.table column scoping)
  mm <- merge(nw[, .(parameter, est_new = est)], old, by = "parameter"); sv <- screen(nw)
  d <- abs(mm$est_new - mm$est_old) / mm$sd_old
  data.table(m = mi, max_diff_sd = if (all(is.na(d))) NA_real_ else max(d, na.rm = TRUE), refit_valid = sv$valid, reason = sv$reason) }))
fwrite(cmp, file.path(OUT_DIR, sprintf("bsregen_refit_compare_%s.csv", RUN_TAG)))
cat("\n==== full run vs stored sweep fit (max difference over all 16 parameters, in SD of the stored fit) ====\n"); print(cmp)
use_m <- cmp[refit_valid == TRUE, m]; drop_m <- setdiff(todo, use_m)
cat(sprintf("\nValid under the frozen screen in the regular full run: %d of %d. Not valid: %s\n", length(use_m), length(todo),
            if (length(drop_m)) paste(sprintf("m%d (%s)", drop_m, cmp[match(drop_m, m), reason]), collapse = ", ") else "none"))
# The regular full-run fit is the definitive fit of an imputation (as for every other imputation); the sweep (07b) was a prefilter. The comparison above is informational.
covs <- c(lapply(old_valid, function(m) { f <- readRDS(ex_path_old(m))$fits[[LABEL]]$cov; stopifnot(is.null(f$error), !is.null(f$cov)); f }),
          lapply(use_m, function(m) { f <- readRDS(ex_path(m))$fits[[LABEL]]$cov; stopifnot(is.null(f$error), !is.null(f$cov)); f }))
pc <- pool_cov(covs); saveRDS(pc, file.path(OUT_DIR, sprintf("pooled_%s_bsregen_cov.rds", RUN_TAG)))
fwrite(as.data.table(pc$T, keep.rownames = "parameter"), file.path(OUT_DIR, sprintf("pooled_%s_bsregen_cov.csv", RUN_TAG)))
# estimates pooled from the same fits (regular pipeline fits only): old valid fits + valid full runs
fin <- rbind(allr[m %in% old_valid], full[m %in% use_m], fill = TRUE); pf <- pool_all(fin)
fwrite(pf, file.path(OUT_DIR, sprintf("pooled_%s_bsregen_regular.csv", RUN_TAG)))
low <- screen(fin)[ratio < 0.5, m]; pl <- pool_all(fin[!m %in% low]); fwrite(pl, file.path(OUT_DIR, sprintf("pooled_%s_bsregen_regular_nolowratio.csv", RUN_TAG)))
cat(sprintf("\nPooled %s from %d regular-pipeline fits (covariance from the same %d). Sensitivity without fits of alpha-SD ratio < 0.5 (m %s): %d fits.\n",
            LABEL, uniqueN(fin$m), pc$M, paste(low, collapse = ","), uniqueN(fin[!m %in% low]$m)))
cmpx <- merge(pf[, .(parameter, m_eff, est = round(pooled_est, 4), se = round(se, 4), fmi = round(fmi, 2))], pl[, .(parameter, m_nolow = m_eff, est_nolow = round(pooled_est, 4), se_nolow = round(se, 4))], by = "parameter", sort = FALSE)
print(cmpx)
chk <- data.table(parameter = colnames(pc$T), se_from_cov = sqrt(diag(pc$T)))[pf[, .(parameter, se_pooled_table = se)], on = "parameter"]; chk[, ratio := se_from_cov / se_pooled_table]
cat("\nDiagonal of the pooled covariance vs the pooled SE table (ratio should be near 1):\n"); print(chk)
