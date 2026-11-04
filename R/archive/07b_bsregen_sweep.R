# ARCHIVE (2026-10-22): used for the 'int' run only. For the 'gs' run and later use 07_redraw.R.
# ==============================================================================
# 07b_bsregen_sweep.R   RUN ON THE MACHINE WITH INLA (resumable; batches of N_WORKERS)
# Purpose: finish black_spruce->regen (the only transition without 30 valid fits) by fitting ONLY that transition on fresh candidate imputations.
#   Same data build as the full MI loop (same functions, same seed rule set.seed(SEED * 1000 + m), same stage-2 models); only transition k is fitted
#   (fit_inla_models is a plain lapply over transitions, so the INLA call for k is identical). Posterior samples for the other 90 transitions are not drawn
#   (they are not needed here). The other 90 transitions are NOT touched: they already have 30 valid fits (see 07_redraw_failed.R).
# FROZEN SCREEN (identical rule for every transition; constants for this transition fixed in advance, from the 114 existing fits):
#   a fit is VALID if: no INLA error, est/sd not NA, alpha SD ratio  sd / (0.78 * alpha / sqrt(events)) >= 0.1,
#   and |intercept - ICPT_REF| / ICPT_SD_REF <= Z_MAX   (intercept stuck/shifted = collapse signature; healthy fits are within 0.62 of this scale, failures >= 11).
# STOPPING RULE (set before the run): stop when TARGET_VALID valid fits exist among m = 1..max (original 30 + redraw m31-114 + sweep); pool ALL valid fits.
# CHECK BEFORE SWEEP: the rebuilt data of one imputation with a saved dataset (baddata_int_mXX.rds) must be identical to the saved data (else stop).
# Outputs (data/processed/fit_int/): mi_raw_sweep_int.csv (every candidate, black_spruce->regen rows), sweep_screen_int.csv (screen result of ALL candidates m=1..max),
#   pooled_int_bsregen_final.csv (Rubin, all valid), pooled_int_bsregen_drop3.csv (sensitivity: m30, m99, m109 removed)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(INLA) })
SIM_DIR <- "/Users/lukas/Desktop/RProjects/Quebec-forest-succession/R-scripts/Sim4stateExact/p1p2_v2_test"
SEED <- 2026L; N_WORKERS <- 6L; M_MI <- 30L; SIM_M <- 200L; RUN_TAG <- "int"; SD_MIN <- 1e-3; N_COV_SAMPLES <- 0L
TARGET_VALID <- 33L; MAX_CAND <- 400L; LAST_M_EXISTING <- 114L
RATIO_MIN <- 0.1; ICPT_REF <- -4.985135842; ICPT_SD_REF <- 0.06137571849; Z_MAX <- 4          # frozen screen constants
LABEL <- "black_spruce->regen"; SUSPECT_M <- c(30L, 99L, 109L)                                  # passed every screen; sensitivity only
OUT_DIR <- here("data", "processed", "fit_int"); options(width = 200)
source(file.path(SIM_DIR, "model_fitting.R")); source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED); fit_fun <- fit_inla_models
inp <- prepare_fit_inputs(here("data", "processed"), subsample_n = NULL, seed = SEED)
s2 <- readRDS(file.path(OUT_DIR, "stage2.rds"))
K <- match(LABEL, inp$td$label); stopifnot(!is.na(K))
inp_k <- list(td = inp$td[K], mask = inp$mask[K, , drop = FALSE])                                # for tidy_fit() on this one transition
cols <- c("from", "to", "status", "Tstop", "truncation", inp$terms)

build_data_k <- function(m) {                                                                     # identical to run_one_imputation() up to fit_rows()
  set.seed(SEED * 1000L + m)
  th <- draw_theta_m(s2$models, inp, s2$theta)
  p1 <- p1_entries(th, inp, SIM_M, "draw", 1L)
  dg <- new.env()
  b <- build_rows(inp$iv, setNames(p1$entry, p1$id_pe), inp$states, inp$tr, p2_fun = make_p2_moment(th, inp, dg), keep_cols = inp$keep_cols)
  d <- b$rows[kept == TRUE]; idx <- match(paste(d$from, d$to, sep = "->"), inp$td$label); stopifnot(!anyNA(idx))
  as.data.frame(d[idx == K, ..cols])
}
fit_one <- function(m) {
  dk <- build_data_k(m)
  mod <- fit_fun(list(dk), rhs_formula = inp$rhs[K], exact = TRUE)
  res <- tidy_fit(mod, inp_k); res[, `:=`(k = K, m = m, n_events = sum(dk$status == 1L))]; res[] }

screen <- function(raw) {                                                                         # one row per m for the transition: valid TRUE/FALSE + reason
  a <- raw[parameter == "alpha", .(m, alpha = est, a_sd = sd, n_events, fit_error)]
  a[, ratio := a_sd / (0.78 * alpha / sqrt(n_events))]
  i <- raw[parameter == "(Intercept)", .(m, icpt = est)]
  x <- merge(a, i, by = "m", all.x = TRUE); x[, z_icpt := (icpt - ICPT_REF) / ICPT_SD_REF]
  x[, reason := fifelse(fit_error %in% TRUE | is.na(alpha) | is.na(a_sd) | is.na(icpt), "error/NA",
                fifelse(ratio < RATIO_MIN, "alpha_sd_ratio", fifelse(abs(z_icpt) > Z_MAX, "intercept", "ok")))]
  x[, valid := reason == "ok"]; x[order(m)] }

# ---- existing fits of this transition (original 30 + redraw m31-114) ----------------------------------------------------------------------------------
old <- rbind(fread(file.path(OUT_DIR, sprintf("mi_raw_%s.csv", RUN_TAG))), fread(file.path(OUT_DIR, sprintf("mi_raw_redraw_%s.csv", RUN_TAG))), fill = TRUE)[label == LABEL]
sweep_path <- file.path(OUT_DIR, sprintf("mi_raw_sweep_%s.csv", RUN_TAG))
sw <- function() if (file.exists(sweep_path)) fread(sweep_path) else old[0]
n_valid <- function() sum(screen(rbind(old, sw(), fill = TRUE))$valid)
cat(sprintf("Existing fits of %s: %d candidates, valid under the frozen screen: %d (target %d)\n", LABEL, uniqueN(old$m), n_valid(), TARGET_VALID))
print(screen(old)[valid == FALSE, .(m, alpha = round(alpha, 3), ratio = round(ratio, 3), icpt = round(icpt, 2), z_icpt = round(z_icpt, 1), reason)][order(m)][1:min(.N, 100)])

# ---- data identity check (once) ----------------------------------------------------------------------------------------------------------------------------
if (!file.exists(sweep_path)) {
  bf <- list.files(OUT_DIR, "^baddata_int_m[0-9]+\\.rds$", full.names = TRUE); chk <- NULL
  for (f in bf) { b <- readRDS(f); if (LABEL %in% names(b)) { chk <- list(m = as.integer(sub(".*_m([0-9]+)\\.rds", "\\1", f)), d = b[[LABEL]]$data); break } }
  if (is.null(chk)) stop("no saved dataset of ", LABEL, " found for the identity check")
  rb <- build_data_k(chk$m)
  ok <- isTRUE(all.equal(unname(as.list(rb)), unname(as.list(chk$d)), tolerance = 1e-10, check.attributes = FALSE))
  cat(sprintf("Identity check m=%d: rebuilt %d rows vs saved %d rows -> %s\n", chk$m, nrow(rb), nrow(chk$d), if (ok) "IDENTICAL" else "MISMATCH"))
  if (!ok) stop("rebuilt data differ from the saved fitting data; sweep not started")
}

# ---- sweep ---------------------------------------------------------------------------------------------------------------------------------------------
done <- sort(unique(sw()$m)); t0 <- proc.time()["elapsed"]
while (n_valid() < TARGET_VALID && length(done) < MAX_CAND) {
  b0 <- max(c(LAST_M_EXISTING, done)) + 1L; batch <- b0:(b0 + N_WORKERS - 1L)
  out <- mclapply(batch, function(m) tryCatch(fit_one(m), error = function(e) structure(list(msg = conditionMessage(e), m = m), class = "sweep_err")), mc.cores = N_WORKERS)
  for (r in out) { if (inherits(r, "sweep_err") || !is.data.table(r)) { cat(sprintf("m=%s FAILED to run: %s\n", if (is.list(r)) r$m else "?", if (is.list(r)) r$msg else "unknown")); next }
    fwrite(r, sweep_path, append = file.exists(sweep_path)) }
  done <- sort(unique(sw()$m))
  cat(sprintf("  batch to m=%d (%.0f s) | sweep candidates %d | valid in total %d (target %d)\n", max(batch), proc.time()["elapsed"] - t0, length(done), n_valid(), TARGET_VALID))
}
if (n_valid() < TARGET_VALID) cat("!! target not reached within MAX_CAND candidates\n")

# ---- record + pooling -----------------------------------------------------------------------------------------------------------------------------------
allr <- rbind(old, sw(), fill = TRUE); sc <- screen(allr); fwrite(sc, file.path(OUT_DIR, sprintf("sweep_screen_%s.csv", RUN_TAG)))
cat("\nScreen result, all candidates m=1..max: "); print(sc[, .N, by = reason])
fin <- allr[m %in% sc[valid == TRUE, m]]
pf <- pool_all(fin); fwrite(pf, file.path(OUT_DIR, sprintf("pooled_%s_bsregen_final.csv", RUN_TAG)))
pd <- pool_all(fin[!m %in% SUSPECT_M]); fwrite(pd, file.path(OUT_DIR, sprintf("pooled_%s_bsregen_drop3.csv", RUN_TAG)))
cat(sprintf("\nPooled %s: %d valid fits (all) | %d (m30, m99, m109 removed)\n", LABEL, uniqueN(fin$m), uniqueN(fin[!m %in% SUSPECT_M]$m)))
cmp <- merge(pf[, .(parameter, est_all = round(pooled_est, 4), se_all = round(se, 4), m_all = m_eff, fmi_all = round(fmi, 2))],
             pd[, .(parameter, est_drop3 = round(pooled_est, 4), se_drop3 = round(se, 4), m_drop3 = m_eff)], by = "parameter", sort = FALSE)
print(cmp)
