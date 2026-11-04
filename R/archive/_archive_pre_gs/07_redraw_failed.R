# ==============================================================================
# 07_redraw_failed.R   RUN ON THE MACHINE WITH INLA (long: full imputation draws, batches of 6; resumable)
# PER-TRANSITION replacement (default, MODE = "transition"): only the transitions with a numerical failure get replacement fits; all other transitions are untouched.
#   failure  = INLA error / NA alpha, or alpha SD ratio  sd / (0.78 * alpha / sqrt(events)) < 0.1   (empty gap 0.009-0.265 in the 2727 fits)
#   original: black_spruce->regen 8 failing fits, other_temperate->white_birch 3 (computed from mi_raw_diag.csv, not hard-coded)
#   candidates = fresh imputations m = 31, 32, ... (same code and seed rule as the MI loop: set.seed(SEED * 1000 + m)); the loop stops when every failing
#   transition has as many valid NEW fits as it lost. For each such transition the first valid candidates (increasing m) replace its failing fits.
#   (MODE = "imputation": whole imputations are replaced for all 91 transitions; the earlier version is kept as 07_redraw_failed_wholeimputation.R)
# The original 30 (mi_raw_diag.csv etc.) are NOT modified. Outputs (data/processed/fit/):
#   mi_raw_redraw_diag.csv, mi_diag_redraw_diag.csv, redraw_screen_diag.csv (every candidate x transition that fails), redraw_accepted_diag.csv,
#   pooled_diag_redraw.csv (final: 30 fits per transition), redraw_vs_asis_diag.csv (shift per parameter in SE units; only the replaced transitions can differ)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(INLA) })
SIM_DIR <- "/Users/lukas/Desktop/RProjects/Quebec-forest-succession/R-scripts/Sim4stateExact/p1p2_v2_test"
SEED <- 2026L; N_WORKERS <- 6L; M_MI <- 30L; SIM_M <- 200L; RUN_TAG <- "int"; SD_MIN <- 1e-3; N_COV_SAMPLES <- 1000L; MAX_CAND <- 80L; RATIO_MIN <- 0.1
OUT_DIR <- here("data", "processed", "fit_int"); options(width = 200)
source(file.path(SIM_DIR, "model_fitting.R")); source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED); fit_fun <- fit_inla_models
inp <- prepare_fit_inputs(here("data", "processed"), subsample_n = NULL, seed = SEED)
s2 <- readRDS(file.path(OUT_DIR, "stage2.rds"))
# pre-specified screen
fail_table <- function(raw) {                       # one row per (m, transition) that fails
  a <- raw[parameter == "alpha"]
  a[, ratio := sd / (0.78 * est / sqrt(n_events))]
  a[fit_error == TRUE | is.na(est) | is.na(sd) | is.na(ratio) | ratio < RATIO_MIN, .(m, k, label, ratio, fit_error)] }
raw0 <- fread(file.path(OUT_DIR, sprintf("mi_raw_%s.csv", RUN_TAG)))
f0 <- fail_table(raw0); need <- f0[, .(n_need = .N, failed_m = paste(sort(m), collapse = ",")), by = .(k, label)]
cat("Original 30 fits: failing fits per transition (these, and only these, get replaced):\n"); print(need)
raw_path <- file.path(OUT_DIR, sprintf("mi_raw_redraw_%s.csv", RUN_TAG)); diag_path <- file.path(OUT_DIR, sprintf("mi_diag_redraw_%s.csv", RUN_TAG))
scr_path <- file.path(OUT_DIR, sprintf("redraw_screen_%s.csv", RUN_TAG)); acc_path <- file.path(OUT_DIR, sprintf("redraw_accepted_%s.csv", RUN_TAG))
done <- if (file.exists(diag_path)) fread(diag_path)$m else integer(0)
valid_per_k <- function() {                      # valid NEW fits per failing transition
  if (!file.exists(raw_path)) return(setNames(rep(0L, nrow(need)), need$k))
  r <- fread(raw_path); bad <- fail_table(r)
  vapply(need$k, function(kk) { mm <- unique(r[k == kk, m]); sum(!(mm %in% bad[k == kk, m])) }, 1L) |> setNames(need$k) }
enough <- function() all(valid_per_k() >= need$n_need)
t0 <- proc.time()["elapsed"]
while (!enough() && length(done) < MAX_CAND) {
  b0 <- max(c(M_MI, done)) + 1L; batch <- b0:(b0 + N_WORKERS - 1L)
  out <- mclapply(batch, function(m) tryCatch({ set.seed(SEED * 1000L + m)
      run_one_imputation(m, inp, s2$models, s2$theta, fit_fun, M = SIM_M, keep_models = FALSE) }, error = function(e) list(error = conditionMessage(e), m = m)), mc.cores = N_WORKERS)
  for (r in out) {
    if (is.null(r$res)) { cat(sprintf("m=%s FAILED to run: %s\n", r$m, r$error)); next }
    if (!is.null(r$bad_data)) saveRDS(r$bad_data, file.path(OUT_DIR, sprintf("baddata_%s_m%02d.rds", RUN_TAG, r$diag$m)))
    saveRDS(r$extras, file.path(OUT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, r$diag$m)))
    fwrite(r$res, raw_path, append = file.exists(raw_path)); fwrite(r$diag, diag_path, append = file.exists(diag_path)) }
  done <- fread(diag_path)$m
  cat(sprintf("  batch to m=%d done (%.0f s) | candidates %d | valid new fits: %s (needed %s)\n", max(batch), proc.time()["elapsed"] - t0, length(done),
      paste(valid_per_k(), collapse = "/"), paste(need$n_need, collapse = "/")))
}
rawc <- fread(raw_path); fc <- fail_table(rawc)
# screening record: every candidate x transition (all 91) that fails, plus candidate failure rate per transition
fwrite(fc[, .(m, k, label, ratio, fit_error)], scr_path)
cat("\n==== candidate fits failing the screen, per transition (all transitions) ====\n"); print(fc[, .(n_fail = .N, of = uniqueN(rawc$m)), by = label])
# accepted: per failing transition, first n_need valid candidates in increasing m
acc <- rbindlist(lapply(seq_len(nrow(need)), function(i) {
  kk <- need$k[i]; ok <- sort(setdiff(unique(rawc[k == kk, m]), fc[k == kk, m]))
  data.table(k = kk, label = need$label[i], m_new = head(ok, need$n_need[i]), m_replaced = head(sort(f0[k == kk, m]), length(head(ok, need$n_need[i])))) }))
fwrite(acc, acc_path); cat("\n==== replacements (new candidate m replaces failed original m, per transition) ====\n"); print(acc)
if (!enough()) cat("!! not all replacements found within MAX_CAND candidates\n")
# final: 30 fits per failing transition; untouched for all others
drop <- f0[, .(m, k)]; keep0 <- raw0[!drop, on = .(m, k)]
add <- merge(rawc, acc[, .(k, m = m_new)], by = c("k", "m"))
fin <- rbind(keep0, add, fill = TRUE); cat(sprintf("Final: fits per transition min %d max %d\n", min(fin[parameter == "alpha", .N, by = k]$N), max(fin[parameter == "alpha", .N, by = k]$N)))
pf <- pool_all(fin); fwrite(pf, file.path(OUT_DIR, sprintf("pooled_%s_redraw.csv", RUN_TAG)))
pa <- fread(file.path(OUT_DIR, sprintf("pooled_%s.csv", RUN_TAG))); px <- fread(file.path(OUT_DIR, sprintf("pooled_%s_excl.csv", RUN_TAG)))
cmp <- merge(pa[, .(k, label, parameter, est_asis = pooled_est, se_asis = se, m_asis = m_eff)], pf[, .(k, label, parameter, est_redraw = pooled_est, se_redraw = se, m_redraw = m_eff)], by = c("k", "label", "parameter"))
cmp <- merge(cmp, px[, .(k, label, parameter, est_excl = pooled_est, se_excl = se)], by = c("k", "label", "parameter"), all.x = TRUE)
cmp[, `:=`(shift_se = (est_redraw - est_asis) / se_asis, shift_vs_excl_se = (est_redraw - est_excl) / se_excl)]; fwrite(cmp, file.path(OUT_DIR, sprintf("redraw_vs_asis_%s.csv", RUN_TAG)))
cat("\n==== black_spruce->regen alpha and the other_temperate->white_birch alpha ====\n")
print(cmp[parameter == "alpha" & label %in% c("black_spruce->regen", "other_temperate->white_birch"), .(label, est_asis = round(est_asis, 3), se_asis = round(se_asis, 4), m_asis, est_excl = round(est_excl, 3), se_excl = round(se_excl, 4), est_redraw = round(est_redraw, 3), se_redraw = round(se_redraw, 4), m_redraw)])
cat("\n==== shift of the pooled estimate, redraw vs as-is, in SE units, all parameters ====\n")
print(cmp[, .(n = .N, median_abs = round(median(abs(shift_se), na.rm = TRUE), 3), p95_abs = round(quantile(abs(shift_se), .95, na.rm = TRUE), 3), max_abs = round(max(abs(shift_se), na.rm = TRUE), 3))])
print(cmp[order(-abs(shift_se))][1:6, .(label, parameter, shift_se = round(shift_se, 3))])
cat("\nSE ratio redraw / as-is: median", round(median(cmp$se_redraw / cmp$se_asis, na.rm = TRUE), 3), "\n")
