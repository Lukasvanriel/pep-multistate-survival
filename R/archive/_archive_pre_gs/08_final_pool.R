# ==============================================================================
# 08_final_pool.R   (no INLA needed; seconds to run; run from the project root: Rscript R/08_final_pool.R)
# Purpose: final pooled results for ALL 91 transitions with the SAME rule for every transition.
#   1. Pull every regular-pipeline fit: mi_raw_int.csv (m1-30), mi_raw_redraw_int.csv (m31-114), mi_raw_sweepfull_int.csv (regular full runs of the sweep-valid imputations).
#      (The 07b sweep fits are NOT used: they were a prefilter; the regular full-run fit decides validity.)
#   2. Apply the frozen screen to every fit of every transition. VALID if: no INLA error/NA; alpha-SD ratio sd/(0.78*alpha/sqrt(events)) >= 0.1;
#      |intercept - ref| / ref_sd <= 4, where ref / ref_sd = median intercept / median intercept SD of the transition's ratio-passing fits in m1-114.
#      (For black_spruce->regen this reproduces the constants frozen in 07b; checked below.)
#   3. Use the first TARGET (30) valid fits per transition in increasing m (selection by index only, never by results). Stop if a transition has fewer.
#   4. Rubin-pool estimates (pool_all) and the joint covariance (Rubin on the stored covariances) over exactly those fits.
# Outputs (data/processed/fit_int/): fit_validity_int.csv (every fit: screen result, used yes/no), pooled_int_final.csv, pooled_int_final_cov.rds (named list by transition)
# ==============================================================================
suppressPackageStartupMessages(library(data.table)); options(width = 200)
RUN_TAG <- "int"; TARGET <- 30L; RATIO_MIN <- 0.1; Z_MAX <- 4; REF_M_MAX <- 114L; SD_MIN <- 1e-3
BS_LABEL <- "black_spruce->regen"; BS_ICPT_REF <- -4.985135842; BS_ICPT_SD_REF <- 0.06137571849   # only used to CHECK the reference reproduces the frozen constants
OUT_DIR <- file.path("data", "processed", "fit_int")
ex <- parse(file.path("R", "05a_fit_funs.R"))                                                       # take pool_rubin / flag_fits / pool_all only (no INLA needed)
for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]])[1] %in% c("pool_rubin", "flag_fits", "pool_all")) eval(e)
pool_cov <- function(covs) {
  M <- length(covs); nm <- colnames(covs[[1]]$cov); stopifnot(all(vapply(covs, function(x) identical(colnames(x$cov), nm), logical(1))))
  Q <- do.call(rbind, lapply(covs, function(x) x$mean[nm])); Wbar <- Reduce(`+`, lapply(covs, function(x) x$cov)) / M
  B <- stats::cov(Q); Tm <- Wbar + (1 + 1 / M) * B; dimnames(Tm) <- list(nm, nm); list(M = M, mean = colMeans(Q), W = Wbar, B = B, T = Tm) }

# ---- 1. all regular-pipeline fits ---------------------------------------------------------------------------------------------------------------------------
files <- sprintf(c("mi_raw_%s.csv", "mi_raw_redraw_%s.csv", "mi_raw_sweepfull_%s.csv"), RUN_TAG)
raw <- rbindlist(lapply(files, function(f) fread(file.path(OUT_DIR, f))), fill = TRUE)
stopifnot(!anyDuplicated(raw[, .(k, m, parameter)]))
cat(sprintf("Fits pulled: %d (imputation x transition) from %d imputations, %d transitions\n", uniqueN(raw[, .(k, m)]), uniqueN(raw$m), uniqueN(raw$k)))

# ---- 2. frozen screen on every fit ----------------------------------------------------------------------------------------------------------------------
a <- raw[parameter == "alpha", .(k, label, m, alpha = est, a_sd = sd, n_events, fit_error)]
a[, ratio := a_sd / (0.78 * alpha / sqrt(n_events))]
i <- raw[parameter == "(Intercept)", .(k, m, icpt = est, i_sd = sd)]
x <- merge(a, i, by = c("k", "m"), all.x = TRUE)
x[, ref := m <= REF_M_MAX & !(fit_error %in% TRUE) & !is.na(ratio) & ratio >= RATIO_MIN]
x[, `:=`(icpt_ref = median(icpt[ref]), icpt_sd_ref = median(i_sd[ref])), by = k]
chk <- unique(x[label == BS_LABEL, .(icpt_ref, icpt_sd_ref)])
cat(sprintf("Reference check, %s: intercept ref %.9f (frozen %.9f), SD ref %.9f (frozen %.9f)\n", BS_LABEL, chk$icpt_ref, BS_ICPT_REF, chk$icpt_sd_ref, BS_ICPT_SD_REF))
stopifnot(abs(chk$icpt_ref - BS_ICPT_REF) < 1e-6, abs(chk$icpt_sd_ref - BS_ICPT_SD_REF) < 1e-6)
x[, z_icpt := (icpt - icpt_ref) / icpt_sd_ref]
x[, reason := fifelse(fit_error %in% TRUE | is.na(alpha) | is.na(a_sd) | is.na(icpt), "error/NA",
              fifelse(ratio < RATIO_MIN, "alpha_sd_ratio", fifelse(abs(z_icpt) > Z_MAX, "intercept", "ok")))]
x[, valid := reason == "ok"]

# ---- 3. first TARGET valid fits per transition, by increasing m -------------------------------------------------------------------------------------------
setorder(x, k, m)
x[, rank_valid := ifelse(valid, cumsum(valid), NA_integer_), by = k]
x[, used := valid & rank_valid <= TARGET]
summ <- x[, .(label = label[1], candidates = .N, valid = sum(valid), used = sum(used), max_m_used = max(m[used])), by = k]
cat("\nValid fits and fits used per transition (summary):\n"); print(summ[, .(transitions = .N), by = .(valid_ge_target = valid >= TARGET, used)])
cat(sprintf("Transitions whose %d used fits include replacements (max m used > %d): %d\n", TARGET, TARGET, summ[max_m_used > TARGET, .N])); print(summ[max_m_used > TARGET][order(-max_m_used)][1:min(.N, 40)])
if (any(summ$used < TARGET)) { print(summ[used < TARGET]); stop("some transition has fewer than ", TARGET, " valid fits; nothing pooled") }
fwrite(x[, .(k, label, m, alpha, ratio, icpt, z_icpt, reason, valid, used)], file.path(OUT_DIR, sprintf("fit_validity_%s.csv", RUN_TAG)))

# ---- 4. pool estimates and joint covariance over exactly the used fits ---------------------------------------------------------------------------------------
fin <- raw[x[used == TRUE, .(k, m)], on = .(k, m)]
pf <- pool_all(fin); stopifnot(all(pf$m_eff == TARGET))
fwrite(pf, file.path(OUT_DIR, sprintf("pooled_%s_final.csv", RUN_TAG)))
cat(sprintf("\nPooled: %d rows (%d transitions x parameters), m_eff = %d everywhere\n", nrow(pf), uniqueN(pf$k), TARGET))
old <- fread(file.path(OUT_DIR, sprintf("pooled_%s_redraw.csv", RUN_TAG)))
cm <- merge(pf[, .(k, label, parameter, est = pooled_est, se)], old[, .(k, parameter, est_old = pooled_est, se_old = se, m_old = m_eff)], by = c("k", "parameter"))
cm[, d_se := abs(est - est_old) / se]
cat("Comparison with the earlier pooled_int_redraw.csv (SE units; transitions with identical fit sets should show 0):\n")
print(cm[, .(max_shift_in_SE = round(max(d_se, na.rm = TRUE), 4), m_old = m_old[1]), by = label][max_shift_in_SE > 1e-8 | m_old != TARGET])

# ---- joint covariance: correct alpha's spread, then Rubin-pool -----------------------------------------------------------------------------------------------
# Why: inla.posterior.sample() draws the hyperparameter (alpha) only from INLA's integration grid (step dz = 0.75 SD, stopping when the log density has fallen by
#   diff.logdens = 6; r-inla.org docs), whereas the intercept/betas are drawn from their exact Gaussian conditional on alpha. So the betas are fine, but alpha's draws
#   can miss its long right tail (sparse transitions with alpha near 1). INLA's own posterior summary of alpha (summary.hyperpar: mean, sd) comes from the full marginal.
# Fix (law of total variance, per fit): the sample covariance is Cov = E[Sigma(alpha)] + b b' Var_s(alpha), with b = Cov_s(beta, alpha) / Var_s(alpha) the slope of the
#   conditional mean on alpha. Keep E[Sigma] = Cov_s(beta) - b b' Var_s(alpha) and the slope b, and put in INLA's own alpha mean and variance:
#   Var(alpha) = Var_INLA(alpha), mean(alpha) = mean_INLA(alpha), Cov(beta, alpha) = b Var_INLA(alpha). The beta/intercept block and means are left as sampled (see validation).
#   Assumes the conditional mean is ~linear in alpha and the conditional covariance ~constant over alpha's range. Validation below: corrected beta SDs vs INLA's
#   own summary.fixed SDs (which are not sampled), before and after the correction.
correct_alpha <- function(cv, sh, sf) {
  a <- ncol(cv$cov); nm <- colnames(cv$cov); stopifnot(nm[a] == "alpha")
  hr <- sh[grepl("alpha", rownames(sh)), ][1, ]; ma <- hr[["mean"]]; va <- hr[["sd"]]^2; s2 <- cv$cov[a, a]
  if (!is.finite(va) || !is.finite(s2) || s2 <= 0) return(NULL)
  bb <- cv$cov[-a, a] / s2; mu <- cv$mean; C <- cv$cov                                              # only alpha's mean/variance and its covariances with the betas are changed;
  mu[a] <- ma; C[-a, a] <- bb * va; C[a, -a] <- bb * va; C[a, a] <- va                                # the beta/intercept block is kept as sampled (exact Gaussian given alpha)
  list(mean = mu, cov = C, n = cv$n, psd = min(eigen(C, symmetric = TRUE, only.values = TRUE)$values) > 0) }
psd_fail <- logical(0); store <- list(); stor0 <- list(); val <- list(); no_cov <- list(); need <- x[used == TRUE, .(k, label, m)]
for (mm in sort(unique(need$m))) {
  f <- file.path(OUT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, mm)); if (!file.exists(f)) { no_cov[[length(no_cov) + 1L]] <- need[m == mm]; next }
  fits <- readRDS(f)$fits
  for (lb in need[m == mm, label]) { fx <- fits[[lb]]; cv <- fx$cov
    if (is.null(cv) || !is.null(cv$error)) { no_cov[[length(no_cov) + 1L]] <- need[m == mm & label == lb]; next }
    cc <- correct_alpha(cv, as.data.frame(fx$summary_hyper), fx$summary_fixed); if (is.null(cc)) { no_cov[[length(no_cov) + 1L]] <- need[m == mm & label == lb]; next }
    stor0[[lb]][[as.character(mm)]] <- cv; store[[lb]][[as.character(mm)]] <- cc; psd_fail <- c(psd_fail, !cc$psd)
    sd_i <- fx$summary_fixed[, "sd"]; names(sd_i) <- rownames(fx$summary_fixed); b <- setdiff(colnames(cv$cov), "alpha"); b <- b[b %in% names(sd_i)]
    val[[length(val) + 1L]] <- data.table(label = lb, m = mm, parameter = b, r_before = sqrt(diag(cv$cov)[b]) / sd_i[b], r_after = sqrt(diag(cc$cov)[b]) / sd_i[b],
                                          a_before = sqrt(cv$cov["alpha", "alpha"]) / (va_ <- as.data.frame(fx$summary_hyper)[grepl("alpha", rownames(fx$summary_hyper)), "sd"][1]),
                                          a_after = sqrt(cc$cov["alpha", "alpha"]) / va_) } }
cov_list <- list(); cov_list0 <- list()
for (lb in unique(need$label)) if (length(store[[lb]]) == TARGET) { cov_list[[lb]] <- pool_cov(unname(store[[lb]])); cov_list0[[lb]] <- pool_cov(unname(stor0[[lb]])) }
saveRDS(cov_list, file.path(OUT_DIR, sprintf("pooled_%s_final_cov.rds", RUN_TAG))); saveRDS(cov_list0, file.path(OUT_DIR, sprintf("pooled_%s_final_cov_uncorrected.rds", RUN_TAG)))
cat(sprintf("\nJoint covariance pooled for %d of %d transitions (each over %d fits): corrected version saved as pooled_%s_final_cov.rds, uncorrected as ..._uncorrected.rds\n", length(cov_list), uniqueN(need$label), TARGET, RUN_TAG))
if (length(no_cov)) { nc <- rbindlist(no_cov); cat("!! used fits without a usable stored covariance:\n"); print(nc[order(label, m)]) }
cat(sprintf("\nCorrected per-fit matrices not positive definite: %d of %d\n", sum(psd_fail), length(psd_fail)))
val <- rbindlist(val)
fwrite(val, file.path(OUT_DIR, sprintf("cov_correction_check_%s.csv", RUN_TAG)))
cat("\nValidation per fit, ratio to INLA's own SD (1 = agrees). Betas/intercept (INLA summary.fixed) and alpha (summary.hyperpar):\n")
print(rbind(data.table(what = "beta  before", t(round(quantile(val$r_before, c(0, .01, .5, .99, 1), na.rm = TRUE), 3))), data.table(what = "beta  after ", t(round(quantile(val$r_after, c(0, .01, .5, .99, 1), na.rm = TRUE), 3))),
            data.table(what = "alpha before", t(round(quantile(unique(val[, .(label, m, a_before)])$a_before, c(0, .01, .5, .99, 1), na.rm = TRUE), 3))),
            data.table(what = "alpha after ", t(round(quantile(unique(val[, .(label, m, a_after)])$a_after, c(0, .01, .5, .99, 1), na.rm = TRUE), 3)))))
d <- rbindlist(lapply(names(cov_list), function(lb) { p <- cov_list[[lb]]; p0 <- cov_list0[[lb]]; s <- pf[label == lb & parameter %in% colnames(p$T), .(parameter, se)]
  data.table(label = lb, parameter = s$parameter, r0 = sqrt(diag(p0$T))[s$parameter] / s$se, r = sqrt(diag(p$T))[s$parameter] / s$se) }))
cat("\nPooled covariance diagonal / pooled SE table, all transitions and parameters (before -> after correction):\n")
print(rbind(data.table(what = "alpha before", t(round(quantile(d[parameter == "alpha", r0], c(0, .5, 1)), 3))), data.table(what = "alpha after ", t(round(quantile(d[parameter == "alpha", r], c(0, .5, 1)), 3))),
            data.table(what = "other before", t(round(quantile(d[parameter != "alpha", r0], c(0, .5, 1)), 3))), data.table(what = "other after ", t(round(quantile(d[parameter != "alpha", r], c(0, .5, 1)), 3)))))
cat("\nbsregen pooled (check):\n"); print(pf[label == BS_LABEL, .(parameter, m_eff, est = round(pooled_est, 4), se = round(se, 4), fmi = round(fmi, 2))])
