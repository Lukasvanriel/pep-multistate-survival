# ==============================================================================
# 08_final_pool.R   (no INLA needed; seconds to run; run from the project root: Rscript R/08_final_pool.R)
# Purpose: final pooled results for ALL kept transitions with the SAME rule for every transition.
#   1. Pull every regular-pipeline fit: mi_raw_<tag>.csv (m1-30, 05_fit_inla.R) and mi_raw_redraw_<tag>.csv (candidates, 07_redraw.R).
#   2. Apply the screen (user decision 2026-11-05; replaces the old rule 3 "intercept vs. screen_ref median", z <= 4):
#      VALID if (1) no INLA error and est/sd not NA; (2) alpha-SD ratio sd/(0.78*alpha/sqrt(events)) >= 0.1;
#      (3') |z_map_alpha| <= 3 and |z_map_intercept| <= 3 (B6 columns of mi_raw / mi_raw_redraw: INLA mode vs. direct maximisation of the same log-posterior, fit_map_v0).
#      The screen needs no reference table (screen_ref_<tag>.csv is no longer read).
#   3. Use the first TARGET (30) valid fits per transition in increasing m (selection by index only, never by results). Stop if a transition has fewer.
#   4. Rubin-pool estimates (pool_all) and the joint covariance (Rubin on the stored per-fit covariances, as sampled by inla.posterior.sample;
#      no post-hoc correction of alpha's spread: the correction was tried for the "int" run and rejected, see PIPELINE_NOTES 2026-10-20).
#   5. (WP3, B3) intercept on the YEARS scale per used fit, b0 - alpha*log(10), Rubin-pooled -> intercept_years_<tag>.csv (reporting only).
#   6. (WP3, B6) cap_check_<tag>.csv: per transition, share of used fits with maxH_map / maxH_inla > 22.85 and median z_map_max / z_map_alpha / z_map_intercept (diagnostic).
# Reads FIT_DIR (data-output/<tag>/fits/); outputs in RESULTS_DIR (data-output/<tag>/results/): fit_validity_<tag>.csv (every fit: screen result, used yes/no), pooled_<tag>_final.csv,
#   pooled_<tag>_final_cov.rds (named list by transition), cov_check_<tag>.csv (sampled SD / INLA SD per fit and parameter).
# ==============================================================================
suppressPackageStartupMessages(library(data.table)); options(width = 200)
source(here::here("R", "config.R"))   # RUN_TAG, FIT_DIR, RESULTS_DIR
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
TARGET <- 30L; RATIO_MIN <- 0.1; Z_MAP_MAX <- 3
OUT_DIR <- FIT_DIR
ex <- parse(file.path("R", "05a_fit_funs.R"))                                                       # take pool_rubin / flag_fits / pool_all only (no INLA needed)
for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]])[1] %in% c("pool_rubin", "flag_fits", "pool_all")) eval(e)
pool_cov <- function(covs) {
  M <- length(covs); nm <- colnames(covs[[1]]$cov); stopifnot(all(vapply(covs, function(x) identical(colnames(x$cov), nm), logical(1))))
  Q <- do.call(rbind, lapply(covs, function(x) x$mean[nm])); Wbar <- Reduce(`+`, lapply(covs, function(x) x$cov)) / M
  B <- stats::cov(Q); Tm <- Wbar + (1 + 1 / M) * B; dimnames(Tm) <- list(nm, nm); list(M = M, mean = colMeans(Q), W = Wbar, B = B, T = Tm) }

# ---- 1. all regular-pipeline fits ---------------------------------------------------------------------------------------------------------------------------
files <- file.path(OUT_DIR, sprintf(c("mi_raw_%s.csv", "mi_raw_redraw_%s.csv"), RUN_TAG)); files <- files[file.exists(files)]
raw <- rbindlist(lapply(files, fread), fill = TRUE)
stopifnot(!anyDuplicated(raw[, .(k, m, parameter)]))
stopifnot("mi_raw lacks the B6 cap columns (run 05/07 of the final pipeline)" = all(c("maxH_inla", "n_above_inla", "maxH_map", "n_above_map", "z_map_max", "z_map_alpha", "z_map_intercept") %in% names(raw)))   # before any computation or output
cat(sprintf("Fits pulled: %d (imputation x transition) from %d imputations, %d transitions\n", uniqueN(raw[, .(k, m)]), uniqueN(raw$m), uniqueN(raw$k)))

# ---- 2. screen on every fit ----------------------------------------------------------------------------------------------------------------------
a <- raw[parameter == "alpha", .(k, label, m, alpha = est, a_sd = sd, n_events, fit_error, z_map_alpha, z_map_intercept)]   # B6 columns repeat over a fit's rows
i <- raw[parameter == "(Intercept)", .(k, m, icpt = est, i_sd = sd)]
x <- merge(a, i, by = c("k", "m"), all.x = TRUE); x[, ratio := a_sd / (0.78 * alpha / sqrt(n_events))]
x[, reason := fifelse(fit_error %in% TRUE | is.na(alpha) | is.na(a_sd) | is.na(icpt), "error/NA",
              fifelse(ratio < RATIO_MIN, "alpha_sd_ratio",
              fifelse(is.na(z_map_alpha) | is.na(z_map_intercept) | abs(z_map_alpha) > Z_MAP_MAX | abs(z_map_intercept) > Z_MAP_MAX, "z_map", "ok")))]
x[, valid := reason == "ok"]

# ---- 3. first TARGET valid fits per transition, by increasing m -------------------------------------------------------------------------------------------
setorder(x, k, m)
x[, rank_valid := ifelse(valid, cumsum(valid), NA_integer_), by = k]
x[, used := valid & rank_valid <= TARGET]
summ <- x[, .(label = label[1], candidates = .N, valid = sum(valid), used = sum(used), max_m_used = max(c(0L, m[used]))), by = k]
cat(sprintf("\nTransitions whose %d used fits include replacements (max m used > 30): %d\n", TARGET, summ[max_m_used > 30, .N])); print(summ[max_m_used > 30][order(-max_m_used)])
cat("\nInvalid fits by transition and reason:\n"); inv <- x[valid == FALSE, .N, by = .(label, reason)]
if (nrow(inv)) print(dcast(inv, label ~ reason, value.var = "N", fill = 0)) else cat("  none\n")
fwrite(x[, .(k, label, m, alpha, ratio, icpt, z_map_alpha, z_map_intercept, reason, valid, used)], file.path(RESULTS_DIR, sprintf("fit_validity_%s.csv", RUN_TAG)))
if (any(summ$used < TARGET)) { print(summ[used < TARGET]); stop("some transition has fewer than ", TARGET, " valid fits; nothing pooled (run 07_redraw.R further)") }

# ---- 4. pool estimates and joint covariance over exactly the used fits -------------------------------------------------------------------------------------
fin <- raw[x[used == TRUE, .(k, m)], on = .(k, m)]
pf <- pool_all(fin); stopifnot(all(pf$m_eff == TARGET))
fwrite(pf, file.path(RESULTS_DIR, sprintf("pooled_%s_final.csv", RUN_TAG)))
cat(sprintf("\nPooled: %d rows (%d transitions x parameters), m_eff = %d everywhere\n", nrow(pf), uniqueN(pf$k), TARGET))

need <- x[used == TRUE, .(k, label, m)]; store <- list(); chk <- list(); no_cov <- list()
for (mm in sort(unique(need$m))) {
  f <- file.path(OUT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, mm)); if (!file.exists(f)) { no_cov[[length(no_cov) + 1L]] <- need[m == mm]; next }
  fits <- readRDS(f)$fits
  for (lb in need[m == mm, label]) { fx <- fits[[lb]]; cv <- fx$cov
    if (is.null(cv) || !is.null(cv$error)) { no_cov[[length(no_cov) + 1L]] <- need[m == mm & label == lb]; next }
    store[[lb]][[as.character(mm)]] <- cv
    sh <- as.data.frame(fx$summary_hyper); sd_a <- sh[grepl("alpha", rownames(sh)), "sd"][1]
    sd_i <- fx$summary_fixed[, "sd"]; names(sd_i) <- rownames(fx$summary_fixed); b <- intersect(setdiff(colnames(cv$cov), "alpha"), names(sd_i))
    chk[[length(chk) + 1L]] <- data.table(label = lb, m = mm, parameter = c(b, "alpha"), ratio_sampled_to_inla = c(sqrt(diag(cv$cov)[b]) / sd_i[b], sqrt(cv$cov["alpha", "alpha"]) / sd_a)) } }
if (length(no_cov)) { cat("!! used fits without a usable stored covariance:\n"); print(rbindlist(no_cov)[order(label, m)]) }
# ---- 5. B3 (reporting only): intercept on the YEARS scale, b0 - alpha*log(10), per used fit, then Rubin-pooled -----------------------------------------------
# b0 (decade scale, variant 0) and alpha = INLA means from mi_raw; variance = var(b0) + log(10)^2 var(alpha) - 2 log(10) cov(b0, alpha) from the fit's stored
# sampled covariance. Fits without a usable stored covariance have NA variance and drop out of the pooling (m_eff shows it). Pooling and screen stay on the decade scale.
L10 <- log(10)
iy <- merge(fin[parameter == "(Intercept)", .(k, label, m, b0 = est)], fin[parameter == "alpha", .(k, m, alpha = est)], by = c("k", "m"))
iy[, est_years := b0 - alpha * L10]
iy[, var_years := mapply(function(lb, mm) { cv <- store[[lb]][[as.character(mm)]]
  if (is.null(cv)) return(NA_real_); C <- cv$cov; C["(Intercept)", "(Intercept)"] + L10^2 * C["alpha", "alpha"] - 2 * L10 * C["(Intercept)", "alpha"] }, label, m)]
iy_pool <- iy[, { p <- pool_rubin(est_years, var_years); c(list(label = label[1]), as.list(p), list(se = sqrt(p$total_var))) }, by = k][order(k)]
fwrite(iy_pool, file.path(RESULTS_DIR, sprintf("intercept_years_%s.csv", RUN_TAG)))
cat(sprintf("\nIntercept on the years scale pooled for %d transitions -> intercept_years_%s.csv (m_eff < %d for %d)\n", nrow(iy_pool), RUN_TAG, TARGET, iy_pool[m_eff < TARGET, .N]))

# ---- 6. B6 cap check: per transition, over the used fits (diagnostic, not a validity gate) --------------------------------------------------------------------
CAP <- -log(0.001 * 2^-23)    # 22.85, as in 05a
cap_check <- fin[parameter == "alpha", .(label = label[1], n_used = .N, n_map = sum(!is.na(maxH_map)),
                 share_maxH_map_gt_cap = mean(maxH_map > CAP, na.rm = TRUE), share_maxH_inla_gt_cap = mean(maxH_inla > CAP, na.rm = TRUE),
                 median_maxH_map = median(maxH_map, na.rm = TRUE), median_z_map_max = median(z_map_max, na.rm = TRUE),
                 median_z_map_alpha = median(z_map_alpha, na.rm = TRUE), median_z_map_intercept = median(z_map_intercept, na.rm = TRUE)), by = k][order(k)]
fwrite(cap_check, file.path(RESULTS_DIR, sprintf("cap_check_%s.csv", RUN_TAG)))
cat(sprintf("Cap check -> cap_check_%s.csv: %d transitions with share_maxH_map_gt_cap > 0\n", RUN_TAG, cap_check[share_maxH_map_gt_cap > 0, .N])); print(cap_check[share_maxH_map_gt_cap > 0])

cov_list <- list(); for (lb in unique(need$label)) if (length(store[[lb]]) == TARGET) cov_list[[lb]] <- pool_cov(unname(store[[lb]]))
saveRDS(cov_list, file.path(RESULTS_DIR, sprintf("pooled_%s_final_cov.rds", RUN_TAG)))
chk <- rbindlist(chk); fwrite(chk, file.path(RESULTS_DIR, sprintf("cov_check_%s.csv", RUN_TAG)))
cat(sprintf("\nJoint covariance pooled for %d of %d transitions. Sampled SD / INLA SD per fit (1 = agrees):\n", length(cov_list), uniqueN(need$label)))
print(chk[, .(n = .N, min = round(min(ratio_sampled_to_inla), 3), p01 = round(quantile(ratio_sampled_to_inla, .01), 3), median = round(median(ratio_sampled_to_inla), 3), max = round(max(ratio_sampled_to_inla), 3)),
          by = .(type = fifelse(parameter == "alpha", "alpha", "intercept/betas"))])
