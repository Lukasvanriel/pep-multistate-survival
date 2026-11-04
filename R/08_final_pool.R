# ==============================================================================
# 08_final_pool.R   (no INLA needed; seconds to run; run from the project root: Rscript R/08_final_pool.R)
# Purpose: final pooled results for ALL kept transitions with the SAME rule for every transition.
#   1. Pull every regular-pipeline fit: mi_raw_<tag>.csv (m1-30, 05_fit_inla.R) and mi_raw_redraw_<tag>.csv (candidates, 07_redraw.R).
#   2. Apply the frozen screen (references from screen_ref_<tag>.csv, written by 07_redraw.R before any candidate was fitted):
#      VALID if no INLA error/NA; alpha-SD ratio sd/(0.78*alpha/sqrt(events)) >= 0.1; |intercept - ref| / ref_sd <= 4.
#   3. Use the first TARGET (30) valid fits per transition in increasing m (selection by index only, never by results). Stop if a transition has fewer.
#   4. Rubin-pool estimates (pool_all) and the joint covariance (Rubin on the stored per-fit covariances, as sampled by inla.posterior.sample;
#      no post-hoc correction of alpha's spread: the correction was tried for the "int" run and rejected, see PIPELINE_NOTES 2026-10-20).
# Reads FIT_DIR (data-output/<tag>/fits/); outputs in RESULTS_DIR (data-output/<tag>/results/): fit_validity_<tag>.csv (every fit: screen result, used yes/no), pooled_<tag>_final.csv,
#   pooled_<tag>_final_cov.rds (named list by transition), cov_check_<tag>.csv (sampled SD / INLA SD per fit and parameter).
# ==============================================================================
suppressPackageStartupMessages(library(data.table)); options(width = 200)
source(here::here("R", "config.R"))   # RUN_TAG, FIT_DIR, RESULTS_DIR
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
TARGET <- 30L; RATIO_MIN <- 0.1; Z_MAX <- 4
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
cat(sprintf("Fits pulled: %d (imputation x transition) from %d imputations, %d transitions\n", uniqueN(raw[, .(k, m)]), uniqueN(raw$m), uniqueN(raw$k)))

# ---- 2. frozen screen on every fit ----------------------------------------------------------------------------------------------------------------------
ref <- fread(file.path(OUT_DIR, sprintf("screen_ref_%s.csv", RUN_TAG)))                              # frozen by 07_redraw.R
a <- raw[parameter == "alpha", .(k, label, m, alpha = est, a_sd = sd, n_events, fit_error)]
i <- raw[parameter == "(Intercept)", .(k, m, icpt = est, i_sd = sd)]
x <- merge(a, i, by = c("k", "m"), all.x = TRUE); x[, ratio := a_sd / (0.78 * alpha / sqrt(n_events))]
x <- merge(x, ref, by = "k", all.x = TRUE); x[, z_icpt := (icpt - icpt_ref) / icpt_sd_ref]
x[, reason := fifelse(fit_error %in% TRUE | is.na(alpha) | is.na(a_sd) | is.na(icpt), "error/NA",
              fifelse(ratio < RATIO_MIN, "alpha_sd_ratio", fifelse(is.na(z_icpt) | abs(z_icpt) > Z_MAX, "intercept", "ok")))]
x[, valid := reason == "ok"]

# ---- 3. first TARGET valid fits per transition, by increasing m -------------------------------------------------------------------------------------------
setorder(x, k, m)
x[, rank_valid := ifelse(valid, cumsum(valid), NA_integer_), by = k]
x[, used := valid & rank_valid <= TARGET]
summ <- x[, .(label = label[1], candidates = .N, valid = sum(valid), used = sum(used), max_m_used = max(c(0L, m[used]))), by = k]
cat(sprintf("\nTransitions whose %d used fits include replacements (max m used > 30): %d\n", TARGET, summ[max_m_used > 30, .N])); print(summ[max_m_used > 30][order(-max_m_used)])
cat("\nInvalid fits by transition and reason:\n"); print(dcast(x[valid == FALSE, .N, by = .(label, reason)], label ~ reason, value.var = "N", fill = 0))
fwrite(x[, .(k, label, m, alpha, ratio, icpt, z_icpt, reason, valid, used)], file.path(RESULTS_DIR, sprintf("fit_validity_%s.csv", RUN_TAG)))
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
cov_list <- list(); for (lb in unique(need$label)) if (length(store[[lb]]) == TARGET) cov_list[[lb]] <- pool_cov(unname(store[[lb]]))
saveRDS(cov_list, file.path(RESULTS_DIR, sprintf("pooled_%s_final_cov.rds", RUN_TAG)))
chk <- rbindlist(chk); fwrite(chk, file.path(RESULTS_DIR, sprintf("cov_check_%s.csv", RUN_TAG)))
cat(sprintf("\nJoint covariance pooled for %d of %d transitions. Sampled SD / INLA SD per fit (1 = agrees):\n", length(cov_list), uniqueN(need$label)))
print(chk[, .(n = .N, min = round(min(ratio_sampled_to_inla), 3), p01 = round(quantile(ratio_sampled_to_inla, .01), 3), median = round(median(ratio_sampled_to_inla), 3), max = round(max(ratio_sampled_to_inla), 3)),
          by = .(type = fifelse(parameter == "alpha", "alpha", "intercept/betas"))])
