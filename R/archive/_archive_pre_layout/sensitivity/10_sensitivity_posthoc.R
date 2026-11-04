# ==============================================================================
# 10_sensitivity_posthoc.R   (no INLA; run from the project root: Rscript R/sensitivity/10_sensitivity_posthoc.R)   needs 08_final_pool.R and 09_hazard_ratios.R run first
# Purpose: how much do the main results depend on analysis choices that can be changed WITHOUT new INLA fits. Main analysis = first 30 valid fits per transition (pooled_int_final.csv).
#   S1  Interactions with enough disturbed events. Main analysis keeps every GDD x disturbance interaction whose disturbance term exists (>= 15 disturbed events); interactions with 15-19
#       disturbed events are lower confidence (imputation attenuates interactions toward 0). Compare: all interactions vs those with >= MIN_EVENTS_S1 (20) disturbed events -
#       number, share with CI excluding 1, sign agreement and the median ratio-of-HR. Answers: do conclusions about climate x disturbance rest on the sparse interactions?
#   S2  Drop fits with a very small alpha SD RELATIVE to their own transition (alpha-SD ratio < REL_RATIO_S2 = 0.5 x the median ratio of the transition's used fits; they pass the frozen screen but are
#       suspicious; absolute ratios are naturally < 1 in sparse transitions, so a fixed cut-off would drop a third of all fits). Re-pool every transition without them; shift in SE units per parameter.
#   S3  Monte Carlo error of using 30 imputations: pool ALL valid fits (up to ~128 imputations per transition) instead of the first 30; shift in SE units, SE ratio.
# Outputs (data/processed/fit_<tag>/): sens_pooled_compare_<tag>.csv (S2, S3: one row per scenario x transition x parameter), sens_interactions_int.csv (S1: one row per interaction), summary printed.
# NOT covered here (need new INLA runs, see PIPELINE_NOTES 2026-10-20): TSD upper bound (imputed_prev), start-visit flags/TSD.
# ==============================================================================
suppressPackageStartupMessages(library(data.table)); options(width = 200)
source(here::here("R", "config.R"))   # RUN_TAG, TABLE_DIR, FIT_DIR (T0)
MIN_EVENTS_S1 <- 20L; REL_RATIO_S2 <- 0.5; MIN_FITS_S2 <- 20L
OUT_DIR <- FIT_DIR; SD_MIN <- 1e-3
ex <- parse(file.path("R", "05a_fit_funs.R"))
for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]])[1] %in% c("pool_rubin", "flag_fits", "pool_all")) eval(e)
rf <- file.path(OUT_DIR, sprintf(c("mi_raw_%s.csv", "mi_raw_redraw_%s.csv", "mi_raw_sweepfull_%s.csv"), RUN_TAG)); raw <- rbindlist(lapply(rf[file.exists(rf)], fread), fill = TRUE)   # sweepfull: int run only
v <- fread(file.path(OUT_DIR, sprintf("fit_validity_%s.csv", RUN_TAG))); pf <- fread(file.path(OUT_DIR, sprintf("pooled_%s_final.csv", RUN_TAG)))
key <- function(d) d[, .(k, label, parameter, est = pooled_est, se, m_eff)]
base <- key(pf)
compare <- function(alt, name) { a <- key(alt); x <- merge(base, a, by = c("k", "label", "parameter"), suffixes = c("", "_alt"))
  x[, `:=`(scenario = name, shift_se = (est_alt - est) / se, se_ratio = se_alt / se)][, .(scenario, k, label, parameter, est, se, m_eff, est_alt, se_alt, m_eff_alt, shift_se, se_ratio)] }

# ---- S2: without low-ratio fits ----------------------------------------------------------------------------------------------------------------------------------
used <- v[used == TRUE, .(k, m)]; v[used == TRUE, rel_ratio := ratio / median(ratio), by = k]; low <- v[used == TRUE & rel_ratio < REL_RATIO_S2, .(k, m)]
cat(sprintf("S2: used fits with alpha-SD ratio < %.2f x their transition's median: %d of %d, in %d transitions\n", REL_RATIO_S2, nrow(low), nrow(used), uniqueN(low$k)))
keep <- used[!low, on = .(k, m)]; cnt <- keep[, .N, by = k]; ok_k <- cnt[N >= MIN_FITS_S2, k]
if (nrow(cnt[N < MIN_FITS_S2])) cat("  transitions with < ", MIN_FITS_S2, " fits left (skipped): ", paste(v[k %in% cnt[N < MIN_FITS_S2, k], unique(label)], collapse = ", "), "\n")
s2 <- compare(pool_all(raw[keep[k %in% ok_k], on = .(k, m)]), "S2_drop_low_ratio")
# ---- S3: all valid fits instead of the first 30 -----------------------------------------------------------------------------------------------------------------
allv <- v[valid == TRUE, .(k, m)]; cat(sprintf("S3: all valid fits: %d (per transition min %d, median %d, max %d)\n", nrow(allv), allv[, .N, by = k][, min(N)], allv[, .N, by = k][, as.integer(median(N))], allv[, .N, by = k][, max(N)]))
s3 <- compare(pool_all(raw[allv, on = .(k, m)]), "S3_all_valid_fits")
sens <- rbind(s2, s3); fwrite(sens, file.path(OUT_DIR, sprintf("sens_pooled_compare_%s.csv", RUN_TAG)))
cat("\nShift of the pooled estimate in units of the main SE, and SE ratio (alt / main), over all transitions x parameters:\n")
print(sens[, .(rows = as.numeric(.N), m_eff_alt_median = as.numeric(median(m_eff_alt)), shift_med_abs = round(median(abs(shift_se)), 3), shift_p95_abs = round(quantile(abs(shift_se), .95), 3), shift_max_abs = round(max(abs(shift_se)), 3),
               n_shift_gt_0.5SE = as.numeric(sum(abs(shift_se) > 0.5)), se_ratio_med = round(median(se_ratio), 3), se_ratio_min = round(min(se_ratio), 3), se_ratio_max = round(max(se_ratio), 3)), by = scenario])
cat("\nLargest shifts (|shift| in SE units):\n"); print(sens[order(-abs(shift_se))][1:12, .(scenario, label, parameter, est = round(est, 3), est_alt = round(est_alt, 3), se = round(se, 3), shift_se = round(shift_se, 2), m_eff_alt)])
cat("\nBy parameter type (alpha / intercept / beta), share of rows with |shift| > 0.5 SE:\n")
print(sens[, .(rows = as.numeric(.N), gt0.5 = as.numeric(sum(abs(shift_se) > 0.5)), max = round(max(abs(shift_se)), 2)), by = .(scenario, type = fifelse(parameter == "alpha", "alpha", fifelse(parameter == "(Intercept)", "intercept", "beta")))])

# ---- S1: interactions with enough disturbed events ---------------------------------------------------------------------------------------------------------------
hr <- fread(file.path(OUT_DIR, sprintf("hazard_ratios_%s.csv", RUN_TAG)))[grepl("_gdd_interaction$", contrast)]
hr[, `:=`(excl1 = hr_lo > 1 | hr_hi < 1, subset = n_dist >= MIN_EVENTS_S1)]
fwrite(hr[, .(label, contrast, dist, n_events, n_dist, subset, low_events, hr, hr_lo, hr_hi, p, fmi, alpha_cov_low_share)][order(n_dist)], file.path(OUT_DIR, sprintf("sens_interactions_%s.csv", RUN_TAG)))
summ <- function(d) d[, .(n = .N, ci_excl_1 = sum(excl1), share = round(mean(excl1), 3), hr_lt_1 = sum(hr < 1), median_hr = round(median(hr), 3), median_log_hr = round(median(log_hr), 3))]
cat(sprintf("\nS1: GDD x disturbance interactions (ratio of the +1 SD GDD hazard ratio, disturbed vs undisturbed): %d in the main analysis; %d with >= %d disturbed events; %d with %d-%d (lower confidence)\n",
            nrow(hr), hr[subset == TRUE, .N], MIN_EVENTS_S1, hr[low_events == TRUE, .N], 15L, MIN_EVENTS_S1 - 1L))
print(rbind(data.table(set = "all", summ(hr)), data.table(set = sprintf(">= %d disturbed events", MIN_EVENTS_S1), summ(hr[subset == TRUE])), data.table(set = sprintf("15-%d (lower confidence)", MIN_EVENTS_S1 - 1L), summ(hr[subset == FALSE]))))
cat("\nBy disturbance type, share of interactions with CI excluding 1: all vs >= threshold:\n")
print(hr[, .(n_all = .N, excl_all = sum(excl1), n_sub = sum(subset), excl_sub = sum(excl1 & subset), median_hr_all = round(median(hr), 2), median_hr_sub = round(median(hr[subset]), 2)), by = dist])
cat("\nInteractions with < ", MIN_EVENTS_S1, " disturbed events (the ones dropped in S1):\n"); print(hr[subset == FALSE][order(n_dist), .(label, contrast, n_dist, hr = round(hr, 2), lo = round(hr_lo, 2), hi = round(hr_hi, 2), p = signif(p, 2))])
