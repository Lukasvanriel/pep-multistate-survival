# ==============================================================================
# 09_hazard_ratios.R   (no INLA; run from the project root: Rscript R/09_hazard_ratios.R)   needs 08_final_pool.R to have been run
# Hazard ratios for every kept transition (WP7, D1-D4, D8). Model (weibullsurv variant 0): h(t|x) = alpha t^(alpha-1) exp(eta), eta = Intercept + sum_j beta_j x_j.
#   => for two covariate profiles x1, x0 the log hazard ratio is  c'beta  with c = x1 - x0 (same at every t: proportional hazards).  HR = exp(c'beta); alpha does not enter.
# Per imputation m (one of the 30 used fits of the transition):  q_m = c'beta_m (INLA posterior means, as in mi_raw),  u_m = c' V_bb,m c,
#   V_bb,m = the beta block of the stored joint covariance of that fit (the rows/columns of the intercept and alpha are not used).
# Rubin on the LOG scale, Barnard-Rubin df with complete-data df = events - parameters, t interval; then exp().
# Covariate coding: z_gdd_3, z_cmi_gs, z_soil_pc1 are z-scores (effects per +1 SD); disturbance flags are 0/1 (vs undisturbed); tsd_* are f(TSD in years) with f from TSD_FORM
#   (log(TSD + 1) in the main build; 0 where undisturbed). The disturbance effect at GDD z and time since disturbance tau (years) is
#   log HR(d vs none | z, tau) = beta_d + beta_gddxd * z + beta_tsd * f(tau).
# Contrasts (a term that is not in the transition's model gives no row):
#   per_sd_<cov>                    +1 SD of GDD / CMI / soil, undisturbed plots (no per_sd_gdd_3 where z_gdd_3 is not in the model, A11)
#   <d>_vs_none z=<z> tau=<t>y      disturbance d (fire, pest, harv_light, harv_heavy) vs undisturbed, z in Z_GRID, tau in TAU_GRID_Y (years). Where z_gdd_3 is not in the model
#                                   (A11) only z = 0 (D4 point 5, D8). A level with only its main term (A6; no GDD x, no TSD term) gives ONE row "<d>_vs_none" (main term alone).
#   <d>_gdd_per_sd                  +1 SD GDD within disturbed plots (beta_gdd + beta_gddxd; beta_gddxd alone where z_gdd_3 is not in the model)
#   <d>_gdd_interaction             ratio of the +1 SD GDD hazard ratio, disturbed vs undisturbed (beta_gddxd only; not where z_gdd_3 is not in the model)
#   <d>_tsd_coef                    the TSD coefficient itself (per unit of f(TSD)); "hr" is exp of it and has no direct reading; input of 12's block D
#   <d>_tsd_1_vs_10y, <d>_tsd_10_vs_30y   log HR at 10 y minus log HR at 1 y (30 y minus 10 y): c * (f(10) - f(1)), c * (f(30) - f(10)). Both are fixed multiples of the TSD
#                                   coefficient (same sign and p-value as _tsd_coef): not independent evidence.
# One-sided bounds (D4): for a "<d>_vs_none" row whose reference has 0 events (harvest: n_harv0 == 0; fire / pest: n_events - n_fire / n_pest == 0) `bound` = "lower" and
#   hr_bound = exp(log_hr - qt(0.95, df) se); for a level with 0 events (n_dist == 0) `bound` = "upper" and hr_bound = exp(log_hr + qt(0.95, df) se). In these rows hr_lo, hr_hi and p are NA.
# Ranges (D2): tsd_risk_* = min / 5th / 95th percentile / max of the raw TSD (years, end-visit value of the interval) over the intervals that start in the origin state with the
#   level's flag on at the end visit (at risk); tsd_ev_* = the same over the intervals that end in the destination state (events). Repeated on every row of that transition x level.
# Output: RESULTS_DIR (data-output/<tag>/results/)/hazard_ratios_<tag>.csv (one row per transition x contrast) and alpha_report_<tag>.csv (D3, one row per transition).
#   n_dist = disturbed events of that level (transition_set); low_events = "<15" (below the event rule on its own: n_dist < MIN_EVENTS; the level is in the model only because light + heavy
#   harvest together reach MIN_EVENTS) or "15-19" (n_dist <= MIN_EVENTS + 4), else NA.
# ==============================================================================
suppressPackageStartupMessages(library(data.table)); options(width = 200)
source(here::here("R", "config.R"))   # RUN_TAG, TABLE_DIR, FIT_DIR, RESULTS_DIR, MIN_EVENTS, TSD_FORM
source(here::here("R", "covariate_spec.R"))   # tsd_fun()
Z_GRID <- c(-1, 0, 1); TAU_GRID_Y <- c(1, 10, 30, 60); LOW_EVENTS_MAX <- MIN_EVENTS + 4L
f_tsd <- tsd_fun(TSD_FORM)
OUT_DIR <- FIT_DIR
DIST <- list(fire       = c(main = "fire_f",     gdd = "gdd_x_fire",       tsd = "tsd_fire",       n = "n_fire"),
             pest       = c(main = "pest_f",     gdd = "gdd_x_pest",       tsd = "tsd_pest",       n = "n_pest"),
             harv_light = c(main = "harv_light", gdd = "gdd_x_harv_light", tsd = "tsd_harv_light", n = "n_harv1"),
             harv_heavy = c(main = "harv_heavy", gdd = "gdd_x_harv_heavy", tsd = "tsd_harv_heavy", n = "n_harv2"))

val <- fread(file.path(RESULTS_DIR, sprintf("fit_validity_%s.csv", RUN_TAG)))[used == TRUE, .(k, label, m)]
rf <- file.path(OUT_DIR, sprintf(c("mi_raw_%s.csv", "mi_raw_redraw_%s.csv"), RUN_TAG)); raw <- rbindlist(lapply(rf[file.exists(rf)], fread), fill = TRUE)
raw <- raw[val, on = .(k, m)]; stopifnot(nrow(val) == 90 * 30L)
ts <- readRDS(file.path(TABLE_DIR, "transition_set.rds")); setDT(ts); ts[, label := paste(from, to, sep = "->")]

# ---- contrasts as named coefficient vectors, restricted to the terms present in a transition ------------------------------------------------------------------
make_contrasts <- function(terms) {
  has <- function(...) all(c(...) %in% terms); L <- list()
  has_g <- "z_gdd_3" %in% terms
  for (v in c("z_gdd_3", "z_cmi_gs", "z_soil_pc1")) if (has(v)) L[[paste0("per_sd_", sub("^z_", "", v))]] <- setNames(1, v)
  for (d in names(DIST)) { s <- DIST[[d]]; if (!has(s[["main"]])) next
    if (!has(s[["gdd"]]) && !has(s[["tsd"]])) { L[[paste0(d, "_vs_none")]] <- setNames(1, s[["main"]]); next }   # level with the main term only (A6, D4 point 4)
    stopifnot("GDD x and TSD terms of a disturbance level must come together" = has(s[["gdd"]], s[["tsd"]]))
    for (z in if (has_g) Z_GRID else 0) for (tau in TAU_GRID_Y) L[[sprintf("%s_vs_none z=%g tau=%gy", d, z, tau)]] <- setNames(c(1, z, f_tsd(tau)), s[c("main", "gdd", "tsd")])
    L[[paste0(d, "_gdd_per_sd")]] <- if (has_g) setNames(c(1, 1), c("z_gdd_3", s[["gdd"]])) else setNames(1, s[["gdd"]])
    if (has_g) L[[paste0(d, "_gdd_interaction")]] <- setNames(1, s[["gdd"]])
    L[[paste0(d, "_tsd_coef")]]       <- setNames(1, s[["tsd"]])
    L[[paste0(d, "_tsd_1_vs_10y")]]   <- setNames(f_tsd(10) - f_tsd(1), s[["tsd"]])
    L[[paste0(d, "_tsd_10_vs_30y")]]  <- setNames(f_tsd(30) - f_tsd(10), s[["tsd"]]) }
  L }
rubin <- function(q, u, nu_com) {                                              # q, u: per-imputation log HR and variance
  M <- length(q); Qb <- mean(q); W <- mean(u); B <- var(q); Tt <- W + (1 + 1 / M) * B
  r <- if (W > 0) (1 + 1 / M) * B / W else NA_real_; lam <- r / (1 + r)
  nu_old <- if (is.finite(lam) && lam > 0) (M - 1) / lam^2 else Inf; nu_obs <- (nu_com + 1) / (nu_com + 3) * nu_com * (1 - lam)
  df <- if (is.finite(nu_old)) 1 / (1 / nu_old + 1 / nu_obs) else nu_obs
  se <- sqrt(Tt); list(log_hr = Qb, se = se, df = df, fmi = (r + 2 / (df + 3)) / (r + 1), M = M) }

# ---- per fit: contrast estimates and variances (beta block of the stored covariance) -----------------------------------------------------------------------------
res <- list()
for (mm in sort(unique(val$m))) {
  fits <- readRDS(file.path(OUT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, mm)))$fits
  for (lb in val[m == mm, label]) {
    cv <- fits[[lb]]$cov; stopifnot(!is.null(cv), is.null(cv$error)); r <- raw[label == lb & m == mm]
    est <- setNames(r$est, r$parameter); V <- cv$cov; CC <- make_contrasts(setdiff(colnames(V), c("(Intercept)", "alpha")))
    miss <- setdiff(unlist(lapply(CC, names)), names(est))
    if (length(miss)) stop(sprintf("09: %s (m = %d): contrast term(s) missing from the fit's estimates: %s", lb, mm, paste(unique(miss), collapse = ", ")), call. = FALSE)   # no silent NA
    for (cn in names(CC)) { cc <- CC[[cn]]
      res[[length(res) + 1L]] <- data.table(label = lb, m = mm, contrast = cn, q = sum(cc * est[names(cc)]), u = as.numeric(t(cc) %*% V[names(cc), names(cc), drop = FALSE] %*% cc)) } } }
res <- rbindlist(res)
nu <- raw[, .(nu_com = max(mean(n_events[parameter == "alpha"]) - uniqueN(parameter), 1)), by = label]
out <- merge(res[, c(rubin(q, u, nu$nu_com[match(label[1], nu$label)])), by = .(label, contrast)], nu[, .(label)], by = "label")
out[, `:=`(hr = exp(log_hr), hr_lo = exp(log_hr - qt(0.975, df) * se), hr_hi = exp(log_hr + qt(0.975, df) * se), p = 2 * pt(-abs(log_hr / se), df))]

# ---- disturbed-event counts, low_events label, one-sided bounds (D4) ----------------------------------------------------------------------------------------------
out[, dist := fifelse(grepl("^(fire|pest|harv_light|harv_heavy)_", contrast), sub("_(vs_none.*|gdd_per_sd|gdd_interaction|tsd_.*)$", "", contrast), NA_character_)]
out <- merge(out, ts[, .(label, n_events, n_fire, n_pest, n_harv0, n_harv1, n_harv2)], by = "label", all.x = TRUE)
out[, n_dist := as.integer(fcase(dist == "fire", n_fire, dist == "pest", n_pest, dist == "harv_light", n_harv1, dist == "harv_heavy", n_harv2))]
out[, n_ref  := as.integer(fcase(dist == "fire", n_events - n_fire, dist == "pest", n_events - n_pest, dist %in% c("harv_light", "harv_heavy"), n_harv0))]
out[, low_events := fifelse(is.na(n_dist), NA_character_, fifelse(n_dist < MIN_EVENTS, sprintf("<%d", MIN_EVENTS), fifelse(n_dist <= LOW_EVENTS_MAX, sprintf("%d-%d", MIN_EVENTS, LOW_EVENTS_MAX), NA_character_)))]
vs <- grepl("_vs_none", out$contrast) & !is.na(out$dist)
stopifnot("a row whose reference and level both have 0 events" = !any(vs & out$n_ref == 0L & out$n_dist == 0L, na.rm = TRUE))
out[, bound := fifelse(vs & n_ref == 0L, "lower", fifelse(vs & n_dist == 0L, "upper", NA_character_))]
out[, hr_bound := fcase(bound == "lower", exp(log_hr - qt(0.95, df) * se), bound == "upper", exp(log_hr + qt(0.95, df) * se), default = NA_real_)]
out[!is.na(bound), `:=`(hr_lo = NA_real_, hr_hi = NA_real_, p = NA_real_)]

# ---- TSD ranges (D2): raw TSD in years over the at-risk and the event intervals of the origin state -------------------------------------------------------------
iv <- readRDS(file.path(TABLE_DIR, "visit_intervals.rds")); setDT(iv)
q4 <- function(v) if (!length(v)) setNames(rep(NA_real_, 4), c("min", "p05", "p95", "max")) else c(min = min(v), p05 = unname(quantile(v, 0.05)), p95 = unname(quantile(v, 0.95)), max = max(v))
rng <- rbindlist(lapply(which(ts$kept), function(i) {
  a <- iv[state_start == ts$from[i]]
  rbindlist(lapply(names(DIST), function(d) {
    s <- switch(d, fire = a[fire == 1L], pest = a[pest == 1L], harv_light = a[harv == 1L], harv_heavy = a[harv == 2L])
    cl <- if (d == "fire") "tsd_fire_y" else if (d == "pest") "tsd_pest_y" else "tsd_harv_y"; e <- s[state_end == ts$to[i] & is_change == TRUE]
    r <- q4(s[[cl]]); v <- q4(e[[cl]])
    data.table(label = ts$label[i], dist = d, n_risk = nrow(s), tsd_risk_min = r[["min"]], tsd_risk_p05 = r[["p05"]], tsd_risk_p95 = r[["p95"]], tsd_risk_max = r[["max"]],
               n_ev_int = nrow(e), tsd_ev_min = v[["min"]], tsd_ev_p05 = v[["p05"]], tsd_ev_p95 = v[["p95"]], tsd_ev_max = v[["max"]]) })) }))
out <- merge(out, rng, by = c("label", "dist"), all.x = TRUE)
out[, c("n_fire", "n_pest", "n_harv0", "n_harv1", "n_harv2", "n_ref") := NULL]
setcolorder(out, c("label", "contrast", "dist", "hr", "hr_lo", "hr_hi", "p", "log_hr", "se", "fmi", "df", "M", "n_events", "n_dist", "low_events", "bound", "hr_bound"))
fwrite(out[order(label, contrast)], file.path(RESULTS_DIR, sprintf("hazard_ratios_%s.csv", RUN_TAG)), na = "NA")   # na = "NA": fread reads the character columns bound / low_events with real NA
cat(sprintf("Hazard ratios: %d rows, %d transitions, %d contrast types; M = %s\n", nrow(out), uniqueN(out$label), uniqueN(sub(" z=.*", "", out$contrast)), paste(unique(out$M), collapse = "/")))

# ---- D3: alpha report ----------------------------------------------------------------------------------------------------------------------------------------------
pf <- fread(file.path(RESULTS_DIR, sprintf("pooled_%s_final.csv", RUN_TAG)))
cf <- file.path(RESULTS_DIR, sprintf("cap_check_%s.csv", RUN_TAG)); stopifnot("cap_check file missing (run 08)" = file.exists(cf)); cap <- fread(cf)
ar <- pf[parameter == "alpha", .(k, label, alpha = pooled_est, alpha_se = se, alpha_log_pooled = est_from_log, alpha_ci_lo_log = ci_lo_log, alpha_ci_hi_log = ci_hi_log)]
ar <- merge(ar, ts[, .(label, events = n_events)], by = "label")
ar[, class := fifelse(sub("->.*", "", label) == "regen", "regen_exit", fifelse(sub(".*->", "", label) == "regen", "to_regen", "between"))][, starts_at_disturbance := class == "regen_exit"]
tsdw <- dcast(pf[parameter %in% c("tsd_fire", "tsd_pest", "tsd_harv_light", "tsd_harv_heavy"), .(label, parameter, est = pooled_est, se)], label ~ parameter, value.var = c("est", "se"))
setnames(tsdw, names(tsdw), sub("^est_(.*)$", "\\1", sub("^se_(.*)$", "\\1_se", names(tsdw))))
ar <- merge(ar, tsdw, by = "label", all.x = TRUE)
ar <- merge(ar, cap[, .(label, n_used, share_maxH_map_gt_cap, share_maxH_inla_gt_cap, cap_flag = share_maxH_map_gt_cap > 0, median_maxH_map, median_z_map_max, median_z_map_alpha, median_z_map_intercept)], by = "label", all.x = TRUE)
setorder(ar, k); fwrite(ar, file.path(RESULTS_DIR, sprintf("alpha_report_%s.csv", RUN_TAG)))
cat(sprintf("\nAlpha report: %d rows (%d regen_exit = starts at the disturbance, %d to_regen, %d between); cap_flag TRUE in %d\n", nrow(ar), ar[class == "regen_exit", .N], ar[class == "to_regen", .N], ar[class == "between", .N], ar[cap_flag == TRUE, .N]))
chk_a <- merge(ar[, .(label, alpha)], pf[parameter == "alpha", .(label, pooled_est)], by = "label")   # alpha_report's pooled alpha vs the alpha row of pooled_<tag>_final.csv, by label
stopifnot(nrow(ar) == 90L, nrow(chk_a) == 90L, all(abs(chk_a$alpha - chk_a$pooled_est) < 1e-12))

# ---- checks and prints ---------------------------------------------------------------------------------------------------------------------------------------------
cat("\nD1 check: pooled log HR per SD CMI vs pooled beta of z_cmi_gs (same estimates: identical up to floating point):\n")
ck <- merge(out[contrast == "per_sd_cmi_gs", .(label, log_hr, se)], pf[parameter == "z_cmi_gs", .(label, pooled_est, se_inla = se)], by = "label")
cat(sprintf("  max |difference| = %.2e over %d transitions; SE (sampled covariance) / SE (INLA sd): median %.3f, range %.3f-%.3f\n", max(abs(ck$log_hr - ck$pooled_est)), nrow(ck), median(ck$se / ck$se_inla), min(ck$se / ck$se_inla), max(ck$se / ck$se_inla)))
cat("\nD2 check: transition x level pairs with the four tau rows (z = 0) and filled range columns:\n")
tt <- out[grepl("_vs_none z=0 tau=", contrast), .(n_tau = .N, risk_filled = all(!is.na(tsd_risk_max)), ev_filled = all(!is.na(tsd_ev_max))), by = .(label, dist)]
print(tt[, .(pairs = .N, all_four_tau = sum(n_tau == length(TAU_GRID_Y)), risk_range_filled = sum(risk_filled), event_range_filled = sum(ev_filled)), by = dist])
stopifnot(all(tt$n_tau == length(TAU_GRID_Y)), all(tt$risk_filled), all(tt$ev_filled))
cat("\nD4: rows with one-sided bounds (expected on the main build: regen->sugar_maple lower, other_temperate->yellow_birch upper):\n")
print(out[!is.na(bound), .(label, contrast, bound, hr = round(hr, 2), hr_bound = round(hr_bound, 2), n_dist)][order(label, contrast)], nrows = 100)
if (uniqueN(out[!is.na(bound), label]) != 2L) warning("the number of transitions with a one-sided bound is not 2")
cat("\nD8: transitions without a GDD main term (no per_sd_gdd_3 row):\n")
nog <- setdiff(out$label, out[contrast == "per_sd_gdd_3", label]); print(nog)
for (lb in nog) {
  bad <- out[label == lb & grepl("_gdd_interaction$|_vs_none z=(-1|1) ", contrast), .N]
  d8 <- out[label == lb & grepl("_gdd_per_sd$", contrast), .(contrast, log_hr)][, parameter := paste0("gdd_x_", sub("_gdd_per_sd$", "", contrast))][pf[label == lb, .(parameter, pooled_est)], on = "parameter", nomatch = 0]
  cat(sprintf("  %s: excluded rows present = %d; gdd_per_sd == pooled gdd_x_ coefficient:\n", lb, bad)); print(d8[, .(contrast, log_hr = round(log_hr, 6), pooled_gdd_x = round(pooled_est, 6), equal = abs(log_hr - pooled_est) < 1e-9)])
  stopifnot(bad == 0L, all(abs(d8$log_hr - d8$pooled_est) < 1e-9)) }
cat("\nblack_spruce->regen (reference, the transition with the most events):\n")
print(out[label == "black_spruce->regen" & (!grepl("z=", contrast) | grepl("z=0 tau=(1|10)y$", contrast)), .(contrast, hr = round(hr, 2), lo = round(hr_lo, 2), hi = round(hr_hi, 2), p = signif(p, 2), fmi = round(fmi, 2))][order(contrast)])
