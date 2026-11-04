# ==============================================================================
# 09_hazard_ratios.R   (no INLA; run from the project root: Rscript R/09_hazard_ratios.R)   needs 08_final_pool.R to have been run
# Hazard ratios for every kept transition. Model (weibullsurv variant 1): h(t|x) = alpha t^(alpha-1) exp(alpha * eta), eta = Intercept + sum_j beta_j x_j.
#   => for two covariate profiles x1, x0 the log hazard ratio is  alpha * c'beta  with c = x1 - x0 (same at every t: proportional hazards).  HR = exp(alpha * c'beta), NOT exp(beta).
# Per imputation m (one of the 30 used fits of the transition):  est_m = alpha_m * c'beta_m   (INLA posterior means, as in mi_raw / pooled_int_final.csv)
#   var_m = g' V_m g  (delta method), g = (d/dbeta = alpha_m * c, d/dalpha = c'beta_m), V_m = stored joint covariance of (betas, alpha) of that fit.
# Rubin on the LOG scale (log HR is closer to normal): Qbar, Wbar, B, T = Wbar + (1+1/M) B, Barnard-Rubin df with complete-data df = events - parameters, t interval; then exp().
# Covariate coding: z_gdd_3, z_cmi_gs, z_soil_pc1 are z-scores (effects are per +1 SD); disturbance flags are 0/1 (vs undisturbed); tsd_* are RAW decades since the disturbance
#   (0 = at the disturbance, structurally 0 when undisturbed). Interactions are GDD x flag, so a disturbance effect depends on GDD (z) and on time since (tau, decades):
#   log HR(disturbance vs none | z, tau) = alpha * (beta_dist + beta_gddxdist * z + beta_tsd * tau).    Main effect = effect at mean GDD (z = 0) and tau = 0.
# Contrasts (constants below can be changed):
#   per_sd_<cov>              +1 SD of GDD / CMI / soil, undisturbed (reference) plots
#   <d>_vs_none z=.. tau=..   disturbance d (fire, pest, harv_light, harv_heavy, harv_any) vs undisturbed at GDD z and tau decades since
#   <d>_gdd_per_sd            +1 SD GDD within disturbed plots  (beta_gdd + beta_gddxd)
#   <d>_gdd_interaction       ratio of the +1 SD GDD hazard ratio, disturbed vs undisturbed (beta_gddxd only)
#   <d>_tsd_per_decade        +1 decade since disturbance, among disturbed plots (beta_tsd)
# Output: data/processed/fit_<tag>/hazard_ratios_<tag>.csv  (one row per transition x contrast)
#   columns: alpha_cov_low_share (share of the 30 used fits whose sampled alpha SD is < 0.9 x INLA's), n_dist (disturbed events of that type from transition_set.rds), low_events (n_dist 15-19 = lower confidence), alpha_cov_ratio_med (median over used fits of sampled alpha SD / INLA alpha SD;
#   alpha_cov_ratio_med < 0.9 or a high low_share means the sampled alpha spread, hence var_m of contrasts involving alpha, is somewhat too small; see PIPELINE_NOTES 2026-10-20).
# ==============================================================================
suppressPackageStartupMessages(library(data.table)); options(width = 200)
source(here::here("R", "config.R"))   # RUN_TAG, TABLE_DIR, FIT_DIR (T0)
Z_GRID <- c(-1, 0, 1); TAU_GRID <- c(0, 1, 3); LOW_EVENTS_MAX <- 19L
OUT_DIR <- FIT_DIR
DIST <- list(fire = c(main = "fire_f", gdd = "gdd_x_fire", tsd = "tsd_fire", n = "n_fire"),
             pest = c(main = "pest_f", gdd = "gdd_x_pest", tsd = "tsd_pest", n = "n_pest"),
             harv_light = c(main = "harv_light", gdd = "gdd_x_harv_light", tsd = "tsd_harv", n = "n_harv1"),
             harv_heavy = c(main = "harv_heavy", gdd = "gdd_x_harv_heavy", tsd = "tsd_harv", n = "n_harv2"),
             harv_any = c(main = "harv_any", gdd = "gdd_x_harv_any", tsd = "tsd_harv", n = "n_harv_any"))

val <- fread(file.path(OUT_DIR, sprintf("fit_validity_%s.csv", RUN_TAG)))[used == TRUE, .(k, label, m)]
rf <- file.path(OUT_DIR, sprintf(c("mi_raw_%s.csv", "mi_raw_redraw_%s.csv", "mi_raw_sweepfull_%s.csv"), RUN_TAG)); raw <- rbindlist(lapply(rf[file.exists(rf)], fread), fill = TRUE)   # sweepfull: int run only
raw <- raw[val, on = .(k, m)]; stopifnot(nrow(val) == 91 * 30L)
ts <- readRDS(file.path(TABLE_DIR, "transition_set.rds")); ts[, `:=`(label = paste(from, to, sep = "->"), n_harv_any = n_harv1 + n_harv2)]

# ---- contrasts as named coefficient vectors, restricted to the terms present in a transition ------------------------------------------------------------------
make_contrasts <- function(terms) {
  has <- function(...) all(c(...) %in% terms); L <- list()
  for (v in c("z_gdd_3", "z_cmi_gs", "z_soil_pc1")) if (has(v)) L[[paste0("per_sd_", sub("^z_", "", v))]] <- setNames(1, v)
  for (d in names(DIST)) { s <- DIST[[d]]; if (!has(s[["main"]], s[["gdd"]], s[["tsd"]])) next
    for (z in Z_GRID) for (tau in TAU_GRID) L[[sprintf("%s_vs_none z=%g tau=%g", d, z, tau)]] <- setNames(c(1, z, tau), s[c("main", "gdd", "tsd")])
    L[[paste0(d, "_gdd_per_sd")]] <- setNames(c(1, 1), c("z_gdd_3", s[["gdd"]]))
    L[[paste0(d, "_gdd_interaction")]] <- setNames(1, s[["gdd"]])
    L[[paste0(d, "_tsd_per_decade")]] <- setNames(1, s[["tsd"]]) }
  L }
rubin <- function(q, u, nu_com) {                                              # q, u: per-imputation log HR and variance
  M <- length(q); Qb <- mean(q); W <- mean(u); B <- var(q); Tt <- W + (1 + 1 / M) * B
  r <- if (W > 0) (1 + 1 / M) * B / W else NA_real_; lam <- r / (1 + r)
  nu_old <- if (is.finite(lam) && lam > 0) (M - 1) / lam^2 else Inf; nu_obs <- (nu_com + 1) / (nu_com + 3) * nu_com * (1 - lam)
  df <- if (is.finite(nu_old)) 1 / (1 / nu_old + 1 / nu_obs) else nu_obs
  se <- sqrt(Tt); list(log_hr = Qb, se = se, df = df, fmi = (r + 2 / (df + 3)) / (r + 1), M = M) }

# ---- per fit: contrast estimates and delta-method variances ---------------------------------------------------------------------------------------------------
res <- list(); acr <- list()
for (mm in sort(unique(val$m))) {
  fits <- readRDS(file.path(OUT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, mm)))$fits
  for (lb in val[m == mm, label]) {
    cv <- fits[[lb]]$cov; stopifnot(!is.null(cv), is.null(cv$error)); r <- raw[label == lb & m == mm]
    est <- setNames(r$est, r$parameter); al <- est[["alpha"]]; V <- cv$cov; nm <- colnames(V)
    asd <- r[parameter == "alpha", sd]; acr[[length(acr) + 1L]] <- data.table(label = lb, m = mm, ratio = sqrt(V["alpha", "alpha"]) / asd)
    for (cn in names(make_contrasts(setdiff(nm, c("(Intercept)", "alpha"))))) {
      cc <- make_contrasts(setdiff(nm, c("(Intercept)", "alpha")))[[cn]]; cb <- sum(cc * est[names(cc)])
      g <- setNames(numeric(length(nm)), nm); g[names(cc)] <- al * cc; g["alpha"] <- cb
      res[[length(res) + 1L]] <- data.table(label = lb, m = mm, contrast = cn, q = al * cb, u = as.numeric(t(g) %*% V %*% g)) } } }
res <- rbindlist(res); acr <- rbindlist(acr)[, .(alpha_cov_ratio_med = median(ratio), alpha_cov_low_share = mean(ratio < 0.9)), by = label]
nu <- raw[, .(nu_com = max(mean(n_events[parameter == "alpha"]) - uniqueN(parameter), 1)), by = label]
out <- merge(res[, c(rubin(q, u, nu$nu_com[match(label[1], nu$label)])), by = .(label, contrast)], nu[, .(label)], by = "label")
out[, `:=`(hr = exp(log_hr), hr_lo = exp(log_hr - qt(0.975, df) * se), hr_hi = exp(log_hr + qt(0.975, df) * se), p = 2 * pt(-abs(log_hr / se), df))]
# disturbed-event counts and flags
out[, dist := fifelse(grepl("^(fire|pest|harv_light|harv_heavy|harv_any)", contrast), sub("(_vs_none.*|_gdd_per_sd|_gdd_interaction|_tsd_per_decade)$", "", contrast), NA_character_)]
out <- merge(out, ts[, .(label, n_events, n_fire, n_pest, n_harv1, n_harv2, n_harv_any)], by = "label", all.x = TRUE)
out[, n_dist := mapply(function(d, i) if (is.na(d)) NA_integer_ else as.integer(c(fire = n_fire[i], pest = n_pest[i], harv_light = n_harv1[i], harv_heavy = n_harv2[i], harv_any = n_harv_any[i])[[d]]), dist, seq_len(.N))]
out[, low_events := !is.na(n_dist) & n_dist <= LOW_EVENTS_MAX]
out <- merge(out, acr, by = "label"); out[, c("n_fire", "n_pest", "n_harv1", "n_harv2", "n_harv_any") := NULL]
setcolorder(out, c("label", "contrast", "dist", "hr", "hr_lo", "hr_hi", "p", "log_hr", "se", "fmi", "df", "M", "n_events", "n_dist", "low_events", "alpha_cov_ratio_med", "alpha_cov_low_share"))
fwrite(out[order(label, contrast)], file.path(OUT_DIR, sprintf("hazard_ratios_%s.csv", RUN_TAG)))
cat(sprintf("Hazard ratios: %d rows, %d transitions, %d contrast types; M = %s\n", nrow(out), uniqueN(out$label), uniqueN(sub(" z=.*", "", out$contrast)), paste(unique(out$M), collapse = "/")))
# ---- sanity checks -----------------------------------------------------------------------------------------------------------------------------------------------
pf <- fread(file.path(OUT_DIR, sprintf("pooled_%s_final.csv", RUN_TAG)))
ck <- merge(out[contrast == "per_sd_cmi_gs", .(label, log_hr, se)], dcast(pf[parameter %in% c("alpha", "z_cmi_gs"), .(label, parameter, pooled_est)], label ~ parameter, value.var = "pooled_est"), by = "label")
cat(sprintf("\nCheck: pooled log HR per SD CMI vs pooled alpha x pooled beta (equal up to the alpha-beta covariance across imputations): max |difference| = %.4f, median %.4f\n", max(abs(ck$log_hr - ck$alpha * ck$z_cmi_gs)), median(abs(ck$log_hr - ck$alpha * ck$z_cmi_gs))))
cat("\nblack_spruce->regen (reference, the transition with the most events):\n")
print(out[label == "black_spruce->regen" & (!grepl("z=", contrast) | grepl("z=0 tau=(0|1)$", contrast)), .(contrast, hr = round(hr, 2), lo = round(hr_lo, 2), hi = round(hr_hi, 2), p = signif(p, 2), fmi = round(fmi, 2))][order(contrast)])
cat("\nAlpha-sampling flag: transitions where >= 20% of the used fits have sampled alpha SD < 0.9 x INLA's:", out[alpha_cov_low_share >= 0.2, uniqueN(label)], "\n"); print(unique(out[alpha_cov_low_share >= 0.2, .(label, alpha_cov_low_share = round(alpha_cov_low_share, 2))])[order(-alpha_cov_low_share)])
