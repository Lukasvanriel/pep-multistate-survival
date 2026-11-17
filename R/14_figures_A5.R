# ==============================================================================
# 14_figures_A5.R   (no INLA; run from the project root: Rscript R/14_figures_A5.R)
# Exploratory figure set A, part 5: COMBINATIONS of coefficients (model output evaluated at chosen covariate values). NOT manuscript figures.
# Needs: 08 (pooled_<tag>_final.csv + pooled_<tag>_final_cov.rds), 09 (hazard_ratios_<tag>.csv), 04 (visit_intervals.rds).
# Writes: data-output/<tag>/figures/A/A5*  (and A5_by_transition.pdf, one page per transition with any disturbance term)
#
# Disturbed vs undisturbed log HR at GDD z and time since disturbance tau (years):  beta_d + beta_gddxd * z + beta_tsd_d * log(tau + 1)
#   (only the terms present in that transition). Same contrast as 09's "<d>_vs_none z= tau=" rows; checked against them below.
# Point estimates: pooled INLA means (pooled_<tag>_final.csv). Uncertainty: Rubin W + (1 + 1/M) B, W from the pooled covariance file (08), B from the INLA means of the used fits (as 09), df as in 09 (Barnard-Rubin, nu_com = events - number of parameters),
#   95% t interval. 09's one-sided-bound cases (reference or level with 0 events) are NOT drawn (not identified by the data).
#
#   A5a  disturbance vs none over tau (0 to the max TSD at risk), GDD z = 0
#   A5b  disturbance vs none over GDD z (-2.5 to 2.5), at tau = 10 y and 30 y
#   A5c  GDD (+1 SD) HR: undisturbed vs within each disturbance (09 rows per_sd_gdd_3, <d>_gdd_per_sd)
#   A5d  heavy vs light harvest over tau, GDD z = 0 (transitions with both harvest levels identified)
#   A5e/A5f  disturbance vs none on 09's tau grid (1, 10, 30, 60 y), z = 0: heatmaps / forest plots (moved from 13, formerly A4c/A4d)
# Display (user, 2026-11-08: "make both so we can choose later"): curves solid inside the observed range (5-95th pct), dashed outside;
#   observed TSD range = TSD at the EVENT (09 tsd_ev_p05..p95); observed GDD range = z of the at-risk disturbed intervals of the origin state
#   (visit_intervals, state_start = origin, flag on; 5-95th pct). Both an overview per disturbance type and one page per transition.
# ==============================================================================
source(here::here("R", "13_fig_helpers.R"))
source(here::here("R", "covariate_spec.R")); f_tsd <- tsd_fun(TSD_FORM)
cvl <- readRDS(file.path(RESULTS_DIR, sprintf("pooled_%s_final_cov.rds", RUN_TAG)))
stopifnot(setequal(names(cvl), modelled))
DIST  <- c("fire", "pest", "harv_light", "harv_heavy"); DL <- c(fire = "fire", pest = "pest", harv_light = "light harvest", harv_heavy = "heavy harvest")
DCOL  <- c(fire = "#D95F02", pest = "#7570B3", harv_light = "#66A61E", harv_heavy = "#1B6E3C")
TERMS <- list(fire = c("fire_f", "gdd_x_fire", "tsd_fire"), pest = c("pest_f", "gdd_x_pest", "tsd_pest"),
              harv_light = c("harv_light", "gdd_x_harv_light", "tsd_harv_light"), harv_heavy = c("harv_heavy", "gdd_x_harv_heavy", "tsd_harv_heavy"))
nu_com <- sapply(modelled, function(lb) max(ar[label == lb, events] - length(cvl[[lb]]$mean), 1))

# linear combinations C (rows = grid points, named columns = terms) for one transition -> est, se, df, lo, hi (as 09's rubin())
# point estimates: pooled INLA means (pooled_<tag>_final.csv, as 09); the cov file's $mean is a posterior-sample mean (MC noise, differs by up to ~0.02), not used
PE <- split(pooled[, .(label, parameter, pooled_est)], by = "label", keep.by = FALSE); PE <- lapply(PE, function(x) setNames(x$pooled_est, x$parameter))
# B (between imputations) from the INLA means of the 30 used fits (mi_raw + mi_raw_redraw, fit_validity used = TRUE), as 09 does (var of the per-fit
#   estimates); the cov file's $B is built from posterior-sample means (08 pool_cov) and is not used. W = the cov file's $W (mean per-fit covariance), as 09.
fv  <- fread(file.path(RESULTS_DIR, sprintf("fit_validity_%s.csv", RUN_TAG)))[used == TRUE, .(label, m)]
mir <- rbindlist(lapply(file.path(FIT_DIR, sprintf(c("mi_raw_%s.csv", "mi_raw_redraw_%s.csv"), RUN_TAG)), function(f) if (file.exists(f)) fread(f)[, .(label, m, parameter, est)]))
mir <- unique(mir[fv, on = c("label", "m"), nomatch = NULL], by = c("label", "m", "parameter"))
stopifnot(all(mir[, uniqueN(m), by = label]$V1 == 30L))
BI <- lapply(split(mir, by = "label"), function(x) { Q <- as.matrix(dcast(x, m ~ parameter, value.var = "est")[, -1]); stats::cov(Q) })
lincomb <- function(lb, C) {
  m <- cvl[[lb]]; C <- C[, colnames(C) %in% names(m$mean), drop = FALSE]; nm <- colnames(C)
  est <- drop(C %*% PE[[lb]][nm]); W <- rowSums((C %*% m$W[nm, nm, drop = FALSE]) * C); B <- rowSums((C %*% BI[[lb]][nm, nm, drop = FALSE]) * C)
  M <- m$M; Tt <- W + (1 + 1 / M) * B; r <- ifelse(W > 0, (1 + 1 / M) * B / W, NA_real_); lam <- r / (1 + r)
  nu_old <- ifelse(is.finite(lam) & lam > 0, (M - 1) / lam^2, Inf); nc <- nu_com[[lb]]; nu_obs <- (nc + 1) / (nc + 3) * nc * (1 - lam)
  df <- ifelse(is.finite(nu_old), 1 / (1 / nu_old + 1 / nu_obs), nu_obs); se <- sqrt(Tt); q <- qt(0.975, df)
  data.table(est = est, se = se, df = df, lo = est - q * se, hi = est + q * se)
}
dist_C <- function(d, z, tau) { tm <- TERMS[[d]]; C <- cbind(1, z, f_tsd(tau)); colnames(C) <- tm; C }
has_level <- function(lb, d) TERMS[[d]][1] %in% names(cvl[[lb]]$mean)
bound_pairs <- unique(hr[!is.na(bound) & grepl("_vs_none", contrast), .(label, dist)])

# ---- check against 09 (tau grid x z grid) -------------------------------------------------------------------------------------------------------------------
chk <- hr[grepl("^(fire|pest|harv_light|harv_heavy)_vs_none z=", contrast) & is.na(bound)]
chk[, `:=`(z = as.numeric(sub(".* z=(-?[0-9]+) .*", "\\1", contrast)), tau = as.numeric(sub(".*tau=([0-9]+)y", "\\1", contrast)))]
mine <- chk[, { r <- lincomb(label[1], dist_C(dist[1], z, tau)); .(z, tau, est = r$est, lo = r$lo, hi = r$hi, log_hr, lhl = log(hr_lo), lhh = log(hr_hi)) }, by = .(label, dist)]
dmax <- mine[, .(est = max(abs(est - log_hr)), lo = max(abs(lo - lhl)), hi = max(abs(hi - lhh)))]
cat(sprintf("check vs 09 (%d contrasts): max |diff| log HR %.2e, lower %.2e, upper %.2e\n", nrow(mine), dmax$est, dmax$lo, dmax$hi))
stopifnot(dmax$est < 1e-6, dmax$lo < 1e-3, dmax$hi < 1e-3)

# ---- observed ranges ------------------------------------------------------------------------------------------------------------------------------------------
sc <- readRDS(file.path(FIT_DIR, "scaling.rds")); gm <- sc[var == "gdd_3", mean]; gs <- sc[var == "gdd_3", sd]
flag_on <- list(fire = quote(fire == 1), pest = quote(pest == 1), harv_light = quote(harv == 1), harv_heavy = quote(harv == 2))
zr <- rbindlist(lapply(DIST, function(d) vi[eval(flag_on[[d]]), .(dist = d, z_p05 = quantile((gdd_3 - gm) / gs, .05, na.rm = TRUE),
                                                                   z_p95 = quantile((gdd_3 - gm) / gs, .95, na.rm = TRUE), n_risk_int = .N), by = .(from = state_start)]))
tr <- unique(hr[!is.na(dist) & !is.na(tsd_risk_max), .(label, dist, tsd_risk_max, tsd_ev_p05, tsd_ev_p95)])

# ---- A5a curves over tau, A5b curves over GDD ----------------------------------------------------------------------------------------------------------------
pairs <- rbindlist(lapply(modelled, function(lb) data.table(label = lb, dist = DIST[sapply(DIST, has_level, lb = lb)])))
pairs <- pairs[!bound_pairs, on = c("label", "dist")]
pairs <- merge(merge(pairs, tr, by = c("label", "dist"), all.x = TRUE), unique(trow[, .(label, from, to)]), by = "label")
pairs <- merge(pairs, zr, by = c("from", "dist"), all.x = TRUE)
cat(sprintf("A5: %d transition x disturbance pairs drawn (%d one-sided-bound pairs left out)\n", nrow(pairs), nrow(bound_pairs)))
cA <- pairs[, { tau <- seq(0, max(tsd_risk_max, 60, na.rm = TRUE), by = 0.5); r <- lincomb(label, dist_C(dist, 0, tau))
               .(tau = tau, est = r$est, lo = r$lo, hi = r$hi, inside = !is.na(tsd_ev_p05) & tau >= tsd_ev_p05 & tau <= tsd_ev_p95) }, by = .(label, from, to, dist)]
cB <- pairs[, rbindlist(lapply(c(10, 30), function(tt) { z <- seq(-2.5, 2.5, by = 0.05); r <- lincomb(label, dist_C(dist, z, tt))
               data.table(tau = tt, z = z, est = r$est, lo = r$lo, hi = r$hi, inside = !is.na(z_p05) & z >= z_p05 & z <= z_p95) })), by = .(label, from, to, dist)]
for (x in list(cA, cB)) x[, `:=`(hr = exp(est), hr_lo = exp(lo), hr_hi = exp(hi), hr_in = fifelse(inside, exp(est), NA_real_))]
cA[, hr_pt := fifelse(inside & sum(inside) == 1L, hr, NA_real_), by = .(label, dist)]   # observed range = one point (e.g. a single event interval): drawn as a dot
cB[, hr_pt := fifelse(inside & sum(inside) == 1L, hr, NA_real_), by = .(label, dist, tau)]
fwrite(cA[, !c("hr_in", "hr_pt")], file.path(OUT, "A5a_curves_tau.csv")); fwrite(cB[, !c("hr_in", "hr_pt")], file.path(OUT, "A5b_curves_gdd.csv"))

# ---- A5d heavy vs light harvest over tau ------------------------------------------------------------------------------------------------------------------------
hl <- pairs[, .(n = .N), by = label][label %in% pairs[dist == "harv_light", label] & label %in% pairs[dist == "harv_heavy", label], label]
cD <- rbindlist(lapply(hl, function(lb) { tau <- seq(0, 60, by = 0.5)
  C <- cbind(1, -1, f_tsd(tau), -f_tsd(tau)); colnames(C) <- c("harv_heavy", "harv_light", "tsd_harv_heavy", "tsd_harv_light"); r <- lincomb(lb, C)
  rng <- tr[label == lb & dist %in% c("harv_light", "harv_heavy"), .(lo = max(tsd_ev_p05), hi = min(tsd_ev_p95))]   # inside = within BOTH levels' event range
  data.table(label = lb, tau = tau, est = r$est, lo = r$lo, hi = r$hi, inside = tau >= rng$lo & tau <= rng$hi) }))
cD[, `:=`(hr = exp(est), hr_lo = exp(lo), hr_hi = exp(hi), hr_in = fifelse(inside, exp(est), NA_real_))]
cD[, hr_pt := fifelse(inside & sum(inside) == 1L, hr, NA_real_), by = label]
cD <- merge(cD, unique(trow[, .(label, from, to)]), by = "label"); fwrite(cD[, !c("hr_in", "hr_pt")], file.path(OUT, "A5d_heavy_vs_light.csv"))

# ---- overview figures: one per disturbance type, all transitions as small panels ------------------------------------------------------------------------------
ord_lab <- function(x) factor(x, levels = trow[order(match(from, ORDER), match(to, ORDER)), label])
plab <- function(x) gsub("_", " ", sub("->", " → ", x))
curve_panel <- function(d, xvar, xlab_txt, colour_by = NULL) {
  d <- copy(d); d[, label := ord_lab(label)]
  g <- ggplot(d, aes(.data[[xvar]], hr)) + geom_hline(yintercept = 1, linetype = 2, colour = "grey50", linewidth = 0.3)
  if (is.null(colour_by)) g <- g + geom_ribbon(aes(ymin = hr_lo, ymax = hr_hi), fill = "grey80") +
      geom_line(linetype = 2, linewidth = 0.3) + geom_line(aes(y = hr_in), linewidth = 0.6, na.rm = TRUE) + geom_point(aes(y = hr_pt), size = 1, na.rm = TRUE)
  else g <- g + geom_ribbon(aes(ymin = hr_lo, ymax = hr_hi, fill = .data[[colour_by]]), alpha = 0.15) +
      geom_line(aes(colour = .data[[colour_by]]), linetype = 2, linewidth = 0.3) + geom_line(aes(y = hr_in, colour = .data[[colour_by]]), linewidth = 0.6, na.rm = TRUE) + geom_point(aes(y = hr_pt, colour = .data[[colour_by]]), size = 1, na.rm = TRUE)
  g + scale_y_log10(labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)) +
    facet_wrap(~ label, scales = "free_y", labeller = as_labeller(plab), ncol = 8) + labs(x = xlab_txt, y = "hazard ratio vs undisturbed (log scale)") +
    theme(strip.text = element_text(size = 6), axis.text = element_text(size = 6))
}
for (d in DIST) {
  a <- cA[dist == d]; if (!nrow(a)) next; n <- uniqueN(a$label); hgt <- 1.6 * ceiling(n / 8) + 1.5
  save_fig(curve_panel(a, "tau", "time since disturbance (years)") +
             labs(title = sprintf("A5a  %s vs undisturbed over time since disturbance", DL[d]),
                  subtitle = "band = 95% CI; dashed = outside the observed range"), paste0("A5a_tau_", d), 16, hgt, NULL)
  b <- cB[dist == d][, tau_lab := factor(sprintf("%d y", tau), levels = c("10 y", "30 y"))]
  save_fig(curve_panel(b, "z", "GDD (SD from the mean)", "tau_lab") + scale_colour_manual(values = c(`10 y` = "#1B9E77", `30 y` = "#7570B3"), name = "time since\ndisturbance") +
             scale_fill_manual(values = c(`10 y` = "#1B9E77", `30 y` = "#7570B3"), name = "time since\ndisturbance") +
             labs(title = sprintf("A5b  %s vs undisturbed along GDD", DL[d]),
                  subtitle = "band = 95% CI; dashed = outside the observed range"),
           paste0("A5b_gdd_", d), 16, hgt, NULL)
}
save_fig(curve_panel(cD, "tau", "time since harvest (years)") + labs(y = "hazard ratio, heavy / light harvest (log scale)", title = "A5d  Heavy vs light harvest over time since harvest",
           subtitle = "band = 95% CI; dashed = outside the observed range"),
         "A5d_heavy_vs_light", 16, 1.6 * ceiling(uniqueN(cD$label) / 8) + 1.5, NULL)

# ---- A5c GDD effect undisturbed vs within each disturbance (09 rows) ------------------------------------------------------------------------------------------
gc <- hr[contrast == "per_sd_gdd_3" | grepl("_gdd_per_sd$", contrast), .(label, contrast, hr, hr_lo, hr_hi)]
gc[, cond := fifelse(contrast == "per_sd_gdd_3", "undisturbed", DL[sub("_gdd_per_sd$", "", contrast)])]
gc[, `:=`(sig = hr_lo > 1 | hr_hi < 1, cond = factor(cond, levels = rev(c("undisturbed", DL))))]
fwrite(gc, file.path(OUT, "A5c_gdd_by_disturbance.csv"))

# ---- one page per transition (A5_by_transition.pdf) ------------------------------------------------------------------------------------------------------------
fcol <- c(`TRUE` = "#B2182B", `FALSE` = "#2166AC")
page <- function(lb) {
  ttl <- sprintf("%s  (%d events)", plab(lb), trow[label == lb, events]); P <- list()
  a <- cA[label == lb]; b <- cB[label == lb]; g5 <- gc[label == lb]; dd <- cD[label == lb]
  ylog <- scale_y_log10(labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE))
  if (nrow(a)) P$a <- ggplot(a, aes(tau, hr, colour = dist, fill = dist)) + geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
    geom_ribbon(aes(ymin = hr_lo, ymax = hr_hi), alpha = 0.12, colour = NA) + geom_line(linetype = 2, linewidth = 0.3) +
    geom_line(aes(y = hr_in), linewidth = 0.8, na.rm = TRUE) + geom_point(aes(y = hr_pt), size = 1.5, na.rm = TRUE) + ylog + scale_colour_manual(values = DCOL, labels = DL, name = NULL) +
    scale_fill_manual(values = DCOL, labels = DL, name = NULL) + labs(x = "time since disturbance (years)", y = "hazard ratio vs undisturbed", subtitle = "A5a  over time since disturbance")
  if (nrow(b)) P$b <- ggplot(b, aes(z, hr, colour = dist, fill = dist)) + geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
    geom_ribbon(aes(ymin = hr_lo, ymax = hr_hi), alpha = 0.12, colour = NA) + geom_line(linetype = 2, linewidth = 0.3) +
    geom_line(aes(y = hr_in), linewidth = 0.8, na.rm = TRUE) + geom_point(aes(y = hr_pt), size = 1.5, na.rm = TRUE) + ylog + facet_wrap(~ sprintf("%d y since disturbance", tau)) +
    scale_colour_manual(values = DCOL, labels = DL, name = NULL) + scale_fill_manual(values = DCOL, labels = DL, name = NULL) +
    labs(x = "GDD (SD from the mean)", y = "hazard ratio vs undisturbed", subtitle = "A5b  along GDD")
  if (nrow(g5)) P$c <- ggplot(g5, aes(hr, cond, colour = hr > 1)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_errorbar(aes(xmin = hr_lo, xmax = hr_hi), width = 0, orientation = "y") + geom_point(aes(shape = sig), size = 2, fill = "white", stroke = 0.8) +
    geom_text(aes(x = max(hr_hi) * 1.6, label = sprintf("%s [%s, %s]", fmt_hr(hr), fmt_hr(hr_lo), fmt_hr(hr_hi)), fontface = ifelse(sig, "bold", "plain")), hjust = 0, size = 2.4) +
    scale_x_log10(expand = expansion(mult = c(0.05, 0.9))) + scale_colour_manual(values = fcol, guide = "none") + scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), guide = "none") +
    labs(x = "GDD hazard ratio (log scale)", y = NULL, subtitle = "A5c  GDD effect: undisturbed vs within each disturbance")
  if (nrow(dd)) P$d <- ggplot(dd, aes(tau, hr)) + geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
    geom_ribbon(aes(ymin = hr_lo, ymax = hr_hi), fill = "grey80") + geom_line(linetype = 2, linewidth = 0.3) + geom_line(aes(y = hr_in), linewidth = 0.8, na.rm = TRUE) + geom_point(aes(y = hr_pt), size = 1.5, na.rm = TRUE) +
    ylog + labs(x = "time since harvest (years)", y = "hazard ratio, heavy / light", subtitle = "A5d  heavy vs light harvest")
  if (!length(P)) return(NULL)
  wrap_plots(P, ncol = 2, guides = "collect") + plot_annotation(title = ttl,
    caption = "band = 95% CI; dashed = outside the observed range")
}
n_pg <- 0L
cairo_pdf(file.path(OUT, "A5_by_transition.pdf"), width = 11, height = 7.5, onefile = TRUE)
for (lb in levels(ord_lab(modelled))) { pg <- page(lb); if (!is.null(pg)) { print(pg); n_pg <- n_pg + 1L } }
invisible(dev.off()); cat(sprintf("  wrote A5_by_transition.pdf (%d pages)\n", n_pg))

# ---- A5e / A5f (moved from 13): disturbance vs none on the tau grid of 09 ----------------------------------------------------------------------------------------
TAUS <- c(1, 10, 30, 60)
# A5e heatmaps / A5f forest plots: disturbance vs none at GDD z = 0, at the four tau of 09
dv <- rbindlist(lapply(DIST, function(d) rbindlist(lapply(TAUS, function(tau) {
  x <- hr[contrast == sprintf("%s_vs_none z=0 tau=%dy", d, tau)]
  m <- hr[contrast == sprintf("%s_vs_none", d) & !(label %in% x$label)]         # level with only its main term: no tau (09, Q1)
  rbind(x[, main_only := FALSE], m[, main_only := TRUE], fill = TRUE)[, `:=`(dist = d, tau = tau)]
}))))
rng <- unique(hr[!is.na(tsd_ev_p05), .(label, dist, tsd_ev_p05, tsd_ev_p95)])
dv <- merge(dv[, !c("tsd_ev_p05", "tsd_ev_p95")], rng, by = c("label", "dist"), all.x = TRUE)
dv[, `:=`(is_bound = !is.na(bound), outside = !main_only & !is.na(tsd_ev_p05) & (tau < tsd_ev_p05 | tau > tsd_ev_p95))]
dv[, sig := !is_bound & (hr_lo > 1 | hr_hi < 1)]
limd <- cap_sym(log(dv[is_bound == FALSE, hr]))
for (d in DIST) {
  gg <- rbindlist(lapply(TAUS, function(tau) {
    tau_i <- tau; g <- merge(full_grid(), dv[dist == d & tau == tau_i], by = c("from", "to"), all.x = TRUE)
    g[, lab2 := paste0(from, "->", to)]
    g[, status := fifelse(!(lab2 %in% modelled), "notmod", fifelse(is.na(hr), "noterm", fifelse(is_bound, "bound", "value")))]
    g[, `:=`(fillv = log(hr), bold = sig %in% TRUE, tau_lab = sprintf("%d y since disturbance", tau_i))]
    g[status == "value", txt := paste0(fmt_hr(hr), fifelse(outside, "+", ""), fifelse(main_only, "*", ""))]
    g[status == "bound", txt := paste0(fifelse(bound == "lower", "≥", "≤"), fmt_hr(hr_bound))]
    g
  }))
  gg[, tau_lab := factor(tau_lab, levels = sprintf("%d y since disturbance", TAUS))]
  p <- heat(gg, "hazard\nratio", limits = c(-limd, limd), trans_lab = function(x) fmt_hr(exp(x)), breaks = hr_breaks(limd),
            title = sprintf("A5e  %s vs undisturbed", lab(d)),
            subtitle = "bold = 95% CI excludes 1; + = outside the observed range; * = no time-since-disturbance term; ≥ / ≤ = one-sided bound") +
    facet_wrap(~ tau_lab, ncol = 2)
  save_fig(p, paste0("A5e_dist_vs_none_", d), 15, 14, gg[status %in% c("value", "bound"), .(label = lab2, tau, hr, hr_lo, hr_hi, bound, hr_bound, outside, main_only, sig)])
  fd <- merge(dv[dist == d & is_bound == FALSE], trow[, .(label, ylab)], by = "label")
  fd[, `:=`(from = factor(from, levels = ORDER), tau = factor(tau, levels = TAUS))]
  p <- ggplot(fd, aes(hr, ylab, colour = tau)) + geom_vline(xintercept = 1, linetype = 2) +
    geom_errorbar(aes(xmin = hr_lo, xmax = hr_hi, group = tau), width = 0, position = position_dodge(width = 0.7), linewidth = 0.3, orientation = "y") +
    geom_point(aes(shape = outside, group = tau), position = position_dodge(width = 0.7), size = 1.2) +   # group = tau in both layers: shape must not change the dodge (user 2026-11-10)
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), labels = c(`FALSE` = "within", `TRUE` = "outside"), name = "observed range") +
    scale_colour_viridis_d(end = 0.85, name = "years since\ndisturbance") + scale_x_log10(labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)) +
    facet_grid(from ~ ., scales = "free_y", space = "free_y", labeller = as_labeller(lab), switch = "y") +
    labs(x = "hazard ratio vs undisturbed (log scale)", y = "destination (events)", title = sprintf("A5f  %s vs undisturbed by origin", lab(d))) +
    theme(strip.text.y.left = element_text(angle = 0), strip.placement = "outside", panel.spacing = unit(1, "pt"), axis.text.y = element_text(size = 6))
  save_fig(p, paste0("A5f_forest_dist_", d), 8, max(4, 0.12 * nrow(fd) / 4 * 4 + 2), NULL)
}

cat(sprintf("\nDone (part 5): %d A5 files in %s\n", length(list.files(OUT, pattern = "^A5")), OUT))
