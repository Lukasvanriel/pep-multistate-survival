# ==============================================================================
# 13_figures_A.R   (no INLA; ~1 min; run from the project root: Rscript R/13_figures_A.R)
# Exploratory figure set A, parts 1-4, for the final run (part 5 = 14_figures_A5.R). NOT manuscript figures.
#   A1 data and model structure   A2 baseline rate   A3 Weibull shape alpha   A4 hazard ratios
# Needs: 04 (data/tables/<tag>/), 08, 09 (data-output/<tag>/results/).
# Writes: data-output/<tag>/figures/A/  -> <name>.png (300 dpi), <name>.pdf, <name>.csv (the plotted values)
#
# Model facts used (05a, 08, 09): Weibull variant 0. Years scale: S(t) = exp(-exp(b0_y) t^alpha), b0_y = b0 - alpha*log(10)
#   (pooled in intercept_years_<tag>.csv). Hazard h(t) = alpha t^(alpha-1) exp(b0_y + beta'x). HR = exp(c'beta) (alpha does not enter).
#
# Display decisions (user, 2026-11-08):
#   - baseline (A2) = hazard at t = 10 y in the state; A2b = single-route median time t50 (S = 0.5; user 2026-11-10, replaced t_c); undisturbed, GDD/CMI/soil z = 0; single-route quantity (competing exits ignored)
#   - bold text = 95% CI excludes 1 (HRs, alpha); no extra flag for low-event transitions
#   - HRs: heatmaps + one forest plot per origin; h(t) for every transition, faceted by origin
#   - two neutral fills: "transition not modelled" (light grey) vs "term not in this transition's model" (mid grey, "-")
# Display choices made here, to be confirmed (marked CONFIRM below):
#   - colour scales on log(HR), symmetric, capped at the 98th percentile of |log HR| of the figure (values beyond are squished)
#   - A3c and A2 use log axes / log colour scales
#   - A4c: disturbance vs none at GDD z = 0 for the four tau of 09 (1, 10, 30, 60 y); "+" marks a tau outside the 5-95th percentile
#     of TSD at the event (tsd_ev_p05..p95, from 09); 09's one-sided bounds are printed as ">= x" / "<= x" on a white cell
# ==============================================================================
source(here::here("R", "13_fig_helpers.R"))   # helpers + inputs (ORDER, PAL, heat(), save_fig(), pooled, hr, ar, iy, ts, vi, trow)

# ==============================================================================
# A1 DATA AND MODEL STRUCTURE
# ==============================================================================
cat("A1 data and structure\n")
# A1a observed state changes, origin x destination (all observed pairs; black outline = modelled, >= MIN_EVENTS)
g <- merge(full_grid(), ev[, .(from, to, n_events, kept)], by = c("from", "to"), all.x = TRUE)
g[is.na(n_events), `:=`(n_events = 0L, kept = FALSE)]
g[, `:=`(status = fifelse(n_events > 0, "value", "notmod"), fillv = log10(pmax(n_events, 1)), txt = as.character(n_events), bold = FALSE)]
p <- heat(g, "events", diverging = FALSE, trans_lab = function(x) round(10^x),
          title = "A1a  Observed state changes", subtitle = sprintf("black outline = modelled transition (≥ %d events); grey = none observed", MIN_EVENTS)) +
  geom_tile(data = g[kept == TRUE][, `:=`(from = factor(from, levels = rev(ORDER)), to = factor(to, levels = ORDER))],
            fill = NA, colour = "black", linewidth = 0.4)
save_fig(p, "A1a_observed_events", 7.5, 7, g[, .(from, to, n_events, modelled = kept)])

# A1b terms in each modelled transition
terms <- setdiff(unique(pooled$parameter), c("alpha", "(Intercept)"))
TERM_ORDER <- c("z_gdd_3", "z_cmi_gs", "z_soil_pc1", "fire_f", "pest_f", "harv_light", "harv_heavy",
                "gdd_x_fire", "gdd_x_pest", "gdd_x_harv_light", "gdd_x_harv_heavy", "tsd_fire", "tsd_pest", "tsd_harv_light", "tsd_harv_heavy")
stopifnot(setequal(terms, TERM_ORDER))
tm <- CJ(label = modelled, parameter = TERM_ORDER)
tm[, present := paste(label, parameter) %in% pooled[, paste(label, parameter)]]
tm <- merge(tm, trow[, .(label, from, ylab)], by = "label")
tm[, `:=`(from = factor(from, levels = ORDER), parameter = factor(parameter, levels = TERM_ORDER))]
p <- ggplot(tm, aes(parameter, ylab, fill = present)) + geom_tile(colour = "white") +
  scale_fill_manual(values = c(`TRUE` = "grey25", `FALSE` = "grey92"), labels = c(`TRUE` = "in model", `FALSE` = "not in model")) +
  facet_grid(from ~ ., scales = "free_y", space = "free_y", labeller = as_labeller(lab), switch = "y") +
  labs(x = NULL, y = "destination (events), grouped by origin", fill = NULL, title = "A1b  Terms in each transition's model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text.y.left = element_text(angle = 0), strip.placement = "outside",
        panel.grid = element_blank(), panel.spacing = unit(1, "pt"))
save_fig(p, "A1b_terms_by_transition", 7, 13, tm[, .(label, parameter, present)])

# A1c covariates by state (one value per visit interval, raw units; plots with more intervals weigh more)
cv <- melt(vi[, .(state = state_start, `GDD` = gdd_3, `CMI (growing season)` = cmi_gs, `soil PC1` = soil_pc1)],
           id.vars = "state", variable.name = "covariate")
cv[, state := factor(state, levels = ORDER)]
p <- ggplot(cv, aes(state, value, fill = state)) + geom_boxplot(outlier.size = 0.2, linewidth = 0.3) +
  scale_fill_manual(values = PAL, guide = "none") + scale_x_discrete(labels = lab) +
  facet_wrap(~ covariate, ncol = 1, scales = "free_y") + labs(x = "state at interval start", y = NULL,
  title = "A1c  Covariates by state", subtitle = sprintf("one value per visit interval (n = %s intervals, %s plots)", format(nrow(vi), big.mark = ","), format(uniqueN(vi$id_pe), big.mark = ","))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p, "A1c_covariates_by_state", 7, 8, cv[, .(n = .N, p05 = quantile(value, .05, na.rm = TRUE), median = median(value, na.rm = TRUE),
                                                    p95 = quantile(value, .95, na.rm = TRUE), n_na = sum(is.na(value))), by = .(covariate, state)])

# A1d TSD ranges per transition x disturbance type (from 09: at risk and at the event)
tr <- unique(hr[!is.na(dist) & !is.na(tsd_risk_p05), .(label, from, to, dist, tsd_risk_min, tsd_risk_p05, tsd_risk_p95, tsd_risk_max,
                                                       tsd_ev_min, tsd_ev_p05, tsd_ev_p95, tsd_ev_max, n_ev_int)])
stopifnot(!anyDuplicated(tr[, .(label, dist)]))
tr <- merge(tr, trow[, .(label, ylab)], by = "label")
DIST_L <- c(fire = "fire", pest = "pest", harv_light = "light harvest", harv_heavy = "heavy harvest")
tr[, `:=`(from = factor(from, levels = ORDER), dist = factor(dist, levels = c("fire", "pest", "harv_light", "harv_heavy")))]
p <- ggplot(tr, aes(y = ylab)) +
  geom_linerange(aes(xmin = tsd_risk_p05, xmax = tsd_risk_p95), colour = "grey70", linewidth = 2.2) +
  geom_linerange(aes(xmin = tsd_ev_min, xmax = tsd_ev_max), colour = "black", linewidth = 0.3) +
  geom_linerange(aes(xmin = tsd_ev_p05, xmax = tsd_ev_p95), colour = "black", linewidth = 1) +
  facet_grid(from ~ dist, scales = "free_y", space = "free_y", labeller = labeller(from = as_labeller(lab), dist = as_labeller(DIST_L)), switch = "y") +
  geom_vline(xintercept = c(10, 30, 60), linetype = 3, linewidth = 0.3) +
  labs(x = "time since disturbance (years)", y = NULL, title = "A1d  Time since disturbance covered by the data",
       subtitle = "grey = stands at risk (5–95th percentile); black = at the event (thick 5–95th percentile, thin range)") +
  theme(strip.text.y.left = element_text(angle = 0), strip.placement = "outside", panel.spacing = unit(1, "pt"), axis.text.y = element_text(size = 5.5))
save_fig(p, "A1d_tsd_ranges", 9, 13, tr)

# A1e composition of measured plots by calendar period (network grew over time: NOT a successional trajectory)
vis <- unique(rbind(vi[, .(id_pe, year = year_start, state = state_start)],
                    vi[vi[, .I[which.max(interval)], by = id_pe]$V1, .(id_pe, year = year_end, state = state_end)]))
vis[, period := paste0(floor(year / 5) * 5, "s")]
cmp <- vis[, .N, by = .(period, state)][, share := N / sum(N), by = period][, n_period := sum(N), by = period]
cmp[, state := factor(state, levels = ORDER)]
p <- ggplot(cmp, aes(period, share, fill = state)) + geom_col(width = 0.85) + scale_fill_manual(values = PAL, labels = lab) +
  geom_text(data = unique(cmp[, .(period, n_period)]), aes(period, 1.02, label = n_period), inherit.aes = FALSE, size = 2.2) +
  labs(x = "5-year period of the visit", y = "share of visits", fill = NULL, title = "A1e  State of measured plots by period",
       subtitle = "numbers = visits; plots were added over time, so this is not a trajectory") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p, "A1e_composition_by_period", 8, 5, cmp)

# A1f disturbance flags by state (share of visit intervals with the flag on)
df <- vi[, .(fire = mean(fire == 1), pest = mean(pest == 1), harv_light = mean(harv == 1), harv_heavy = mean(harv == 2), n = .N), by = .(state = state_start)]
stopifnot(all(vi$harv %in% 0:2))
dfl <- melt(df, id.vars = c("state", "n"), variable.name = "disturbance", value.name = "share")
dfl[, state := factor(state, levels = ORDER)]
p <- ggplot(dfl, aes(state, share, fill = state)) + geom_col() + facet_wrap(~ disturbance, nrow = 1, labeller = as_labeller(DIST_L)) +
  scale_fill_manual(values = PAL, guide = "none") + scale_x_discrete(labels = lab) + coord_flip() +
  labs(x = NULL, y = "share of visit intervals flagged", title = "A1f  Disturbance by state at interval start")
save_fig(p, "A1f_disturbance_by_state", 9, 4.5, dfl)

# A1g map of the plots coloured by each covariate's z-score (model scaling: data-output/<tag>/fits/scaling.rds, mean/SD over all model
#     intervals). One value per plot = mean of its visit intervals (raw), then z = (x - mean) / sd. White = model mean (z = 0), i.e. where
#     the baseline of A2/A3 sits geographically. Coordinates: PEP.gpkg layer "placette" (latitude, longitude).
con <- DBI::dbConnect(RSQLite::SQLite(), here::here("data", "raw", "pep", "PEP.gpkg"), flags = RSQLite::SQLITE_RO)
xy <- as.data.table(DBI::dbGetQuery(con, "SELECT id_pe, latitude, longitude FROM placette")); DBI::dbDisconnect(con)
sc <- readRDS(file.path(FIT_DIR, "scaling.rds"))
pc <- vi[, .(gdd_3 = mean(gdd_3, na.rm = TRUE), cmi_gs = mean(cmi_gs, na.rm = TRUE), soil_pc1 = mean(soil_pc1, na.rm = TRUE)), by = id_pe]
pc <- merge(pc, xy, by = "id_pe", all.x = TRUE)
cat(sprintf("  A1g: %d plots, %d without coordinates\n", nrow(pc), pc[is.na(latitude) | is.na(longitude), .N]))
for (v in sc$var) pc[, (paste0("z_", v)) := (get(v) - sc[var == v, mean]) / sc[var == v, sd]]
COVN <- c(z_gdd_3 = sprintf("GDD (mean %.0f, SD %.0f)", sc[var == "gdd_3", mean], sc[var == "gdd_3", sd]),
          z_cmi_gs = sprintf("CMI (mean %.1f, SD %.1f)", sc[var == "cmi_gs", mean], sc[var == "cmi_gs", sd]),
          z_soil_pc1 = sprintf("soil PC1 (mean %.2f, SD %.2f)", sc[var == "soil_pc1", mean], sc[var == "soil_pc1", sd]))
pm <- melt(pc[!is.na(latitude)], id.vars = c("id_pe", "latitude", "longitude"), measure.vars = names(COVN), variable.name = "cov", value.name = "z")
pm[, cov := factor(COVN[as.character(cov)], levels = COVN)]
setorder(pm, cov, -z)[, ord := abs(z)]; setorder(pm, cov, -ord)                 # near-mean plots drawn last (on top)
bb <- pm[, .(x = range(longitude) + c(-1, 1), y = range(latitude) + c(-0.7, 0.7))]
wm <- as.data.table(ggplot2::map_data("world", region = c("Canada", "USA")))
p <- ggplot() + geom_polygon(data = wm, aes(long, lat, group = group), fill = "grey97", colour = "grey60", linewidth = 0.2) +
  geom_point(data = pm, aes(longitude, latitude, colour = z), size = 0.25) +
  scale_colour_gradient2(low = "#2166AC", mid = "grey95", high = "#B2182B", midpoint = 0, limits = c(-2, 2), oob = scales::squish,
                         name = "SD from\nthe mean") +
  facet_wrap(~ cov, nrow = 1) + coord_quickmap(xlim = bb$x, ylim = bb$y) +
  labs(x = NULL, y = NULL, title = "A1g  Climate and soil of the plots relative to the mean",
       subtitle = sprintf("%s plots; light = near the mean (the baseline of A2–A3); colours capped at ± 2 SD", format(uniqueN(pm$id_pe), big.mark = ","))) +
  theme(panel.background = element_rect(fill = "#EAF2F8"), panel.grid = element_blank())
save_fig(p, "A1g_map_covariate_z", 15, 4.4, pc)

# ==============================================================================
# A2 BASELINE RATE: h(10 y), undisturbed, z = 0
# ==============================================================================
cat("A2 baseline\n")
b <- merge(iy[, .(label, from, to, b0_y = pooled_est)], ar[, .(label, alpha)], by = "label")
b[, h10 := alpha * 10^(alpha - 1) * exp(b0_y)]                       # per year
g <- merge(full_grid(), b, by = c("from", "to"), all.x = TRUE)
g[, `:=`(status = fifelse(is.na(h10), "notmod", "value"), fillv = h10 * 1e3, txt = formatC(signif(h10 * 1e3, 2), format = "fg", digits = 2), bold = FALSE)]
# linear colour scale in units of 1e-3 per year (user 2026-11-10): spreads the 1e-3 range, values below ~1e-4 all near the lightest colour
p <- heat(g, expression(atop(h(10), (10^{-3}~year^{-1}))), diverging = FALSE, limits = c(0, ceiling(max(b$h10) * 1e3)), breaks = 0:8,
          title = "A2  Baseline hazard at 10 y in the state, h(10)",
          subtitle = "undisturbed, mean climate and soil; each route on its own")
save_fig(p, "A2_baseline_h10", 7.5, 7, b[, .(label, from, to, b0_years = b0_y, alpha, h10_per_year = h10)])
# A2b single-route median time t50 (user 2026-11-10, replaces the characteristic time t_c = exp(-b0_y/alpha) with H = 1):
#     H(t50) = exp(b0_y) t50^alpha = ln 2  ->  t50 = exp((log(log(2)) - b0_y) / alpha); S = 0.5 on that route alone, i.e. half of the stands would
#     have taken this route if it were the only exit. Another view of the same hazards (ranks routes as t_c does); NOT a residence time (B2).
#     Undisturbed, z = 0, competing exits ignored; values > 55 y extrapolate the Weibull form.
b[, t50 := exp((log(log(2)) - b0_y) / alpha)]
g <- merge(full_grid(), b, by = c("from", "to"), all.x = TRUE)
# linear colour scale 0-1000 y (user 2026-11-10): spreads the fast routes; everything >= 1000 y gets the lightest colour (squished); cells = plain years
g[, `:=`(status = fifelse(is.na(t50), "notmod", "value"), fillv = t50, bold = FALSE,
         txt = formatC(signif(t50, 2), format = "d", big.mark = ","))]
p <- heat(g, "t50 (years)", diverging = FALSE, seq_dir = 1, limits = c(0, 1000), breaks = seq(0, 1000, 250),
          trans_lab = function(x) ifelse(x >= 1000, "\u2265 1,000", formatC(x, format = "d", big.mark = ",")),
          title = "A2b  Single-route median time, t50",
          subtitle = "years until half the stands would have taken this route if it were the only exit; undisturbed, mean climate and soil")
save_fig(p, "A2b_single_route_median_time", 7.5, 7, b[, .(label, from, to, b0_years = b0_y, alpha, t50_years = t50)])

# ==============================================================================
# A3 WEIBULL SHAPE ALPHA
# ==============================================================================
cat("A3 alpha\n")
a <- ar[, .(label, from, to, class, events, alpha, lo = alpha_ci_lo_log, hi = alpha_ci_hi_log, cap_flag)]
a[, sig := lo > 1 | hi < 1]
# A3a forest, faceted by origin
a2 <- merge(a, trow[, .(label, ylab)], by = "label")[, from := factor(from, levels = ORDER)]
a2[, ylab := factor(ylab)]
p <- ggplot(a2, aes(alpha, reorder(ylab, alpha))) + geom_vline(xintercept = 1, linetype = 2) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, colour = "grey40", orientation = "y") +
  geom_point(aes(shape = cap_flag, colour = class), size = 1.6) +
  geom_text(aes(x = hi, label = sprintf("%.2f", alpha), fontface = ifelse(sig, "bold", "plain")), hjust = -0.2, size = 2) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 2), labels = c(`FALSE` = "no", `TRUE` = "yes"), name = "numerical cap") +
  scale_x_log10(expand = expansion(mult = c(0.02, 0.15))) +
  facet_grid(from ~ ., scales = "free_y", space = "free_y", labeller = as_labeller(lab), switch = "y") +
  labs(x = "Weibull shape α (log scale)", y = "destination (events)", colour = "class",
       title = "A3a  Weibull shape α by transition", subtitle = "α > 1 = hazard rises with time in state; bold = 95% CI excludes 1") +
  theme(strip.text.y.left = element_text(angle = 0), strip.placement = "outside", panel.spacing = unit(1, "pt"), axis.text.y = element_text(size = 6))
save_fig(p, "A3a_alpha_forest", 7, 13, a)

# A3b alpha heatmap
g <- merge(full_grid(), a, by = c("from", "to"), all.x = TRUE)
g[, `:=`(status = fifelse(is.na(alpha), "notmod", "value"), fillv = log(alpha), txt = sprintf("%.2f", alpha), bold = sig %in% TRUE)]
lim <- cap_sym(g$fillv)
p <- heat(g, "α", limits = c(-lim, lim), trans_lab = function(x) sprintf("%.2f", exp(x)), breaks = log(c(0.5, 0.75, 1, 1.5, 2, 3)),
          title = "A3b  Weibull shape α", subtitle = "red = hazard rises with time in state; blue = falls; bold = 95% CI excludes 1")
save_fig(p, "A3b_alpha_heatmap", 7.5, 7, NULL)

# A3c h(t) for every transition, faceted by origin (undisturbed, z = 0)
tt <- seq(1, 50, by = 0.5)
hc <- b[, .(t = tt, h = alpha * tt^(alpha - 1) * exp(b0_y) * 1000), by = .(label, from, to)]   # h in units of 1e-3 per year
hc[, `:=`(from = factor(from, levels = ORDER), to = factor(to, levels = ORDER))]
p <- ggplot(hc, aes(t, h, colour = to)) + geom_line(linewidth = 0.5) + scale_y_log10(labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)) +
  scale_colour_manual(values = PAL, labels = lab, drop = FALSE) +
  facet_wrap(~ from, ncol = 4, scales = "free_y", labeller = as_labeller(lab)) +
  labs(x = "time in state (years)", y = expression(h(t)~~(10^{-3}~year^{-1})*",  log scale"), colour = "destination",
       title = "A3c  Baseline hazard over time in state", subtitle = "undisturbed, mean climate and soil; each route on its own")
save_fig(p, "A3c_hazard_curves", 10, 8, hc[t %in% c(1, 5, 10, 20, 30, 40, 50)])

# ==============================================================================
# A4 HAZARD RATIOS
# ==============================================================================
cat("A4 hazard ratios\n")
COVS <- c(per_sd_gdd_3 = "GDD", per_sd_cmi_gs = "CMI", per_sd_soil_pc1 = "soil PC1")   # per +1 SD (standardised); not written on the figures (user 2026-11-10)
cc <- hr[contrast %in% names(COVS), .(label, from, to, contrast, hr, hr_lo, hr_hi)]
cc[, sig := hr_lo > 1 | hr_hi < 1]
lim <- cap_sym(log(cc$hr))
# A4a heatmaps per covariate (one file each + one combined)
pl <- list()
for (k in names(COVS)) {
  g <- merge(full_grid(), cc[contrast == k], by = c("from", "to"), all.x = TRUE)
  g[, status := fifelse(label %in% modelled | paste0(from, "->", to) %in% modelled, fifelse(is.na(hr), "noterm", "value"), "notmod")]
  g[, `:=`(fillv = log(hr), txt = fmt_hr(hr), bold = sig %in% TRUE)]
  pl[[k]] <- heat(g, "hazard\nratio", limits = c(-lim, lim), trans_lab = function(x) fmt_hr(exp(x)), breaks = hr_breaks(lim),
                  title = paste("A4a  Hazard ratio:", COVS[[k]]), subtitle = "bold = 95% CI excludes 1; '–' = term not in this transition")
  save_fig(pl[[k]], paste0("A4a_hr_", sub("per_sd_", "", k)), 7.5, 7, NULL)
}
save_fig(wrap_plots(pl, ncol = 3, guides = "collect"), "A4a_hr_climate_soil_combined", 20, 7, cc)

# A4b forest per origin: climate and soil
f <- merge(cc, trow[, .(label, ylab)], by = "label")
COVS_B <- COVS
f[, `:=`(from = factor(from, levels = ORDER), covariate = factor(COVS_B[contrast], levels = COVS_B))]
# colour = direction (red HR > 1, blue HR < 1), filled = CI excludes 1; no legend (user 2026-11-10)
p <- ggplot(f, aes(hr, ylab)) + geom_vline(xintercept = 1, linetype = 2) +
  geom_errorbar(aes(xmin = hr_lo, xmax = hr_hi, colour = hr > 1), width = 0, orientation = "y") +
  geom_point(aes(colour = hr > 1, fill = ifelse(sig, ifelse(hr > 1, "up", "down"), "ns")), shape = 21, size = 1.5, stroke = 0.6) +
  scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"), guide = "none") +
  scale_fill_manual(values = c(up = "#B2182B", down = "#2166AC", ns = "white"), guide = "none") +
  scale_x_log10() + facet_grid(from ~ covariate, scales = "free_y", space = "free_y", labeller = labeller(from = lab), switch = "y") +
  labs(x = "hazard ratio (log scale)", y = "destination (events)", title = "A4b  Climate and soil hazard ratios by origin") +
  theme(strip.text.y.left = element_text(angle = 0), strip.placement = "outside", panel.spacing.y = unit(1, "pt"), axis.text.y = element_text(size = 6))
save_fig(p, "A4b_forest_climate_soil", 10, 13, NULL)

# A4c / A4d (disturbance vs none on the tau grid) moved to 14_figures_A5.R as A5e / A5f (user, 2026-11-08: combinations belong to part 5)

# A4e raw coefficients per transition: exp(beta) with the pooled 95% CI (08), one row per model term, nothing combined.
#     Reading of each row (as fitted): climate/soil = per +1 SD; disturbance main term = HR vs none at TSD = 0 (log(TSD+1) = 0) and GDD at the mean;
#     gdd_x_<d> = ratio of the +1 SD GDD HRs, disturbed vs undisturbed; tsd_<d> = per unit log(TSD + 1). alpha and intercept excluded (A2, A3).
#     Output: one multi-page PDF (one page per transition) + one PNG per transition in A4e_forest_by_transition/.
TERM_LAB <- c(z_gdd_3 = "GDD", z_cmi_gs = "CMI", z_soil_pc1 = "soil PC1",   # "(+1 SD)" dropped from row names (user 2026-11-10); stated in the caption
              fire_f = "fire", pest_f = "pest", harv_light = "light harvest", harv_heavy = "heavy harvest",
              gdd_x_fire = "GDD × fire", gdd_x_pest = "GDD × pest", gdd_x_harv_light = "GDD × light harvest", gdd_x_harv_heavy = "GDD × heavy harvest",
              tsd_fire = "TSD fire", tsd_pest = "TSD pest", tsd_harv_light = "TSD light harvest", tsd_harv_heavy = "TSD heavy harvest")
TERM_GRP <- c(rep("climate & soil", 3), rep("disturbance", 4), rep("GDD × disturbance", 4), rep("time since disturbance", 4))
names(TERM_GRP) <- names(TERM_LAB); stopifnot(setequal(names(TERM_LAB), TERM_ORDER))
rc <- pooled[parameter %in% names(TERM_LAB), .(label, from, to, parameter, beta = pooled_est, lo = ci_lo, hi = ci_hi)]
rc[, `:=`(hr = exp(beta), hr_lo = exp(lo), hr_hi = exp(hi), sig = lo > 0 | hi < 0)]
rc[, term := TERM_LAB[parameter]]
no_gdd <- setdiff(unique(rc$label), rc[parameter == "z_gdd_3", label])          # A11: no GDD main term -> gdd_x_<d> = GDD slope WITHIN level d
rc[label %in% no_gdd & grepl("^gdd_x_", parameter), term := sub("^GDD × (.*)$", "GDD within \\1 *", term)]
cat(sprintf("  A4e: transitions without a GDD main term: %s\n", paste(no_gdd, collapse = ", ")))
rc[, `:=`(term = factor(term, levels = rev(unique(c(TERM_LAB, sort(unique(term)))))), grp = factor(TERM_GRP[parameter], levels = unique(TERM_GRP)))]
rc[, txt := sprintf("%s [%s, %s]", fmt_hr(hr), fmt_hr(hr_lo), fmt_hr(hr_hi))]
fwrite(rc[, .(label, parameter, beta, ci_lo = lo, ci_hi = hi, hr, hr_lo, hr_hi, sig)], file.path(OUT, "A4e_raw_coefficients.csv"))
OUT_E <- file.path(OUT, "A4e_forest_by_transition"); dir.create(OUT_E, showWarnings = FALSE)
labs_e <- trow[order(match(from, ORDER), match(to, ORDER)), label]
forest_raw <- function(lb) {
  d <- rc[label == lb]; ev <- trow[label == lb, events]
  xr <- range(c(d$hr_lo, d$hr_hi, 1), finite = TRUE); xr <- c(xr[1] / 1.2, xr[2] * (xr[2] / xr[1])^0.45)   # room for the value column
  ggplot(d, aes(hr, term)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
    geom_errorbar(aes(xmin = hr_lo, xmax = hr_hi, colour = hr > 1), width = 0, orientation = "y") +
    geom_point(aes(colour = hr > 1, shape = sig), size = 2, fill = "white", stroke = 0.8) +
    geom_text(aes(x = xr[2], label = txt, fontface = ifelse(sig, "bold", "plain"), colour = hr > 1), hjust = 1, size = 2.6) +
    scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"), guide = "none") +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), guide = "none") +
    scale_x_log10(limits = xr, labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE, big.mark = ",", trim = TRUE)) +
    facet_grid(grp ~ ., scales = "free_y", space = "free_y", labeller = label_wrap_gen(14)) +
    labs(x = "hazard ratio (log scale)", y = NULL, title = sprintf("%s  (%d events)", gsub("_", " ", sub("->", " → ", lb)), ev),
         subtitle = "hazard ratios with 95% CI; filled and bold = CI excludes 1",
         caption = if (lb %in% no_gdd) "* no GDD main term in this transition (no events without this disturbance): the row is the GDD effect within it" else NULL) +
    theme(strip.text.y = element_text(angle = 0), plot.caption = element_text(size = 6.5, hjust = 0), plot.title.position = "plot")
}
cairo_pdf(file.path(OUT, "A4e_forest_by_transition.pdf"), width = 7, height = 5, onefile = TRUE)
for (lb in labs_e) print(forest_raw(lb))
invisible(dev.off())
for (i in seq_along(labs_e)) { lb <- labs_e[i]; n <- rc[label == lb, .N]
  ggsave(file.path(OUT_E, sprintf("%02d_%s.png", i, gsub("->", "_to_", lb))), forest_raw(lb), width = 7, height = 1.6 + 0.22 * n + 0.25 * uniqueN(rc[label == lb, grp]), dpi = 200) }
cat(sprintf("  wrote A4e_forest_by_transition.pdf (%d pages) + %d PNGs + A4e_raw_coefficients.csv\n", length(labs_e), length(list.files(OUT_E))))
cat(sprintf("\nDone: %d files in %s\n", length(list.files(OUT)), OUT))
