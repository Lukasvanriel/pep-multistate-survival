# ==============================================================================
# 15_figures_A6.R   (no INLA; seconds; run from the project root: Rscript R/15_figures_A6.R)
# Exploratory figure set A, part 6: the cross-transition synthesis of 12 (meta-analysis / meta-regression across transitions). NOT manuscript figures.
# Needs: 12 run (synthesis_inputs_<tag>.csv = the per-transition estimates entering each model; synthesis_<tag>.csv = the fitted models).
# Writes: data-output/<tag>/figures/A/A6*  (png/pdf; csv where the plotted values are not already in 12's outputs)
# Every figure shows the per-transition estimates (points; size ~ 1/se) and 12's fitted model for the "all" subset (since 2026-11-10; the "robust"
#   subset, >= 25 events, >= 50 for alpha, no low_events label, is no longer drawn: see SHOW_SUB below). Lines of meta-regressions are drawn from 12's coefficients
#   (no confidence band: 12 does not store the intercept-slope covariance); slope, CI and p are printed in each panel.
# Block B: inputs whose tau lies outside the 5-95th pct of TSD at the event (09 tsd_ev_p05..p95) are drawn with a cross (x).
# Circularity: STI/SMI/SSI are computed from the same plots (12 header); a positive climate slope along dIndex is partly expected by construction.
# ==============================================================================
source(here::here("R", "13_fig_helpers.R"))
inp <- fread(file.path(RESULTS_DIR, sprintf("synthesis_inputs_%s.csv", RUN_TAG)))
res <- fread(file.path(RESULTS_DIR, sprintf("synthesis_%s.csv", RUN_TAG)))
# 2026-11-10 (user): figures show the "all" subset only. The robust subset (>= 25 events, >= 50 for alpha, no low_events label) gave the same
#   conclusions everywhere (precision weighting already down-weights sparse transitions); it stays in synthesis_<tag>.csv and MUST be reported
#   in the text as a sensitivity check ("results unchanged when restricted to transitions with >= 25 events"). Set SHOW_SUB to c("all", "robust") to draw both.
SHOW_SUB <- c("all"); res_both <- res; res <- res[subset %in% SHOW_SUB]
inp[, `:=`(hr = exp(log_hr), w = 1 / se, rob = fifelse(robust, "robust", "not robust"))]
SUB_SH <- c(all = 16, robust = 21); SUB_LT <- c(all = "solid", robust = "dashed")
CLS <- c(regen_exit = "regen → X", between = "X → Y", to_regen = "X → regen")
DL <- c(fire = "fire", pest = "pest", harv_light = "light harvest", harv_heavy = "heavy harvest")
DCOL <- c(fire = "#D95F02", pest = "#7570B3", harv_light = "#66A61E", harv_heavy = "#1B6E3C")
COVL <- c(gdd = "GDD", cmi = "CMI", soil = "soil PC1")   # "(+1 SD)" dropped from headers (user 2026-11-10); unit on the y axis
IDXL <- c(gdd = "ΔSTI (thermal)", cmi = "ΔSMI (moisture)", soil = "ΔSSI (soil)")
ylog <- function(...) scale_y_log10(labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE), ...)
mods <- function(mod) { r <- res[model == mod & !is.na(est)]; r[, `:=`(m = exp(est), m_lo = exp(ci_lo), m_hi = exp(ci_hi))]; r[] }
slope_txt <- function(mod, trm) { r <- res[model == mod & term == trm]
  paste(sapply(SHOW_SUB, function(s) { x <- r[subset == s]; if (!nrow(x) || is.na(x$est)) return(paste(s, ": NA"))
    sprintf("%s×%.2f (%.2f–%.2f), %s", if (length(SHOW_SUB) > 1) paste0(s, ": ") else "", exp(x$est), exp(x$ci_lo), exp(x$ci_hi), ptxt(x$p)) }), collapse = "\n") }
ptxt <- function(p) if (p < 0.001) "p < 0.001" else paste("p =", formatC(signif(p, 2), format = "fg", digits = 2))
# per-transition points + model means (both subsets) on a categorical x
cell_plot <- function(pts, xvar, mdl, xlabs = waiver(), colour = NULL) {
  g <- ggplot(pts, aes(.data[[xvar]], hr)) + geom_hline(yintercept = 1, linetype = 2, colour = "grey50")
  g <- if (is.null(colour)) g + geom_point(aes(size = w), colour = "grey60", alpha = 0.5, position = position_jitter(width = 0.15, height = 0, seed = 1))
       else g + geom_point(aes(size = w, colour = .data[[colour]]), alpha = 0.45, position = position_jitter(width = 0.15, height = 0, seed = 1))
  g + geom_errorbar(data = mdl, aes(x = .data[[xvar]], ymin = m_lo, ymax = m_hi, group = subset), inherit.aes = FALSE, width = 0, linewidth = 0.7, position = position_dodge(width = 0.5)) +
    geom_point(data = mdl, aes(x = .data[[xvar]], y = m, shape = subset), inherit.aes = FALSE, size = 2.6, fill = "white", stroke = 1, position = position_dodge(width = 0.5)) +
    scale_shape_manual(values = SUB_SH, name = "model (12)", guide = if (length(SHOW_SUB) > 1) "legend" else "none") + scale_size_continuous(range = c(0.4, 3), guide = "none") + scale_x_discrete(labels = xlabs) + ylog()
}

# ---- A6.0 state indices ------------------------------------------------------------------------------------------------------------------------------------
sc <- readRDS(file.path(FIT_DIR, "scaling.rds")); zf <- function(x, v) (x - sc[var == v, mean]) / sc[var == v, sd]
idx <- vi[, .(STI = mean(zf(gdd_3, "gdd_3")), SMI = mean(zf(cmi_gs, "cmi_gs")), SSI = mean(zf(soil_pc1, "soil_pc1")), n = .N), by = .(state = state_end)]   # as 12
vz <- vi[, .(state = state_end, STI = zf(gdd_3, "gdd_3"), SMI = zf(cmi_gs, "cmi_gs"), SSI = zf(soil_pc1, "soil_pc1"))]
il <- melt(vz, id.vars = "state", variable.name = "index")[, .(mean = mean(value, na.rm = TRUE), q05 = quantile(value, .05, na.rm = TRUE), q25 = quantile(value, .25, na.rm = TRUE),
          q75 = quantile(value, .75, na.rm = TRUE), q95 = quantile(value, .95, na.rm = TRUE), n = .N), by = .(state, index)]
stopifnot(all.equal(il[index == "STI"][order(state), mean], idx[order(state), STI]))
il[, state := factor(state, levels = idx[order(STI), state])]
p <- ggplot(il, aes(y = state, colour = state)) + geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_linerange(aes(xmin = q05, xmax = q95), linewidth = 0.4) + geom_linerange(aes(xmin = q25, xmax = q75), linewidth = 2, alpha = 0.6) +
  geom_point(aes(x = mean), size = 2.4, colour = "black") +
  scale_colour_manual(values = PAL, guide = "none") + scale_y_discrete(labels = lab) +
  facet_wrap(~ index, nrow = 1, labeller = as_labeller(c(STI = "STI (GDD)", SMI = "SMI (CMI)", SSI = "SSI (soil PC1)"))) +
  labs(x = "z-score", y = NULL, title = "A6.0  Climate and soil index of each state",
       subtitle = "dot = mean (the index); thick bar = interquartile range; thin line = 5–95th percentile")
save_fig(p, "A6_0_state_indices", 10, 5, il)

# ---- A6A climate ------------------------------------------------------------------------------------------------------------------------------------------
cl <- inp[block == "A_climate"][, `:=`(cov = factor(sub, levels = names(COVL)), cls = factor(class, levels = names(CLS)))]
mA1 <- rbindlist(lapply(names(COVL), function(v) mods(sprintf("A1_%s_class_means", v))[, cov := v]))
mA1[, `:=`(cls = factor(sub("^class", "", term), levels = names(CLS)), cov = factor(cov, levels = names(COVL)))]
p <- cell_plot(cl, "cls", mA1, CLS) + facet_wrap(~ cov, nrow = 1, labeller = as_labeller(COVL)) +
  labs(x = NULL, y = "hazard ratio (log scale)", title = "A6A1  Climate and soil effects by transition class",
       subtitle = "grey = transitions (size ~ precision); black = class mean")
save_fig(p, "A6A1_climate_class_means", 11, 4.5, NULL)
reg_plot <- function(pts, xvar, mod_of, trm, xl) {
  ln <- rbindlist(lapply(names(COVL), function(v) { r <- res[model == mod_of(v) & !is.na(est)]; if (!nrow(r)) return(NULL)
    rbindlist(lapply(SHOW_SUB, function(s) { b <- r[subset == s]; if (!nrow(b)) return(NULL)
      xs <- range(pts[cov == v, get(xvar)], na.rm = TRUE); xx <- seq(xs[1], xs[2], length.out = 50)
      data.table(cov = v, subset = s, x = xx, y = exp(b[term == "intrcpt", est] + b[term == trm, est] * xx)) })) }))
  ln[, cov := factor(cov, levels = names(COVL))]
  tx <- data.table(cov = factor(names(COVL), levels = names(COVL)), lab = sapply(names(COVL), function(v) slope_txt(mod_of(v), trm)))
  ggplot(pts, aes(.data[[xvar]], hr)) + geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
    geom_point(aes(size = w, shape = rob), colour = "grey30", fill = "white", alpha = 0.7) +
    geom_line(data = ln, aes(x, y, linetype = subset), linewidth = 0.8, colour = "#B2182B") +
    geom_text(data = tx, aes(label = lab), x = -Inf, y = Inf, hjust = -0.03, vjust = 1.2, size = 2.3, inherit.aes = FALSE) +
    scale_shape_manual(values = c(robust = 16, `not robust` = 21), name = "transition") + scale_linetype_manual(values = SUB_LT, name = "model (12)", guide = if (length(SHOW_SUB) > 1) "legend" else "none") +
    scale_size_continuous(range = c(0.4, 3), guide = "none") + ylog(expand = expansion(mult = c(0.05, 0.35))) +
    facet_wrap(~ cov, nrow = 1, scales = "free_x", labeller = as_labeller(COVL)) + labs(x = xl, y = "hazard ratio (log scale)")
}
p <- reg_plot(cl[class == "between"], "dIndex", function(v) sprintf("A2_%s_between_vs_dIndex", v), "dIndex",
              "Δ index, destination − origin (STI / SMI / SSI)") +
  labs(title = "A6A2  Between-state transitions: climate effect vs direction of the move",
       subtitle = "x > 0 = move toward a warmer / wetter / higher-soil-index state")
save_fig(p, "A6A2_climate_vs_dIndex", 13, 5, NULL)
p <- reg_plot(cl[class == "regen_exit"], "Index_to", function(v) sprintf("A3_%s_regenexit_vs_destIndex", v), "Index_to", "index of the destination state") +
  labs(title = "A6A3  Regen exits: climate effect vs the destination's index")
save_fig(p, "A6A3_regen_exits_vs_destIndex", 13, 5, NULL)

# ---- A6B disturbance vs none (tau grid) -------------------------------------------------------------------------------------------------------------------
db <- inp[block == "B_disturbance"][, `:=`(dist = sub, cls = factor(class, levels = names(CLS)))]
rng <- unique(hr[!is.na(tsd_ev_p05), .(label, dist, tsd_ev_p05, tsd_ev_p95)])
db <- merge(db, rng, by = c("label", "dist"), all.x = TRUE)
db[, outside := !is.na(tsd_ev_p05) & (tau_y < tsd_ev_p05 | tau_y > tsd_ev_p95)]
db[, `:=`(tau_lab = factor(sprintf("%d y since disturbance", tau_y), levels = sprintf("%d y since disturbance", c(1, 10, 30, 60))), dist = factor(dist, levels = names(DL)))]
fwrite(db[, .(label, dist, tau_y, class, hr, se, robust, dIndex, outside, tsd_ev_p05, tsd_ev_p95)], file.path(OUT, "A6B_inputs_with_range_flag.csv"))
mB1 <- rbindlist(lapply(c(1, 10, 30, 60), function(tt) mods(sprintf("B1_disturbance_cell_means_tau%dy", tt))[, tau_y := tt]))
mB1[, c("dist", "cls") := tstrsplit(sub("^cell", "", term), ":", fixed = TRUE)]
mB1[, `:=`(tau_lab = factor(sprintf("%d y since disturbance", tau_y), levels = levels(db$tau_lab)), cls = factor(cls, levels = names(CLS)), dist = factor(dist, levels = names(DL)))]
p <- ggplot(db, aes(cls, hr)) + geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_point(aes(size = w, colour = dist, shape = outside), alpha = 0.5, position = position_jitter(width = 0.15, height = 0, seed = 1)) +
  geom_errorbar(data = mB1, aes(x = cls, ymin = m_lo, ymax = m_hi, group = subset), inherit.aes = FALSE, width = 0, linewidth = 0.7, position = position_dodge(width = 0.5)) +
  geom_point(data = mB1, aes(x = cls, y = m, fill = subset, group = subset), inherit.aes = FALSE, shape = 21, size = 2.4, stroke = 1, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = DCOL, labels = DL, guide = "none") + scale_fill_manual(values = c(all = "black", robust = "white"), name = "model (12)", guide = if (length(SHOW_SUB) > 1) "legend" else "none") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), labels = c(`FALSE` = "within", `TRUE` = "outside"), name = "observed TSD range") +
  scale_size_continuous(range = c(0.4, 3), guide = "none") + scale_x_discrete(labels = CLS) + ylog(breaks = 10^(-2:3)) +   # label 1 and 100 too (user 2026-11-10)
  facet_grid(dist ~ tau_lab, labeller = labeller(dist = DL)) + theme(strip.text.y = element_text(angle = 0)) +
  labs(x = NULL, y = "hazard ratio vs undisturbed (log scale)", title = "A6B1  Disturbance effect by transition class",
       subtitle = "points = transitions; black = cell mean")
save_fig(p, "A6B1_disturbance_class_means", 14, 11, NULL)
lnB <- rbindlist(lapply(c(1, 10, 30, 60), function(tt) { r <- res[model == sprintf("B2_disturbance_between_vs_dSTI_tau%dy", tt) & !is.na(est)]
  rbindlist(lapply(SHOW_SUB, function(s) { b <- r[subset == s]; if (!nrow(b)) return(NULL); bb <- setNames(b$est, b$term)
    rbindlist(lapply(names(DL), function(d) { off <- if (paste0("dist", d) %in% names(bb)) bb[[paste0("dist", d)]] else if (d == "fire") 0 else NA_real_
      if (is.na(off) || !nrow(db[tau_y == tt & dist == d & class == "between"])) return(NULL)
      xs <- range(db[tau_y == tt & dist == d & class == "between", dIndex]); xx <- seq(xs[1], xs[2], length.out = 30)
      data.table(tau_y = tt, subset = s, dist = d, x = xx, y = exp(bb[["intrcpt"]] + off + bb[["dSTI"]] * xx)) })) })) }))
lnB[, `:=`(tau_lab = factor(sprintf("%d y since disturbance", tau_y), levels = levels(db$tau_lab)), dist = factor(dist, levels = names(DL)))]
txB <- data.table(tau_lab = factor(levels(db$tau_lab), levels = levels(db$tau_lab)),
                  lab = sapply(c(1, 10, 30, 60), function(tt) slope_txt(sprintf("B2_disturbance_between_vs_dSTI_tau%dy", tt), "dSTI")))
p <- ggplot(db[class == "between"], aes(dIndex, hr)) + geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_point(aes(size = w, colour = dist, shape = outside), alpha = 0.55) +
  geom_line(data = lnB, aes(x, y, colour = dist, linetype = subset), linewidth = 0.7) +
  geom_text(data = txB, aes(label = lab), x = -Inf, y = Inf, hjust = -0.03, vjust = 1.2, size = 2.1, inherit.aes = FALSE) +
  scale_colour_manual(values = DCOL, labels = DL, name = "disturbance") + scale_linetype_manual(values = SUB_LT, name = "model (12)", guide = if (length(SHOW_SUB) > 1) "legend" else "none") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), labels = c(`FALSE` = "within", `TRUE` = "outside"), name = "observed TSD range") +
  scale_size_continuous(range = c(0.4, 3), guide = "none") + ylog(breaks = 10^(-2:3), expand = expansion(mult = c(0.05, 0.4))) + facet_wrap(~ tau_lab, nrow = 1) +
  labs(x = "ΔSTI, destination − origin", y = "hazard ratio vs undisturbed (log scale)",
       title = "A6B2  Between-state transitions: disturbance effect vs ΔSTI")
save_fig(p, "A6B2_disturbance_vs_dSTI", 16, 5.5, NULL)

# ---- A6C interactions ---------------------------------------------------------------------------------------------------------------------------------------
ci <- inp[block == "C_interaction"][, `:=`(dist = factor(sub, levels = names(DL)), cls = factor(class, levels = names(CLS)))]
m1 <- mods("C1_interaction_overall")[, `:=`(grp = "overall", lev = "all")]
m2 <- mods("C2_interaction_by_dist")[, `:=`(grp = "by disturbance", lev = sub("^cell", "", term))]
m3 <- mods("C3_interaction_by_class")[, `:=`(grp = "by class", lev = sub("^cell", "", term))]
mC <- rbind(m1, m2, m3, fill = TRUE)
LEVL <- c(all = "all transitions", DL, CLS)
mC[, `:=`(grp = factor(grp, levels = c("overall", "by disturbance", "by class")), lev = factor(LEVL[lev], levels = rev(unique(LEVL))))]
pts <- rbind(ci[, .(grp = "overall", lev = "all", hr, w)], ci[, .(grp = "by disturbance", lev = as.character(dist), hr, w)], ci[, .(grp = "by class", lev = as.character(cls), hr, w)])
pts[, `:=`(grp = factor(grp, levels = levels(mC$grp)), lev = factor(LEVL[lev], levels = levels(mC$lev)))]
pC1 <- ggplot(pts, aes(hr, lev)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_point(aes(size = w), colour = "grey60", alpha = 0.4, position = position_jitter(height = 0.15, width = 0, seed = 1)) +
  geom_errorbar(data = mC, aes(y = lev, xmin = m_lo, xmax = m_hi, group = subset), inherit.aes = FALSE, width = 0, linewidth = 0.7, orientation = "y", position = position_dodge(width = 0.5)) +
  geom_point(data = mC, aes(x = m, y = lev, shape = subset), inherit.aes = FALSE, size = 2.6, fill = "white", stroke = 1, position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = SUB_SH, name = "model (12)", guide = if (length(SHOW_SUB) > 1) "legend" else "none") + scale_size_continuous(range = c(0.4, 3), guide = "none") +
  scale_x_log10(labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)) +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y") + labs(x = "GDD hazard ratio, disturbed / undisturbed (log scale)", y = NULL,
  subtitle = "< 1 = GDD effect weaker in disturbed stands") + theme(strip.text.y = element_text(angle = 0))
ln4 <- rbindlist(lapply(SHOW_SUB, function(s) { b <- res[model == "C4_interaction_between_vs_dSTI" & subset == s & !is.na(est)]; if (!nrow(b)) return(NULL)
  xs <- range(ci[class == "between", dIndex]); xx <- seq(xs[1], xs[2], length.out = 30); data.table(subset = s, x = xx, y = exp(b[term == "intrcpt", est] + b[term == "dSTI", est] * xx)) }))
pC4 <- ggplot(ci[class == "between"], aes(dIndex, hr)) + geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_point(aes(size = w, colour = dist), alpha = 0.55) + geom_line(data = ln4, aes(x, y, linetype = subset), linewidth = 0.8) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.03, vjust = 1.2, size = 2.3, label = slope_txt("C4_interaction_between_vs_dSTI", "dSTI")) +
  scale_colour_manual(values = DCOL, labels = DL, name = "disturbance") + scale_linetype_manual(values = SUB_LT, name = "model (12)", guide = if (length(SHOW_SUB) > 1) "legend" else "none") +
  scale_size_continuous(range = c(0.4, 3), guide = "none") + ylog(expand = expansion(mult = c(0.05, 0.3))) +
  labs(x = "ΔSTI, destination − origin", y = "GDD hazard ratio, disturbed / undisturbed (log scale)", subtitle = "between-state transitions, vs ΔSTI")
save_fig(pC1 + pC4 + plot_layout(widths = c(1.1, 1)) + plot_annotation(title = "A6C  GDD × disturbance interaction",
         caption = "points = transitions (size ~ precision); black point and line = meta-regression estimate"),
         "A6C_interactions", 14, 6.5, NULL)

# ---- A6D TSD coefficient ------------------------------------------------------------------------------------------------------------------------------------
dt <- inp[block == "D_tsd"][, `:=`(dist = factor(sub, levels = names(DL)), cls = factor(class, levels = names(CLS)))]
mD <- mods("D1_tsd_cell_means"); mD[, c("dist", "cls") := tstrsplit(sub("^cell", "", term), ":", fixed = TRUE)]
mD[, `:=`(dist = factor(dist, levels = names(DL)), cls = factor(cls, levels = names(CLS)))]
p <- cell_plot(dt, "cls", mD, CLS) + facet_wrap(~ dist, nrow = 1, labeller = as_labeller(DL)) +
  labs(x = NULL, y = "TSD hazard ratio (log scale)", title = "A6D  Time-since-disturbance effect by disturbance type and class",
       subtitle = "< 1 = disturbance effect fades with time; grey = transitions; black = cell mean")
save_fig(p, "A6D_tsd_cell_means", 13, 4.5, NULL)

# ---- A6E alpha ----------------------------------------------------------------------------------------------------------------------------------------------
al <- inp[block == "E_alpha"][, alpha := exp(log_hr)]
GRP <- c("regen", "intolerant", "intermediate", "tolerant")
ae <- al[!is.na(tol_group) & tol_group != ""][, grpf := factor(tol_group, levels = GRP)]
mE1 <- mods("E1_logalpha_tolgroup_means")[, grpf := factor(sub("^tol_group", "", term), levels = GRP)]
ae[, hr := alpha]
pE1 <- cell_plot(ae, "grpf", mE1) + geom_hline(yintercept = 1, linetype = 1, colour = "grey30", linewidth = 0.3) +
  labs(x = "shade tolerance of the origin state", y = "Weibull shape α (log scale)", subtitle = sprintf("difference between groups: %s",
       ptxt(res[model == "E1_logalpha_tolgroup_diff" & subset == "all", QMp][1])))
mE2 <- mods("E2_logalpha_origin_means")[, origin := sub("^cell", "", term)]
ordE <- mE2[subset == "all"][order(m), origin]
mE2[, origin := factor(origin, levels = ordE)]; a2 <- al[from %in% ordE][, `:=`(origin = factor(from, levels = ordE), hr = alpha)]
pE2 <- cell_plot(a2, "origin", mE2, lab) + coord_flip() + labs(x = "origin state", y = "Weibull shape α (log scale)")
save_fig(pE1 + pE2 + plot_layout(widths = c(1, 1.1), guides = "collect") +
         plot_annotation(title = "A6E  Weibull shape α across transitions", caption = "grey = transitions (size ~ precision); black = mean; α = 1: hazard constant with time in state"),
         "A6E_alpha_groups", 13, 5.5, NULL)
cat(sprintf("\nDone (part 6): %d A6 files in %s\n", length(list.files(OUT, pattern = "^A6")), OUT))
