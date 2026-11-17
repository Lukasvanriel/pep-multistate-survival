# ==============================================================================
# 16b_figures_B_plots.R   (seconds; run after 16_figures_B.R; from the project root)
# Draws figure set B from data-output/<tag>/figures/B/B_results.rds. Definitions and scenarios: header of 16_figures_B.R. NOT manuscript figures.
# ==============================================================================
source(here::here("R", "13_fig_helpers.R"))
OUT <- file.path(FIG_DIR, "B"); R <- readRDS(file.path(OUT, "B_results.rds")); res <- R$res; share <- R$share
DL <- c(none = "undisturbed", fire = "fire", pest = "pest", harv_light = "light harvest", harv_heavy = "heavy harvest")
DCOL <- c(none = "grey30", fire = "#D95F02", pest = "#7570B3", harv_light = "#66A61E", harv_heavy = "#1B6E3C")
olab <- function(x) gsub("_", " ", x); fo <- function(x) factor(x, levels = ORDER)
flag_lab <- function(o, d) { s <- share[origin == o & dist == d]; if (!nrow(s)) return(""); sprintf("%.0f%% assumed%s", 100 * s$share_assumed, if (s$flag) " *" else "") }
pct <- function(x) sprintf("%.0f", 100 * x)

# ---- B0 model check ---------------------------------------------------------------------------------------------------------------------------------------------
b0 <- copy(R$B0o)[, origin := fo(origin)]
tot <- b0[, .(n = sum(n), o = sum(obs_change), p = sum(pred_change)), by = grp][, lab := sprintf("%s: n = %d, observed %d, predicted %.1f (ratio %.2f)", grp, n, o, p, o / p)]
p <- ggplot(b0, aes(y = origin)) + geom_linerange(aes(xmin = obs_lo, xmax = obs_hi), colour = "grey50") +
  geom_point(aes(x = obs_p, shape = "observed (95% CI)"), size = 2) + geom_point(aes(x = pred_p, shape = "predicted by the model"), size = 2.4, colour = "#B2182B") +
  geom_text(aes(x = 1, label = sprintf("n = %d", n)), hjust = 1, size = 2.3, colour = "grey40") +
  scale_shape_manual(values = c(`observed (95% CI)` = 16, `predicted by the model` = 4), name = NULL) + scale_y_discrete(labels = olab, limits = rev) +
  scale_x_continuous(limits = c(0, 1)) + facet_wrap(~ grp) +
  labs(x = "share of intervals with a state change", y = "origin state", title = "B0  Model check: observed vs predicted state changes",
       subtitle = paste(c("undisturbed intervals", tot$lab), collapse = "\n"))
save_fig(p, "B0_model_check", 11, 6, NULL)

# ---- B1 where stands go and how fast ---------------------------------------------------------------------------------------------------------------------
cv_all <- melt(R$curves, id.vars = c("t", "origin", "dist", "zg"), variable.name = "dest", value.name = "p", na.rm = TRUE)
cv <- cv_all[zg == 0]   # B1a / B1c: GDD 0
cv[, dest := factor(fifelse(dest == "stay", "stays", as.character(dest)), levels = c("stays", ORDER))]
PALS <- c(stays = "grey88", PAL)
# one dashed line at the end of the record (55 y) replaces the 55-100 y shading and the 50-y line (user 2026-11-10)
area <- function(d, ttl, sub) ggplot(d, aes(t, p, fill = dest)) +
  geom_area(colour = "white", linewidth = 0.1) + geom_vline(xintercept = 55, linetype = 2, linewidth = 0.4) + scale_fill_manual(values = PALS, labels = olab, name = "state") +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) + labs(x = "years in the origin state", y = "probability", title = ttl, subtitle = sub)
p <- area(cv[dist == "none"][, origin := fo(origin)], "B1a  First move out of each state (undisturbed)",
          "dashed = end of record (55 y)") +
  facet_wrap(~ origin, ncol = 5, labeller = as_labeller(olab))
save_fig(p, "B1a_first_move_curves", 16, 9.5, NULL)
for (d in names(DL)[-1]) { x <- cv[dist == d]; if (!nrow(x)) next; x[, origin := fo(origin)]
  fl <- unique(x[, .(origin)])[, lab := paste0(olab(origin), "\n", mapply(flag_lab, as.character(origin), d))]
  p <- area(x, sprintf("B1c  First move out of each state after %s", DL[[d]]),
            "panel: share of the exit rate assumed unaffected by the disturbance (* = > 25%); dashed = end of record (55 y)") +
    facet_wrap(~ origin, ncol = 5, labeller = as_labeller(setNames(fl$lab, as.character(fl$origin))))
  save_fig(p, paste0("B1c_first_move_curves_", d), 16, 2.2 * ceiling(nrow(fl) / 5) + 1.8, NULL) }
# B1d (user 2026-11-10): one figure per origin, all scenarios side by side: undisturbed GDD -1 / 0 / +1, then each disturbance at entry (GDD 0).
#   Destinations whose exit has no term for that disturbance (kept at the undisturbed rate, option c) are drawn faded; a disturbance with no term
#   in any exit of the origin gets an empty panel saying so. One fill legend per figure.
TERM1 <- c(fire = "fire_f", pest = "pest_f", harv_light = "harv_light", harv_heavy = "harv_heavy")
has_t <- function(o, j, d) paste0(o, "->", j) %in% pooled[parameter == TERM1[[d]], label]
SCN <- c(`none_-1` = "undisturbed\nGDD -1 SD", none_0 = "undisturbed\nGDD 0", none_1 = "undisturbed\nGDD +1 SD",
         fire_0 = "fire", pest_0 = "pest", harv_light_0 = "light harvest", harv_heavy_0 = "heavy harvest")
for (o in ORDER) { x <- cv_all[origin == o][, scen := factor(SCN[paste0(dist, "_", zg)], levels = SCN)]
  x[, dest := factor(fifelse(dest == "stay", "stays", as.character(dest)), levels = c("stays", ORDER))]
  x[, faded := FALSE][dist != "none" & dest != "stays", faded := !mapply(has_t, o, as.character(dest), dist)]
  miss <- data.table(scen = factor(SCN[setdiff(names(SCN), unique(paste0(x$dist, "_", x$zg)))], levels = SCN))
  miss[, lab := paste0("no ", sub("\\n.*", "", as.character(scen)), " term")]
  p <- ggplot(x, aes(t, p, fill = dest, alpha = faded)) + geom_area(colour = "white", linewidth = 0.1) +
    geom_vline(data = unique(x[, .(scen)]), aes(xintercept = 55), linetype = 2, linewidth = 0.4) +   # no line on empty panels (user 2026-11-10)
    scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.3), guide = "none") + scale_fill_manual(values = PALS, labels = olab, name = "state") +
    scale_y_continuous(labels = function(v) paste0(100 * v, "%")) + facet_wrap(~ scen, nrow = 1, drop = FALSE) +
    labs(x = "years in the starting state", y = "probability", title = sprintf("B1  First move out of %s", olab(o)),
         subtitle = "disturbance at entry, GDD 0; faded = exit with no term for that disturbance (kept at the undisturbed rate); dashed = end of record (55 y)")
  if (nrow(miss)) p <- p + geom_text(data = miss, aes(x = 50, y = 0.5, label = lab), inherit.aes = FALSE, size = 3, colour = "grey40")
  save_fig(p, paste0("B1d_first_move_", o), 17, 3.6, NULL) }
heat_p <- function(d, ttl, sub, xl = "state at 50 y (diagonal = stayed)") { g <- CJ(origin = ORDER, dest = ORDER); g <- merge(g, d, by = c("origin", "dest"), all.x = TRUE)
  g[, `:=`(origin = factor(origin, levels = rev(ORDER)), dest = factor(dest, levels = ORDER))]
  g[, txt := fifelse(is.na(est), "", sprintf("%s\n[%s-%s]", pct(est), pct(lo), pct(hi)))]
  g[, tc := fifelse(!is.na(est) & est > 0.55, "white", "black")]
  ggplot(g, aes(dest, origin)) + geom_tile(aes(fill = est), colour = "white") + geom_text(aes(label = txt, colour = tc), size = 1.9, lineheight = 0.85) +
    scale_colour_identity() + scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 1), na.value = "grey95", labels = function(x) paste0(100 * x, "%"), name = "probability") +
    scale_x_discrete(labels = olab, position = "top") + scale_y_discrete(labels = olab) + coord_equal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 0), panel.grid = element_blank()) + labs(x = xl, y = "origin", title = ttl, subtitle = sub) }
b1 <- rbind(res[scen == "none" & zg == 0 & quantity == "cif50", .(origin, dest, est, lo, hi)], res[scen == "none" & zg == 0 & quantity == "stay50", .(origin, dest = origin, est, lo, hi)])
save_fig(heat_p(b1, "B1b  State 50 years after entering a state (undisturbed)", "% [95% CI]; rows sum to 100%"),
         "B1b_probability_at_50y", 9, 8.5, b1)

# ---- B2 residence time -------------------------------------------------------------------------------------------------------------------------------------------
# (user 2026-11-10) "TSD fixed at 10 y" dropped from the figure (differs from "at entry" by <= 4 y; stays in B_quantities.csv); undisturbed GDD 0 in black;
#   climate and disturbance rows as two stacked plots, each with its own legend; 50-y line only on the mean panel (cap), 55-y line only on the median panel.
QL <- c(rmst50 = "years spent in the state, out of the first 50")   # median dropped; RMST50 renamed (user 2026-11-10)
r2 <- rbind(res[scen == "none" & quantity %in% names(QL), .(origin, scenario = gsd(zg), group = "climate", quantity, est, lo, hi, flag = FALSE)],
            res[scen == "dist" & mode == "entry" & zg == 0 & quantity %in% names(QL), .(origin, scenario = DL[dist], group = "disturbance", quantity, est, lo, hi, flag)])
r2[, `:=`(origin = factor(origin, levels = rev(ORDER)), quantity = factor(QL[quantity], levels = QL))]
xr <- r2[, .(x = range(c(lo, hi, est), na.rm = TRUE)), by = quantity]                       # same x range per column in both rows
vl <- data.table(quantity = factor(QL, levels = QL), x = 50, lt = "dotted")
b2p <- function(d, cols, nm, strips) ggplot(d, aes(est, origin, colour = scenario)) +
  geom_vline(data = vl, aes(xintercept = x, linetype = lt), colour = "grey50", show.legend = FALSE) + scale_linetype_identity() +
  geom_blank(data = xr, aes(x = x), inherit.aes = FALSE) +
  geom_errorbar(aes(xmin = lo, xmax = hi, group = scenario), width = 0, orientation = "y", position = position_dodge(width = 0.75), na.rm = TRUE) +
  geom_point(aes(shape = flag, group = scenario), position = position_dodge(width = 0.75), size = 1.6, fill = "white", na.rm = TRUE) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), breaks = "TRUE", labels = c(`TRUE` = "> 25% of exit rate\nassumed unaffected"), name = NULL, drop = FALSE) +   # triangle = flag (user 2026-11-10)
  scale_colour_manual(values = cols, name = nm) + facet_wrap(~ quantity, scales = "free_x") + scale_y_discrete(labels = olab) + labs(x = NULL, y = NULL) +
  (if (!strips) theme(strip.text = element_blank(), strip.background = element_blank()) else theme())
CLIM_COL <- setNames(c("#2166AC", "black", "#B2182B"), gsd(c(-1, 0, 1)))
d1 <- r2[group == "climate"][, scenario := factor(scenario, levels = names(CLIM_COL))]
d2 <- r2[group == "disturbance"][, scenario := factor(scenario, levels = DL[-1])]
p <- b2p(d1, CLIM_COL, "undisturbed", TRUE) / (b2p(d2, setNames(DCOL[-1], DL[-1]), "disturbance at entry\n(GDD 0)", FALSE) + labs(x = "years")) +
  plot_annotation(title = "B2  How long stands stay in a state",
                  subtitle = "from entry into the state; dotted = 50 (stayed the whole time)")
save_fig(p, "B2_residence_time", 9, 11, r2)

# ---- B3 next state -------------------------------------------------------------------------------------------------------------------------------------------------
for (zz in c(-1, 0, 1)) { x <- res[scen == "none" & zg == zz & quantity == "next50", .(origin, dest, est, lo, hi)]
  save_fig(heat_p(x, sprintf("B3  Next state, for stands that moved within 50 y (undisturbed, %s)", gsd(zz)), "% [95% CI]; rows sum to 100%", xl = "next state"),
           sprintf("B3_next_state_gdd%+d", zz), 9, 8.5, x) }

# ---- B4 along gradients -----------------------------------------------------------------------------------------------------------------------------------------
g4 <- res[scen == "gdd_gradient" & quantity == "next50"][, `:=`(origin = fo(origin), dest = factor(dest, levels = ORDER))]
p <- ggplot(g4, aes(zg, est, fill = dest)) + geom_area(colour = "white", linewidth = 0.1) + facet_wrap(~ origin, ncol = 5, labeller = as_labeller(olab)) +
  scale_fill_manual(values = PAL, labels = olab, name = "next state") + scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(x = "GDD (SD from the mean)", y = "probability of the next state", title = "B4a  Next state along GDD (undisturbed)")
save_fig(p, "B4a_next_state_vs_gdd", 16, 9.5, g4[, .(origin, zg, dest, est, lo, hi)])
s4 <- res[scen == "gdd_gradient" & quantity == "shift50"][, origin := fo(origin)]
p <- ggplot(s4, aes(zg, est)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey50") + geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#F4A582", alpha = 0.6) + geom_line(colour = "#B2182B") +
  facet_wrap(~ origin, ncol = 5, labeller = as_labeller(olab)) +
  labs(x = "GDD (SD from the mean)", y = "expected ΔSTI of the next state", title = "B4b  Thermal direction of the first move along GDD (undisturbed)",
       subtitle = "> 0 = next state on average warmer-affinity; band = 95% CI")
save_fig(p, "B4b_thermal_shift_vs_gdd", 16, 9.5, s4[, .(origin, zg, est, lo, hi)])

# ---- B5 under disturbance ----------------------------------------------------------------------------------------------------------------------------------------------
for (d in names(DL)[-1]) { x <- res[scen == "dist" & dist == d & mode == "entry" & zg == 0 & quantity == "next50", .(origin, dest, est, lo, hi)]; if (!nrow(x)) next
  fl <- share[dist == d & flag == TRUE & any_term == TRUE, origin]
  save_fig(heat_p(x, sprintf("B5  Next state after %s", DL[[d]]),
                  sprintf("%% [95%% CI]; > 25%% of the exit rate assumed unaffected: %s", if (length(fl)) paste(olab(fl), collapse = ", ") else "none"), xl = "next state"),
           paste0("B5_next_state_", d), 9, 8.5, x) }

# shapes (user 2026-11-10): filled / open = 95% CI excludes / includes 0; triangle = > 25% of the exit rate assumed unaffected (circle otherwise)
SHK <- c(`FALSE_TRUE` = 16, `FALSE_FALSE` = 21, `TRUE_TRUE` = 17, `TRUE_FALSE` = 24)
SHL <- c(`FALSE_TRUE` = "CI excludes 0", `FALSE_FALSE` = "CI includes 0", `TRUE_TRUE` = "> 25% assumed unaffected,\nCI excludes 0", `TRUE_FALSE` = "> 25% assumed unaffected,\nCI includes 0")
shk <- function(d) d[, sk := factor(paste(flag %in% TRUE, lo > 0 | hi < 0, sep = "_"), levels = names(SHK))]
# ---- B6 climate vs disturbance on a common scale ---------------------------------------------------------------------------------------------------------------
# (user 2026-11-10, option b) climate block on top, disturbance block below, each with its own legend; alternating grey bands per state;
#   same x range per column in both blocks; quantity labels aligned with B2.
b6 <- copy(R$B6)[, origin := factor(origin, levels = rev(ORDER))]
Q6 <- c(rmst = "change in years in the state, out of the first 50", stay = "change in P(still in the state at 50 y)", shift = "change in thermal direction of the first move (ΔSTI)")
b6[, quantity := factor(Q6[quantity], levels = Q6)]
CLC <- c(`GDD +1 SD` = "#B2182B", `CMI +1 SD` = "#2166AC", `soil PC1 +1 SD` = "#8C510A")
b6[, grp := fifelse(contrast %in% names(CLC), "climate", "disturbance")]
x6 <- b6[, .(x = range(c(lo, hi, est), na.rm = TRUE)), by = quantity]
bands <- data.table(origin = factor(levels(b6$origin)[seq(2, nlevels(b6$origin), by = 2)], levels = levels(b6$origin)))
b6p <- function(d, cols, nm, strips) ggplot(d, aes(est, origin, colour = contrast)) +
  geom_tile(data = bands, aes(x = 0, y = origin), width = Inf, height = 1, fill = "grey92", colour = NA, inherit.aes = FALSE) +
  geom_blank(data = x6, aes(x = x), inherit.aes = FALSE) + geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi, group = contrast), width = 0, orientation = "y", position = position_dodge(width = 0.7)) +
  geom_point(aes(shape = sk, group = contrast), position = position_dodge(width = 0.7), size = 1.7, fill = "white") +
  scale_shape_manual(values = SHK, labels = SHL, name = NULL, drop = FALSE) +
  scale_colour_manual(values = cols, name = nm) + scale_y_discrete(labels = olab, drop = FALSE) + facet_wrap(~ quantity, scales = "free_x") + labs(x = NULL, y = NULL) +
  theme(panel.grid.major.y = element_blank()) + (if (!strips) theme(strip.text = element_blank(), strip.background = element_blank()) else theme())
shk(b6)
d1 <- b6[grp == "climate"][, contrast := factor(contrast, levels = names(CLC))]
d2 <- b6[grp == "disturbance"][, contrast := factor(contrast, levels = DL[-1])]
p <- b6p(d1, CLC, "climate", TRUE) / b6p(d2, setNames(DCOL[-1], DL[-1]), "disturbance\nat entry", FALSE) +
  plot_annotation(title = "B6  Climate vs disturbance on a common scale",
                  subtitle = "change vs an undisturbed stand at mean climate; 95% CI")
save_fig(p, "B6_climate_vs_disturbance", 15, 12, NULL)

# ---- B7 interaction on the probability scale ------------------------------------------------------------------------------------------------------------------
b7 <- copy(R$B7)[, `:=`(origin = factor(origin, levels = rev(ORDER)), dist = factor(dist, levels = names(DL)[-1]))]
b7[quantity == "stay", `:=`(est = -est, lo2 = -hi, hi2 = -lo)][quantity == "stay", `:=`(lo = lo2, hi = hi2)][, c("lo2", "hi2") := NULL]
shk(b7)
Q7 <- c(stay = "P(left the state by 50 y)", rmst = "years in the state, out of the first 50", shift = "thermal direction of the first move (ΔSTI)")
b7[, quantity := factor(Q7[quantity], levels = Q7)]
p <- ggplot(b7, aes(est, origin, colour = dist)) +
  geom_tile(data = bands, aes(x = 0, y = origin), width = Inf, height = 1, fill = "grey92", colour = NA, inherit.aes = FALSE) +   # bands as B6 (user 2026-11-10)
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi, group = dist), width = 0, orientation = "y", position = position_dodge(width = 0.7)) +
  geom_point(aes(shape = sk, group = dist), position = position_dodge(width = 0.7), size = 1.7, fill = "white") +
  scale_shape_manual(values = SHK, labels = SHL, name = NULL, drop = FALSE) +
  scale_colour_manual(values = DCOL[-1], labels = DL[-1], name = "disturbance") + scale_y_discrete(labels = olab, drop = FALSE) + facet_wrap(~ quantity, scales = "free_x") +
  theme(panel.grid.major.y = element_blank()) +
  labs(x = "effect of +1 SD GDD, disturbed − undisturbed", y = NULL,
       title = "B7  GDD × disturbance interaction on the outcome scale",
       subtitle = "0 = same warming effect with and without disturbance; 95% CI")
save_fig(p, "B7_interaction_outcome_scale", 15, 8, NULL)
cat(sprintf("\nDone: %d B files in %s\n", length(list.files(OUT, pattern = "^B[0-9].*(png)$")), OUT))
