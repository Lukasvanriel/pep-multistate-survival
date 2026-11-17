# ==============================================================================
# 23c_extended_figures.R   (no INLA; after 23b; run from the project root). Figures E1-E7 (exploratory, NOT manuscript figures).
# Median return time is NOT shown: inside a 100-y window it is set by the window (returns accumulate up to 100 y), not by the dynamics.
# ==============================================================================
source(here::here("R", "13_fig_helpers.R")); OUTE <- file.path(FIG_DIR, "E"); OUT <- OUTE
x <- readRDS(file.path(OUTE, "E_summaries.rds")); m <- readRDS(file.path(OUTE, "cache", "meta.rds")); S <- m$S
SL <- c(none_zm1 = "undisturbed, GDD -1 SD", none_z0 = "undisturbed, GDD 0", none_zp1 = "undisturbed, GDD +1 SD", fire = "fire", pest = "pest", harv_light = "light harvest", harv_heavy = "heavy harvest")
PL <- c(gdd_p1 = "GDD +1 SD", gdd_m1 = "GDD -1 SD", fire = "fire", pest = "pest", harv_light = "light harvest", harv_heavy = "heavy harvest")
lab_s <- function(v) gsub("_", " ", v)
# E1
E1 <- fread(file.path(OUTE, "E1_next_state_by_exit_time.csv"))
d <- E1[!grepl("^<", dest) & scen == "none_z0"]; d[, dest := factor(dest, levels = ORDER)]; d[, origin := factor(origin, levels = ORDER)]
p1 <- ggplot(d, aes(t, est, fill = dest)) + geom_area(colour = "white", linewidth = 0.1) + geom_vline(xintercept = 55, linetype = 2, linewidth = 0.4) +
  facet_wrap(~ origin, labeller = as_labeller(lab_s), ncol = 5) + scale_fill_manual(values = PAL, labels = lab_s, name = "next state") +
  labs(x = "time in state at exit (years)", y = "probability of the next state", title = "E1a  Where a stand goes depends on when it leaves (undisturbed, GDD 0)", subtitle = "dashed = end of record (55 y)")
save_fig(p1, "E1a_next_state_by_exit_time", 15, 9, NULL)
d2 <- E1[dest == "<expected dSTI>" & scen %in% c("none_zm1", "none_z0", "none_zp1")]; d2[, origin := factor(origin, levels = ORDER)]
p2 <- ggplot(d2, aes(t, est, colour = scen, fill = scen)) + geom_hline(yintercept = 0, colour = "grey60") + geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line() +
  facet_wrap(~ origin, labeller = as_labeller(lab_s), ncol = 5, scales = "free_y") + scale_colour_manual(values = c(none_zm1 = "#2166AC", none_z0 = "grey30", none_zp1 = "#B2182B"), labels = SL, name = NULL) +
  scale_fill_manual(values = c(none_zm1 = "#2166AC", none_z0 = "grey30", none_zp1 = "#B2182B"), labels = SL, name = NULL) + geom_vline(xintercept = 55, linetype = 2, linewidth = 0.4) +
  labs(x = "time in state at exit (years)", y = "expected ΔSTI of the next state", title = "E1b  Thermal direction of the next move by exit time (undisturbed)", subtitle = "band = 95% CI; dashed = end of record (55 y)")
save_fig(p2, "E1b_thermal_shift_by_exit_time", 15, 9, NULL)
# E2
e2 <- x$E2[outcome == "thermal index at 50 y (multi-step)" & start == "observed composition"]; e2[, `:=`(pert = factor(PL[pert], PL), part = factor(part, c("holding", "jump", "total")))]
p3 <- ggplot(e2, aes(est, pert, colour = part)) + geom_vline(xintercept = 0, colour = "grey60") + geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", position = position_dodge(width = 0.6)) +
  geom_point(position = position_dodge(width = 0.6), size = 2) + scale_colour_manual(values = c(holding = "#E08214", jump = "#542788", total = "black"), labels = c(holding = "how long (holding)", jump = "which way (jump)", total = "total"), name = NULL) +
  labs(x = "change in the thermal index of the composition at 50 y", y = NULL, title = "E2a  How long vs which way: split of each driver's effect", subtitle = "from the observed first-visit composition; 95% CI")
e2s <- x$E2[outcome == "thermal index at 50 y (multi-step)" & start != "observed composition" & part != "total"]; e2s[, `:=`(pert = factor(PL[pert], PL), start = factor(start, rev(ORDER)))]
p4 <- ggplot(e2s, aes(est, start, colour = part)) + geom_vline(xintercept = 0, colour = "grey60") + geom_point(position = position_dodge(width = 0.5)) + facet_wrap(~ pert, nrow = 1) +
  scale_colour_manual(values = c(holding = "#E08214", jump = "#542788"), labels = c(holding = "how long", jump = "which way"), name = NULL) + scale_y_discrete(labels = lab_s) + labs(x = "change in the thermal index at 50 y", y = "start state", subtitle = "by start state")
save_fig(p3 / p4 + plot_layout(heights = c(1, 2)), "E2_holding_vs_jump", 15, 10, NULL)
# E3
e3 <- x$E3[outcome %in% c("thermal index", "share boreal", "share temperate", "share regen")]
top <- e3[outcome == "thermal index", .(s = max(abs(central_per_10pct))), by = label][order(-s)][1:25, label]
e3t <- e3[label %in% top]; e3t[, label := factor(label, rev(top))]
e3t[, param := factor(fifelse(param == "b0", "baseline hazard ±10%", "Weibull shape α ±10%"), c("baseline hazard ±10%", "Weibull shape α ±10%"))]
p5 <- ggplot(e3t, aes(outcome, label, fill = central_per_10pct)) + geom_tile() + facet_wrap(~ param) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", name = "effect of\n+10%") + labs(x = NULL, y = NULL,
  title = "E3  Which transitions control the composition at 50 y", subtitle = "25 most influential for the thermal index; compare within a panel, not across")
save_fig(p5, "E3_elasticity", 13, 9, NULL)
# E4
e4 <- x$E4[horizon == 50 & scen == "none_z0"]; e4[, `:=`(start = factor(start, rev(ORDER)), state = factor(state, ORDER))]
p6 <- ggplot(e4, aes(state, start, fill = est)) + geom_tile() + geom_text(aes(label = ifelse(est >= 0.5, sprintf("%.0f", est), "")), size = 2.6) +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "years") + scale_x_discrete(labels = lab_s) + scale_y_discrete(labels = lab_s) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "state", y = "start state", title = "E4  Expected years in each state within 50 y (undisturbed, GDD 0)", subtitle = "returns included")
save_fig(p6, "E4_expected_years", 11, 8, NULL)
# E5
e5 <- x$E5[start != "observed composition" & quantity %in% c("state changes within 100 y", "P(return to start within 100 y)")]; e5[, `:=`(start = factor(start, rev(ORDER)), scen = factor(SL[scen], SL))]
p7 <- ggplot(e5, aes(est, start, colour = scen)) + geom_point(position = position_dodge(width = 0.7), size = 1.3) + facet_wrap(~ quantity, scales = "free_x") + scale_y_discrete(labels = lab_s) +
  scale_colour_manual(values = c("#2166AC", "grey30", "#B2182B", "#D95F02", "#7570B3", "#66A61E", "#1B6E3C"), name = NULL) + labs(x = NULL, y = NULL, title = "E5  Turnover and returns within 100 y")
save_fig(p7, "E5_turnover_returns", 13, 7, NULL)
# E6
e6 <- fread(file.path(OUTE, "E6_distance_to_long_run.csv")); e6[, scen := factor(SL[scen], SL[1:3])]
p8 <- ggplot(e6[t <= 500], aes(t, tvd, colour = scen)) + geom_hline(yintercept = c(0.05, 0.10), linetype = 3) + geom_vline(xintercept = 55, linetype = 2, colour = "grey60") + geom_line() +
  scale_colour_manual(values = c("#2166AC", "grey30", "#B2182B"), name = NULL) + labs(x = "years", y = "distance to the long-run composition", title = "E6  Convergence from the observed first-visit composition (undisturbed)",
  subtitle = "±1 SD = adjustment of today's composition to a shifted climate; dashed = end of record (55 y)")
save_fig(p8, "E6_convergence", 10, 5.5, NULL)
# E7
e7 <- x$E7[horizon == 50 & quantity %in% c("scored moves to MORE tolerant", "scored moves to LESS tolerant", "moves into regen")]; e7[, scen := factor(SL[scen], rev(SL))]
p9 <- ggplot(e7, aes(est, scen, colour = quantity)) + geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", position = position_dodge(width = 0.6)) + geom_point(position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = c("#1B7837", "#C51B7D", "grey40"), name = NULL) + labs(x = "expected moves per stand within 50 y", y = NULL,
  title = "E7  Successional direction (from the observed first-visit composition)", subtitle = "shade tolerance from Niinemets & Valladares (2006); moves involving regen, other pioneer or other temperate unscored")
save_fig(p9, "E7_successional_direction", 11, 5, NULL)
cat("23c figures done\n")
