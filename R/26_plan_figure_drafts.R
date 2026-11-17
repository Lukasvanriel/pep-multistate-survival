# ==============================================================================
# 26_plan_figure_drafts.R   (no INLA; <1 min; run from the project root)
# Draft panels for the chapter plan that do not exist yet: F4 (null model with class intercepts), F7 (hindcast slopes with both placebo ranges),
#   F8 (two equilibria + no-harvest). Other plan figures are montages of existing panels (made with ImageMagick in the plan step). Output figures/plan.
# ==============================================================================
source(here::here("R", "25_I_common.R")); OUTP <- file.path(FIG_DIR, "plan"); dir.create(OUTP, showWarnings = FALSE); OUT <- OUTP
# ---- F4: null model ----
src <- readLines(here::here("R", "24_null_model_gdd_directionality.R")); e <- new.env(); eval(parse(text = src[1:(grep("^t2 <- lm", src) - 1)]), envir = e)
obs <- e$obs; arc <- fread(file.path(RESULTS_DIR, sprintf("alpha_report_%s.csv", RUN_TAG)))[, .(label, class)]
obs <- merge(obs, arc, by = "label"); CL <- c(between = "between established states", to_regen = "into regen", regen_exit = "out of regen")
obs[, cls := factor(CL[class], CL)]; m2 <- lm(log_hr ~ 0 + cls + cls:null, obs, weights = 1 / se^2)
gr <- obs[, .(null = seq(min(null), max(null), length.out = 50)), by = cls]; gr[, fit := predict(m2, gr)]
pa <- e$pr; stopifnot(!is.null(pa), "se_D" %in% names(pa))
# panel a colour (user 2026-11-12): red / blue = observed difference significantly above / below detailed balance (|z| > 1.96, z = (D_obs - D_null) / se_D; se of the observed difference only, as 24)
SIDE_COL <- c("warm-ward excess" = "#B2182B", "cold-ward excess" = "#2166AC", "consistent with detailed balance" = "grey55")
pa[, side := factor(fifelse((D_obs - D_null) / se_D > 1.96, "warm-ward excess", fifelse((D_obs - D_null) / se_D < -1.96, "cold-ward excess", "consistent with detailed balance")), names(SIDE_COL))]
CC <- c(`between established states` = "grey35", `into regen` = "#2166AC", `out of regen` = "#B2182B")
p4a <- ggplot(pa, aes(D_null, D_obs)) + geom_abline(linetype = 2, colour = "grey50") + geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_errorbar(aes(ymin = D_obs - 1.96 * se_D, ymax = D_obs + 1.96 * se_D), width = 0, colour = "grey70") + geom_point(aes(colour = side), size = 1.8) +
  scale_colour_manual(values = SIDE_COL, name = "deviation from detailed balance (95%)", drop = FALSE) + theme(legend.position = "bottom") +
  labs(x = "GDD effect on the number ratio, warmer ÷ colder state (log)", y = "GDD effect on the rate ratio, move to warmer ÷ move to colder (log)", subtitle = "a  pairs with moves both ways: dashed = rates hold today's distribution (detailed balance)")
p4b <- ggplot(obs, aes(null, log_hr, colour = cls)) + geom_hline(yintercept = 0, colour = "grey85") + geom_point(aes(size = 1 / se), alpha = 0.6) +
  geom_line(data = gr, aes(y = fit), linewidth = 0.9) + scale_colour_manual(values = CC, name = NULL) + scale_size_continuous(range = c(0.6, 3), guide = "none") +
  labs(x = "GDD effect on the number ratio, destination ÷ origin (log)", y = "GDD effect on the transition rate (log hazard ratio)", subtitle = "b  all transitions: height at zero = warmth speeding the move whatever its direction")
save_fig(p4a + p4b + plot_annotation(title = "F4  Climate direction is sorting; climate pace is regen closure"), "F4_null_model_draft", 12, 5, NULL)
# ---- F7: hindcast with placebos ----
H <- fread(file.path(FIG_DIR, "H", "H_within_plot_slopes_controls.csv"))[period_FE == TRUE]
P1 <- unique(fread(file.path(FIG_DIR, "I", "I9bc_placebo_controls.csv"))[, .(origin, metric, w_lo = placebo_min, w_hi = placebo_max)])
P2 <- fread(file.path(FIG_DIR, "I", "I9d_placebo_crossplot.csv"))[, .(origin, metric, c_lo = placebo_min, c_hi = placebo_max)]
X <- merge(merge(H, P1, by = c("origin", "metric")), P2, by = c("origin", "metric"))
ML <- c(chg = "any state change", dsti = "thermal shift", tmp = "end in a temperate state", bor = "end in a boreal state")
X[, `:=`(metric = factor(ML[metric], ML), origin = factor(origin, rev(c("all", "boreal", "pioneer", "regen", "temperate"))))]; X[, y := as.numeric(origin)]
p7 <- ggplot(X) + geom_vline(xintercept = 0, linetype = 2, colour = "grey50") + geom_vline(xintercept = 1, colour = "grey60") +
  geom_rect(aes(xmin = c_lo, xmax = c_hi, ymin = y - 0.32, ymax = y + 0.32, fill = "another plot's warming (same region and period)")) +
  geom_rect(aes(xmin = w_lo, xmax = w_hi, ymin = y - 0.16, ymax = y + 0.16, fill = "the plot's own warming in random order")) +
  geom_errorbarh(aes(xmin = lo, xmax = hi, y = y), height = 0, colour = "black") + geom_point(aes(slope, y), size = 2) +
  scale_fill_manual(values = c("another plot's warming (same region and period)" = "#F4C7A1", "the plot's own warming in random order" = "#9ECAE1"), name = "range of slopes when each interval gets") +
  scale_y_continuous(breaks = 1:5, labels = levels(X$origin)) + facet_wrap(~ metric, nrow = 1) +
  labs(x = "within-plot slope (1 = spatial effects predict the plot's change as it warmed)", y = "origin group",
       title = "F7  Do climate effects learned across space predict change in time?", subtitle = "points = real slope with 95% CI; bands = range over 20 random reassignments of the warming") +
  theme(legend.position = "bottom")
save_fig(p7, "F7_hindcast_placebo_draft", 13, 4.8, NULL)
# ---- F8: two equilibria ----
G <- fread(file.path(FIG_DIR, "I", "I13_group_shares.csv"))
GL <- melt(G, id.vars = c("bin", "zg", "grp"), variable.name = "src")[src %in% c("obs", "I8", "SC", "SC_noharv", "SC_harvonly")]   # light-only / heavy-only versions (added to I13 later) not shown
SL <- c(obs = "observed now", I8 = "equilibrium, no disturbance", SC = "equilibrium, observed disturbance regime", SC_noharv = "equilibrium, fire and pest only", SC_harvonly = "equilibrium, harvest only")
GL[, `:=`(src = factor(SL[as.character(src)], SL), xlab = sprintf("%s\n(GDD %+.1f)", bin, zg))]
GL[, grp := factor(grp, c("boreal", "pioneer", "temperate", "regen"))]   # pioneer added (user 2026-11-16)
p8 <- ggplot(GL, aes(xlab, value, colour = src, group = src)) + geom_line() + geom_point(size = 1.5) +
  facet_wrap(~ grp, nrow = 1) + scale_colour_manual(values = c("black", "#B2182B", "#1B6E3C", "#7FBF7B", "#E08214"), name = NULL) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(x = "plots binned by mean GDD (quintiles)", y = "share of stands", title = "F8  Two equilibria and the landscape between them",
       subtitle = "1000-y composition implied by the rates; a direction, not a forecast") + theme(legend.position = "bottom", legend.text = element_text(size = 11), axis.text = element_text(size = 10), axis.text.x = element_text(size = 9), axis.title = element_text(size = 12), strip.text = element_text(size = 12), plot.title = element_text(size = 13), plot.subtitle = element_text(size = 11))
save_fig(p8, "F8_two_equilibria_draft", 17, 5.2, NULL)
cat("26 done\n")
