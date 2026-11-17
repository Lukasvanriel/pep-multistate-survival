# 20b_D3a_shaded.R (no INLA; 2026-11-13, user request) -- redraws D3a from D3_averaged_effects.csv (no recomputation) with alternating row shading,
#   then rebuilds the F6 draft montage (D3a + D2a + I14, stacked at 2400 px width) for the proposed-route folder.
source(here::here("R", "13_fig_helpers.R"))
OUTD <- file.path(FIG_DIR, "D"); OUT <- OUTD; SUM <- fread(file.path(OUTD, "D3_averaged_effects.csv"))
PERT_LV <- c("GDD +1 SD", "CMI +1 SD", "soil PC1 +1 SD", "fire", "pest", "light harvest", "heavy harvest")
LV <- c("all plots", rev(ORDER)); SUM[, `:=`(state = factor(state, levels = rev(LV)), pert = factor(pert, levels = PERT_LV))]
BANDS <- data.table(i = seq(2, length(LV), by = 2))
PC <- c(`GDD +1 SD` = "#B2182B", `CMI +1 SD` = "#2166AC", `soil PC1 +1 SD` = "#8C510A", fire = "#D95F02", pest = "#7570B3", `light harvest` = "#66A61E", `heavy harvest` = "#1B6E3C")
mk <- function(v, lo, hi, xl, yl = NULL) ggplot(SUM, aes(.data[[v]], state, colour = pert)) + geom_blank() +
  geom_rect(data = BANDS, aes(ymin = i - 0.5, ymax = i + 0.5), xmin = -Inf, xmax = Inf, fill = "grey92", colour = NA, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = .data[[lo]], xmax = .data[[hi]]), width = 0, orientation = "y", position = position_dodge(width = 0.75)) +
  geom_point(size = 1.5, position = position_dodge(width = 0.75)) + scale_colour_manual(values = PC, name = NULL) +
  scale_y_discrete(labels = function(x) gsub("_", " ", x)) + labs(x = xl, y = yl) + theme(panel.grid.major.y = element_blank())
save_fig(mk("d_pleave", "pl_lo", "pl_hi", "change in P(leaving the state within 50 y)", "starting state (state at the last visit)") + mk("d_shift", "sh_lo", "sh_hi", "change in expected ΔSTI of the next state") +
         plot_layout(guides = "collect") + plot_annotation(title = "D3a  Effects averaged over the real plots",
         subtitle = "each plot with its own climate and soil, state at its last visit; disturbance at entry; 95% CI"),
         "D3a_averaged_effects", 15, 8, NULL)
f <- function(x) shQuote(x); P <- function(...) file.path(FIG_DIR, ...)
out <- here::here("docs", "chapter2_proposed_route", "figures", "main_figures", "F6_pace_direction_draft.png")
system(paste("convert", "\\(", f(P("D", "D3a_averaged_effects.png")), "-resize 2400x \\)", "\\( -size 2400x150 xc:white \\)",
             "\\(", f(P("D", "D2a_group_direction.png")), "-resize 2400x \\)", "\\( -size 2400x150 xc:white \\)",
             "\\(", f(P("I", "I14_first_passage_temperate.png")), "-resize 2400x \\)", "-background white -append", f(out)))
cat("done\n")
