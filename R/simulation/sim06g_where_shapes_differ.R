# R/simulation/sim06g_where_shapes_differ.R  (no INLA; 2026-11-12) -- for four transitions: observed rate ratio + fitted shapes (top),
#   Exposure in the rows is in decades (TIME_SCALE); bars show stand-years (x 10).
#   and how much data each band of years since disturbance has: stand-years watched (bars) and observed changes (numbers) (bottom).
source(here::here("R", "simulation", "sim06f_observed_vs_fitted.R"))
suppressPackageStartupMessages(library(patchwork))
P4 <- levels(OBS$panel)[c(1, 4, 2, 5)]
o <- OBS[panel %in% P4]; o[, panel := factor(panel, P4)]; cu <- CUR[panel %in% P4]; cu[, panel := factor(panel, P4)]
o[, `:=`(x0 = BR[as.integer(bin)], x1 = pmin(BR[as.integer(bin) + 1], 75))]
top <- ggplot() + geom_hline(yintercept = 1, colour = "grey50") +
  geom_errorbar(data = o[events > 0], aes(x = t_mid, ymin = pmax(lo, 1e-3), ymax = hi), width = 0, colour = "grey40") +
  geom_point(data = o[events > 0], aes(t_mid, rr), colour = "grey20", size = 2) +
  geom_line(data = cu, aes(t, hr, colour = shapeL), linewidth = 1) + scale_y_log10() + coord_cartesian(xlim = c(0, 70)) +
  facet_wrap(~ panel, nrow = 1, scales = "free_y") + scale_colour_manual(values = SC, name = NULL) +
  labs(x = NULL, y = "rate ratio vs undisturbed\n(log scale)") + theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top", strip.background = element_rect(fill = "grey95"))
bot <- ggplot(o) + geom_rect(aes(xmin = x0, xmax = x1, ymin = 0, ymax = 10 * expo), fill = "grey75", colour = "white") +
  geom_text(aes(x = (x0 + x1) / 2, y = 10 * expo, label = events), vjust = -0.3, size = 3) + coord_cartesian(xlim = c(0, 70)) +
  facet_wrap(~ panel, nrow = 1, scales = "free_y") + scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "years since disturbance", y = "stand-years watched\n(number = changes seen)") + theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), strip.text = element_blank(), strip.background = element_blank())
ggsave(file.path(S2, "sim06g_where_shapes_differ.png"), top / bot + plot_layout(heights = c(2, 1)), width = 13, height = 6.5, dpi = 180)
print(o[, .(panel, bin, stand_years = round(10 * expo), events, rr = signif(rr, 2))])
