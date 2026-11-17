# 25_I13f_tsd_shape_figure.R (no INLA; after 25_I13e; 2026-11-12) -- figure: shift caused by the disturbance regime, per group and GDD bin, log vs exp A vs exp B.
source(here::here("R", "25_I_common.R"))
G <- fread(file.path(OUT, "I13e_tsd_shape_B_group_shares.csv"))
G[, shape := factor(form, c("log_B", "exp_A", "exp_B"), c("log (current)", "exponential, one τ per type (A)", "exponential, τ per type × class (B)"))]
G[, grp := factor(grp, c("boreal", "temperate", "pioneer", "regen"))]
p <- ggplot(G, aes(bin, 100 * shift, colour = shape, group = shape)) + geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(linewidth = 0.9, position = position_dodge(0.25)) + geom_point(size = 2, position = position_dodge(0.25)) +
  facet_wrap(~ grp, nrow = 1) + scale_colour_manual(values = c("black", "#F4A582", "#B2182B"), name = NULL) +
  labs(x = "plots binned by mean GDD (cold Q1 → warm Q5)", y = "change in share of stands\n(percentage points)",
       title = "Effect of the disturbance regime on equilibrium composition, by TSD shape",
       subtitle = "equilibrium with the observed disturbance regime minus equilibrium without disturbances; same 5 imputations and random numbers for all three shapes") +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), legend.position = "bottom", strip.background = element_rect(fill = "grey95"))
ggsave(file.path(OUT, "I13f_tsd_shape.png"), p, width = 11, height = 4.2, dpi = 180)
