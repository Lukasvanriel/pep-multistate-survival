# ==============================================================================
# R/simulation/sim06b_figures.R   (no INLA; after sim06_tsd_tau.R; run from the project root)
# Figure 1: WAIC difference vs log across tau, per disturbance type x transition class, one line per imputation (summed over transitions),
#   linear form as a point at the right. Figure 2: disturbance effect (hazard ratio) over TSD, log vs exp at the class-best tau (out of regen,
#   between states; into regen left out: unreliable fits), thin = transitions (mean coefficients over imputations), thick = median.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(here) })
D <- here("data-output", "simulation", "sim06_tsd_tau"); TAU <- c(2, 3.5, 6, 10, 15, 25, 40, 60, 100, 160)
MAIN <- c(tsd_fire = "fire_f", tsd_pest = "pest_f", tsd_harv_light = "harv_light", tsd_harv_heavy = "harv_heavy")
TL <- c(tsd_fire = "fire", tsd_pest = "pest", tsd_harv_light = "light harvest", tsd_harv_heavy = "heavy harvest")
CL <- c(to_regen = "into regen", regen_exit = "out of regen", between = "between states")
theme_set(theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95")))
W <- fread(file.path(D, "sim06_waic.csv"))
E <- W[grepl("^exp_", form), .(dwaic = sum(dwaic, na.rm = TRUE)), by = .(term, class, tau, m)]
E[, `:=`(term = factor(TL[term], TL), class = factor(CL[class], CL))]
M <- E[, .(dwaic = mean(dwaic)), by = .(term, class, tau)]
L <- W[grepl("^lin_", form), .(dwaic = sum(dwaic, na.rm = TRUE) / uniqueN(m)), by = .(term, class)][, `:=`(term = factor(TL[term], TL), class = factor(CL[class], CL))]
p1 <- ggplot(E, aes(tau, dwaic)) + geom_hline(yintercept = 0, colour = "grey50") +
  geom_line(aes(group = m), colour = "grey70", linewidth = 0.3) + geom_line(data = M, colour = "#B2182B", linewidth = 0.9) + geom_point(data = M, colour = "#B2182B", size = 1) +
  geom_point(data = L, aes(x = 260, y = dwaic), shape = 17, colour = "#2166AC", size = 2) +
  scale_x_log10(breaks = c(2, 6, 15, 40, 100, 260), labels = c("2", "6", "15", "40", "100", "linear")) +
  facet_grid(class ~ term, scales = "free_y") +
  labs(x = "fade time τ of the exponential form (years, log scale); triangle = linear form", y = "WAIC minus WAIC of the log form (summed over transitions; < 0 = better than log)",
       title = "Which time-since-disturbance form fits best?", subtitle = "grey = each of the 5 imputations, red = mean; into regen: erratic jumps = failed fits (no screen)")
ggsave(file.path(D, "sim06_fig1_profiles.png"), p1, width = 11, height = 7, dpi = 200)
# figure 2
ar <- fread(here("data-output", "final", "results", "alpha_report_final.csv"))[, .(label, class)]
BEST <- data.table(class = rep(c("regen_exit", "between"), each = 4), term = rep(c("tsd_harv_heavy", "tsd_harv_light", "tsd_pest", "tsd_fire"), 2), tau = c(10, 3.5, 10, 15, 10, 40, 15, 160))
R <- rbindlist(lapply(list.files(file.path(D, "fits"), "^p1_", full.names = TRUE), function(f) { r <- readRDS(f); cl <- ar[label == r$label, class]
  if (!cl %in% BEST$class) return(NULL)
  rbindlist(lapply(intersect(names(MAIN), r$terms), function(T) { b <- BEST[class == cl & term == T]; k <- match(b$tau, TAU)
    g <- function(fit) if (is.na(fit$error)) fit$coef[match(c(MAIN[[T]], T), parameter), mean] else c(NA, NA)
    lg <- g(r$fits$log_all); ex <- g(r$fits[[sprintf("exp_%s_%03d", T, k)]])
    data.table(label = r$label, m = r$m, class = cl, term = T, tau = b$tau, lg_main = lg[1], lg_b = lg[2], ex_main = ex[1], ex_b = ex[2]) })) }))
C <- R[, lapply(.SD, mean, na.rm = TRUE), by = .(label, class, term, tau), .SDcols = c("lg_main", "lg_b", "ex_main", "ex_b")]
tg <- c(0, 0.5, 1, 2, 3, 5, 7.5, 10, 15, 20, 25, 30, 40, 50, 60)
H <- rbindlist(lapply(tg, function(t) C[, .(t = t, `log(TSD + 1)` = exp(lg_main + lg_b * log1p(t)), `exponential (best τ)` = exp(ex_main + ex_b * exp(-t / tau))), by = .(label, class, term, tau)]))
H <- melt(H, id.vars = c("label", "class", "term", "tau", "t"), variable.name = "form", value.name = "hr")
H[, `:=`(panel = sprintf("%s, %s (τ = %g y)", CL[class], TL[term], tau))]
H[, panel := factor(panel, unique(H[order(match(class, names(CL)), match(term, names(TL))), panel]))]
Hm <- H[, .(hr = median(hr)), by = .(panel, form, t)]
p2 <- ggplot(H, aes(t, hr, colour = form)) + geom_hline(yintercept = 1, colour = "grey50") +
  geom_line(aes(group = interaction(label, form)), alpha = 0.15, linewidth = 0.3) + geom_line(data = Hm, linewidth = 1) +
  scale_y_log10() + coord_cartesian(ylim = c(0.05, 30)) + facet_wrap(~ panel, ncol = 4) +
  scale_colour_manual(values = c(`log(TSD + 1)` = "#2166AC", `exponential (best τ)` = "#B2182B"), name = NULL) +
  labs(x = "years since disturbance", y = "hazard ratio, disturbed vs undisturbed (log scale)", title = "Disturbance effect over time: log form vs bounded exponential",
       subtitle = "thin = transitions (mean over 5 imputations), thick = median; into regen not shown (unreliable fits)") + theme(legend.position = "bottom")
ggsave(file.path(D, "sim06_fig2_effect_curves.png"), p2, width = 12, height = 6.5, dpi = 200)
cat("figures done\n")
