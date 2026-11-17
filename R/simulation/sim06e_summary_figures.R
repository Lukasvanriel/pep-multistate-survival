# ==============================================================================
# R/simulation/sim06e_summary_figures.R   (no INLA; after sim06 and sim06d; run from the project root)
# Combines sim06 (out of regen, between states; no screen needed: smooth profiles) and sim06d (into regen; main-run screen, valid pairs only).
# Fig 4: fit of linear and of the exponential across tau vs log, per disturbance type x class.  Fig 5: class-median effect over TSD, 3 shapes
# (log, linear, exponential at the class-best tau).  Fig 6: individual transitions, 3 shapes.  Table: hazard ratio at 1/10/30/60 y.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(here); library(patchwork) })
S1 <- here("data-output", "simulation", "sim06_tsd_tau"); S2 <- here("data-output", "simulation", "sim06d_regen_screen"); OUT <- S2
TAU <- c(2, 3.5, 6, 10, 15, 25, 40, 60, 100, 160)
MAIN <- c(tsd_fire = "fire_f", tsd_pest = "pest_f", tsd_harv_light = "harv_light", tsd_harv_heavy = "harv_heavy")
TL <- c(tsd_fire = "fire", tsd_pest = "pest", tsd_harv_light = "light harvest", tsd_harv_heavy = "heavy harvest")
CL <- c(to_regen = "into regen", regen_exit = "out of regen", between = "between states")
theme_set(theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95")))
ar <- fread(here("data-output", "final", "results", "alpha_report_final.csv"))[, .(label, class)]
# ---- fit profiles: sim06 (non-regen-entry classes) + sim06d (into regen, valid pairs) ----
W1 <- fread(file.path(S1, "sim06_waic.csv"))[class != "to_regen" & grepl("^(exp|lin)_", form)]
P1 <- W1[, .(d = sum(dwaic, na.rm = TRUE) / uniqueN(m)), by = .(term, class, form, tau)]
W2 <- merge(fread(file.path(S2, "sim06d_waic.csv")), ar, by = "label")
P2 <- W2[both_valid == TRUE, .(d = mean(dwaic)), by = .(label, term, class, form, tau)][, .(d = sum(d)), by = .(term, class, form, tau)]
P <- rbind(P1, P2); P[, shape := fifelse(grepl("^lin_", form), "linear", "exponential")]
BT <- P[shape == "exponential", .SD[which.min(d)], by = .(term, class)][, .(term, class, tau_best = tau, d_exp = d)]
BT <- merge(BT, P[shape == "linear", .(term, class, d_lin = d)], by = c("term", "class")); fwrite(BT, file.path(OUT, "sim06e_best_tau_by_class.csv"))
P[, `:=`(termL = factor(TL[term], TL), classL = factor(CL[class], CL))]
pf <- ggplot(P[shape == "exponential"], aes(tau, d)) + geom_hline(yintercept = 0, colour = "black") + geom_line(colour = "#B2182B", linewidth = 0.8) + geom_point(colour = "#B2182B", size = 1.2) +
  geom_point(data = P[shape == "linear"], aes(x = 260, y = d), shape = 17, colour = "#2166AC", size = 2.3) +
  scale_x_log10(breaks = c(2, 6, 15, 40, 100, 260), labels = c("2", "6", "15", "40", "100", "linear")) + facet_grid(classL ~ termL, scales = "free_y") +
  labs(x = "fade time τ of the exponential (years, log scale); blue triangle = linear", y = "WAIC minus WAIC of log, summed over transitions (< 0 = fits better than log)",
       title = "How well each shape fits, compared with the current log form (the black line at 0)",
       subtitle = "red = exponential at each τ; into regen: only fits passing the main run's check, 10 imputations")
ggsave(file.path(OUT, "sim06e_fig4_fit.png"), pf, width = 11, height = 7, dpi = 180)
# ---- coefficients per transition: log, linear, exponential at the class-best tau ----
get_coefs <- function(dir, pat, screened) rbindlist(lapply(list.files(file.path(dir, "fits"), pat, full.names = TRUE), function(f) { r <- readRDS(f)
  cl <- ar[label == r$label, class]
  rbindlist(lapply(intersect(names(MAIN), r$terms), function(T) { tb <- BT[term == T & class == cl]; if (!nrow(tb)) return(NULL)
    fm <- c(log = "log_all", linear = paste0("lin_", T), exponential = sprintf("exp_%s_%03d", T, match(tb$tau_best, TAU)))
    rbindlist(lapply(names(fm), function(s) { x <- r$fits[[fm[[s]]]]; ok <- !is.null(x$coef) && is.na(x$error) && (!screened || isTRUE(x$screen$valid))
      data.table(label = r$label, m = r$m, class = cl, term = T, shape = s, tau = tb$tau_best, ok = ok,
                 main = if (ok) x$coef[parameter == MAIN[[T]], mean] else NA_real_, b = if (ok) x$coef[parameter == T, mean] else NA_real_) })) })) }))
K <- rbind(get_coefs(S1, "^p1_", FALSE)[class != "to_regen"], get_coefs(S2, "^m[0-9]+_", TRUE))
K <- K[ok == TRUE, .(main = mean(main), b = mean(b), n_imp = .N), by = .(label, class, term, shape, tau)]
fcurve <- function(s, t, tau) switch(s, log = log1p(t), linear = t / 10, exponential = exp(-t / tau))
tg <- c(0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7.5, 10, 12.5, 15, 20, 25, 30, 40, 50, 60)
H <- K[, .(t = tg, hr = exp(main + b * fcurve(shape, tg, tau))), by = .(label, class, term, shape, tau)]
TAB <- H[t %in% c(1, 10, 30, 60), .(n = .N, hr_median = signif(median(hr), 3), share_below_1 = round(mean(hr < 1), 2)), by = .(class, term, shape, t)]
fwrite(TAB, file.path(OUT, "sim06e_hr_table.csv"))
SH <- c(log = "log (current)", linear = "linear", exponential = "exponential (best τ)"); SC <- c(`log (current)` = "black", linear = "#2166AC", `exponential (best τ)` = "#B2182B")
Hm <- H[, .(hr = median(hr)), by = .(class, term, shape, t, tau)]
Hm[, `:=`(shapeL = factor(SH[shape], SH), panel = sprintf("%s\n%s (τ = %g y)", CL[class], TL[term], tau))]
Hm[, panel := factor(panel, unique(Hm[order(match(class, names(CL)), match(term, names(TL))), panel]))]
pc <- ggplot(Hm, aes(t, hr, colour = shapeL)) + geom_hline(yintercept = 1, colour = "grey50") + geom_line(linewidth = 0.9) +
  scale_y_log10() + facet_wrap(~ panel, ncol = 4, scales = "free_y") + scale_colour_manual(values = SC, name = NULL) +
  labs(x = "years since disturbance", y = "hazard ratio, disturbed vs undisturbed (median over transitions, log scale)",
       title = "What each shape says about the disturbance effect over time", subtitle = "1 = no effect; above 1 = the move runs faster after the disturbance") + theme(legend.position = "bottom")
ggsave(file.path(OUT, "sim06e_fig5_curves_by_class.png"), pc, width = 12, height = 9, dpi = 180)
# ---- individual transitions ----
SEL <- data.table(label = c("black_spruce->regen", "fir->regen", "jack_pine->regen", "sugar_maple->regen", "regen->fir", "regen->black_spruce", "fir->white_birch", "white_birch->fir"),
                  term = c("tsd_harv_heavy", "tsd_harv_light", "tsd_fire", "tsd_harv_heavy", "tsd_harv_heavy", "tsd_harv_heavy", "tsd_pest", "tsd_pest"))
Hi <- merge(H, SEL, by = c("label", "term"))
Hi[, `:=`(shapeL = factor(SH[shape], SH), panel = sprintf("%s\n%s (%s; τ = %g y)", gsub("_", " ", gsub("->", " → ", label)), TL[term], CL[class], tau))]
Hi[, panel := factor(panel, unique(Hi[match(paste(SEL$label, SEL$term), paste(label, term)), panel]))]
pi <- ggplot(Hi, aes(t, hr, colour = shapeL)) + geom_hline(yintercept = 1, colour = "grey50") + geom_line(linewidth = 0.9) + scale_y_log10() +
  facet_wrap(~ panel, ncol = 4, scales = "free_y") + scale_colour_manual(values = SC, name = NULL) +
  labs(x = "years since disturbance", y = "hazard ratio, disturbed vs undisturbed (log scale)", title = "Individual transitions: the same data fitted with the three shapes",
       subtitle = "mean coefficients over the imputations (into regen: fits passing the check only)") + theme(legend.position = "bottom")
ggsave(file.path(OUT, "sim06e_fig6_individual.png"), pi, width = 12, height = 7, dpi = 180)
options(width = 200); print(BT[order(class, term)])
print(dcast(TAB[class == "to_regen"], term + t ~ shape, value.var = "hr_median"))
print(dcast(TAB[class == "to_regen"], term + t ~ shape, value.var = "share_below_1"))
print(K[, .N, by = .(class, shape)])
