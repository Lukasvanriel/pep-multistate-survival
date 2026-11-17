# sim06c_explain.R  (no INLA): (A) one transition, every form drawn, with its fit; (B) per type x class: log vs linear vs best exponential.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(here); library(patchwork) })
D <- here("data-output", "simulation", "sim06_tsd_tau"); TAU <- c(2, 3.5, 6, 10, 15, 25, 40, 60, 100, 160)
theme_set(theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank()))
LB <- "regen->fir"; T <- "tsd_harv_heavy"; MN <- "harv_heavy"
fl <- list.files(file.path(D, "fits"), sprintf("^p1_m.._%s\\.rds$", gsub("->", "_to_", LB)), full.names = TRUE)
A <- rbindlist(lapply(fl, function(f) { r <- readRDS(f); b <- r$fits$log_all$waic
  fm <- c("log_all", paste0("lin_", T), sprintf("exp_%s_%03d", T, seq_along(TAU)))
  rbindlist(lapply(fm, function(n) { x <- r$fits[[n]]; data.table(form = n, m = r$m, main = x$coef[parameter == MN, mean], b = x$coef[parameter == T, mean], dwaic = x$waic - b) })) }))
A <- A[, .(main = mean(main), b = mean(b), dwaic = mean(dwaic)), by = form]
A[, `:=`(shape = fifelse(form == "log_all", "log", fifelse(grepl("^lin_", form), "linear", "exponential")), tau = fifelse(grepl("^exp_", form), TAU[as.integer(sub(".*_", "", form))], NA_real_))]
tg <- seq(0, 60, 0.5)
CUR <- rbindlist(lapply(seq_len(nrow(A)), function(i) { a <- A[i]; x <- switch(a$shape, log = log1p(tg), linear = tg / 10, exponential = exp(-tg / a$tau))
  data.table(form = a$form, shape = a$shape, tau = a$tau, t = tg, hr = exp(a$main + a$b * x)) }))
A[, lab := fifelse(shape == "exponential", sprintf("exponential, τ = %g y", tau), shape)]
A[, lab := sprintf("%s   (WAIC %+.0f vs log)", lab, dwaic)]
CUR <- merge(CUR, A[, .(form, lab)], by = "form")
show <- A[shape != "exponential" | tau %in% c(2, 6, 10, 40, 160), form]
lv <- A[match(show, form), lab][order(match(show, c("log_all", paste0("lin_", T), sprintf("exp_%s_%03d", T, match(c(2, 6, 10, 40, 160), TAU)))))]
cols <- setNames(c("black", "#2166AC", "#FDAE61", "#F46D43", "#B2182B", "#8C6BB1", "#7A7A7A"), lv)
pA <- ggplot(CUR[form %in% show], aes(t, hr, colour = factor(lab, lv))) + geom_hline(yintercept = 1, colour = "grey60") + geom_line(linewidth = 0.9) +
  scale_colour_manual(values = cols, name = NULL) + scale_y_log10() +
  labs(x = "years since heavy harvest", y = "how much faster regen → fir runs (hazard ratio, log scale)",
       title = "A  One transition (regen → fir after heavy harvest): the same data, fitted with each shape",
       subtitle = "each line = the best curve of that shape; WAIC < 0 = fits better than log (best here: τ = 3.5–6 y)") +
  theme(legend.position = c(0.68, 0.25), legend.background = element_rect(fill = alpha("white", 0.8)))
# B: per type x class
W <- fread(file.path(D, "sim06_waic.csv")); TL <- c(tsd_fire = "fire", tsd_pest = "pest", tsd_harv_light = "light harvest", tsd_harv_heavy = "heavy harvest")
CL <- c(to_regen = "into regen", regen_exit = "out of regen", between = "between states")
S <- W[grepl("^(exp|lin)_", form), .(d = sum(dwaic, na.rm = TRUE) / uniqueN(m)), by = .(term, class, form, tau)]
B <- rbind(S[grepl("^lin_", form), .(term, class, shape = "linear", d)], S[grepl("^exp_", form), .SD[which.min(d)], by = .(term, class)][, .(term, class, shape = sprintf("best exponential (τ %g)", tau), d)], fill = TRUE)
B[, best_tau := sub(".*τ ([0-9.]+).*", "\\1", shape)]; B[, kind := fifelse(shape == "linear", "linear", "best exponential")]
B[, `:=`(term = factor(TL[term], TL), class = factor(CL[class], CL))]
B[, txt := fifelse(kind == "best exponential", paste0("τ ", best_tau), "")]
B[, unreliable := class == "into regen"]
pB <- ggplot(B, aes(d, term, colour = kind)) + geom_vline(xintercept = 0, colour = "black") + geom_point(size = 2.5) +
  geom_text(aes(label = txt), nudge_y = 0.3, size = 2.8, show.legend = FALSE) + facet_wrap(~ class, scales = "free_x") +
  scale_colour_manual(values = c(linear = "#2166AC", `best exponential` = "#B2182B"), name = NULL) +
  labs(x = "WAIC minus WAIC of log (summed over transitions)   ← better than log  |  worse than log →", y = NULL,
       title = "B  All transitions: log (the line at 0) vs linear vs the best exponential, by disturbance type and class",
       subtitle = "into regen: not reliable (failed fits in 4 of 5 imputations); out of regen and between states: reliable") + theme(legend.position = "bottom")
ggsave(file.path(D, "sim06_fig3_explained.png"), pA / pB + plot_layout(heights = c(1.1, 1)), width = 11, height = 11, dpi = 170)
print(A[, .(lab)])
