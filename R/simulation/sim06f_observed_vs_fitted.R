# ==============================================================================
# R/simulation/sim06f_observed_vs_fitted.R   (no INLA; run from the project root)
# Observed vs fitted: for selected transitions x disturbance types, the CRUDE observed rate ratio by years since disturbance
#   (events / exposure in disturbed rows of that type, per TSD bin, divided by events / exposure in that transition's undisturbed rows;
#   imputation 1; Poisson 95% CI), with the three fitted shapes on top (hazard ratio at mean GDD; coefficients as in sim06e).
# Crude = not adjusted for climate, soil, other disturbances or time in state; it shows the data the curves have to follow, not the same quantity.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(here) })
source(here("R", "config.R")); SIM_M <- 200L; M_MI <- 30L; SD_MIN <- 1e-3; N_COV_SAMPLES <- 0L
source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
S2 <- here("data-output", "simulation", "sim06d_regen_screen")
SEL <- data.table(label = c("black_spruce->regen", "fir->regen", "jack_pine->regen", "sugar_maple->regen", "regen->fir", "regen->black_spruce", "fir->white_birch", "white_birch->fir"),
                  term = c("tsd_harv_heavy", "tsd_harv_light", "tsd_fire", "tsd_harv_heavy", "tsd_harv_heavy", "tsd_harv_heavy", "tsd_pest", "tsd_pest"))
FLAG <- c(tsd_fire = "fire_f", tsd_pest = "pest_f", tsd_harv_light = "harv_light", tsd_harv_heavy = "harv_heavy")
TL <- c(tsd_fire = "fire", tsd_pest = "pest", tsd_harv_light = "light harvest", tsd_harv_heavy = "heavy harvest")
inp <- prepare_fit_inputs(TABLE_DIR)
ex <- readRDS(file.path(FIT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, 1L))); p1 <- as.data.table(ex$p1); p2 <- as.data.table(ex$p2); rm(ex)
fe <- setNames(p1$entry, as.character(p1$id_pe)); key_iv <- paste(inp$iv$id_pe, inp$iv$interval)
ex_row <- rep(NA_real_, nrow(inp$iv)); ex_row[match(paste(p2$id_pe, p2$interval), key_iv)] <- p2$exit_abs; b_abs <- inp$iv$t_end
p2_fun <- function(a_dur, b_dur, from, to, r) { E <- b_abs[r] - b_dur; min(max(ex_row[r] - E, a_dur), b_dur) }
rows <- build_rows(inp$iv, fe, inp$states, inp$tr, p2_fun = p2_fun, keep_cols = inp$keep_cols)$rows[kept == TRUE]; rows[, label := paste(from, to, sep = "->")]
rows <- rows[label %in% SEL$label]
BR <- c(0, 2, 5, 10, 15, 20, 30, 45, 70, Inf)
OBS <- rbindlist(lapply(seq_len(nrow(SEL)), function(i) { lb <- SEL$label[i]; T <- SEL$term[i]; d <- rows[label == lb]; d[, expo := Tstop - truncation]
  und <- d[fire_f == 0 & pest_f == 0 & harv_light == 0 & harv_heavy == 0]; r0 <- sum(und$status) / sum(und$expo)
  dd <- d[get(FLAG[[T]]) == 1]; dd[, yrs := expm1(get(T))]; dd[, bin := cut(yrs, BR, right = FALSE)]
  o <- dd[, .(events = sum(status), expo = sum(expo), t_mid = median(yrs)), by = bin]
  o[, `:=`(rr = (events / expo) / r0, lo = (qchisq(0.025, 2 * events) / 2 / expo) / r0, hi = (qchisq(0.975, 2 * events + 2) / 2 / expo) / r0)]
  o[, `:=`(label = lb, term = T, undist_events = sum(und$status))][] }))
# fitted curves from sim06e coefficients (recomputed as there)
ar <- fread(here("data-output", "final", "results", "alpha_report_final.csv"))[, .(label, class)]
BT <- fread(file.path(S2, "sim06e_best_tau_by_class.csv")); TAU <- c(2, 3.5, 6, 10, 15, 25, 40, 60, 100, 160)
get1 <- function(f, screened) { r <- readRDS(f); cl <- ar[label == r$label, class]
  rbindlist(lapply(SEL[label == r$label, term], function(T) { tb <- BT[term == T & class == cl]
    fm <- c(log = "log_all", linear = paste0("lin_", T), exponential = sprintf("exp_%s_%03d", T, match(tb$tau_best, TAU)))
    rbindlist(lapply(names(fm), function(s) { x <- r$fits[[fm[[s]]]]; ok <- !is.null(x$coef) && is.na(x$error) && (!screened || isTRUE(x$screen$valid))
      data.table(label = r$label, term = T, shape = s, tau = tb$tau_best, ok = ok, main = if (ok) x$coef[parameter == FLAG[[T]], mean] else NA_real_, b = if (ok) x$coef[parameter == T, mean] else NA_real_) })) })) }
fl1 <- unlist(lapply(SEL$label[5:8], function(lb) list.files(here("data-output", "simulation", "sim06_tsd_tau", "fits"), sprintf("^p1_m.._%s\\.rds$", gsub("->", "_to_", lb)), full.names = TRUE)))
fl2 <- unlist(lapply(SEL$label[1:4], function(lb) list.files(file.path(S2, "fits"), sprintf("^m.._%s\\.rds$", gsub("->", "_to_", lb)), full.names = TRUE)))
K <- rbind(rbindlist(lapply(fl1, get1, FALSE)), rbindlist(lapply(fl2, get1, TRUE)))[ok == TRUE, .(main = mean(main), b = mean(b)), by = .(label, term, shape, tau)]
tg <- seq(0, 70, 0.5); fc <- function(s, t, tau) switch(s, log = log1p(t), linear = t / 10, exponential = exp(-t / tau))
CUR <- K[, .(t = tg, hr = exp(main + b * fc(shape, tg, tau))), by = .(label, term, shape, tau)]
SH <- c(log = "log (current)", linear = "linear", exponential = "exponential (best τ)"); SC <- c(`log (current)` = "black", linear = "#2166AC", `exponential (best τ)` = "#B2182B")
pan <- function(lb, T) sprintf("%s, %s", gsub("_", " ", gsub("->", " → ", lb)), TL[T])
CUR[, `:=`(shapeL = factor(SH[shape], SH), panel = pan(label, term))]; OBS[, panel := pan(label, term)]
lv <- pan(SEL$label, SEL$term); CUR[, panel := factor(panel, lv)]; OBS[, panel := factor(panel, lv)]
OBS[, lab := sprintf("%d", events)]
p <- ggplot() + geom_hline(yintercept = 1, colour = "grey50") +
  geom_errorbar(data = OBS[events > 0], aes(x = t_mid, ymin = pmax(lo, 1e-3), ymax = hi), width = 0, colour = "grey40") +
  geom_point(data = OBS[events > 0], aes(t_mid, rr, size = events), colour = "grey25") +
  geom_point(data = OBS[events == 0], aes(t_mid, hi), shape = 25, colour = "grey40", fill = "grey80", size = 2) +
  geom_line(data = CUR, aes(t, hr, colour = shapeL), linewidth = 0.9) +
  scale_y_log10() + coord_cartesian(xlim = c(0, 70)) + facet_wrap(~ panel, ncol = 4, scales = "free_y") +
  scale_colour_manual(values = SC, name = NULL) + scale_size_area(max_size = 4, name = "observed changes in the bin") +
  labs(x = "years since disturbance", y = "rate ratio, disturbed vs undisturbed stands (log scale)",
       title = "What the data show (grey points) and what each shape makes of it (lines)",
       subtitle = "points = observed rate in each band of years since disturbance relative to undisturbed stands (95% CI; triangle = no changes, upper bound);\nnot adjusted for climate or time in the state, so points and lines measure similar but not identical things") +
  theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), legend.position = "bottom", strip.background = element_rect(fill = "grey95"))
ggsave(file.path(S2, "sim06f_observed_vs_fitted.png"), p, width = 13, height = 7.5, dpi = 180)
fwrite(OBS, file.path(S2, "sim06f_observed_bins.csv")); print(OBS[, .(panel, bin, events, rr = signif(rr, 2))], nrows = 200)
