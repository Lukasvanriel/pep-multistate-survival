# ==============================================================================
# 25_I2_rates_by_period.R   (no INLA; seconds; run from the project root)
# I2 (user 2026-11-11): model-free test of non-stationarity. Intervals split by the end visit: <= 2000 vs > 2000. Crude rate i -> j = observed
#   changes i -> j / plot-years at risk in i (sum of interval lengths of intervals starting in i), per 1000 plot-years; ignores multiple moves within
#   an interval (as A1a). Rate ratio late / early with a Poisson (Wald, log scale) 95% CI. Routes: all observed; forest plot for the 20 busiest.
# Caveat carried: interval length and disturbance mix differ between periods, and different plots are measured in each (A1e, H).
# ==============================================================================
source(here::here("R", "25_I_common.R"))
d <- vi[, .(id_pe, from = state_start, to = state_end, chg = is_change, len = year_end - year_start,
            per = fifelse(year_end <= 2000, "early", "late"), dist = (fire > 0) | (harv > 0) | (pest > 0))]
cat("intervals by period:\n"); print(d[, .(n = .N, plots = uniqueN(id_pe), mean_len = round(mean(len), 1), share_disturbed = round(mean(dist), 2), share_changed = round(mean(chg), 3)), by = per])
expo <- d[, .(py = sum(len)), by = .(from, per)]
cnt  <- d[chg == TRUE, .(n = .N), by = .(from, to, per)]
g <- CJ(from = ORDER, to = ORDER, per = c("early", "late"))[from != to]
g <- merge(merge(g, cnt, by = c("from", "to", "per"), all.x = TRUE), expo, by = c("from", "per"), all.x = TRUE); g[is.na(n), n := 0L]
g[, rate := 1000 * n / py]
w <- dcast(g, from + to ~ per, value.var = c("n", "py", "rate"))
w[, `:=`(n_tot = n_early + n_late, rr = (n_late / py_late) / (n_early / py_early))]
w[, se := sqrt(1 / n_early + 1 / n_late)]; w[, `:=`(rr_lo = exp(log(rr) - 1.96 * se), rr_hi = exp(log(rr) + 1.96 * se))]
w[n_early == 0 | n_late == 0, c("rr", "rr_lo", "rr_hi") := NA_real_]
fwrite(w, file.path(OUT, "I2_rates_by_period.csv"))
# any exit per origin
o <- d[, .(n = sum(chg), py = sum(len)), by = .(from, per)]; ow <- dcast(o, from ~ per, value.var = c("n", "py"))
ow[, `:=`(rate_early = 1000 * n_early / py_early, rate_late = 1000 * n_late / py_late, rr = (n_late / py_late) / (n_early / py_early), se = sqrt(1 / n_early + 1 / n_late))]
ow[, `:=`(rr_lo = exp(log(rr) - 1.96 * se), rr_hi = exp(log(rr) + 1.96 * se))]; fwrite(ow, file.path(OUT, "I2_exit_rates_by_period.csv"))
tot <- d[, .(n = sum(chg), py = sum(len)), by = per]; cat(sprintf("all exits: early %.1f, late %.1f per 1000 plot-years; ratio %.2f\n",
  1000 * tot[per == "early", n / py], 1000 * tot[per == "late", n / py], tot[per == "late", n / py] / tot[per == "early", n / py]))
# figures
PER <- c(early = "end visit ≤ 2000", late = "end visit > 2000")
lim <- max(g$rate[g$n > 0]); gg <- copy(g)[, `:=`(status = fifelse(paste0(from, "->", to) %in% modelled | n > 0, "value", "notmod"),
  fillv = rate, txt = fifelse(n > 0, formatC(signif(rate, 2), format = "fg", digits = 2), ""), bold = FALSE)]
pa <- lapply(names(PER), function(k) heat(gg[per == k], "rate per 1000\nplot-years", diverging = FALSE, limits = c(0, lim), subtitle = PER[[k]]))
save_fig(pa[[1]] + pa[[2]] + plot_layout(guides = "collect") + plot_annotation(title = "I2a  Observed transition rates by period"), "I2a_rates_by_period", 15, 7.5, NULL)
hb <- copy(w)[, `:=`(status = fifelse(!is.na(rr) & n_tot >= 15, "value", "notmod"), fillv = log(rr), txt = fmt_hr(rr), bold = !is.na(rr_lo) & (rr_lo > 1 | rr_hi < 1))]
lr <- cap_sym(hb[status == "value", fillv])
save_fig(heat(hb, "rate ratio", limits = c(-lr, lr), trans_lab = function(v) fmt_hr(exp(v)), breaks = hr_breaks(lr),
              title = "I2b  Rate ratio, end visit > 2000 vs ≤ 2000", subtitle = "routes with ≥ 15 observed changes and events in both periods; bold = 95% CI excludes 1"),
         "I2b_rate_ratio", 7.5, 7, NULL)
top <- w[order(-n_tot)][1:20][, lbl := sprintf("%s → %s (%d)", lab(from), lab(to), n_tot)]
top[, lbl := factor(lbl, levels = rev(lbl))]
pc <- ggplot(top, aes(rr, lbl)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = rr_lo, xmax = rr_hi), width = 0, orientation = "y") +
  geom_point(aes(shape = rr_lo > 1 | rr_hi < 1), size = 1.8, fill = "white") + scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), guide = "none") +
  scale_x_log10() + labs(x = "rate ratio, end visit > 2000 vs ≤ 2000 (log scale)", y = NULL, title = "I2c  The 20 busiest routes",
                         subtitle = "filled = 95% CI excludes 1; in brackets: observed changes")
save_fig(pc, "I2c_rate_ratio_busiest", 7, 6, top)
cat("I2 done\n")
