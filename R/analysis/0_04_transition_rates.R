# ==============================================================================
# 0_04_transition_rates.R -- Chapter 2 analysis, section 0 (Data, no model), row 4: raw transition rates over calendar time.
# Run from the project root:  Rscript R/analysis/0_04_transition_rates.R
#
# Crude (constant-hazard) rate of the route i -> j = observed changes i -> j / plot-years at risk in i, per 1000 plot-years.
# Event = a visit interval whose end state differs from its start state (interval-censored: hidden paths inside an interval are not counted);
# time at risk in i = sum of the lengths of the intervals that START in i (as old I2, which this row extends).
# Choices (Lukas, 2026-11-17): four periods (<= 1990, 1991-2000, 2001-2010, > 2010) assigned by the interval MIDPOINT year (an interval from
#   1985 to 2005 has midpoint 1995 -> 1991-2000); CIs = plot bootstrap (200 draws, whole plots resampled, seed 26; percentile 2.5-97.5%);
#   resolution = 13 states, 4 thermal groups, 4 tolerance groups (other_temperate = "none"); a group route counts only intervals whose group
#   changes, exposure = plot-years in the origin group; also the exit rate per origin (any change) per period and for the whole network;
#   no p-value (plots are the unit; the bootstrap CI is the uncertainty); figures: heatmap per period (states = the modelled routes, thermal
#   groups = all routes), forest plot of the period rate relative to the pooled rate of the route (1 = time-homogeneous; the 20 busiest state
#   routes, all thermal routes), exit rates per period (states, thermal groups). Tolerance groups: tables only.
# Inputs (read-only): data/tables/<RUN_TAG>/{visit_intervals.rds, transition_set.csv}; R/00_groups.R.
# Checks against old output (read-only): figures/I/I2_rates_by_period.csv (old two periods by end-visit year <= 2000 / > 2000: events and plot-years).
# Outputs: data-output/<RUN_TAG>/analysis/0_data/0_04_*
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(here) })
source(here::here("R", "config.R")); source(here::here("R", "00_groups.R"))
OUT <- file.path(RUN_OUT_DIR, "analysis", "0_data"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
options(width = 200); set.seed(26); N_BOOT <- 200L
ORDER <- c("regen", "jack_pine", "white_spruce", "black_spruce", "fir", "white_birch", "trembling_aspen", "other_pioneer",
           "yellow_birch", "red_maple", "sugar_maple", "cedar", "other_temperate")
GCOL <- c(boreal = "#08306B", pioneer = "#E6B800", regen = "#777777", temperate = "#B2182B")
PERL <- c("≤ 1990", "1991–2000", "2001–2010", "> 2010"); PCOL <- setNames(c("#C6DBEF", "#6BAED6", "#2171B5", "#08306B"), PERL)
stopifnot(setequal(names(GRP), ORDER), setequal(names(TOL_GROUP), ORDER))

iv <- as.data.table(readRDS(file.path(TABLE_DIR, "visit_intervals.rds"))); ts <- fread(file.path(TABLE_DIR, "transition_set.csv"))
stopifnot("states outside the 13" = all(c(iv$state_start, iv$state_end) %in% ORDER), all(iv$year_end > iv$year_start))
d <- iv[, .(id_pe, from = state_start, to = state_end, len = year_end - year_start, mid = (year_start + year_end) / 2, chg = is_change, yend = year_end)]
d[, per := as.character(cut(mid, c(-Inf, 1990, 2000, 2010, Inf), labels = PERL))]
stopifnot(!anyNA(d$per), all(d$chg == (d$from != d$to)))
lv <- function(x, lvl) switch(lvl, state = x, thermal = unname(GRP[x]), tolerance = fifelse(is.na(TOL_GROUP[x]), "none", unname(TOL_GROUP[x])))

# ---- 1. rates with plot bootstrap ---------------------------------------------------------------------------------------------
# dd: id_pe, from, to, len, per (character), chg (logical). Returns observed table + bootstrap percentiles.
boot_rates <- function(dd, cats) {
  P <- unique(dd$id_pe); Pn <- length(P); dd[, pidx := match(id_pe, P)]
  EX <- dd[, .(py = sum(len)), by = .(pidx, from, per)]
  EV <- dd[chg == TRUE, .(n = .N), by = .(pidx, from, to, per)]
  EV <- rbindlist(list(EV, EV[, .(n = sum(n), to = "any exit"), by = .(pidx, from, per)]), use.names = TRUE)
  GR <- CJ(from = cats, to = c(cats, "any exit"), per = c(PERL, "all"))[from != to]
  est <- function(w) {
    ex <- EX[, .(py = sum(w[pidx] * py)), by = .(from, per)];  ev <- EV[, .(n = sum(w[pidx] * n)), by = .(from, to, per)]
    ex <- rbindlist(list(ex, ex[, .(py = sum(py), per = "all"), by = from]), use.names = TRUE)
    ev <- rbindlist(list(ev, ev[, .(n = sum(n), per = "all"), by = .(from, to)]), use.names = TRUE)
    g <- merge(merge(GR, ev, by = c("from", "to", "per"), all.x = TRUE), ex, by = c("from", "per"), all.x = TRUE)
    g[is.na(n), n := 0]; g[, rate := fifelse(is.na(py) | py == 0, NA_real_, 1000 * n / py)]
    g[, ratio := rate / rate[per == "all"], by = .(from, to)]; g[!is.finite(ratio), ratio := NA_real_]
    g
  }
  obs <- est(rep(1, Pn))
  bs <- rbindlist(lapply(seq_len(N_BOOT), function(r) est(tabulate(sample.int(Pn, Pn, replace = TRUE), Pn))[, .(from, to, per, rate, ratio)]))
  ci <- bs[, .(lo = quantile(rate, .025, na.rm = TRUE), hi = quantile(rate, .975, na.rm = TRUE),
               ratio_lo = quantile(ratio, .025, na.rm = TRUE), ratio_hi = quantile(ratio, .975, na.rm = TRUE)), by = .(from, to, per)]
  merge(obs, ci, by = c("from", "to", "per"))
}
mk <- function(lvl) { dd <- d[, .(id_pe, from = lv(from, lvl), to = lv(to, lvl), len, per)][, chg := from != to]
                      boot_rates(dd, if (lvl == "state") ORDER else sort(unique(dd$from))) }
cat("bootstrap (", N_BOOT, "draws x 4 levels) ...\n")
RES <- list(state = mk("state"), thermal = mk("thermal"), tolerance = mk("tolerance"))
nw <- boot_rates(d[, .(id_pe, from = "network", to = fifelse(chg, "change", "stay"), len, per, chg)], "network")[to == "any exit"]
RES <- c(RES, list(network = nw))
ST <- RES$state; tk <- ts[kept == TRUE, paste(from, to)]
res <- rbindlist(RES, idcol = "level")
res[, modelled := level == "state" & paste(from, to) %in% tk]
setnames(res, "per", "period"); res[period == "all", period := "all periods"]
setcolorder(res, c("level", "from", "to", "period", "n", "py", "rate", "lo", "hi", "ratio", "ratio_lo", "ratio_hi", "modelled"))
res <- res[order(match(level, c("state", "thermal", "tolerance", "network")), from, to, match(period, c(PERL, "all periods")))]
fwrite(res, file.path(OUT, "0_04_rates.csv"))
psum <- d[, .(intervals = .N, plots = uniqueN(id_pe), plot_years = round(sum(len)), mean_length_y = round(mean(len), 1), changes = sum(chg), share_changed = round(mean(chg), 3)), by = per]
psum <- psum[match(PERL, per)]; fwrite(psum, file.path(OUT, "0_04_period_summary.csv"))

# ---- 2. checks ---------------------------------------------------------------------------------------------------------------
sa <- ST[per == "all" & to != "any exit"]
stopifnot("events over all routes != changes in the interval table" = sum(sa$n) == sum(d$chg),
          "state-route events != transition_set n_events" = all(merge(ts[, .(from, to, n_events)], sa[, .(from, to, n)], by = c("from", "to"))[, n_events == n]),
          "90 modelled routes expected" = length(tk) == 90, "plot-years over origins != total" = abs(sum(ST[per == "all" & to == "any exit", py]) - sum(d$len)) < 1e-6)
old <- fread(file.path(FIG_DIR, "I", "I2_rates_by_period.csv"))
mine <- d[, .(py = sum(len)), by = .(from, p2 = fifelse(yend <= 2000, "early", "late"))]
mn <- d[chg == TRUE, .(n = .N), by = .(from, to, p2 = fifelse(yend <= 2000, "early", "late"))]
o <- melt(old[, .(from, to, n_early, n_late)], id.vars = c("from", "to"), variable.name = "p2", value.name = "n_old")[, p2 := sub("n_", "", p2)]
o <- merge(o, mn, by = c("from", "to", "p2"), all.x = TRUE); o[is.na(n), n := 0L]
op <- melt(unique(old[, .(from, py_early, py_late)]), id.vars = "from", variable.name = "p2", value.name = "py_old")[, p2 := sub("py_", "", p2)]; op <- merge(op, mine, by = c("from", "p2"))
stopifnot("old I2 events differ" = all(o$n == o$n_old), "old I2 plot-years differ" = all(abs(op$py - op$py_old) < 1e-6))
cat("check vs old I2: events (", nrow(o), "cells) and plot-years (", nrow(op), "cells) identical under the old two-period rule\n")

# ---- 3. figures --------------------------------------------------------------------------------------------------------------
theme_set(theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95")))
lab <- function(x) gsub("_", " ", x)
save2 <- function(p, name, w, h) { ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300); ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf) }
hm <- function(lvl, cats, only_modelled, w, h, tsz, name) {
  g <- res[level == lvl & period != "all periods" & to != "any exit"]
  if (only_modelled) g <- g[modelled == TRUE]
  g <- merge(CJ(from = cats, to = cats, period = PERL)[from != to], g[, .(from, to, period, n, rate)], by = c("from", "to", "period"), all.x = TRUE)
  g[, `:=`(period = factor(period, PERL), from = factor(from, rev(cats)), to = factor(to, cats),
           txt = fifelse(is.na(n), "", fifelse(n == 0, "0", formatC(signif(rate, 2), format = "fg", digits = 2))), fillv = fifelse(is.na(n) | n == 0, NA_real_, rate))]
  p <- ggplot(g, aes(to, from, fill = fillv)) + geom_tile(colour = "white") + geom_text(aes(label = txt), size = tsz) + facet_wrap(~ period, nrow = 2) +
    scale_fill_distiller(palette = "Blues", direction = 1, trans = "log10", na.value = "grey96", name = "rate per 1000\nplot-years") +
    scale_x_discrete(drop = FALSE, labels = lab, position = "top") + scale_y_discrete(drop = FALSE, labels = lab) + labs(x = "to", y = "from") +
    theme(axis.text.x = element_text(angle = 40, hjust = 0), panel.grid = element_blank())
  save2(p, name, w, h)
}
hm("state", ORDER, TRUE, 15, 13, 2.2, "0_04_rates_state_heatmap")
hm("thermal", c("boreal", "pioneer", "regen", "temperate"), FALSE, 9, 7.5, 3, "0_04_rates_thermal_heatmap")
fo <- function(lvl, top, name, w, h) {
  g <- res[level == lvl & to != "any exit" & period != "all periods"]
  nall <- res[level == lvl & to != "any exit" & period == "all periods", .(from, to, n_all = n)]
  tp <- nall[order(-n_all)][seq_len(min(top, .N))][, lbl := sprintf("%s → %s (%d)", lab(from), lab(to), as.integer(n_all))]
  g <- merge(g, tp[, .(from, to, lbl)], by = c("from", "to")); g[, `:=`(lbl = factor(lbl, rev(tp$lbl)), period = factor(period, PERL))]
  p <- ggplot(g, aes(ratio, lbl, colour = period)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_errorbar(aes(xmin = ratio_lo, xmax = ratio_hi), width = 0, orientation = "y", position = position_dodge(width = 0.7)) +
    geom_point(position = position_dodge(width = 0.7), size = 1.5) + scale_colour_manual(values = PCOL, name = "period") +
    scale_x_log10() + labs(x = "rate relative to the route's rate over the whole record (log scale)", y = NULL)
  save2(p, name, w, h)
}
fo("state", 20, "0_04_rates_state_forest", 8.5, 9); fo("thermal", 12, "0_04_rates_thermal_forest", 8.5, 5.5)
ex <- function(lvl, name, ncol, w, h, cols = NULL) {
  g <- res[level == lvl & to == "any exit"]
  pooled <- g[period == "all periods", .(from, pooled = rate)]; g <- merge(g[period != "all periods"], pooled, by = "from")
  g[, `:=`(period = factor(period, PERL), from = factor(from, if (lvl == "state") ORDER else sort(unique(from))))]
  p <- ggplot(g, aes(period, rate)) + geom_hline(aes(yintercept = pooled), linetype = 2, colour = "grey50") + geom_errorbar(aes(ymin = lo, ymax = hi), width = 0) +
    geom_point(aes(colour = from), size = 1.8) + facet_wrap(~ from, ncol = ncol, scales = "free_y", labeller = labeller(from = lab)) + expand_limits(y = 0) +
    labs(x = NULL, y = "exit rate per 1000 plot-years") + theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "none")
  p <- p + scale_colour_manual(values = if (!is.null(cols)) cols else rep("black", 13))
  save2(p, name, w, h)
}
ex("state", "0_04_exit_rates_state", 5, 10, 7.5); ex("thermal", "0_04_exit_rates_thermal", 4, 9, 3.2, GCOL)

# ---- 4. console summary ------------------------------------------------------------------------------------------------------
rd <- function(x) if (is.numeric(x)) round(x, 1) else x; PC <- c("from", PERL, "all periods")
cat("\nPeriods (by interval midpoint):\n"); print(psum, row.names = FALSE)
cat("\nNetwork exit rate per 1000 plot-years (95% plot-bootstrap CI):\n")
print(res[level == "network", .(period, n, plot_years = round(py), rate = round(rate, 1), lo = round(lo, 1), hi = round(hi, 1))], row.names = FALSE)
cat("\nExit rate per origin thermal group:\n"); print(dcast(res[level == "thermal" & to == "any exit"], from ~ period, value.var = "rate")[, PC, with = FALSE][, lapply(.SD, rd)], row.names = FALSE)
cat("\nExit rate per origin state:\n"); print(dcast(res[level == "state" & to == "any exit"], from ~ period, value.var = "rate")[match(ORDER, from), PC, with = FALSE][, lapply(.SD, rd)], row.names = FALSE)
m <- res[modelled == TRUE & period != "all periods"]
cat("\nModelled state routes (90) whose period/pooled ratio CI excludes 1 (descriptive count; no multiplicity adjustment):\n")
print(m[, .(routes = .N, excl1 = sum(ratio_lo > 1 | ratio_hi < 1, na.rm = TRUE), higher = sum(ratio_lo > 1, na.rm = TRUE), lower = sum(ratio_hi < 1, na.rm = TRUE)), by = period][match(PERL, period)], row.names = FALSE)
cat("\nSaved in", OUT, "\n")
