# ==============================================================================
# 0_05_climate_trend.R -- Chapter 2 analysis, section 0 (Data, no model), row 5: climate trend in the plots over the record.
# Run from the project root:  Rscript R/analysis/0_05_climate_trend.R
#
# Series = the YEARLY raw climate of each retained plot, 1960-2022 (NRCan; data/build/climate_by_year_raw.rds): GDD = gdd_3 (growing-season GDD);
#   CMI = cmi_gs, rebuilt per plot-year with the growing-season weighting rule of R/03c_build_climate_covariates.R step 4b (function copied unchanged).
#   Soil PC1 is constant within every plot, so it has no trend and is not analysed.
# Slope per plot = (value in the last year - value in the first year) / (last year - first year) x 10 = change per decade, between the plot's first
#   and last VISIT year (the years the plot was measured; a last visit after 2022 is clipped to 2022 because the series ends there).
#   Main estimate = the mean of the plot slopes over plots whose span is >= 10 years (plot bootstrap, 200 draws, seed 26, percentile 2.5-97.5%).
#   Check next to it: the least-squares slope over the same years of the plot, and the mean over all plots with a span > 0.
# Units: raw (GDD or CMI per decade) and z (SD per decade; model scaling data-output/<RUN_TAG>/fits/scaling.rds).
# Splits (Lukas, 2026-11-17): all plots, period (of the MIDPOINT of the plot's span; cuts as rows 3-4), latitude tertile of the plots (as rows 2-3).
# Also the yearly network mean 1960-2022 (plot bootstrap per year) with its linear trend per decade.
# Inputs (read-only): data/tables/<RUN_TAG>/visit_intervals.rds; data/build/climate_by_year_raw.rds; <FIT_DIR>/scaling.rds; PEP.gpkg (latitude).
# Checks: the rebuilt yearly series averaged over each interval's window reproduces the model covariates gdd_3 and cmi_gs (visit_intervals) to 1e-6.
#   Informational: the old statement in R/22_hindcast_space_for_time.R (median +0.49 SD GDD, first to last visit).
# Outputs: data-output/<RUN_TAG>/analysis/0_data/0_05_*
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(here) })
source(here::here("R", "config.R"))
OUT <- file.path(RUN_OUT_DIR, "analysis", "0_data"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
options(width = 200); set.seed(26); N_BOOT <- 200L; YMIN <- 1960L; YMAX <- 2022L; MIN_SPAN <- 10L
PERL <- c("≤ 1990", "1991–2000", "2001–2010", "> 2010")

iv  <- as.data.table(readRDS(file.path(TABLE_DIR, "visit_intervals.rds")))
scz <- readRDS(file.path(FIT_DIR, "scaling.rds")); SD <- c(gdd = scz[var == "gdd_3", sd], cmi = scz[var == "cmi_gs", sd])
con <- DBI::dbConnect(RSQLite::SQLite(), here("data", "raw", "pep", "PEP.gpkg"), flags = RSQLite::SQLITE_RO)
pl  <- unique(as.data.table(DBI::dbGetQuery(con, "SELECT id_pe, latitude FROM placette")), by = "id_pe"); DBI::dbDisconnect(con)
ids <- unique(iv$id_pe); stopifnot("plot without latitude" = all(ids %in% pl[!is.na(latitude), id_pe]))
q <- quantile(pl[id_pe %in% ids, latitude], c(1/3, 2/3))
BAND_LAB <- c(sprintf("south (≤ %.1f°N)", q[1]), sprintf("middle (%.1f–%.1f°N)", q[1], q[2]), sprintf("north (> %.1f°N)", q[2]))

# ---- 1. yearly series ---------------------------------------------------------------------------------------------------------
# copied unchanged from R/03c_build_climate_covariates.R (step 4b), 2026-11-17
cmi_gs_matched <- function(year, start_gs, end_gs, cmi_mat) {
  stopifnot(ncol(cmi_mat) == 12L, length(year) == nrow(cmi_mat), length(start_gs) == length(year), length(end_gs) == length(year))
  leap  <- (year %% 4L == 0L & year %% 100L != 0L) | year %% 400L == 0L
  md    <- rbind(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31),
                 c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))
  mdays  <- md[leap + 1L, , drop = FALSE]
  mend   <- t(apply(md, 1, cumsum))[leap + 1L, , drop = FALSE]
  mstart <- mend - mdays + 1
  inside <- pmax(pmin(mend, end_gs) - pmax(mstart, start_gs) + 1, 0)
  w <- inside / mdays
  cm <- as.matrix(cmi_mat); cm[!is.na(w) & w == 0] <- 0
  list(cmi_gs = rowSums(cm * w), w = w)
}
cr <- as.data.table(readRDS(file.path(BUILD_DIR, "climate_by_year_raw.rds")))[id_pe %in% ids]
cr[, cmi_gs := cmi_gs_matched(year, start_gs, end_gs, cr[, paste0("cmi_", 1:12), with = FALSE])$cmi_gs]
cr <- cr[, .(id_pe, year, gdd = gdd_3, cmi = cmi_gs)]; setkey(cr, id_pe, year)
stopifnot(all(cr$year >= YMIN & cr$year <= YMAX), uniqueN(cr[, .(id_pe, year)]) == nrow(cr), all(ids %in% cr$id_pe))
cat("retained plots:", length(ids), "| plot-years:", nrow(cr), "| NA gdd:", sum(is.na(cr$gdd)), "| NA cmi:", sum(is.na(cr$cmi)), "\n")

# check: windows of the model covariates
W0 <- iv[, .(year = seq.int(win_start, win_end)), by = .(id_pe, interval)]
W0 <- merge(W0, cr, by = c("id_pe", "year"), all.x = TRUE)
A <- merge(W0[, .(gdd_w = mean(gdd), cmi_w = mean(cmi), ny = .N), by = .(id_pe, interval)], iv[, .(id_pe, interval, gdd_3, cmi_gs, n_years_used)], by = c("id_pe", "interval"))
chk <- data.table(variable = c("gdd_3", "cmi_gs"), intervals = nrow(A), max_abs_diff = c(max(abs(A$gdd_w - A$gdd_3)), max(abs(A$cmi_w - A$cmi_gs))))
stopifnot("window length != n_years_used" = all(A$ny == A$n_years_used), "yearly series does not reproduce the model covariates" = all(chk$max_abs_diff < 1e-6))
fwrite(chk, file.path(OUT, "0_05_series_check.csv")); cat("check: yearly series averaged over each interval's window reproduces gdd_3 and cmi_gs (max abs diff",
                                                           format(max(chk$max_abs_diff), digits = 2), ")\n")

# ---- 2. slope per plot --------------------------------------------------------------------------------------------------------
P <- iv[, .(first_year = min(year_start), last_year = max(year_end), n_intervals = .N), by = id_pe]
P[, `:=`(fy = pmax(first_year, YMIN), ly = pmin(last_year, YMAX))]; P[, `:=`(span = ly - fy, clipped = last_year > YMAX, mid = (fy + ly) / 2)]
P[, `:=`(per = as.character(cut(mid, c(-Inf, 1990, 2000, 2010, Inf), labels = PERL)), latitude = pl$latitude[match(id_pe, pl$id_pe)])]
P[, band := as.character(cut(latitude, c(-Inf, q, Inf), labels = BAND_LAB))]; stopifnot(!anyNA(P$per), !anyNA(P$band))
L <- P[span > 0, .(year = seq.int(fy, ly)), by = id_pe]; L <- merge(L, cr, by = c("id_pe", "year"), all.x = TRUE)
for (v in c("gdd", "cmi")) {
  f <- cr[P[, .(id_pe, year = fy)], on = c("id_pe", "year"), get(v)]; l <- cr[P[, .(id_pe, year = ly)], on = c("id_pe", "year"), get(v)]
  P[, paste0("x_first_", v) := f]; P[, paste0("x_last_", v) := l]
  P[, paste0("sl_two_point_", v) := fifelse(span > 0, 10 * (l - f) / span, NA_real_)]
  o <- L[!is.na(get(v)), .(n = .N, b = 10 * cov(year, get(v)) / var(year)), by = id_pe][n >= 3]
  P[, paste0("sl_ols_", v) := o$b[match(id_pe, o$id_pe)]]
}
cat("plots:", nrow(P), "| span >= ", MIN_SPAN, " y: ", P[span >= MIN_SPAN, .N], " (", round(100 * P[, mean(span >= MIN_SPAN)], 1), "%) | span > 0: ", P[span > 0, .N],
    " | last visit after ", YMAX, " (clipped): ", P[, sum(clipped)], " | NA endpoint (span >= 10): ", P[span >= MIN_SPAN, sum(is.na(sl_two_point_gdd) | is.na(sl_two_point_cmi))], "\n", sep = "")
fwrite(P[, !"mid"], file.path(OUT, "0_05_plot_slopes.csv"))

# ---- 3. mean slope with plot bootstrap ----------------------------------------------------------------------------------------
np <- nrow(P); W <- vapply(seq_len(N_BOOT), function(r) tabulate(sample.int(np, np, replace = TRUE), np), numeric(np))
SPL <- c(list(list("total", "all plots", rep(TRUE, np))), lapply(PERL, function(l) list("period", l, P$per %in% l)), lapply(BAND_LAB, function(l) list("latitude", l, P$band %in% l)))
FILT <- list("span >= 10 y" = P$span >= MIN_SPAN, "all plots with a span" = P$span > 0)
rows <- list()
for (v in c("gdd", "cmi")) for (e in c("two_point", "ols")) { s <- P[[paste0("sl_", e, "_", v)]]
  for (fn in names(FILT)) for (sp in SPL) {
    idx <- which(FILT[[fn]] & sp[[3]] & !is.na(s)); if (!length(idx)) next
    Wi <- W[idx, , drop = FALSE]; bm <- colSums(Wi * s[idx]) / colSums(Wi); qs <- quantile(s[idx], c(.05, .5, .95))
    rows[[length(rows) + 1L]] <- data.table(var = v, estimator = e, filter = fn, split_type = sp[[1]], split = sp[[2]], n_plots = length(idx), mean = mean(s[idx]),
                                            lo = quantile(bm, .025, na.rm = TRUE), hi = quantile(bm, .975, na.rm = TRUE), median = qs[[2]], p05 = qs[[1]], p95 = qs[[3]])
  } }
SUM <- rbindlist(rows); cols <- c("mean", "lo", "hi", "median", "p05", "p95")
SUM[, paste0(cols, "_z") := lapply(.SD, function(x) x / SD[var]), .SDcols = cols]
setnames(SUM, cols, paste0(cols, "_raw")); fwrite(SUM, file.path(OUT, "0_05_slopes_summary.csv"))

# ---- 4. yearly network mean and trend -----------------------------------------------------------------------------------------
yrs <- YMIN:YMAX; ny <- length(yrs); pidx <- match(cr$id_pe, ids); yi <- cr$year - YMIN + 1L
tr <- function(m) { m <- rbind(m); drop((m - rowMeans(m)) %*% (yrs - mean(yrs)) / sum((yrs - mean(yrs))^2) * 10) }
res4 <- lapply(c("gdd", "cmi"), function(v) {
  X <- matrix(0, length(ids), ny); Mk <- matrix(0, length(ids), ny); ok <- !is.na(cr[[v]]); X[cbind(pidx[ok], yi[ok])] <- cr[[v]][ok]; Mk[cbind(pidx[ok], yi[ok])] <- 1
  m0 <- colSums(X) / colSums(Mk); mb <- crossprod(W, X) / crossprod(W, Mk); tb <- tr(mb)
  list(series = data.table(var = v, year = yrs, n_plots = colSums(Mk), mean = m0, lo = apply(mb, 2, quantile, .025, na.rm = TRUE), hi = apply(mb, 2, quantile, .975, na.rm = TRUE)),
       trend = data.table(var = v, years = paste(YMIN, YMAX, sep = "-"), per_decade_raw = tr(m0), lo_raw = quantile(tb, .025), hi_raw = quantile(tb, .975)))
})
NS <- rbindlist(lapply(res4, `[[`, "series")); NT <- rbindlist(lapply(res4, `[[`, "trend"))
NT[, c("per_decade_z", "lo_z", "hi_z") := lapply(.SD, function(x) x / SD[var]), .SDcols = c("per_decade_raw", "lo_raw", "hi_raw")]
fwrite(NS, file.path(OUT, "0_05_network_series.csv")); fwrite(NT, file.path(OUT, "0_05_network_trend.csv"))

# ---- 5. figures ---------------------------------------------------------------------------------------------------------------
theme_set(theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95")))
save2 <- function(p, name, w, h) { ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300); ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf) }
VL <- c(gdd = "GDD", cmi = "CMI")
samp <- sample(ids, 1000)
ln <- melt(cr[id_pe %in% samp], id.vars = c("id_pe", "year"), variable.name = "var", value.name = "value")[!is.na(value)][, var := factor(VL[as.character(var)], VL)]
nsd <- copy(NS)[, var := factor(VL[var], VL)]
fit <- nsd[, .(year = yrs, fit = coef(lm(mean ~ year))[[1]] + coef(lm(mean ~ year))[[2]] * yrs), by = var]
ann <- copy(NT)[, `:=`(var = factor(VL[var], VL), txt = fifelse(var == "gdd", sprintf("%+.0f GDD (%+.2f SD) per decade", per_decade_raw, per_decade_z), sprintf("%+.2f CMI (%+.2f SD) per decade", per_decade_raw, per_decade_z)))]
p1 <- ggplot() + geom_line(data = ln, aes(year, value, group = id_pe), colour = "grey60", linewidth = 0.1, alpha = 0.25) +
  geom_ribbon(data = nsd, aes(year, ymin = lo, ymax = hi), fill = "#08306B", alpha = 0.3) + geom_line(data = nsd, aes(year, mean), colour = "#08306B", linewidth = 0.8) +
  geom_line(data = fit, aes(year, fit), colour = "#B2182B", linetype = 2, linewidth = 0.6) + geom_text(data = ann, aes(x = YMIN, y = Inf, label = txt), hjust = 0, vjust = 1.5, size = 3) +
  facet_wrap(~ var, nrow = 1, scales = "free_y") + labs(x = "year", y = "yearly value at the plot")
save2(p1, "0_05_series_network", 10, 4)
mn <- P[span >= MIN_SPAN]
mkv <- function(unit) rbindlist(lapply(c("gdd", "cmi"), function(v) rbindlist(lapply(SPL, function(sp) {
  i <- which(P$span >= MIN_SPAN & sp[[3]] & !is.na(P[[paste0("sl_two_point_", v)]])); data.table(var = v, split_type = sp[[1]], split = sp[[2]], value = P[[paste0("sl_two_point_", v)]][i] / (if (unit == "z") SD[[v]] else 1)) }))))
for (unit in c("raw", "z")) {
  vd <- mkv(unit); ps <- SUM[estimator == "two_point" & filter == "span >= 10 y"][, .(var, split_type, split, n_plots, mean = get(paste0("mean_", unit)), lo = get(paste0("lo_", unit)), hi = get(paste0("hi_", unit)))]
  SL <- c("all plots", PERL, BAND_LAB); lab_n <- ps[var == "gdd", setNames(sprintf("%s (%s)", split, format(n_plots, big.mark = ",", trim = TRUE)), split)]
  GL <- c(total = "network", period = "timing (period)", latitude = "space (latitude)")
  for (z in list(vd, ps)) z[, `:=`(split = factor(split, rev(SL)), var = factor(VL[var], VL), grp = factor(GL[split_type], GL))]
  p2 <- ggplot(vd, aes(value, split)) + geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_violin(orientation = "y", scale = "width", fill = "grey90", colour = "grey65", linewidth = 0.3) +
    geom_errorbar(data = ps, aes(x = mean, xmin = lo, xmax = hi), width = 0, orientation = "y", colour = "#08306B", linewidth = 0.9) + geom_point(data = ps, aes(x = mean), colour = "#08306B", size = 2) +
    scale_y_discrete(labels = function(x) lab_n[x]) + facet_grid(grp ~ var, scales = "free", space = "free_y") +
    labs(x = if (unit == "raw") "change per decade within the plot (GDD or CMI units)" else "change per decade within the plot (SD)", y = NULL)
  save2(p2, paste0("0_05_slopes_by_split_", unit), 10, 5.5)
}

# ---- 6. console summary -------------------------------------------------------------------------------------------------------
cat("\nMean within-plot change per decade (plots with a span >= ", MIN_SPAN, " y; 95% plot-bootstrap CI):\n", sep = "")
print(SUM[estimator == "two_point" & filter == "span >= 10 y" & split_type == "total", .(var, n_plots, mean_raw = round(mean_raw, 2), lo_raw = round(lo_raw, 2), hi_raw = round(hi_raw, 2), mean_z = round(mean_z, 3), lo_z = round(lo_z, 3), hi_z = round(hi_z, 3), median_raw = round(median_raw, 2), p05_raw = round(p05_raw, 1), p95_raw = round(p95_raw, 1))], row.names = FALSE)
cat("\nChecks next to the main estimate (all plots, per decade, raw):\n")
print(dcast(SUM[split_type == "total", .(var, estimator, filter, n_plots, m = round(mean_raw, 2))], var + estimator ~ filter, value.var = c("m", "n_plots")), row.names = FALSE)
cat("\nBy period and latitude (two-point, span >= 10 y; z per decade):\n")
print(SUM[estimator == "two_point" & filter == "span >= 10 y" & split_type != "total", .(var, split_type, split, n_plots, mean_z = round(mean_z, 3), lo_z = round(lo_z, 3), hi_z = round(hi_z, 3))], row.names = FALSE)
cat("\nNetwork yearly mean, linear trend ", YMIN, "-", YMAX, " per decade:\n", sep = ""); print(NT[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], row.names = FALSE)
f2 <- iv[order(id_pe, interval), .(n = .N, dz = (data.table::last(gdd_3) - data.table::first(gdd_3)) / SD[["gdd"]]), by = id_pe]
cat("\nOld statement (22_hindcast header: median +0.49 SD GDD, first to last visit): model covariate, last vs first interval of the plot: median ", round(f2[n >= 2, median(dz)], 2),
    " SD (plots with >= 2 intervals), ", round(f2[, median(dz)], 2), " SD (all plots)\n", sep = "")
cat("\nSaved in", OUT, "\n")
