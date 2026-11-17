# ==============================================================================
# 0_02_observed_change.R -- Chapter 2 analysis, section 0 (Data, no model), row 2: observed change between visits (thermal index, shade tolerance).
# Run from the project root:  Rscript R/analysis/0_02_observed_change.R
#
# Per visit interval: change in the thermal index (STI) and in shade tolerance between the start and end state, per decade
# (change / interval length x 10); a stay counts as 0. STI of a state = mean z-GDD of the visits ending in that state (as old I1; GDD scaled with
# scaling.rds). Tolerance: TOL_SCORE (R/00_groups.R); an interval with regen or other_temperate at either end has no tolerance change (NA, left out).
# Estimate = mean of the per-interval values (each interval once). 95% CI: plot bootstrap (200 draws, whole plots resampled).
# Choices (Lukas, 2026-11-17): scale = mean per-decade change; groups = thermal group, 13 origin states, shade-tolerance group (of the origin);
#   splits = all intervals, period of the end visit, latitude tertile of the plot; each interval once; observed start -> end occupancy matrix
#   (13 x 13, all intervals, n per origin and median interval length).
# Inputs (read-only): data/tables/<RUN_TAG>/visit_intervals.rds; <FIT_DIR>/scaling.rds; data/raw/pep/PEP.gpkg (plot latitude); R/00_groups.R.
# Checks against old output (read-only): figures/I/I1_*.csv, I1b_*.csv (point estimates), 0_data/0_01_events_by_transition.csv (event counts).
# Outputs: data-output/<RUN_TAG>/analysis/0_data/
#   0_02_observed_change.csv        group_type, group, split_type, split, n, STI and tolerance change per decade with CI
#   0_02_occupancy_matrix.csv       from, to, n, n_from, share, median_len_y (stays on the diagonal)
#   0_02_state_index.csv            STI per state
#   0_02_change_*.(png|pdf), 0_02_occupancy_heatmap.(png|pdf)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(here) })
source(here::here("R", "config.R"))      # RUN_TAG, TABLE_DIR, FIT_DIR, RUN_OUT_DIR, FIG_DIR
source(here::here("R", "00_groups.R"))   # GRP, TOL_GROUP, TOL_SCORE
OUT <- file.path(RUN_OUT_DIR, "analysis", "0_data"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
options(width = 200); set.seed(25)
N_BOOT <- 200L

ORDER <- c("regen", "jack_pine", "white_spruce", "black_spruce", "fir", "white_birch", "trembling_aspen", "other_pioneer",
           "yellow_birch", "red_maple", "sugar_maple", "cedar", "other_temperate")
GCOL <- c(all = "black", boreal = "#08306B", pioneer = "#E6B800", regen = "#777777", temperate = "#B2182B")
TCOL <- c(all = "black", intolerant = "#E6B800", intermediate = "#238B45", tolerant = "#08306B", regen = "#777777", none = "#9E9AC8")
stopifnot("TOL_SCORE not in 00_groups.R" = exists("TOL_SCORE"), setequal(names(TOL_SCORE), ORDER), setequal(names(GRP), ORDER))

iv  <- as.data.table(readRDS(file.path(TABLE_DIR, "visit_intervals.rds")))
scz <- readRDS(file.path(FIT_DIR, "scaling.rds"))
con <- DBI::dbConnect(RSQLite::SQLite(), here("data", "raw", "pep", "PEP.gpkg"), flags = RSQLite::SQLITE_RO)
pl  <- unique(as.data.table(DBI::dbGetQuery(con, "SELECT id_pe, latitude FROM placette")), by = "id_pe"); DBI::dbDisconnect(con)
stopifnot("states outside the 13" = all(c(iv$state_start, iv$state_end) %in% ORDER),
          "plot without latitude" = all(iv$id_pe %in% pl[!is.na(latitude), id_pe]))

# ---- 1. thermal index of each state; change per interval ------------------------------------------------------------------
zg  <- (iv$gdd_3 - scz[var == "gdd_3", mean]) / scz[var == "gdd_3", sd]
STI <- setNames(tapply(zg, iv$state_end, mean), sort(unique(iv$state_end)))[ORDER]
stopifnot("a state is never an end state" = !anyNA(STI))
fwrite(data.table(state = ORDER, sti = unname(STI), n_end_visits = as.integer(table(factor(iv$state_end, ORDER)))), file.path(OUT, "0_02_state_index.csv"))

per_of <- function(y) factor(cut(y, c(-Inf, 1990, 2000, 2010, Inf), labels = c("≤ 1990", "1991–2000", "2001–2010", "> 2010")))
q <- quantile(pl[id_pe %in% iv$id_pe, latitude], c(1/3, 2/3))
BAND_LAB <- c(sprintf("south (≤ %.1f°N)", q[1]), sprintf("middle (%.1f–%.1f°N)", q[1], q[2]), sprintf("north (> %.1f°N)", q[2]))
d <- merge(iv[, .(id_pe, a = state_start, b = state_end, len = year_end - year_start, per = per_of(year_end))], pl, by = "id_pe", sort = FALSE)
d[, band := factor(cut(latitude, c(-Inf, q, Inf), labels = BAND_LAB), levels = BAND_LAB)]
d[, `:=`(dsti = (STI[b] - STI[a]) / len * 10, dtol = (TOL_SCORE[b] - TOL_SCORE[a]) / len * 10,
         thermal = GRP[a], state = a, tolerance = fifelse(is.na(TOL_GROUP[a]), "none", TOL_GROUP[a]))]
stopifnot(nrow(d) == nrow(iv), !anyNA(d$dsti), !anyNA(d$per), !anyNA(d$band), all(d$len > 0))

# ---- 2. cells: (group type x group) x (split type x split) -----------------------------------------------------------------
gdim <- list(all = function(x) rep("all", nrow(x)), thermal = function(x) x$thermal, state = function(x) x$state, tolerance = function(x) x$tolerance)
sdim <- list(total = function(x) rep("all intervals", nrow(x)), period = function(x) as.character(x$per), latitude = function(x) as.character(x$band))
long <- rbindlist(lapply(names(gdim), function(gt) rbindlist(lapply(names(sdim), function(st)
  data.table(group_type = gt, group = gdim[[gt]](d), split_type = st, split = sdim[[st]](d), id_pe = d$id_pe, dsti = d$dsti, dtol = d$dtol)))))
A <- long[, .(n = .N, s_sti = sum(dsti), n_tol = sum(!is.na(dtol)), s_tol = sum(dtol, na.rm = TRUE)), by = .(group_type, group, split_type, split, id_pe)]
A[, cell := .GRP, by = .(group_type, group, split_type, split)]
P <- uniqueN(d$id_pe); A[, pidx := match(id_pe, unique(d$id_pe))]
est <- function(w) { A[, wi := w[pidx]]; A[, .(n = sum(wi * n), sti = sum(wi * s_sti) / sum(wi * n), n_tol = sum(wi * n_tol),
                                                tol = fifelse(sum(wi * n_tol) > 0, sum(wi * s_tol) / sum(wi * n_tol), NA_real_)), by = cell][order(cell)] }
E0 <- est(rep(1, P)); stopifnot(E0[, sum(n)] == nrow(long))
t0 <- Sys.time()
B <- rbindlist(lapply(seq_len(N_BOOT), function(r) est(tabulate(sample.int(P, P, replace = TRUE), P))[, r := r]))
cat(sprintf("bootstrap: %d draws, %.0f s\n", N_BOOT, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
CI <- B[, .(sti_lo = quantile(sti, .025, na.rm = TRUE), sti_hi = quantile(sti, .975, na.rm = TRUE),
            tol_lo = quantile(tol, .025, na.rm = TRUE), tol_hi = quantile(tol, .975, na.rm = TRUE)), by = cell]
key <- unique(A[, .(cell, group_type, group, split_type, split)])
res <- Reduce(function(x, y) merge(x, y, by = "cell"), list(key, E0, CI))
res[, `:=`(n = as.integer(n), n_tol = as.integer(n_tol), cell = NULL)]
res[is.na(tol), c("tol_lo", "tol_hi") := NA_real_]
# natural order (not alphabetical): group types, groups, split types, splits (south-middle-north, periods in time order)
res[, `:=`(group_type = factor(group_type, c("all", "thermal", "state", "tolerance")), split_type = factor(split_type, c("total", "period", "latitude")),
           group = factor(group, unique(c("all", "boreal", "pioneer", "regen", "temperate", "intolerant", "intermediate", "tolerant", "none", ORDER))),
           split = factor(split, c("all intervals", levels(d$per), BAND_LAB)))]
setorder(res, group_type, group, split_type, split)
res[, `:=`(group_type = as.character(group_type), group = as.character(group), split_type = as.character(split_type), split = as.character(split))]
fwrite(res, file.path(OUT, "0_02_observed_change.csv"))

# ---- 3. checks against old I1 / I1b point estimates ----------------------------------------------------------------------
fi <- file.path(FIG_DIR, "I")
o1 <- fread(file.path(fi, "I1_raw_thermophilization.csv")); o2 <- fread(file.path(fi, "I1b_raw_thermophilization_latitude.csv"))
m1 <- merge(o1[, .(group = grp, split = fifelse(per == "all periods", "all intervals", per), n, sti, tol, n_tol)],
            res[group_type %in% c("all", "thermal") & split_type %in% c("total", "period")], by = c("group", "split"), suffixes = c("_old", "_new"))
m2 <- merge(o2[, .(group = grp, rk = rowid(grp), n, sti, tol)][, rk := rowid(group)],
            res[group_type %in% c("all", "thermal") & split_type == "latitude"][, rk := rowid(group)], by = c("group", "rk"), suffixes = c("_old", "_new"))
same <- function(a, b) all(abs(a - b) < 1e-9 | (is.na(a) & is.na(b)))
stopifnot("I1 cells not all found"  = nrow(m1) == nrow(o1),  "I1b cells not all found" = nrow(m2) == nrow(o2),
          "n differs from I1"       = m1[, all(n_old == n_new)],  "STI differs from I1"        = m1[, same(sti_old, sti_new)],
          "tolerance differs from I1" = m1[, same(tol_old, tol_new)], "n_tol differs from I1"   = m1[, all(n_tol_old == n_tol_new)],
          "n differs from I1b"      = m2[, all(n_old == n_new)],  "STI differs from I1b"       = m2[, same(sti_old, sti_new)],
          "tolerance differs from I1b" = m2[, same(tol_old, tol_new)])
cat(sprintf("check vs old I1 (%d cells) and I1b (%d cells): n, STI, tolerance change identical\n", nrow(m1), nrow(m2)))

# ---- 4. observed occupancy: start -> end shares (13 x 13, stays on the diagonal) --------------------------------------------
oc <- iv[, .(n = .N), by = .(from = state_start, to = state_end)]
oc <- merge(CJ(from = ORDER, to = ORDER), oc, by = c("from", "to"), all.x = TRUE); oc[is.na(n), n := 0L]
oc <- merge(oc, iv[, .(n_from = .N, median_len_y = as.numeric(median(year_end - year_start))), by = .(from = state_start)], by = "from")
oc[, share := n / n_from]
r1 <- fread(file.path(OUT, "0_01_events_by_transition.csv"))
stopifnot("matrix total != intervals" = oc[, sum(n)] == nrow(iv), "row sums != n_from" = oc[, sum(n), by = from][oc[, .(n_from = n_from[1]), by = from], on = "from"][, all(V1 == n_from)],
          "off-diagonal counts differ from row 1" = merge(oc[from != to], r1, by = c("from", "to"))[, all(n == n_events)])
fwrite(oc[order(match(from, ORDER), match(to, ORDER))], file.path(OUT, "0_02_occupancy_matrix.csv"))
cat(sprintf("check: matrix total %d intervals; off-diagonal counts = row 1 events (%d)\n", oc[, sum(n)], oc[from != to, sum(n)]))

# ---- 5. figures ---------------------------------------------------------------------------------------------------------
theme_set(theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95")))
lab <- function(x) gsub("_", " ", x)
save2 <- function(p, name, w, h) { ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300); ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf) }
FACET <- c(sti = "thermal index", tol = "shade tolerance")
long_res <- function(r) { L <- melt(r, id.vars = c("group", "split"), measure.vars = list(c("sti", "tol"), c("sti_lo", "tol_lo"), c("sti_hi", "tol_hi")),
                                     variable.name = "q", value.name = c("est", "lo", "hi")); L[, q := factor(FACET[q], levels = FACET)][!is.na(est)] }
fig_cells <- function(gt, st, pal, xlab, name, splits, glev, w = 11, h = 4.6) {
  r <- res[(group_type == gt | group_type == "all") & split_type == st]
  r[, `:=`(group = factor(group, levels = glev), split = factor(split, levels = splits))]
  L <- long_res(r[!is.na(group)])
  p <- ggplot(L, aes(split, est, colour = group, group = group)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, position = position_dodge(width = 0.6)) +
    geom_point(position = position_dodge(width = 0.6), size = 1.8) +
    scale_colour_manual(values = pal, name = NULL, labels = function(x) ifelse(x == "none", "other temperate (no group)", x)) + facet_wrap(~ q, scales = "free_y") + labs(x = xlab, y = "observed change per decade")
  save2(p, name, w, h)
}
PER <- levels(d$per)
fig_cells("thermal",   "period",   GCOL, "period of the end visit", "0_02_change_thermal_by_period",   PER,      names(GCOL))
fig_cells("thermal",   "latitude", GCOL, "latitude band",           "0_02_change_thermal_by_latitude", BAND_LAB, names(GCOL))
fig_cells("tolerance", "period",   TCOL, "period of the end visit", "0_02_change_tolerance_by_period",   PER,      names(TCOL))
fig_cells("tolerance", "latitude", TCOL, "latitude band",           "0_02_change_tolerance_by_latitude", BAND_LAB, names(TCOL))
rs <- res[group_type == "state" & split_type == "total"]; rs[, group := factor(group, levels = rev(ORDER))]
Ls <- long_res(rs); Ls[, thermal := GRP[as.character(group)]]
p <- ggplot(Ls, aes(est, group, colour = thermal)) + geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y") + geom_point(size = 1.8) + scale_colour_manual(values = GCOL, name = "thermal group") +
  scale_y_discrete(labels = lab) + facet_wrap(~ q, scales = "free_x") + labs(x = "observed change per decade", y = "origin state")
save2(p, "0_02_change_by_state", 9, 5)

# occupancy heatmap: share of intervals (%) from origin (row) ending in destination (column); diagonal = stays (grey)
g <- copy(oc)[, pct := 100 * share]
g[, `:=`(ylab = sprintf("%s (%s)", lab(from), format(n_from, big.mark = ",", trim = TRUE)), txt = fifelse(round(pct, 1) >= 10, sprintf("%.0f", pct), fifelse(pct >= 0.1, sprintf("%.1f", pct), "<0.1")))]
YL <- g[, .(ylab = ylab[1]), by = from][match(ORDER, from), ylab]
g[, `:=`(ylab = factor(ylab, levels = rev(YL)), to = factor(to, levels = ORDER), fillv = log10(pmax(pct, 0.01)))]
off <- g[from != to & n > 0]
ph <- ggplot(g, aes(to, ylab)) + geom_tile(data = g[from == to], fill = "grey85", colour = "white") +
  geom_tile(data = g[from != to & n == 0], fill = "grey96", colour = "white") +
  geom_tile(data = off, aes(fill = fillv), colour = "white") +
  geom_text(data = g[from == to], aes(label = txt), size = 2.1, colour = "grey30") +
  geom_text(data = off, aes(label = txt, colour = fillv > 0.65 * max(off$fillv) + 0.35 * min(off$fillv)), size = 2.1) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  scale_fill_viridis_c(option = "magma", direction = -1, breaks = log10(c(0.1, 1, 10)), labels = c(0.1, 1, 10)) +
  scale_x_discrete(labels = lab, position = "top", drop = FALSE) + scale_y_discrete(drop = FALSE) +
  labs(x = "destination", y = "origin (number of intervals)", fill = "% of intervals\n(log scale)", subtitle = "Grey diagonal: stays.") +
  coord_equal() + theme(axis.text.x = element_text(angle = 45, hjust = 0), panel.grid = element_blank())
save2(ph, "0_02_occupancy_heatmap", 8, 7)

# ---- console summary ----------------------------------------------------------------------------------------------------
cat("\nNetwork and thermal groups, all intervals:\n"); print(res[group_type %in% c("all", "thermal") & split_type == "total", .(group, n, sti, sti_lo, sti_hi, n_tol, tol, tol_lo, tol_hi)], row.names = FALSE)
cat("\nShade-tolerance groups, all intervals:\n"); print(res[group_type == "tolerance" & split_type == "total", .(group, n, sti, sti_lo, sti_hi, n_tol, tol, tol_lo, tol_hi)], row.names = FALSE)
cat("\n13 origin states, all intervals:\n"); print(res[group_type == "state" & split_type == "total"][match(ORDER, group), .(group, n, sti, sti_lo, sti_hi, n_tol, tol, tol_lo, tol_hi)], row.names = FALSE)
cat("\nStays (share of intervals):\n"); print(oc[from == to][match(ORDER, from), .(state = from, n_from, stay_pct = round(100 * share, 1), median_len_y)], row.names = FALSE)
cat("\nSaved in", OUT, "\n")
