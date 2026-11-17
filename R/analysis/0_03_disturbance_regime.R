# ==============================================================================
# 0_03_disturbance_regime.R -- Chapter 2 analysis, section 0 (Data, no model), row 3: disturbance frequency, co-occurrence, climate position.
# Run from the project root:  Rscript R/analysis/0_03_disturbance_regime.R
#
# Flags are the model's covariates: per visit interval, the disturbance status at the END visit (fire 0/1, pest 0/1, harvest 0/1 light/2 heavy).
# A partial event stays active until a stand-replacing event ends it (R/03a_disturbance_events.R), so a flag is either
#   NEW    = an event of that type starts at the interval's end visit (event table: not dropped, start_no_mes == visit_end; as old I6b), or
#   LEGACY = the flag is active but its event started earlier (including events active at a plot's first visit).
# Choices (Lukas, 2026-11-17): frequency = share of intervals flagged, split new / legacy, plus the legacy fraction of flags; groups = network,
#   origin thermal group, origin shade-tolerance group, 13 origin states (origin = state at interval start, flag at the end visit);
#   splits = all intervals, period of the end visit, latitude tertile; co-occurrence = active flag combinations (with new/legacy per flag) AND
#   combinations of new events starting in the same interval (counted from the event table); climate position = new events only (as I6b); return intervals as I6c;
#   time since disturbance (TSD) of flagged intervals as median and percentiles, split new / legacy.
# Scope note: climate position and return intervals are network-wide (no group / split); frequency, co-occurrence and TSD carry the splits.
# Inputs (read-only): data/tables/<RUN_TAG>/visit_intervals.rds; data/build/disturbance_events.rds; <FIT_DIR>/scaling.rds; PEP.gpkg (latitude);
#   R/00_groups.R.
# Checks against old output (read-only): figures/A/A1f_*.csv (flag shares by state), figures/I/I6_flag_combinations.csv, I6_climate_position.csv,
#   I6_return_intervals.csv.
# Outputs: data-output/<RUN_TAG>/analysis/0_data/0_03_*
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(here) })
source(here::here("R", "config.R")); source(here::here("R", "00_groups.R"))
OUT <- file.path(RUN_OUT_DIR, "analysis", "0_data"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
options(width = 200); set.seed(26); N_BOOT <- 200L
ORDER <- c("regen", "jack_pine", "white_spruce", "black_spruce", "fir", "white_birch", "trembling_aspen", "other_pioneer",
           "yellow_birch", "red_maple", "sugar_maple", "cedar", "other_temperate")
KINDS <- c("fire", "pest", "light harvest", "heavy harvest"); KS <- c("fire", "pest", "light", "heavy")
stopifnot(setequal(names(GRP), ORDER), setequal(names(TOL_GROUP), ORDER))

iv  <- as.data.table(readRDS(file.path(TABLE_DIR, "visit_intervals.rds")))
scz <- readRDS(file.path(FIT_DIR, "scaling.rds")); zf <- function(x, v) (x - scz[var == v, mean]) / scz[var == v, sd]
con <- DBI::dbConnect(RSQLite::SQLite(), here("data", "raw", "pep", "PEP.gpkg"), flags = RSQLite::SQLITE_RO)
pl  <- unique(as.data.table(DBI::dbGetQuery(con, "SELECT id_pe, latitude FROM placette")), by = "id_pe"); DBI::dbDisconnect(con)
stopifnot("plot without latitude" = all(iv$id_pe %in% pl[!is.na(latitude), id_pe]), all(iv$harv %in% 0:2), all(iv$fire %in% 0:1), all(iv$pest %in% 0:1))
iv[, `:=`(latitude = pl$latitude[match(id_pe, pl$id_pe)], z_gdd = zf(gdd_3, "gdd_3"), z_cmi = zf(cmi_gs, "cmi_gs"), z_soil = zf(soil_pc1, "soil_pc1"))]
per_of <- function(y) factor(cut(y, c(-Inf, 1990, 2000, 2010, Inf), labels = c("≤ 1990", "1991–2000", "2001–2010", "> 2010")))
q <- quantile(pl[id_pe %in% iv$id_pe, latitude], c(1/3, 2/3))
BAND_LAB <- c(sprintf("south (≤ %.1f°N)", q[1]), sprintf("middle (%.1f–%.1f°N)", q[1], q[2]), sprintf("north (> %.1f°N)", q[2]))
iv[, `:=`(per = per_of(year_end), band = factor(cut(latitude, c(-Inf, q, Inf), labels = BAND_LAB), levels = BAND_LAB),
          thermal = GRP[state_start], state = state_start, tolerance = fifelse(is.na(TOL_GROUP[state_start]), "none", TOL_GROUP[state_start]))]
stopifnot(!anyNA(iv$per), !anyNA(iv$band))

# ---- 1. events -> new / legacy flags -------------------------------------------------------------------------------------
ev <- as.data.table(readRDS(file.path(BUILD_DIR, "disturbance_events.rds"))$events)[dropped == FALSE & id_pe %in% iv$id_pe]
ev[, kind := fcase(type == "harv" & level == 1, "light harvest", type == "harv" & level == 2, "heavy harvest", type == "fire", "fire", type == "pest", "pest")]
stopifnot(!anyNA(ev$kind))
fv <- iv[, .(first_visit = min(visit_start)), by = id_pe]; ev <- fv[ev, on = "id_pe"]
ev[, matched := paste(id_pe, start_no_mes) %in% paste(iv$id_pe, iv$visit_end)]
ev[, where := fcase(start_no_mes < first_visit, "before the first visit", start_no_mes == first_visit, "at the first visit",
                    matched, "new: at the end visit of an interval", default = "later visit, no interval ends there")]
acc <- ev[, .N, by = where][order(-N)]; cat("Events in retained plots (not dropped), by where they start:\n"); print(acc, row.names = FALSE)
fwrite(rbind(acc, data.table(where = "total", N = nrow(ev))), file.path(OUT, "0_03_event_accounting.csv"))
EVN <- ev[where == "new: at the end visit of an interval"]
for (i in seq_along(KINDS)) {
  nn <- EVN[kind == KINDS[i], .N, by = .(id_pe, visit_end = start_no_mes)]
  iv[nn, on = c("id_pe", "visit_end"), paste0("nn_", KS[i]) := i.N]; set(iv, which(is.na(iv[[paste0("nn_", KS[i])]])), paste0("nn_", KS[i]), 0L)
}
iv[, `:=`(flag_fire = fire == 1L, flag_pest = pest == 1L, flag_light = harv == 1L, flag_heavy = harv == 2L)]
# A new event can lack an active flag only when another event starts in the same interval and supersedes it (a later stand-replacing event ends the
# older one from its start; R/03a_disturbance_events.R, step 4 and 6). Those events count as new events, but have no flag at the end visit.
NN <- as.matrix(iv[, .(nn_fire, nn_pest, nn_light, nn_heavy)]); FL <- as.matrix(iv[, .(flag_fire, flag_pest, flag_light, flag_heavy)])
sup <- (NN > 0) & !FL; other_new <- rowSums(NN > 0) - (NN > 0) > 0
cat("\nnew events at an end visit without an active flag (superseded in the same interval):", paste(KS, colSums(sup), collapse = "  "), "\n")
stopifnot("a new event has no active flag and no other event starts in its interval: inspect" = all(!sup | other_new))
fwrite(iv[rowSums(sup) > 0, .(id_pe, visit_end, year_end, fire, pest, harv, nn_fire, nn_pest, nn_light, nn_heavy)], file.path(OUT, "0_03_superseded_new_events.csv"))
for (k in KS) { iv[, paste0("new_", k) := get(paste0("flag_", k)) & get(paste0("nn_", k)) > 0]; iv[, paste0("leg_", k) := get(paste0("flag_", k)) & !get(paste0("new_", k))] }
cat(sprintf("flagged intervals: %s\n", paste(sprintf("%s %d (new %d)", KS, sapply(KS, function(k) iv[, sum(get(paste0("flag_", k)))]), sapply(KS, function(k) iv[, sum(get(paste0("new_", k)))])), collapse = "; ")))

# ---- 2. frequency: share of intervals flagged, new / legacy -----------------------------------------------------------------
gdim <- list(all = function(x) rep("all", nrow(x)), thermal = function(x) x$thermal, tolerance = function(x) x$tolerance, state = function(x) x$state)
sdim <- list(total = function(x) rep("all intervals", nrow(x)), period = function(x) as.character(x$per), latitude = function(x) as.character(x$band))
agg <- rbindlist(lapply(names(gdim), function(gt) rbindlist(lapply(names(sdim), function(st)
  iv[, c(list(n_intervals = .N), unlist(lapply(KS, function(k) setNames(list(sum(get(paste0("flag_", k))), sum(get(paste0("new_", k)))), paste0(k, c("_flagged", "_new")))), recursive = FALSE)),
     by = .(group = gdim[[gt]](iv), split = sdim[[st]](iv))][, `:=`(group_type = gt, split_type = st)]))))
fr <- rbindlist(lapply(seq_along(KS), function(i) agg[, .(group_type, group, split_type, split, kind = KINDS[i], n_intervals,
                                                           n_flagged = get(paste0(KS[i], "_flagged")), n_new = get(paste0(KS[i], "_new")))]))
fr[, `:=`(n_legacy = n_flagged - n_new, share_flagged = n_flagged / n_intervals, share_new = n_new / n_intervals,
          share_legacy = (n_flagged - n_new) / n_intervals, legacy_fraction = fifelse(n_flagged > 0, (n_flagged - n_new) / n_flagged, NA_real_))]
fr[, `:=`(group_type = factor(group_type, c("all", "thermal", "tolerance", "state")), split_type = factor(split_type, c("total", "period", "latitude")),
          group = factor(group, unique(c("all", "boreal", "pioneer", "regen", "temperate", "intolerant", "intermediate", "tolerant", "none", ORDER))),
          split = factor(split, c("all intervals", levels(iv$per), BAND_LAB)), kind = factor(kind, KINDS))]
setorder(fr, group_type, group, split_type, split, kind)
fr[, `:=`(group_type = as.character(group_type), group = as.character(group), split_type = as.character(split_type), split = as.character(split), kind = as.character(kind))]
fwrite(fr, file.path(OUT, "0_03_flag_shares.csv"))

# ---- 3. co-occurrence ------------------------------------------------------------------------------------------------------
combo_of <- function(...) apply(cbind(...), 1, function(r) { r <- r[!is.na(r)]; if (length(r)) paste(r, collapse = "+") else "none" })
iv[, combo := combo_of(fifelse(flag_fire, "fire", NA_character_), fifelse(flag_pest, "pest", NA_character_), fifelse(flag_light, "light harvest", NA_character_), fifelse(flag_heavy, "heavy harvest", NA_character_))]
st_lab <- function(flag, new, nm) fifelse(flag, paste0(nm, fifelse(new, "(new)", "(legacy)")), NA_character_)
iv[, combo_status := combo_of(st_lab(flag_fire, new_fire, "fire"), st_lab(flag_pest, new_pest, "pest"), st_lab(flag_light, new_light, "light harvest"), st_lab(flag_heavy, new_heavy, "heavy harvest"))]
# new events starting in the interval, counted from the event table (includes the events superseded within the interval)
iv[, combo_new := combo_of(fifelse(nn_fire > 0, "fire", NA_character_), fifelse(nn_pest > 0, "pest", NA_character_), fifelse(nn_light > 0, "light harvest", NA_character_), fifelse(nn_heavy > 0, "heavy harvest", NA_character_))]
by_split <- function(expr) rbindlist(lapply(names(sdim), function(st) eval(substitute(expr), list(S = quote(sdim[[st]](iv)), st = st))))
co1 <- rbindlist(lapply(names(sdim), function(st) iv[, .(n = .N), by = .(split = sdim[[st]](iv), combo, combo_status)][, `:=`(split_type = st, share = n / sum(n)), by = split]))
co1 <- co1[, .(split_type, split, combo, combo_status, n, share)][order(match(split_type, names(sdim)), match(split, c("all intervals", levels(iv$per), BAND_LAB)), -n)]
fwrite(co1, file.path(OUT, "0_03_cooccurrence_flags.csv"))
co2 <- rbindlist(lapply(names(sdim), function(st) {
  tot <- iv[, .(n_all = .N), by = .(split = sdim[[st]](iv))]
  sub <- iv[combo_new != "none"]
  sub[, .(n = .N), by = .(split = sdim[[st]](sub), combo_new)][tot, on = "split"][, `:=`(split_type = st, share_of_intervals = n / n_all)][
    , share_of_new_event_intervals := n / sum(n), by = split]
}))
co2 <- co2[, .(split_type, split, combo_new, n, share_of_intervals, share_of_new_event_intervals)][order(match(split_type, names(sdim)), match(split, c("all intervals", levels(iv$per), BAND_LAB)), -n)]
fwrite(co2, file.path(OUT, "0_03_cooccurrence_new_events.csv"))

# ---- 4. time since disturbance of flagged intervals, new vs legacy -------------------------------------------------------------
TSDC <- c("tsd_fire_y", "tsd_pest_y", "tsd_harv_y", "tsd_harv_y"); DATC <- c("fire_dating", "pest_dating", "harv_dating", "harv_dating")
tsdt <- rbindlist(lapply(seq_along(KS), function(i) iv[get(paste0("flag_", KS[i])), .(kind = KINDS[i], status = fifelse(get(paste0("new_", KS[i])), "new", "legacy"),
                                                                                   tsd = get(TSDC[i]), dated = get(DATC[i]) == "dated", per, band)]))
stopifnot(!anyNA(tsdt$tsd))
tsdt2 <- rbind(tsdt, tsdt[status == "legacy"][, status := fifelse(dated, "legacy: dated", "legacy: imputed")])
tq <- rbindlist(lapply(names(sdim), function(st) tsdt2[, .(n = .N, n_dated = sum(dated), min = min(tsd), p05 = quantile(tsd, .05), p25 = quantile(tsd, .25), median = median(tsd),
                                                         p75 = quantile(tsd, .75), p95 = quantile(tsd, .95), max = max(tsd)), by = .(kind, status, split = sdim[[st]](tsdt2))][, split_type := st]))
tq <- rbind(tq, tsdt[, .(n = .N, n_dated = sum(dated), min = min(tsd), p05 = quantile(tsd, .05), p25 = quantile(tsd, .25), median = median(tsd), p75 = quantile(tsd, .75),
                         p95 = quantile(tsd, .95), max = max(tsd)), by = .(kind)][, `:=`(status = "all flagged", split = "all intervals", split_type = "total")], fill = TRUE)
tq[, `:=`(kind = factor(kind, KINDS), status = factor(status, c("new", "legacy", "legacy: dated", "legacy: imputed", "all flagged")), split_type = factor(split_type, c("total", "period", "latitude")),
          split = factor(split, c("all intervals", levels(iv$per), BAND_LAB)))]
setorder(tq, split_type, split, kind, status); tq[, `:=`(kind = as.character(kind), status = as.character(status), split_type = as.character(split_type), split = as.character(split))]
setcolorder(tq, c("kind", "status", "split_type", "split")); fwrite(tq, file.path(OUT, "0_03_tsd_flagged.csv"))

# ---- 5. climate position of NEW events (as old I6b) -------------------------------------------------------------------------
ne <- merge(EVN[, .(id_pe, kind, visit_end = start_no_mes)], iv[, .(id_pe, visit_end, z_gdd, z_cmi, z_soil)], by = c("id_pe", "visit_end"))
stopifnot("an event did not find its interval" = nrow(ne) == nrow(EVN))
X <- rbind(ne[, .(id_pe, kind, z_gdd, z_cmi, z_soil)], iv[, .(id_pe, kind = "all intervals", z_gdd, z_cmi, z_soil)])
ids <- unique(iv$id_pe); P <- length(ids)
A <- X[, .(n = .N, sg = sum(z_gdd), sc = sum(z_cmi), ss = sum(z_soil, na.rm = TRUE), ns = sum(!is.na(z_soil))), by = .(kind, id_pe)][, pidx := match(id_pe, ids)]
cp_est <- function(w) { A[, wi := w[pidx]]; A[, .(n = sum(wi * n), gdd = sum(wi * sg) / sum(wi * n), cmi = sum(wi * sc) / sum(wi * n), soil = sum(wi * ss) / sum(wi * ns)), by = kind] }
CP0 <- cp_est(rep(1, P)); CPB <- rbindlist(lapply(seq_len(N_BOOT), function(r) cp_est(tabulate(sample.int(P, P, replace = TRUE), P))))
CP <- melt(CP0, id.vars = c("kind", "n"), variable.name = "cov", value.name = "est")
CP <- merge(CP, melt(CPB, id.vars = c("kind", "n"), variable.name = "cov")[, .(lo = quantile(value, .025, na.rm = TRUE), hi = quantile(value, .975, na.rm = TRUE), lo90 = quantile(value, .05, na.rm = TRUE), hi90 = quantile(value, .95, na.rm = TRUE)), by = .(kind, cov)], by = c("kind", "cov"))
CP[, n := as.integer(CP0$n[match(kind, CP0$kind)])]
CPQ <- melt(X, id.vars = c("id_pe", "kind"), measure.vars = c(gdd = "z_gdd", cmi = "z_cmi", soil = "z_soil"), variable.name = "cov")[!is.na(value), .(n_obs = .N, p05 = quantile(value, .05), p25 = quantile(value, .25),
  median = median(value), p75 = quantile(value, .75), p95 = quantile(value, .95)), by = .(kind, cov)]
CPQ[, cov := factor(c("z_gdd" = "gdd", "z_cmi" = "cmi", "z_soil" = "soil")[as.character(cov)], levels(CP$cov))]
CP <- merge(CP, CPQ, by = c("kind", "cov")); fwrite(CP[order(match(kind, c("all intervals", KINDS)), match(cov, c("gdd", "cmi", "soil")))], file.path(OUT, "0_03_climate_position_new_events.csv"))

# ---- 6. return intervals (dated events; as old I6c) ---------------------------------------------------------------------------
de <- ev[!is.na(yr)][order(id_pe, yr)]
ri_any  <- de[, .(gap = as.numeric(diff(yr)), pair = paste(head(kind, -1), "→", tail(kind, -1))), by = id_pe][gap > 0]
ri_same <- de[, .(gap = as.numeric(diff(yr))), by = .(id_pe, kind)][gap > 0]
sm <- function(x) list(n = length(x), p05 = quantile(x, .05), p25 = quantile(x, .25), median = median(x), p75 = quantile(x, .75), p95 = quantile(x, .95))
ri <- rbind(ri_any[, c(list(what = "any type to any type"), sm(gap))], ri_same[, c(list(what = paste("same type:", kind[1])), sm(gap)), by = kind][, kind := NULL][match(paste("same type:", KINDS), what)])
fwrite(ri, file.path(OUT, "0_03_return_intervals.csv"))

# ---- 7. checks against old outputs ----------------------------------------------------------------------------------------------
same <- function(a, b) all(abs(a - b) < 1e-9 | (is.na(a) & is.na(b)))
a1f <- fread(file.path(FIG_DIR, "A", "A1f_disturbance_by_state.csv"))[, kind := c(fire = "fire", pest = "pest", harv_light = "light harvest", harv_heavy = "heavy harvest")[as.character(disturbance)]]
m <- merge(a1f[, .(group = state, kind, n_old = n, share_old = share)], fr[group_type == "state" & split_type == "total", .(group, kind, n_new = n_intervals, share_new = share_flagged)], by = c("group", "kind"))
stopifnot("A1f cells not all found" = nrow(m) == nrow(a1f), "A1f n differs" = m[, all(n_old == n_new)], "A1f flagged share differs" = m[, same(share_old, share_new)])
cat(sprintf("check vs old A1f (%d cells): n and share flagged identical\n", nrow(m)))
i6a <- fread(file.path(FIG_DIR, "I", "I6_flag_combinations.csv")); mine <- iv[, .(N = .N), by = .(combo = combo)]
mc <- merge(i6a[, .(combo, N_old = N)], mine, by = "combo", all = TRUE)
stopifnot("I6a combos differ" = mc[, !anyNA(N_old) && !anyNA(N) && all(N_old == N)], "status split does not add up" = co1[split_type == "total", sum(n)] == nrow(iv))
cat(sprintf("check vs old I6a (%d combinations): counts identical\n", nrow(mc)))
i6b <- fread(file.path(FIG_DIR, "I", "I6_climate_position.csv")); mb <- merge(i6b[, .(kind, cov, n_old = n, est_old = est)], CP[, .(kind, cov, n_new = n, est_new = est)], by = c("kind", "cov"))
stopifnot("I6b cells not all found" = nrow(mb) == nrow(i6b), "I6b n differs" = mb[, all(n_old == n_new)], "I6b means differ" = mb[, same(est_old, est_new)])
cat(sprintf("check vs old I6b (%d cells): n and means identical\n", nrow(mb)))
i6c <- fread(file.path(FIG_DIR, "I", "I6_return_intervals.csv"), colClasses = c(id_pe = "character"))
stopifnot("I6c rows differ" = nrow(i6c) == nrow(ri_any), "I6c gaps differ" = identical(as.numeric(sort(i6c$gap)), as.numeric(sort(ri_any$gap))),
          "I6c pairs differ" = identical(i6c[, .N, by = pair][order(pair)]$N, ri_any[, .N, by = pair][order(pair)]$N))
cat(sprintf("check vs old I6c (%d return intervals): gaps and pairs identical\n", nrow(i6c)))

# ---- 8. figures ---------------------------------------------------------------------------------------------------------------
theme_set(theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95")))
lab <- function(x) ifelse(x == "none", "other temperate (no group)", gsub("_", " ", x))
save2 <- function(p, name, w, h) { ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300); ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf) }
SC <- c(new = "#08306B", legacy = "#9ECAE1")
stk <- function(x) melt(x, id.vars = c("group_type", "group", "split_type", "split", "kind"), measure.vars = c("share_new", "share_legacy"), variable.name = "status", value.name = "share")[
  , `:=`(status = factor(sub("share_", "", status), c("legacy", "new")), kind = factor(kind, KINDS))]
g1 <- stk(fr[split_type == "total"]); g1[, `:=`(group_type = factor(group_type, c("all", "thermal", "tolerance", "state"), c("network", "thermal group", "tolerance group", "origin state")),
                                               group = factor(group, rev(unique(c("all", "boreal", "pioneer", "regen", "temperate", "intolerant", "intermediate", "tolerant", "none", ORDER)))))]
p1 <- ggplot(g1, aes(share, group, fill = status)) + geom_col(width = 0.75) + facet_grid(group_type ~ kind, scales = "free_y", space = "free_y") +
  scale_fill_manual(values = SC, name = NULL, breaks = c("new", "legacy")) + scale_y_discrete(labels = lab) + labs(x = "share of intervals flagged", y = NULL)
save2(p1, "0_03_flag_shares_by_group", 10, 9)
g2 <- stk(fr[group_type == "all" & split_type != "total"])
mk2 <- function(st, lev, xl) ggplot(g2[split_type == st][, split := factor(split, lev)], aes(split, share, fill = status)) + geom_col(width = 0.75) + facet_wrap(~ kind, nrow = 1) +
  scale_fill_manual(values = SC, name = NULL, breaks = c("new", "legacy")) + coord_cartesian(ylim = c(0, 0.45)) + labs(x = xl, y = "share of intervals flagged") + theme(axis.text.x = element_text(angle = 30, hjust = 1))
p2 <- (mk2("period", levels(iv$per), "period of the end visit") / mk2("latitude", BAND_LAB, "latitude band")) + plot_layout(guides = "collect")
save2(p2, "0_03_flag_shares_by_split", 10, 6.5)
KC <- c(`all intervals` = "black", fire = "#D95F02", pest = "#7570B3", `light harvest` = "#66A61E", `heavy harvest` = "#1B6E3C")
cpl <- copy(CP)[, `:=`(kind = factor(kind, rev(c("all intervals", KINDS))), cov = factor(c(gdd = "GDD", cmi = "CMI", soil = "soil PC1")[as.character(cov)], c("GDD", "CMI", "soil PC1")))]
p3 <- ggplot(cpl, aes(median, kind, colour = kind)) + geom_vline(xintercept = 0, linetype = 2, colour = "grey50") + geom_errorbar(aes(xmin = p05, xmax = p95), width = 0, orientation = "y") +
  geom_point(size = 2) + scale_colour_manual(values = KC, guide = "none") + facet_wrap(~ cov, nrow = 1) +
  scale_y_discrete(labels = function(x) paste0(x, " (", format(CP0$n[match(x, CP0$kind)], big.mark = ",", trim = TRUE), ")")) + labs(x = "z-score at the event (median and 5th to 95th percentile)", y = NULL)
save2(p3, "0_03_climate_position", 10, 3.4)
tp <- tq[split_type == "total" & status %in% c("new", "legacy: dated", "legacy: imputed")][, `:=`(kind = factor(kind, rev(KINDS)), status = factor(status, c("new", "legacy: dated", "legacy: imputed"), c("new", "legacy, dated", "legacy, imputed")))]
p4 <- ggplot(tp, aes(median, kind, colour = status)) + geom_errorbar(aes(xmin = p05, xmax = p95), width = 0, orientation = "y", position = position_dodge(width = 0.7)) +
  geom_point(position = position_dodge(width = 0.7), size = 2) + scale_colour_manual(values = c(new = "#08306B", `legacy, dated` = "#4292C6", `legacy, imputed` = "#9ECAE1"), name = NULL) + labs(x = "time since disturbance (years; median and 5th to 95th percentile)", y = NULL)
save2(p4, "0_03_tsd_flagged", 8, 4.2)

# ---- console summary ------------------------------------------------------------------------------------------------------------
cat("\nNetwork, all intervals: share of intervals flagged, new, legacy (legacy fraction of flags):\n")
print(fr[group_type == "all" & split_type == "total", .(kind, n_intervals, n_flagged, n_new, n_legacy, share_flagged = round(share_flagged, 3), share_new = round(share_new, 3), share_legacy = round(share_legacy, 3), legacy_fraction = round(legacy_fraction, 3))], row.names = FALSE)
cat("\nFlag combinations (all intervals, with new / legacy):\n"); print(co1[split_type == "total"][1:14], row.names = FALSE)
cat("\nCombinations of new events in the same interval (all intervals):\n"); print(co2[split_type == "total"], row.names = FALSE)
cat("\nTSD of flagged intervals (years), all intervals:\n"); print(tq[split_type == "total", .(kind, status, n, n_dated, p05, p25, median, p75, p95, max)], row.names = FALSE)
cat("\nClimate position of new events:\n"); print(CP[order(match(kind, c("all intervals", KINDS)), match(cov, c("gdd", "cmi", "soil")))][, .(kind, cov, n, est = round(est, 3), lo = round(lo, 3), hi = round(hi, 3))], row.names = FALSE)
cat("\nReturn intervals (years, dated events):\n"); print(ri, row.names = FALSE)
cat("\nSaved in", OUT, "\n")
