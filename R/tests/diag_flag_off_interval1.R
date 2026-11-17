# ==============================================================================
# R/tests/diag_flag_off_interval1.R -- READ-ONLY diagnostic (2026-11-05, "Pipeline implementations" project).
# Question: why is a disturbance flag that is ON at a plot's first visit OFF at the end visit of interval 1?
#   (wp4_C1_check.R real: fire 20 / pest 72 / harvest 29 = 121 plots)
# Edits NO pipeline file and NO table. Reads  data/tables/final/{plot_info,visit_intervals}.rds,
#   data/build/{disturbance_events,disturbance_covariates,stand_age}.rds, R/03a_disturbance_events.R (sourced for its functions),
#   and the raw PEP layers pee_ori_sond / station_pe / placette_mes as CSV (data-output/diagnostics/flag_off_raw/, exported from
#   data/raw/pep/PEP.gpkg by R/tests/diag_flag_off_extract_raw.py).
# Writes only to data-output/diagnostics/flag_off/ .
# Usage: Rscript R/tests/diag_flag_off_interval1.R   (from the project root; PEP_ROOT overrides the root)
# ==============================================================================
suppressMessages(library(data.table))
ROOT <- Sys.getenv("PEP_ROOT", getwd()); OUT <- file.path(ROOT, "data-output", "diagnostics", "flag_off"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "R", "03a_disturbance_events.R"))
rd <- function(p) as.data.table(readRDS(file.path(ROOT, p)))
pi <- rd("data/tables/final/plot_info.rds"); iv <- rd("data/tables/final/visit_intervals.rds")
de <- readRDS(file.path(ROOT, "data/build/disturbance_events.rds")); dc <- rd("data/build/disturbance_covariates.rds"); sa <- rd("data/build/stand_age.rds")
V <- as.data.table(de$visits); setorder(V, id_pe, no_mes)
CODES <- list(fire_o = c("BR","CRB"), fire_p = "BRP",   # copied from R/03_build_covariates.R (lists unchanged)
  harv_o = c("CT","CPR","CBA","CBT","CEF","CRS","CPH","CPT","CS","CDV","CPP"),
  harv_p = c("CP","EPC","CJ","CDL","EC","CEA","CAM","CPI","CPF","CJG","CJP","CJT","CPS","CB","CTR","DLD","CPM","CPC","CIP","CA","CPX","ESI","ECE","CE","CPR"),
  pest_o = "ES", pest_p = c("EL","CE"))
FV <- first_visit_partial_ages(V, CODES, 0.5)
ci <- function(t, fn) if (nzchar(t)) cat(t, "\n")

# ---- 1. the 121 ------------------------------------------------------------------------------------------------------------------------------
i1 <- iv[interval == 1L]
m <- merge(pi, i1[, .(id_pe, no1 = visit_start, no2 = visit_end, state_start, state_end, is_change, fire2 = fire, harv2 = harv, pest2 = pest,
                      tsd_fire2 = tsd_fire_y, tsd_pest2 = tsd_pest_y, tsd_harv2 = tsd_harv_y, D_heavy)], by = "id_pe")
stopifnot(nrow(m) == nrow(pi))
cases <- rbind(m[fv_fire == 1L & fire2 == 0L, .(id_pe, T = "fire")], m[fv_pest == 1L & pest2 == 0L, .(id_pe, T = "pest")],
               m[fv_harv  > 0L & harv2 == 0L, .(id_pe, T = "harv")])
cat(sprintf("flagged at first visit: fire %d, pest %d, harvest %d | OFF at the interval-1 end visit: fire %d, pest %d, harvest %d, total %d (distinct plots %d)\n",
            sum(m$fv_fire == 1L), sum(m$fv_pest == 1L), sum(m$fv_harv > 0L), sum(cases$T == "fire"), sum(cases$T == "pest"), sum(cases$T == "harv"), nrow(cases), uniqueN(cases$id_pe)))
cat(sprintf("(not counted as off: harvest heavy -> light between the visits: %d plots)\n", m[fv_harv == 2L & harv2 == 1L, .N]))

# ---- 2. replay of 03a for each case ----------------------------------------------------------------------------------------------------------
raw_pee <- fread(file.path(ROOT, "data-output/diagnostics/flag_off_raw/pee_ori_sond.csv"), select = c("id_pe","no_mes","origine","an_origine","perturb","an_perturb"), colClasses = list(character = c("id_pe","origine","an_origine","perturb","an_perturb")))
raw_sta <- fread(file.path(ROOT, "data-output/diagnostics/flag_off_raw/station_pe.csv"), select = c("id_pe","no_mes","origine","perturb"), colClasses = list(character = c("id_pe","origine","perturb")))
raw_pm  <- fread(file.path(ROOT, "data-output/diagnostics/flag_off_raw/placette_mes.csv"), select = c("id_pe","no_mes","date_sond","statut_mes"), colClasses = list(character = c("id_pe","date_sond","statut_mes")))
codestr <- function(o, ao, p, ap, so, sp) paste0("pee[", ifelse(is.na(o), "-", paste0(o, ifelse(is.na(ao), "", paste0("/", ao)))), " | ", ifelse(is.na(p), "-", paste0(p, ifelse(is.na(ap), "", paste0("/", ap)))), "] sta[", ifelse(is.na(so), "-", so), " | ", ifelse(is.na(sp), "-", sp), "]")
lay <- function(g, v, ty, lv) { # which layer carries the record of type ty / level lv at visit v
  oc <- CODES[[paste0(ty, "_o")]]; pc <- CODES[[paste0(ty, "_p")]]
  inpee <- if (lv == 2L) isTRUE(g$origine[v] %in% oc) else isTRUE(g$perturb[v] %in% pc)
  insta <- if (lv == 2L) isTRUE(g$s_origine[v] %in% oc) else isTRUE(g$s_perturb[v] %in% pc)
  if (inpee && insta) "pee+station" else if (inpee) "pee" else if (insta) "station_only" else "none"
}
tsd_mid <- function(b, T, v) list(fire = b$tsf_imputed_mid, pest = b$sfp_imputed_mid, harv = b$tsh_imputed_mid)[[T]][v]
rows <- list(); trace <- list(); mism_raw <- 0L
for (i in seq_len(nrow(cases))) {
  id <- cases$id_pe[i]; T <- cases$T[i]
  g <- V[id_pe == id]; b <- dc[id_pe == id][match(g$no_mes, no_mes)]; w <- sa[id_pe == id]; W1 <- w$W_j + (g$year[1] - w$year_1); year <- g$year
  # raw layers vs the `visits` table saved by 03 (after the ignore codes): codes must be identical or an ignore code (nulled) in raw
  rp <- raw_pee[id_pe == id][match(g$no_mes, no_mes)]; rs <- raw_sta[id_pe == id][match(g$no_mes, no_mes)]
  mism_raw <- mism_raw + sum(!(is.na(g$origine) | g$origine == rp$origine), !(is.na(g$perturb) | g$perturb == rp$perturb), !(is.na(g$s_origine) | g$s_origine == rs$origine), !(is.na(g$s_perturb) | g$s_perturb == rs$perturb))
  ba_ev <- ba_event_flag(b$total_ba, 5, .25)
  r0 <- make_records(year, g$origine, g$an_origine, g$perturb, g$an_perturb, CODES); nn <- rep(NA_integer_, length(year))
  r0 <- merge_station_records(r0, make_records(year, g$s_origine, nn, g$s_perturb, nn, CODES))$r
  r <- collapse_same_visit(backfill_years(collapse_same_visit(r0, year), year), year)
  ev <- build_events(r, year, ba_ev); ea <- event_activity(ev, year, ba_ev); vv <- visit_values(ev, ea, year, W1, FV)
  stopifnot("replay differs from the stored flags" = identical(as.integer(vv[[T]]$flag), as.integer(switch(T, fire = b$fire, harv = b$harv, pest = b$pest))))
  e <- vv[[T]]$event[1L]; stopifnot(!is.na(e), vv[[T]]$flag[2L] == 0L)
  ks <- which(ev$level == 2L & is.finite(ea$start) & !ea$dropped & seq_along(ev$type) != e & vapply(seq_along(ev$type), function(k) known_older(ea$range[[e]], ea$range[[k]]), logical(1)))
  k <- ks[order(ea$start[ks])][1L]; stopifnot(!is.na(k), ea$start[k] <= 2)
  # (k_lay is computed below; for the "partial raised to stand-origin" case lay() finds no stand-origin code at the first visit of the event)
  rk <- ev$records[ev$records$event == k, , drop = FALSE]
  e_src1 <- lay(g, ev$first[e], ev$type[e], ev$level[e]); k_lay <- lay(g, ev$first[k], ev$type[k], 2L)
  k_ysrc <- if (nrow(rk)) paste(unique(rk$yr_src), collapse = "+") else NA_character_
  ydiff <- if (!is.na(ev$yr[k])) ev$yr[k] - year[1L] else NA_integer_
  kf <- ev$first[k]; ko <- CODES[[paste0(ev$type[k], "_o")]]; kp <- CODES[[paste0(ev$type[k], "_p")]]
  year_origin <- if (is.na(ev$yr[k])) "undated" else if (isTRUE(g$origine[kf] %in% ko) && isTRUE(g$an_origine[kf] == ev$yr[k])) "own year (an_origine at the visit)"
    else if (isTRUE(g$perturb[kf] %in% kp) && isTRUE(g$an_perturb[kf] == ev$yr[k])) "year of a same-visit partial code of that type (pair merge, 03a step 4)"
    else "year borrowed from a LATER visit's dated record (back-fill, 03a step 3)"
  pos <- if (kf == 1L) "recorded at visit 1, start deferred to visit 2 (year >= visit-1 year)" else if (kf > 2L) "first recorded at visit >= 3, dated back into interval 1" else "recorded at visit 2"
  if (k_lay == "none") k_lay <- "partial code only at that visit (no stand-origin code), raised to stand-origin by a later dated record (03a step 4, A2)"
  cause <- paste0(pos, " | ", k_lay, " | ", year_origin)
  cause_grp <- if (grepl("^partial code only", k_lay)) "E ender = partial code at visit 2 raised to stand-origin by a later dated record (A2)" else if (kf != 2L) "D ender recorded at visit 1 or >= 3" else if (year_origin == "undated") "B ender undated at visit 2 (order = visit position only)"
    else if (grepl("^own", year_origin)) "A ender dated by PEP at visit 2 (own year in interval 1)" else "C ender year borrowed by the build (back-fill / pair merge)"
  tc <- function(v) { o <- CODES[[paste0(T, "_o")]]; pp <- CODES[[paste0(T, "_p")]]; isTRUE(g$origine[v] %in% o) || isTRUE(g$perturb[v] %in% pp) || isTRUE(g$s_origine[v] %in% o) || isTRUE(g$s_perturb[v] %in% pp) }
  same_type_code_v2 <- tc(2L)   # the SAME disturbance type is coded again at visit 2 (pee or station), yet the flag is 0
  same_type_dropped_v2 <- any(ev$type == T & ev$first == 2L & ea$dropped)
  flag_later_on <- any(vv[[T]]$flag[-(1:2)] > 0L)
  e_kind <- if (!is.na(ev$yr[e])) "e dated by PEP (year given)" else if (ev$first[e] > 1L) "e undated, first seen after visit 1" else if (ev$level[e] == 2L || T == "fire") "e undated at visit 1: age = stand age W1" else "e undated at visit 1: age = fv_age (median of dated first-visit partial codes)"
  rows[[i]] <- data.table(type = T, id_pe = id, no_mes1 = g$no_mes[1], no_mes2 = g$no_mes[2], year1 = year[1], year2 = year[2], gap_y = year[2] - year[1],
    fire1 = b$fire[1], harv1 = b$harv[1], pest1 = b$pest[1], fire2 = b$fire[2], harv2 = b$harv[2], pest2 = b$pest[2],
    tsd_mid_v1 = tsd_mid(b, T, 1L), tsd_mid_v2 = tsd_mid(b, T, 2L), tsd_in_iv_v2 = m[id_pe == id, c(fire = tsd_fire2, pest = tsd_pest2, harv = tsd_harv2)[[T]]],
    state_start = m[id_pe == id, state_start], state_end = m[id_pe == id, state_end], is_change = m[id_pe == id, is_change],
    e_level = ev$level[e], e_yr = ev$yr[e], e_first_visit = ev$first[e], e_layer_v1 = e_src1, e_kind = e_kind, e_dating_v1 = list(fire = b$fire_dating, harv = b$harv_dating, pest = b$pest_dating)[[T]][1L],
    k_type = ev$type[k], k_level = ev$level[k], k_yr = ev$yr[k], k_yr_minus_year1 = ydiff, k_first_visit = ev$first[k], k_start_visit = ea$start[k], k_layer = k_lay, k_yr_src = k_ysrc, n_enders = length(ks), same_type_code_at_v2 = same_type_code_v2, same_type_event_dropped_by_selection = same_type_dropped_v2, flag_on_again_later = flag_later_on, cause_grp = cause_grp, k_year_origin = year_origin,
    ba1 = b$total_ba[1], ba2 = b$total_ba[2], ba_event2 = ba_ev[2],
    codes_v1 = codestr(g$origine[1], g$an_origine[1], g$perturb[1], g$an_perturb[1], g$s_origine[1], g$s_perturb[1]),
    codes_v2 = codestr(g$origine[2], g$an_origine[2], g$perturb[2], g$an_perturb[2], g$s_origine[2], g$s_perturb[2]), cause = cause)
  trace[[i]] <- data.table(type = T, id_pe = id, cause = cause, visit = g$no_mes, year = year, raw_pee_origine = rp$origine, raw_pee_an_origine = rp$an_origine, raw_pee_perturb = rp$perturb, raw_pee_an_perturb = rp$an_perturb,
    raw_sta_origine = rs$origine, raw_sta_perturb = rs$perturb, total_ba = b$total_ba, ba_event = ba_ev, fire = b$fire, harv = b$harv, pest = b$pest,
    tsf_mid = b$tsf_imputed_mid, tsh_mid = b$tsh_imputed_mid, sfp_mid = b$sfp_imputed_mid)
}
R <- rbindlist(rows); TR <- rbindlist(trace)
cat(sprintf("\nraw layers vs the 03 `visits` table at the visits of the 121 plots: code mismatches (other than ignore codes nulled by 03 step 3): %d\n", mism_raw))
cat("All 121 replays reproduce the stored flags (stopifnot above); in all 121 an ender was found: a live stand-origin event of ANOTHER type, known newer, started at visit 2.\n")
stopifnot(all(R$k_type != R$type))

# ---- 3. counts ---------------------------------------------------------------------------------------------------------------------------------
cat("\n== Cause group (proximate trigger = the other-type stand-origin event that ends the old event), per type ==\n"); print(dcast(R[, .N, by = .(cause_grp, type)], cause_grp ~ type, value.var = "N", fill = 0))
cat("\n== Detail: position | layer | origin of the ender's year ==\n"); print(dcast(R[, .N, by = .(cause, type)], cause ~ type, value.var = "N", fill = 0), nrows = 100)
cat("\n== Old event recorded where? (e_first_visit) ==\n"); print(R[, .N, by = .(type, e_first_visit_is_1 = e_first_visit == 1L)])
cat("\n== Ender order evidence: old event dated by PEP and ender dated by PEP? ==\n"); print(R[, .N, by = .(type, e_dated = !is.na(e_yr), k_dated = !is.na(k_yr))][order(type)])
cat("\n== Type of the old (switched-off) event x type of the ender ==\n"); print(R[, .N, by = .(type, ender = k_type)][order(type, -N)])
cat("\n== How the old event was dated at visit 1 ==\n"); print(dcast(R[, .N, by = .(e_kind, type)], e_kind ~ type, value.var = "N", fill = 0))
cat("\n== Layer of the old event's code at visit 1 ==\n"); print(dcast(R[, .N, by = .(e_layer_v1, type)], e_layer_v1 ~ type, value.var = "N", fill = 0))
cat("\n== Ender year relative to visit 1 (k_yr - year1; NA = undated) ==\n"); print(R[, .N, by = .(type, k_yr_minus_year1)][order(type, k_yr_minus_year1)])
cat("\n== SAME disturbance type coded again at visit 2 (pee or station) although the flag is 0 at visit 2 ==\n"); print(R[, .N, by = .(type, same_type_code_at_v2, dropped_by_selection = same_type_event_dropped_by_selection)][order(type)])
cat("\n== Flag of the old type on again at some later visit (>= 3) ==\n"); print(R[, .N, by = .(type, flag_on_again_later)][order(type)])
cat("\nn_enders > 1 in", R[n_enders > 1, .N], "cases\n")
R[, code1 := sub("\\].*", "]", codes_v1)]
cat("\n== Codes at visit 1 (old event) / visit 2 (ender), most frequent ==\n"); print(R[, .N, by = .(type, codes_v1, codes_v2)][order(-N)][1:25])

# the same machinery on ALL flagged-at-first-visit plots: how many have a stand-origin event of ANOTHER type starting at visit 2 and still keep the flag?
chk <- list()
for (id in pi[fv_fire == 1L | fv_pest == 1L | fv_harv > 0L, id_pe]) {
  g <- V[id_pe == id]; b <- dc[id_pe == id][match(g$no_mes, no_mes)]; w <- sa[id_pe == id]; W1 <- w$W_j + (g$year[1] - w$year_1); year <- g$year
  ba_ev <- ba_event_flag(b$total_ba, 5, .25)
  r0 <- make_records(year, g$origine, g$an_origine, g$perturb, g$an_perturb, CODES); nn <- rep(NA_integer_, length(year))
  r0 <- merge_station_records(r0, make_records(year, g$s_origine, nn, g$s_perturb, nn, CODES))$r
  r <- collapse_same_visit(backfill_years(collapse_same_visit(r0, year), year), year)
  ev <- build_events(r, year, ba_ev); ea <- event_activity(ev, year, ba_ev); vv <- visit_values(ev, ea, year, W1, FV)
  for (T in c("fire","pest","harv")) if (vv[[T]]$flag[1L] > 0L) {
    e <- vv[[T]]$event[1L]
    ks <- which(ev$level == 2L & ev$type != T & is.finite(ea$start) & !ea$dropped & ea$start == 2L)
    if (length(ks)) chk[[length(chk) + 1L]] <- data.table(id_pe = id, type = T, n_other_so_start2 = length(ks), kept = vv[[T]]$flag[2L] > 0L,
      known_older = any(vapply(ks, function(k) known_older(ea$range[[e]], ea$range[[k]]), logical(1))),
      k_yr_eq_year1 = any(!is.na(ev$yr[ks]) & ev$yr[ks] == year[1L]), e_yr_eq_k_yr = any(!is.na(ev$yr[e]) & !is.na(ev$yr[ks]) & ev$yr[ks] == ev$yr[e]))
  }
}
CH <- rbindlist(chk)
cat("\n== ALL flagged-at-first-visit (plot, type) with another-type stand-origin event starting at visit 2: switched off vs kept ==\n"); print(CH[, .N, by = .(type, flag_at_v2 = fifelse(kept, "kept ON", "OFF"), order_known = known_older)][order(type, flag_at_v2)])

# ---- 4. transitions of interval 1 ---------------------------------------------------------------------------------------------------------------
R[, `:=`(heavy_ender_end_visit = (k_type == "harv" & harv2 == 2L) | (k_type == "fire" & fire2 == 1L))]
tr1 <- R[, .(n = .N, n_state_change = sum(is_change)), by = .(type, state_start, state_end)][order(type, -n)]
cat("\n== Interval-1 rows of the 121, by state_start -> state_end ==\n"); print(tr1[n >= 3])
cat(sprintf("\nstate changes among the 121 interval-1 rows: %d of %d (fire %d/%d, pest %d/%d, harvest %d/%d)\n", R[is_change == TRUE, .N], nrow(R),
    R[type == "fire" & is_change, .N], R[type == "fire", .N], R[type == "pest" & is_change, .N], R[type == "pest", .N], R[type == "harv" & is_change, .N], R[type == "harv", .N]))
cat(sprintf("state changes into regen: %d | with D_heavy set (exit rule): %d\n", R[is_change & state_end == "regen", .N], m[id_pe %in% R[is_change & state_end == "regen", id_pe] & !is.na(D_heavy), .N]))
cat(sprintf("intervals with ALL of fire / pest / harvest = 0 at the end visit: %d\n", R[fire2 == 0L & pest2 == 0L & harv2 == 0L, .N]))
fwrite(R, file.path(OUT, "flag_off_cases_121.csv")); fwrite(TR, file.path(OUT, "flag_off_raw_trace_121.csv")); fwrite(tr1, file.path(OUT, "flag_off_by_transition.csv")); fwrite(CH, file.path(OUT, "flag_off_vs_kept_all_flagged.csv"))
fwrite(dcast(R[, .N, by = .(cause_grp, cause, type)], cause_grp + cause ~ type, value.var = "N", fill = 0), file.path(OUT, "flag_off_cause_summary.csv"))
# ---- 5. worked traces: up to 5 examples per type, one per distinct cause where possible ------------------------------------------------------
sink(file.path(OUT, "flag_off_examples_trace.txt"))
for (T in c("fire", "pest", "harv")) {
  pick <- R[type == T][order(cause, id_pe)][, .SD[1L], by = cause][1:min(5L, .N)]
  for (id in pick$id_pe) {
    z <- R[id_pe == id & type == T]; cat(sprintf("\n######## %s | plot %s | %s\n", T, id, z$cause))
    cat(sprintf("interval 1: visit %d (%d) -> visit %d (%d), state %s -> %s | old flag %s TSD(mid) %.1f at visit 1; flag 0 at visit 2 | ender: %s %s level 2, year %s, first recorded at visit %d, start visit %d\n",
        z$no_mes1, z$year1, z$no_mes2, z$year2, z$state_start, z$state_end, T, z$tsd_mid_v1, z$k_type, ifelse(z$k_layer == "none", "?", z$k_layer), z$k_yr, z$k_first_visit, z$k_start_visit))
    print(TR[id_pe == id & type == T, !c("type", "id_pe", "cause")])
    g <- V[id_pe == id]; b <- dc[id_pe == id][match(g$no_mes, no_mes)]; w <- sa[id_pe == id]; year <- g$year; ba_ev <- ba_event_flag(b$total_ba, 5, .25)
    r0 <- make_records(year, g$origine, g$an_origine, g$perturb, g$an_perturb, CODES); nn <- rep(NA_integer_, length(year))
    r0 <- merge_station_records(r0, make_records(year, g$s_origine, nn, g$s_perturb, nn, CODES))$r
    r <- collapse_same_visit(backfill_years(collapse_same_visit(r0, year), year), year); ev <- build_events(r, year, ba_ev); ea <- event_activity(ev, year, ba_ev)
    cat("events (03a replay):\n"); print(data.table(event = seq_along(ev$type), type = ev$type, level = ev$level, yr = ev$yr, first_visit = ev$first, start_visit = ea$start, dropped = ea$dropped,
        range_lo = sapply(ea$range, `[`, 1), range_hi = sapply(ea$range, `[`, 2)))
    cat("records (visit, type, level, yr, yr_src, event, action):\n"); print(ev$records[, intersect(c("v","type","level","yr","yr_src","event","action"), names(ev$records))], row.names = FALSE)
  }
}
sink()
saveRDS(list(R = R, TR = TR), file.path("/tmp/w", "flag_off_RTR.rds"))
