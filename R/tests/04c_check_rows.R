# ==============================================================================
# 04c_check_rows.R   (B8: does the real-data row builder lose or double count anything?)
#
# Builds the Stage-1 dataset (entry = W_j/2, midpoint exits) from the real visit
# intervals and checks it against the visit-interval table and against an
# independent per-sojourn (collapsed) likelihood. Pure checks: stops at the first
# failure, writes nothing.
#
# Input : TABLE_DIR (data/tables/<RUN_TAG>/){visit_intervals,plot_info,transition_set}.rds
# Needs : covariate_spec.R, 04b_build_rows.R (same folder)
# ==============================================================================
suppressPackageStartupMessages(library(data.table))

source(here::here("R", "config.R"))   # TABLE_DIR = current run's 04 tables (T0)
proc_dir <- if (nzchar(Sys.getenv("PEP_PROC_DIR"))) Sys.getenv("PEP_PROC_DIR") else TABLE_DIR
r_dir    <- if (nzchar(Sys.getenv("PEP_R_DIR")))    Sys.getenv("PEP_R_DIR")    else here::here("R")
source(file.path(r_dir, "covariate_spec.R"))
source(file.path(r_dir, "04b_build_rows.R"))

iv <- readRDS(file.path(proc_dir, "visit_intervals.rds")); setDT(iv)
pi <- readRDS(file.path(proc_dir, "plot_info.rds"));       setDT(pi)
tr <- readRDS(file.path(proc_dir, "transition_set.rds"));  setDT(tr)
states <- sort(unique(c(iv$state_start, iv$state_end)))
setorder(iv, id_pe, interval)

# Stage-1 entry: W_j / 2 for every plot (the simulation's naive seed)
first_entry <- setNames(pi$W_j / 2, as.character(pi$id_pe))
t0 <- proc.time()["elapsed"]
b <- build_rows(iv, first_entry, states, tr)
cat(sprintf("Built %d rows from %d segments, %d plots (%.1f s)\n", nrow(b$rows), nrow(b$segs), b$n_plots, proc.time()["elapsed"] - t0))
rows <- b$rows; segs <- b$segs
nd <- length(states) - 1L
eq <- function(a, b) isTRUE(all.equal(a, b))
ok <- function(label, cond) { cat(sprintf("  [%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", label)); if (!isTRUE(cond)) stop("check failed: ", label) }

cat("\nCHECK 1: bookkeeping\n")
ok(sprintf("every segment has %d destination rows", nd), all(rows[, .N, by = seg_id]$N == nd))
ok("status-1 rows == state-change intervals", sum(rows$status) == sum(iv$is_change))
ok("at most one status-1 row per segment", all(rows[, sum(status), by = seg_id]$V1 <= 1L))
ok("status-1 only on event segments, one per event segment", all(rows[, sum(status), by = .(seg_id, type)][, (V1 == 1L) == (type == "event")]))
ok("event segments == state-change intervals", segs[type == "event", .N] == sum(iv$is_change))
ok("stub segments == state-change intervals (one stub after every event, also a final one)", segs[type == "stub", .N] == sum(iv$is_change))
ok("no row with from == to", all(rows$from != rows$to))
ok("realized destination equals the observed state_end", {
  e <- rows[status == 1L, .(id_pe, interval, to)]; m <- merge(e, iv[is_change == TRUE, .(id_pe, interval, state_end)], by = c("id_pe", "interval"))
  nrow(m) == nrow(e) && all(m$to == m$state_end) })
ok("kept-transition events == transition_set", sum(rows$status[rows$kept %in% TRUE]) == sum(tr[kept == TRUE, n_events]))

cat("\nCHECK 2: time bookkeeping\n")
ok("Tstop > 0 and truncation >= 0 and truncation < Tstop (non-stub rows)", all(rows$Tstop > 0) && all(rows$truncation >= 0) && all(rows[type != "stub", truncation < Tstop]))
ok("stub truncation is 0", all(rows[type == "stub", truncation] == 0))
setorder(segs, id_pe, seg_id)
ok("segments are contiguous in time within a plot (start == previous end)",
   segs[, all(abs(start_abs[-1L] - end_abs[-.N]) < 1e-9), by = id_pe][, all(V1)])
span <- segs[, .(risk = sum(end_abs - start_abs)), by = id_pe]
span <- merge(span, iv[, .(span = max(t_end) - min(t_start)), by = id_pe], by = "id_pe")
ok("total time at risk per plot == last visit - first visit (no time lost or doubled)", all(abs(span$risk - span$span) < 1e-8))
ev <- merge(segs[type == "event", .(id_pe, interval, T_abs = end_abs)], iv[, .(id_pe, interval, t_start, t_end)], by = c("id_pe", "interval"))
ok("imputed exits lie inside their visit bracket", all(ev$T_abs >= ev$t_start - 1e-9 & ev$T_abs <= ev$t_end + 1e-9))
ok("first-segment truncation == (W_j - W_j/2)/10 = W_j/20",
   { f <- merge(segs[, .SD[1L], by = id_pe], pi[, .(id_pe, W_j)], by = "id_pe"); all(abs(pmax(f$start_abs - f$E_abs, 0) / TIME_SCALE - f$W_j / 20) < 1e-9) })

cat("\nCHECK 3: covariates follow the right visit\n")
cov_chk <- merge(rows[status == 0L | status == 1L, .(seg_id, id_pe, interval, type, fire_r = fire, harv_r = harv, pest_r = pest)][, .SD[1L], by = seg_id],
                 iv[, .(id_pe, interval, fire, harv, pest)], by = c("id_pe", "interval"))
ok("obs/event segments carry their own interval's end-visit flags", cov_chk[type != "stub", all(fire_r == fire & harv_r == harv & pest_r == pest)])
stb <- merge(cov_chk[type == "stub", .(id_pe, interval_next = interval + 1L, fire_r, harv_r, pest_r)],
             iv[, .(id_pe, interval_next = interval, fire_s = fire_start, harv_s = harv_start, pest_s = pest_start)],
             by = c("id_pe", "interval_next"))
ok("stub flags == flags AT the first visit of the new state (= start flags of the next interval)",
   nrow(stb) == nrow(cov_chk[type == "stub"][iv[, .(last_int = max(interval)), by = id_pe], on = "id_pe"][interval < last_int]) && stb[, all(fire_r == fire_s & harv_r == harv_s & pest_r == pest_s)])
ok("no missing flags or continuous covariates on any row", !anyNA(rows[, c("fire", "harv", "pest", COV_CONT), with = FALSE]))

cat("\nCHECK 4: row likelihood == independent per-sojourn likelihood (fixed Weibull parameters)\n")
pr <- tr[, .(from, to)]; pr <- pr[order(from, to)]
pr[, `:=`(alpha = 0.8 + 0.2 * (.I %% 6), lam = 0.03 + 0.05 * (.I %% 5))]
rr <- merge(rows, pr, by = c("from", "to"))
rr[, ll := status * (log(alpha) + alpha * log(lam) + (alpha - 1) * log(Tstop)) - (lam * Tstop)^alpha + (lam * pmax(truncation, 0))^alpha]
# stub rows have truncation 0 -> (lam*0)^alpha = 0
ll_rows <- sum(rr$ll)
# collapsed side: one term per SOJOURN per destination, built from segs (first start / last end / event)
sj <- segs[, .(E = E_abs[1L], first_start = start_abs[1L], last_end = end_abs[.N], ev_to = event_to[.N], from = from[1L], first_type = type[1L]), by = .(id_pe, sojourn)]
sj[, `:=`(tr0 = fifelse(sojourn == 1L, pmax(first_start - E, 0), 0) / TIME_SCALE, T1 = (last_end - E) / TIME_SCALE)]
sj <- merge(sj, CJ(from = states, to = states)[from != to], by = "from", allow.cartesian = TRUE)
sj <- merge(sj, pr, by = c("from", "to"))
sj[, ll := (!is.na(ev_to) & to == ev_to) * (log(alpha) + alpha * log(lam) + (alpha - 1) * log(T1)) - (lam * T1)^alpha + (lam * tr0)^alpha]
ok(sprintf("log-lik rows %.8f vs collapsed %.8f (diff %.2e)", ll_rows, sum(sj$ll), ll_rows - sum(sj$ll)), abs(ll_rows - sum(sj$ll)) < 1e-6 * max(1, abs(ll_rows)))
ok("number of sojourns == changes + 1 summed over plots", nrow(sj[, .N, by = .(id_pe, sojourn)]) == sum(iv$is_change) + uniqueN(iv$id_pe))

cat("\nCHECK 5 (B7): recount of events on the final rows\n")
rc <- recount_events(rows, tr)
cat(sprintf("  kept transitions %d, events %d, all counts equal transition_set and all thresholds pass\n", nrow(rc), sum(rc$n_events)))

cat("\nCHECK 6 (WP2: A2-A11, B7): term set, TSD columns, heavy-disturbance dates\n")
m <- term_mask(tr)
kt <- tr[kept == TRUE]; lab <- rownames(m)
has_dist <- kt$fire_term | kt$pest_term | kt$harv_term == "3level"
ok("B7: switches defined in config.R (WEIBULL_VARIANT, COEF_PRIOR_SD, TSD_FORM, MIN_EVENTS, XREGEN_EXIT_RULE)",
   all(vapply(c("WEIBULL_VARIANT", "COEF_PRIOR_SD", "TSD_FORM", "MIN_EVENTS", "XREGEN_EXIT_RULE"), exists, logical(1))) && TSD_FORM %in% c("log1p", "linear"))
ok("A2: harv_term is only 3level or none; 3level <=> kept and light + heavy >= MIN_EVENTS (no 'any')",
   all(tr$harv_term %in% c("3level", "none")) && all((tr$harv_term == "3level") == (tr$kept & tr$n_harv1 + tr$n_harv2 >= MIN_EVENTS)))
ok("A2: no harv_any / gdd_x_harv_any / tsd_harv term anywhere in the spec", !any(c("harv_any", "gdd_x_harv_any", "tsd_harv") %in% all_terms()) && !any(grepl("harv_any", colnames(m))))
ok("A3: fire / pest terms <=> kept and >= MIN_EVENTS events with the flag", all(tr$fire_term == (tr$kept & tr$n_fire >= MIN_EVENTS)) && all(tr$pest_term == (tr$kept & tr$n_pest >= MIN_EVENTS)))
ok("A6: harvest level terms: main for 3level; GDD x and TSD only with >= 1 event of the level",
   all(m[, "harv_light"] == (kt$harv_term == "3level")) && all(m[, "harv_heavy"] == (kt$harv_term == "3level")) &&
   all(m[, "gdd_x_harv_light"] == (kt$harv_term == "3level" & kt$n_harv1 > 0)) && all(m[, "tsd_harv_light"] == m[, "gdd_x_harv_light"]) &&
   all(m[, "gdd_x_harv_heavy"] == (kt$harv_term == "3level" & kt$n_harv2 > 0)) && all(m[, "tsd_harv_heavy"] == m[, "gdd_x_harv_heavy"]))
cat("  A6: 3level transitions with a 0-event level (main term only):\n")
print(kt[harv_term == "3level" & (n_harv1 == 0L | n_harv2 == 0L), .(from, to, n_harv0, n_harv1, n_harv2)])
if ("other_temperate->yellow_birch" %in% lab) {
  r <- m["other_temperate->yellow_birch", ]
  ok("A7: other_temperate->yellow_birch has harv_term 3level (heavy main term present)", kt[from == "other_temperate" & to == "yellow_birch", harv_term] == "3level" && r["harv_heavy"] && r["harv_light"])
  cat("  A6/A7 terms of other_temperate->yellow_birch: ", paste(colnames(m)[r], collapse = " + "), "\n", sep = "")
}
ok("A11: z_gdd_3 absent <=> n_undist == 0 and >= 1 disturbance term; every other base term present",
   all(m[, "z_gdd_3"] == !(kt$n_undist == 0L & has_dist)) && all(m[, c("z_cmi_gs", "z_soil_pc1")]))
cat("  A11: transitions without z_gdd_3: ", paste(lab[!m[, "z_gdd_3"]], collapse = ", "), "\n", sep = "")
if ("regen->sugar_maple" %in% lab) ok("A11: regen->sugar_maple has no z_gdd_3 but keeps gdd_x_harv_light / heavy", !m["regen->sugar_maple", "z_gdd_3"] && all(m["regen->sugar_maple", c("gdd_x_harv_light", "gdd_x_harv_heavy")]))
ok("A11: n_undist <= n_events for every transition", all(tr$n_undist <= tr$n_events))

sc <- make_scaling(iv)
dm <- add_model_covariates(iv, sc); dl <- add_model_covariates(iv, sc, tsd_form = "linear")
ok("A4: at most one of tsd_harv_light / tsd_harv_heavy is non-zero per interval", all(dm$tsd_harv_light == 0 | dm$tsd_harv_heavy == 0))
ok("A4: each harvest TSD column is 0 wherever its level flag is 0", all(dm$tsd_harv_light[dm$harv_light == 0L] == 0) && all(dm$tsd_harv_heavy[dm$harv_heavy == 0L] == 0))
ok("A5: model TSD columns = flag * f(TSD years) for TSD_FORM (log1p main)",
   eq(dm$tsd_fire, dm$fire_f * log1p(dm$tsd_fire_y)) && eq(dm$tsd_harv_heavy, dm$harv_heavy * log1p(dm$tsd_harv_y)) && eq(dm$tsd_pest, dm$pest_f * log1p(dm$tsd_pest_y)))
ok("A5: TSD_FORM linear gives TSD in decades (S24)", eq(dl$tsd_harv_light, dl$harv_light * dl$tsd_harv_y / 10) && eq(dl$tsd_fire, dl$fire_f * dl$tsd_fire_y / 10))
cat(sprintf("  A5: max TSD years fire %.1f, pest %.1f, harvest %.1f (A8: no cap) | max model column (log1p) %.2f (plan: log(176) ~ 5.2 on the full table)\n",
            max(iv$tsd_fire_y), max(iv$tsd_pest_y), max(iv$tsd_harv_y), max(dm[, .(tsd_fire, tsd_pest, tsd_harv_light, tsd_harv_heavy)])))

ok("A9: D_harv_heavy / D_fire / D_heavy exist; D_heavy is the earlier of the two", all(c("D_harv_heavy", "D_fire", "D_heavy") %in% names(iv)) &&
   identical(is.na(iv$D_heavy), is.na(iv$D_harv_heavy) & is.na(iv$D_fire)) && all(iv[!is.na(D_heavy), D_heavy == pmin(D_harv_heavy, D_fire, na.rm = TRUE)]))
ok("A9: D_harv_heavy only where the end visit shows heavy harvest, D_fire only where it shows fire", all(iv[!is.na(D_harv_heavy), harv == 2L]) && all(iv[!is.na(D_fire), fire == 1L]))
ok("A9: every date lies in (year_start, year_end]", all(iv[!is.na(D_heavy), D_heavy > year_start & D_heavy <= year_end]))
xr <- iv[is_change == TRUE & state_end == "regen"]
cat(sprintf("  A9: X->regen change intervals %d, with D_heavy %d (heavy harvest %d, fire %d, both %d) -- full table: 1,537 (1,223 / 315)\n",
            nrow(xr), sum(!is.na(xr$D_heavy)), sum(!is.na(xr$D_harv_heavy)), sum(!is.na(xr$D_fire)), sum(!is.na(xr$D_harv_heavy) & !is.na(xr$D_fire))))

cat(sprintf("\nB1 spec preview (placeholder continuous set: %s)\n", paste(COV_CONT, collapse = ", ")))
cat(sprintf("  kept transitions %d | formulas with fire %d, pest %d, harvest light+heavy %d\n",
            nrow(m), sum(m[, "fire_f"]), sum(m[, "pest_f"]), sum(m[, "harv_light"])))
print(head(rhs_formulas(m), 3))
print(sc)
cat("\nALL CHECKS PASSED\n")
