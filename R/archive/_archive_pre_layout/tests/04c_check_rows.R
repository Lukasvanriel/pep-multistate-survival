# ==============================================================================
# 04c_check_rows.R   (B8: does the real-data row builder lose or double count anything?)
#
# Builds the Stage-1 dataset (entry = W_j/2, midpoint exits) from the real visit
# intervals and checks it against the visit-interval table and against an
# independent per-sojourn (collapsed) likelihood. Pure checks: stops at the first
# failure, writes nothing.
#
# Input : data/processed/{visit_intervals,plot_info,transition_set}.rds
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

cat("\nB1 spec preview (placeholder continuous set: ", paste(COV_CONT, collapse = ", "), ")\n", sep = "")
m <- term_mask(tr)
cat(sprintf("  kept transitions %d | formulas with fire %d, pest %d, harv 3level %d, harv any %d\n",
            nrow(m), sum(m[, "fire_f"]), sum(m[, "pest_f"]), sum(m[, "harv_light"]), sum(m[, "harv_any"])))
print(head(rhs_formulas(m), 3))
sc <- make_scaling(iv); print(sc)
cat("\nALL CHECKS PASSED\n")
