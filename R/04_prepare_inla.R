# ==============================================================================
# 04_prepare_inla.R
#
# PART A of the INLA pipeline: data preparation, NO imputation (P1/P2 come later).
#   1. W_j per plot (stand age at the FIRST visit; used only by P1)
#   2. Visit-interval table (one row per consecutive visit pair per plot)
#   3. Transition set and per-transition covariate terms (>= MIN_EVENTS rules)
#
# Input : data/build/combined_covariates.rds   (from 03z)
#         data/raw/pep/PEP.gpkg                (cl_age; origine/perturb codes)
# Output (TABLE_DIR = data/tables/<RUN_TAG>/, set in R/config.R):
#         wj_table.rds (copied from data/build/stand_age.rds, 02b), visit_intervals.rds, plot_info.rds,
#         transition_set.rds (+ transition_set.csv)
#
# Decisions implemented here are documented in docs/PIPELINE_NOTES.md.
# ==============================================================================

library(data.table)
library(sf)
library(here)
source(here::here("R", "config.R"))   # RUN_TAG, BUILD_DIR, TABLE_DIR
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

# MIN_EVENTS (transition kept with >= this many state changes; per-term minimum for fire / pest / harvest) is set in R/config.R (B7).

# ==============================================================================
# 0. LOAD
# ==============================================================================
dt <- as.data.table(readRDS(file.path(BUILD_DIR, "combined_covariates.rds")))
dt[, no_mes := as.integer(no_mes)]
setorder(dt, id_pe, no_mes)
cat(sprintf("Combined covariates: %d visits, %d plots\n\n", nrow(dt), uniqueN(dt$id_pe)))

# ==============================================================================
# 1. W_j PER PLOT -- computed in R/02b_stand_age.R since 2026-10-26 (logic unchanged; it reproduces the W_j table of
# the inttest run exactly). Here: restricted to the plots of the combined table, with a check that 02b's first visit
# and first state are the combined table's.
# ==============================================================================
wj <- as.data.table(readRDS(file.path(BUILD_DIR, "stand_age.rds")))[id_pe %in% dt$id_pe]
chk <- dt[, .(year_first = year[1L], state_first = state[1L]), by = id_pe][wj, on = "id_pe"]
stopifnot("02b first visit / first state differ from the combined table" =
            all(chk$year_first == chk$year_1) && all(chk$state_first == chk$state_1))
cat("W_j rule breakdown:\n"); print(wj[, .N, by = .(rule, include)][order(rule)])
cat(sprintf("\nIncluded: %d / %d plots (%.2f%%)\n", sum(wj$include), nrow(wj), 100 * mean(wj$include)))
saveRDS(wj, file.path(TABLE_DIR, "wj_table.rds"))

# ==============================================================================
# 2. VISIT-INTERVAL TABLE
# One row per consecutive visit pair of an included plot. Time t_abs = W_j +
# (year - first year), so the first visit sits at W_j. Covariates of an
# interval = values at its END visit (all columns of combined_covariates);
# the three disturbance flags at the START visit are carried too (*_start),
# so the start-vs-end choice can be made later without re-running this.
# ==============================================================================
dt2 <- dt[id_pe %in% wj[include == TRUE, id_pe]]
dt2[wj, W_j := i.W_j, on = "id_pe"]
dt2[, t_abs := W_j + (year - year[1L]), by = id_pe]

# ------------------------------------------------------------------------------
# 2a. TIME SINCE DISTURBANCE (TSD) per disturbance type, in YEARS (A5: the model columns are built in covariate_spec.R from TSD_FORM)
# Since 2026-10-26 every flagged visit has a TSD from 03 (dated events: real; undated events: lower / upper / mid;
# undated at the first visit: stand age or first-visit age). TSD_SOURCE is set in R/config.R: "mid" (main),
# "lower" / "upper" (S11). The former W_j trick and the exclusion of plots with undatable disturbances are gone
# (both now in 03); a guard stops if any flagged visit is left without a TSD.
# ------------------------------------------------------------------------------
stopifnot("TSD_SOURCE in R/config.R must be mid, lower or upper" = TSD_SOURCE %in% c("mid", "lower", "upper"))
sfx <- c(lower = "_imputed", upper = "_imputed_prev", mid = "_imputed_mid")[[TSD_SOURCE]]
dt2[, `:=`(tsd_fire_y = get(paste0("tsf", sfx)), tsd_pest_y = get(paste0("sfp", sfx)),
           tsd_harv_y = get(paste0("tsh", sfx)))]
n_na <- dt2[(fire == 1L & is.na(tsd_fire_y)) | (pest == 1L & is.na(tsd_pest_y)) | (harv > 0L & is.na(tsd_harv_y)), .N]
cat(sprintf("\nTSD source '%s' (columns *%s); flagged visits without TSD: %d\n", TSD_SOURCE, sfx, n_na))
if (n_na > 0) stop("flagged visits without a TSD -- 03 should have dated every flagged visit; inspect")

cov_cols <- c(setdiff(names(dt), c("id_pe", "no_mes", "state", "year")), "tsd_fire_y", "tsd_pest_y", "tsd_harv_y")

st <- dt2[, .SD[-.N],  by = id_pe]     # start visits: all but the last of each plot
en <- dt2[, .SD[-1L],  by = id_pe]     # end visits:   all but the first of each plot
stopifnot(nrow(st) == nrow(en), identical(st$id_pe, en$id_pe))

iv <- data.table(
  id_pe = st$id_pe,
  interval = st[, seq_len(.N), by = id_pe]$V1,
  visit_start = st$no_mes, visit_end = en$no_mes,
  year_start = st$year, year_end = en$year,
  t_start = st$t_abs, t_end = en$t_abs,
  state_start = st$state, state_end = en$state,
  W_j = st$W_j
)
iv[, is_change := state_start != state_end]
iv[, `:=`(fire_start = st$fire, harv_start = st$harv, pest_start = st$pest)]
iv <- cbind(iv, en[, ..cov_cols])
# TSD in YEARS, raw (not centred, not capped: A8); structural zero where the flag is 0 (A5)
iv[, `:=`(tsd_fire_y = fifelse(fire == 1L, tsd_fire_y, 0),
          tsd_pest_y = fifelse(pest == 1L, tsd_pest_y, 0),
          tsd_harv_y = fifelse(harv >  0L, tsd_harv_y, 0))]
stopifnot(!anyNA(iv[, .(tsd_fire_y, tsd_pest_y, tsd_harv_y)]), all(iv[, .(tsd_fire_y, tsd_pest_y, tsd_harv_y)] >= 0))
cat(sprintf("TSD (years, end-visit value of flagged intervals) max: fire %.1f, pest %.1f, harvest %.1f (no cap)\n",
            max(iv$tsd_fire_y), max(iv$tsd_pest_y), max(iv$tsd_harv_y)))

# A9 (2026-11-04): date of a HEAVY disturbance inside the interval (input for C2, WP5). TSD follows TSD_SOURCE, so S11 moves the date too.
#   D_harv_heavy = year_end - tsd_harv_y  if the end visit shows heavy harvest and its date lies inside the interval
#   D_fire       = year_end - tsd_fire_y  if the end visit shows fire and its date lies inside the interval
#   D_heavy      = the earlier of the two (C2 decision)
iv[, `:=`(D_harv_heavy = fifelse(harv == 2L & tsd_harv_y < year_end - year_start, year_end - tsd_harv_y, NA_real_),
          D_fire       = fifelse(fire == 1L & tsd_fire_y < year_end - year_start, year_end - tsd_fire_y, NA_real_))]
iv[, D_heavy := pmin(D_harv_heavy, D_fire, na.rm = TRUE)]
stopifnot("A9: a date outside (year_start, year_end]" =
            all(iv[!is.na(D_heavy), D_heavy > year_start & D_heavy <= year_end]))
xr <- iv[is_change == TRUE & state_end == "regen"]
cat(sprintf("A9: X->regen change intervals %d | with a heavy-disturbance date (D_heavy) %d (heavy harvest %d, fire %d, both %d; expected 1,537 / 1,223 / 315 on the full table)\n",
            nrow(xr), sum(!is.na(xr$D_heavy)), sum(!is.na(xr$D_harv_heavy)), sum(!is.na(xr$D_fire)),
            sum(!is.na(xr$D_harv_heavy) & !is.na(xr$D_fire))))
# A9 edge case C(i), logged not solved: a heavy harvest followed by a light one inside the same interval shows only the light one at the
# end visit, so D_harv_heavy is NA. Count the intervals with harv_start != 2, light harvest at the end visit and a (non-dropped) heavy harvest
# record in data/build/disturbance_events.rds that belongs to the interval: a DATED record with year in (year_start, year_end], or an UNDATED
# record first seen at the end visit (first_no_mes == visit_end; its imputed date then lies inside the interval). User decision 2026-11-04.
hv <- as.data.table(readRDS(file.path(BUILD_DIR, "disturbance_events.rds"))$events)[type == "harv" & level == 2L & !(dropped %in% TRUE),
                                                                                    .(id_pe, yr, first_no_mes)]
ci <- merge(iv[harv_start != 2L & harv == 1L, .(id_pe, interval, visit_end, year_start, year_end, is_change, state_end)],
            hv, by = "id_pe", allow.cartesian = TRUE)
ci <- unique(ci[(!is.na(yr) & yr > year_start & yr <= year_end) | (is.na(yr) & first_no_mes == visit_end),
                .(id_pe, interval, is_change, state_end)])
cat(sprintf("A9 C(i) edge case: %d intervals with a missed heavy harvest (light at the end visit; dated heavy record inside the interval or undated heavy record first seen at the end visit); %d of them are X->regen change intervals (the exit rule would miss them)\n",
            nrow(ci), ci[is_change == TRUE & state_end == "regen", .N]))

# A1 (2026-11-04): no NA climate / soil covariate may reach INLA (it would be read as 0 = the mean).
#   iv        : end-visit values of every interval
#   first row : first visit of each plot (never an end visit; needed by P1, C1)
cat(sprintf("\nNA check: intervals gdd_3=%d cmi_gs=%d soil_pc1=%d | first visits gdd_3=%d cmi_gs=%d soil_pc1=%d\n",
            sum(is.na(iv$gdd_3)), sum(is.na(iv$cmi_gs)), sum(is.na(iv$soil_pc1)),
            sum(is.na(dt2[, .SD[1L], by = id_pe]$gdd_3)), sum(is.na(dt2[, .SD[1L], by = id_pe]$cmi_gs)),
            sum(is.na(dt2[, .SD[1L], by = id_pe]$soil_pc1))))
stopifnot("NA in gdd_3 / cmi_gs / soil_pc1 of the visit intervals (A1)" = !anyNA(iv[, .(gdd_3, cmi_gs, soil_pc1)]),
          "NA in gdd_3 / cmi_gs / soil_pc1 of a plot's first visit (A1.5)" =
            !anyNA(dt2[, .SD[1L], by = id_pe][, .(gdd_3, cmi_gs, soil_pc1)]))

stopifnot(all(iv$t_end - iv$t_start == iv$year_end - iv$year_start))
cat(sprintf("\nVisit-interval table: %d intervals, %d plots, %d with a state change\n",
            nrow(iv), uniqueN(iv$id_pe), sum(iv$is_change)))
cat("  (expected = the pre-A1 table, unchanged by A1: 11,286 plots, 36,251 intervals, 8,770 state changes)\n")
cat("  gap between visits (years):\n"); print(summary(iv$t_end - iv$t_start))
cat("  plots with no interval (single visit):",
    sum(!wj[include == TRUE, id_pe] %in% iv$id_pe), "\n")

# WP4 / C1.1: first-visit covariates of every plot, read by P1 (R/05a_fit_funs.R, prepare_fit_inputs). Raw values (the interval scaling is applied
# there). TSD in YEARS per TSD_SOURCE, 0 where the flag is 0 (same rule as the columns of `iv`). The flags and TSD at the first visit are
# never an end visit of an interval, so they are not in `iv`.
plot_info <- dt2[, .(W_j = W_j[1L], first_state = state[1L], first_year = year[1L],
                     n_visits = .N, n_changes = sum(state[-.N] != state[-1L]),
                     fv_gdd_3 = gdd_3[1L], fv_cmi_gs = cmi_gs[1L], fv_soil_pc1 = soil_pc1[1L],
                     fv_fire = fire[1L], fv_harv = harv[1L], fv_pest = pest[1L],
                     fv_tsd_fire_y = fifelse(fire[1L] == 1L, tsd_fire_y[1L], 0),
                     fv_tsd_pest_y = fifelse(pest[1L] == 1L, tsd_pest_y[1L], 0),
                     fv_tsd_harv_y = fifelse(harv[1L] >  0L, tsd_harv_y[1L], 0)), by = id_pe]
# A1.5 on what P1 reads (the WP1 check above runs on dt2; this one on the written columns) + coding checks
fvc <- grep("^fv_", names(plot_info), value = TRUE)
cat(sprintf("\nplot_info first-visit columns: NA counts %s\n", paste(sprintf("%s=%d", fvc, vapply(fvc, function(v) sum(is.na(plot_info[[v]])), 0L)), collapse = " ")))
stopifnot("A1.5: NA in a first-visit column of plot_info" = !anyNA(plot_info[, ..fvc]),
          "plot_info: fv_fire / fv_pest must be 0/1 and fv_harv 0/1/2" =
            all(plot_info$fv_fire %in% 0:1, plot_info$fv_pest %in% 0:1, plot_info$fv_harv %in% 0:2),
          "plot_info: negative first-visit TSD" = all(plot_info[, .(fv_tsd_fire_y, fv_tsd_pest_y, fv_tsd_harv_y)] >= 0),
          "plot_info: first-visit TSD must be 0 where the flag is 0" =
            all(plot_info[fv_fire == 0L, fv_tsd_fire_y] == 0, plot_info[fv_pest == 0L, fv_tsd_pest_y] == 0, plot_info[fv_harv == 0L, fv_tsd_harv_y] == 0))
cat(sprintf("plot_info: %d plots; flagged at the first visit: fire %d, pest %d, harvest light %d / heavy %d; any flag %d, none %d (plan, main build: 6,330 with none)\n",
            nrow(plot_info), sum(plot_info$fv_fire == 1L), sum(plot_info$fv_pest == 1L), sum(plot_info$fv_harv == 1L), sum(plot_info$fv_harv == 2L),
            plot_info[, sum(fv_fire == 1L | fv_pest == 1L | fv_harv > 0L)], plot_info[, sum(fv_fire == 0L & fv_pest == 0L & fv_harv == 0L)]))
saveRDS(iv, file.path(TABLE_DIR, "visit_intervals.rds"))
saveRDS(plot_info, file.path(TABLE_DIR, "plot_info.rds"))

# ==============================================================================
# 3. TRANSITION SET AND COVARIATE TERMS
# Events = state-change intervals; flag counts use the EVENT interval's END-visit
# flags (the *_start counts are reported alongside for comparison).
#   transition kept        : n_events >= MIN_EVENTS
#   fire / pest term       : >= MIN_EVENTS events with the flag = 1
#   harvest (0 / 1 light / 2 heavy), A2 (2026-11-04): "3level" if light + heavy >= MIN_EVENTS events (the no-harvest reference has
#       no minimum); else "none" (no harvest term). "3level" = light and heavy each vs none; there is no "any" level.
#       Level terms: main + GDD x + TSD for a level with events; a level with 0 events keeps its main term only (A6, covariate_spec.R).
#   n_undist (A11): events with none of the transition's disturbance terms (fire = 0 if it has a fire term, pest = 0 if a pest term,
#       harvest = 0 if a harvest term). With 0 the GDD main term is dropped (covariate_spec.R term_mask()).
# ==============================================================================
states <- sort(unique(dt$state))
ev <- iv[is_change == TRUE]
tr <- ev[, .(n_events = .N, n_plots = uniqueN(id_pe),
             n_fire = sum(fire == 1L, na.rm = TRUE), n_pest = sum(pest == 1L, na.rm = TRUE),
             n_harv0 = sum(harv == 0L, na.rm = TRUE), n_harv1 = sum(harv == 1L, na.rm = TRUE),
             n_harv2 = sum(harv == 2L, na.rm = TRUE),
             n_fire_start = sum(fire_start == 1L, na.rm = TRUE), n_pest_start = sum(pest_start == 1L, na.rm = TRUE),
             n_harv1_start = sum(harv_start == 1L, na.rm = TRUE), n_harv2_start = sum(harv_start == 2L, na.rm = TRUE)),
         by = .(from = state_start, to = state_end)]
allp <- CJ(from = states, to = states)[from != to]
tr <- merge(allp, tr, by = c("from", "to"), all.x = TRUE)
num_cols <- setdiff(names(tr), c("from", "to"))
for (cc in num_cols) set(tr, which(is.na(tr[[cc]])), cc, 0L)

tr[, kept := n_events >= MIN_EVENTS]
tr[, fire_term := kept & n_fire >= MIN_EVENTS]
tr[, pest_term := kept & n_pest >= MIN_EVENTS]
tr[, harv_term := fcase(!kept, "none", (n_harv1 + n_harv2) >= MIN_EVENTS, "3level", default = "none")]
# A11: n_undist, computed AFTER the term set is known (end-visit flags of the event intervals, as for all counts: A10)
und <- merge(ev[, .(from = state_start, to = state_end, fire, pest, harv)], tr[, .(from, to, fire_term, pest_term, harv_term)], by = c("from", "to"))
und <- und[, .(n_undist = sum(!((fire_term & fire %in% 1L) | (pest_term & pest %in% 1L) | (harv_term == "3level" & harv %in% 1:2)))),
           by = .(from, to)]
tr <- merge(tr, und, by = c("from", "to"), all.x = TRUE); tr[is.na(n_undist), n_undist := 0L]
setorder(tr, -n_events)

cat(sprintf("\nTransitions: %d possible | %d observed | %d kept (>= %d events)\n",
            nrow(tr), sum(tr$n_events > 0), sum(tr$kept), MIN_EVENTS))
cat(sprintf("Events in kept transitions: %d | in dropped: %d (%.1f%%), touching %d plots\n",
            sum(tr[kept == TRUE, n_events]), sum(tr[kept == FALSE, n_events]),
            100 * sum(tr[kept == FALSE, n_events]) / sum(tr$n_events),
            uniqueN(ev[paste(state_start, state_end) %in% tr[kept == FALSE, paste(from, to)], id_pe])))
cat("Covariate terms among kept transitions:\n")
print(tr[kept == TRUE, .(n_kept = .N, fire = sum(fire_term), pest = sum(pest_term)), ])
print(tr[kept == TRUE, .N, by = harv_term])
cat(sprintf("Transitions with a harvest term (light + heavy vs none): %d\n", tr[kept == TRUE, sum(harv_term == "3level")]))
cat("A11: kept transitions with n_undist == 0 and a disturbance term (GDD main term dropped):\n")
print(tr[kept == TRUE & n_undist == 0L & (fire_term | pest_term | harv_term == "3level"),
         .(from, to, n_events, n_undist, n_harv0, n_harv1, n_harv2, fire_term, pest_term)])

# reachability from regen through kept transitions
reach <- "regen"; repeat { nx <- unique(c(reach, tr[kept == TRUE & from %in% reach, to])); if (length(nx) == length(reach)) break; reach <- nx }
cat("States reachable from regen via kept transitions:", length(reach), "of", length(states), "\n")
if (length(reach) < length(states)) cat("  NOT reachable:", setdiff(states, reach), "\n")

saveRDS(tr, file.path(TABLE_DIR, "transition_set.rds"))
fwrite(tr, file.path(TABLE_DIR, "transition_set.csv"))
cat("\nSaved: wj_table.rds, visit_intervals.rds, plot_info.rds, transition_set.rds (+ csv)\n")

