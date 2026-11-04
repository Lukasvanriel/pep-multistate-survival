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

MIN_EVENTS <- 15L     # transition kept only with >= this many observed state changes;
                      # also the per-level minimum for fire / pest / harvest terms

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
# 2a. TIME SINCE DISTURBANCE (TSD) per disturbance type, in YEARS here (decades below)
# Since 2026-10-26 every flagged visit has a TSD from 03 (dated events: real; undated events: lower / upper / mid;
# undated at the first visit: stand age or first-visit age). TSD_SOURCE is set in R/config.R: "mid" (main),
# "lower" / "upper" (S11). The former W_j trick and the exclusion of plots with undatable disturbances are gone
# (both now in 03); a guard stops if any flagged visit is left without a TSD.
# ------------------------------------------------------------------------------
stopifnot("TSD_SOURCE in R/config.R must be mid, lower or upper" = TSD_SOURCE %in% c("mid", "lower", "upper"))
sfx <- c(lower = "_imputed", upper = "_imputed_prev", mid = "_imputed_mid")[[TSD_SOURCE]]
dt2[, `:=`(tsd_fire = get(paste0("tsf", sfx)), tsd_pest = get(paste0("sfp", sfx)),
           tsd_harv = get(paste0("tsh", sfx)))]
n_na <- dt2[(fire == 1L & is.na(tsd_fire)) | (pest == 1L & is.na(tsd_pest)) | (harv > 0L & is.na(tsd_harv)), .N]
cat(sprintf("\nTSD source '%s' (columns *%s); flagged visits without TSD: %d\n", TSD_SOURCE, sfx, n_na))
if (n_na > 0) stop("flagged visits without a TSD -- 03 should have dated every flagged visit; inspect")

cov_cols <- c(setdiff(names(dt), c("id_pe", "no_mes", "state", "year")), "tsd_fire", "tsd_pest", "tsd_harv")

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
# TSD in DECADES, raw (not centred, not capped); structural zero where the flag is 0
iv[, `:=`(tsd_fire = fifelse(fire == 1L, tsd_fire / 10, 0),
          tsd_pest = fifelse(pest == 1L, tsd_pest / 10, 0),
          tsd_harv = fifelse(harv >  0L, tsd_harv / 10, 0))]
stopifnot(!anyNA(iv[, .(tsd_fire, tsd_pest, tsd_harv)]), all(iv[, .(tsd_fire, tsd_pest, tsd_harv)] >= 0))

stopifnot(all(iv$t_end - iv$t_start == iv$year_end - iv$year_start))
cat(sprintf("\nVisit-interval table: %d intervals, %d plots, %d with a state change\n",
            nrow(iv), uniqueN(iv$id_pe), sum(iv$is_change)))
cat("  (expected: 8881 state changes over all 11,393 plots; 8880 after the 47 excluded plots; 36,508 intervals, 11,346 plots)\n")
cat("  gap between visits (years):\n"); print(summary(iv$t_end - iv$t_start))
cat("  plots with no interval (single visit):",
    sum(!wj[include == TRUE, id_pe] %in% iv$id_pe), "\n")

plot_info <- dt2[, .(W_j = W_j[1L], first_state = state[1L], first_year = year[1L],
                     n_visits = .N, n_changes = sum(state[-.N] != state[-1L])), by = id_pe]
saveRDS(iv, file.path(TABLE_DIR, "visit_intervals.rds"))
saveRDS(plot_info, file.path(TABLE_DIR, "plot_info.rds"))

# ==============================================================================
# 3. TRANSITION SET AND COVARIATE TERMS
# Events = state-change intervals; flag counts use the EVENT interval's END-visit
# flags (the *_start counts are reported alongside for comparison).
#   transition kept        : n_events >= MIN_EVENTS
#   fire / pest term       : >= MIN_EVENTS events with the flag = 1
#   harvest (0 / 1 light / 2 heavy):
#       "3level" if none, light and heavy each have >= MIN_EVENTS events;
#       else "any" (any harvest vs none) if none >= MIN_EVENTS and light+heavy >= MIN_EVENTS;
#       else "none" (harvest term dropped)
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
tr[, harv_term := fcase(
  !kept, "none",
  n_harv0 >= MIN_EVENTS & n_harv1 >= MIN_EVENTS & n_harv2 >= MIN_EVENTS, "3level",
  n_harv0 >= MIN_EVENTS & (n_harv1 + n_harv2) >= MIN_EVENTS, "any",
  default = "none")]
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

# reachability from regen through kept transitions
reach <- "regen"; repeat { nx <- unique(c(reach, tr[kept == TRUE & from %in% reach, to])); if (length(nx) == length(reach)) break; reach <- nx }
cat("States reachable from regen via kept transitions:", length(reach), "of", length(states), "\n")
if (length(reach) < length(states)) cat("  NOT reachable:", setdiff(states, reach), "\n")

saveRDS(tr, file.path(TABLE_DIR, "transition_set.rds"))
fwrite(tr, file.path(TABLE_DIR, "transition_set.csv"))
cat("\nSaved: wj_table.rds, visit_intervals.rds, plot_info.rds, transition_set.rds (+ csv)\n")

