# ==============================================================================
# 04_prepare_inla.R
#
# PART A of the INLA pipeline: data preparation, NO imputation (P1/P2 come later).
#   1. W_j per plot (stand age at the FIRST visit; used only by P1)
#   2. Visit-interval table (one row per consecutive visit pair per plot)
#   3. Transition set and per-transition covariate terms (>= MIN_EVENTS rules)
#
# Input : data/processed/combined_covariates.rds   (from 03z)
#         data/raw/pep/PEP.gpkg                    (cl_age; origine/perturb codes)
# Output: data/processed/wj_table.rds
#         data/processed/visit_intervals.rds
#         data/processed/plot_info.rds
#         data/processed/transition_set.rds   (+ data-output/data/transition_set.csv)
#
# Decisions implemented here are documented in docs/PIPELINE_NOTES.md.
# ==============================================================================

library(data.table)
library(sf)
library(here)
source(here::here("R", "config.R"))   # RUN_TAG, TABLE_DIR, FIT_DIR (T0)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

MIN_EVENTS <- 15L     # transition kept only with >= this many observed state changes;
                      # also the per-level minimum for fire / pest / harvest terms

# ==============================================================================
# 0. LOAD
# ==============================================================================
dt <- as.data.table(readRDS(here("data", "processed", "combined_covariates.rds")))
dt[, no_mes := as.integer(no_mes)]
setorder(dt, id_pe, no_mes)
cat(sprintf("Combined covariates: %d visits, %d plots\n\n", nrow(dt), uniqueN(dt$id_pe)))

# ==============================================================================
# 1. W_j PER PLOT
# cl_age source: pee_ori_sond first, gaps filled from station_pe (per visit).
# Rules: 1 first visit | 2 later visit, back-calculated (never 0 unless first state is regen) | 3 time since
# fire/harvest at the first visit | 4 first state regen AND later exits regen
# -> 0 | otherwise EXCLUDED.
# ==============================================================================
PEP_PATH <- here("data", "raw", "pep", "PEP.gpkg")
pee <- as.data.table(st_read(PEP_PATH, layer = "pee_ori_sond", quiet = TRUE))
sta <- as.data.table(st_read(PEP_PATH, layer = "station_pe",   quiet = TRUE))
pee[, no_mes := as.integer(no_mes)]
sta[, no_mes := as.integer(no_mes)]

sw <- function(x, p) !is.na(x) & startsWith(x, p)
cl_age_to_num <- function(code) {
  code <- as.character(code)
  fcase(
    sw(code, "VIN") | sw(code, "VIR"), 100,
    sw(code, "JIN") | sw(code, "JIR"), 50,
    sw(code, "120"), 120, sw(code, "90"), 90, sw(code, "70"), 70,
    sw(code, "50"),  50,  sw(code, "30"), 30, sw(code, "10"), 10,
    default = NA_real_
  )
}

src <- dt[, .(id_pe, no_mes, year, state)]
src[pee, `:=`(cl_pee = i.cl_age, origine = i.origine, an_origine = i.an_origine,
              perturb = i.perturb, an_perturb = i.an_perturb), on = c("id_pe", "no_mes")]
src[sta, cl_sta := i.cl_age, on = c("id_pe", "no_mes")]
src[, `:=`(num_pee = cl_age_to_num(cl_pee), num_sta = cl_age_to_num(cl_sta))]
src[, cl_num := fifelse(!is.na(num_pee), num_pee, num_sta)]

cat("cl_age availability over the visits in the dataset:\n")
cat(sprintf("  pee_ori_sond numeric: %d / %d | station_pe numeric: %d | pee then station: %d\n\n",
            sum(!is.na(src$num_pee)), nrow(src), sum(!is.na(src$num_sta)), sum(!is.na(src$cl_num))))

first <- src[, .SD[1L], by = id_pe]
wj <- first[, .(id_pe, year_1 = year, state_1 = state, W_j = cl_num,
                rule = fifelse(!is.na(cl_num), 1L, NA_integer_))]

# Rule 2: first later visit with a value, minus elapsed years. If that is <= 0:
#   first state regen     -> W_j = 0
#   first state not regen -> W_j = the later class value itself (not subtracted),
#                            since a non-regen stand cannot be at age 0 (P1 would find no path).
miss <- wj[is.na(W_j), id_pe]
r2 <- src[id_pe %in% miss & !is.na(cl_num), .SD[1L], by = id_pe][, .(id_pe, cl_later = cl_num, year_later = year)]
wj[, `:=`(W_r2 = NA_real_, cl_later = NA_real_)]
wj[r2, on = "id_pe", `:=`(W_r2 = i.cl_later - (i.year_later - year_1), cl_later = i.cl_later)]
n_nonpos <- wj[!is.na(W_r2) & W_r2 <= 0, .(regen = sum(state_1 == "regen"), non_regen = sum(state_1 != "regen"))]
wj[!is.na(W_r2) & W_r2 <= 0, W_r2 := fifelse(state_1 == "regen", 0, cl_later)]
cat(sprintf("Rule 2 back-calculation <= 0: %d plots with first state regen -> W_j 0; %d non-regen -> later class value\n",
            n_nonpos$regen, n_nonpos$non_regen))
wj[is.na(W_j) & !is.na(W_r2), `:=`(W_j = W_r2, rule = 2L)]

# Rule 3: min(time since fire, time since harvest) at the first visit (raw codes; lists from 03)
fire_o <- c("BR", "CRB"); fire_p <- "BRP"
harv_o <- c("CT", "CPR", "CBA", "CBT", "CEF", "CRS", "CPH", "CPT", "CS", "CDV", "CPP")
harv_p <- c("CP", "EPC", "CJ", "CDL", "EC", "CEA", "CAM", "CPI", "CPF", "CJG", "CJP", "CJT", "CPS", "CB",
            "CTR", "DLD", "CPM", "CPC", "CIP", "CA", "CPX", "ESI", "ECE", "CE", "CPR")
f <- src[, .SD[1L], by = id_pe]
f[, `:=`(an_o = as.integer(as.character(an_origine)), an_p = as.integer(as.character(an_perturb)))]
f[, fire_yr := fcase(origine %in% fire_o, an_o, perturb %in% fire_p, an_p)]
f[, harv_yr := fcase(origine %in% harv_o, an_o, perturb %in% harv_p, an_p)]
f[!is.na(fire_yr) & fire_yr > year, fire_yr := NA_integer_]
f[!is.na(harv_yr) & harv_yr > year, harv_yr := NA_integer_]
f[, W_r3 := suppressWarnings(pmin(year - fire_yr, year - harv_yr, na.rm = TRUE))]
f[is.infinite(W_r3), W_r3 := NA_real_]
wj[f, on = "id_pe", W_r3 := i.W_r3]
wj[is.na(W_j) & !is.na(W_r3), `:=`(W_j = W_r3, rule = 3L)]

# Rule 4: still NA, first state regen AND the plot later exits regen -> 0
ex <- dt[, .(any_exit = any(state[-.N] == "regen" & state[-1L] != "regen")), by = id_pe]
wj[ex, on = "id_pe", any_exit := i.any_exit]
wj[is.na(W_j) & state_1 == "regen" & any_exit %in% TRUE, `:=`(W_j = 0, rule = 4L)]

# Rule 5: exclude the rest
wj[, include := !is.na(W_j)]
wj[, c("W_r2", "cl_later", "W_r3") := NULL]

cat("W_j rule breakdown:\n"); print(wj[, .N, by = .(rule, include)][order(rule)])
cat(sprintf("\nIncluded: %d / %d plots (%.2f%%)\n", sum(wj$include), nrow(wj), 100 * mean(wj$include)))
cat("First state of the excluded plots:\n"); print(wj[include == FALSE, .N, by = state_1][order(-N)])
cat("\nW_j distribution (included plots):\n"); print(wj[include == TRUE, .N, by = W_j][order(W_j)])
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
# OPEN DECISION A: which TSD column is used. "imputed" = minimum-plausible TSD (lower
# bound); "imputed_prev" = maximum bound. Set it before running; no default.
# ------------------------------------------------------------------------------
TSD_SOURCE <- "imputed"            # decision A: "imputed" = minimum-plausible TSD (lower bound)
stopifnot("Set TSD_SOURCE (decision A) in 04_prepare_inla.R" = TSD_SOURCE %in% c("imputed", "imputed_prev"))
sfx <- paste0("_", TSD_SOURCE)
dt2[, `:=`(tsd_fire = get(paste0("tsf", sfx)), tsd_pest = get(paste0("sfp", sfx)),
           tsd_harv = get(paste0("tsh", sfx)))]

# Stand-origin (origine) disturbance already present at the plot's FIRST visit without a
# year: TSD = W_j + (visit year - first visit year) = t_abs (the W_j trick). Only fills
# the missing cells of the matching type (fire / harvest).
orig1 <- src[, .SD[1L], by = id_pe][!is.na(origine), id_pe]
n_fill <- c(fire = dt2[id_pe %in% orig1 & fire == 1L & is.na(tsd_fire), .N],
            harv = dt2[id_pe %in% orig1 & harv >  0L & is.na(tsd_harv), .N])
dt2[id_pe %in% orig1 & fire == 1L & is.na(tsd_fire), tsd_fire := t_abs]
dt2[id_pe %in% orig1 & harv >  0L & is.na(tsd_harv), tsd_harv := t_abs]

# Any plot that still has a flagged interval (END-visit row) without a TSD is EXCLUDED
# (partial 'perturb' disturbances present at the first visit, no year, no usable bound).
endr <- dt2[, .SD[-1L], by = id_pe]
miss <- endr[(fire == 1L & is.na(tsd_fire)) | (pest == 1L & is.na(tsd_pest)) | (harv > 0L & is.na(tsd_harv)),
             .(n_intervals = .N), by = id_pe]
tsd_excl <- merge(miss, wj[, .(id_pe, W_j, state_1, year_1)], by = "id_pe")
cat(sprintf("\nTSD (source '%s'): origine W_j-trick fills: %d fire cells, %d harvest cells (%d origine-at-first-visit plots)\n",
            TSD_SOURCE, n_fill["fire"], n_fill["harv"], length(orig1)))
cat(sprintf("TSD exclusions: %d plots (%d flagged intervals without TSD)\n", nrow(tsd_excl), sum(miss$n_intervals)))
saveRDS(tsd_excl, file.path(TABLE_DIR, "tsd_excluded_plots.rds"))
dt2 <- dt2[!id_pe %in% miss$id_pe]

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
fwrite(tr, here("data-output", "data", sprintf("transition_set_%s.csv", RUN_TAG)))
cat("\nSaved: wj_table.rds, visit_intervals.rds, plot_info.rds, transition_set.rds (+ csv)\n")
