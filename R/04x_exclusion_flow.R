# ==============================================================================
# 04x_exclusion_flow.R  (T8, 2026-10-26) -- exclusion-flow table: what each rule removes, in pipeline order.
# Run from the project root AFTER 04_prepare_inla.R:  Rscript R/04x_exclusion_flow.R
#
# Every plot-visit is followed from the raw PEP visits to the INLA input. A visit is attributed to the FIRST step that removes it.
# For every step: plots / visits / intervals / state changes still in the data, and what the step removed. Intervals = consecutive
# remaining visits of a plot (as 04 builds them); state change = the 02 state differs between the two visits. Before step 1 (02) there
# are no states, so only plots and visits are counted.
# Steps: 1 02 state classification (its own exclusions: see the 02 log) | 2 03 drop-plot codes (pee, station only) | 3 03 CHT/DT
#   truncation (pee / station; plot dropped / visits from the code removed) | 4 03 other (not in the disturbance table) | 5 03z inner join
#   (no soil / no climate covariates) | 6 04 no stand age (W_j) | 7 04 old regen (T4 decided 2026-10-26: NO exclusion in
#   the main run; sensitivity S8 excludes first-state-regen plots with W_j >= 70 and writes TABLE_DIR/old_regen_excluded.rds) | 8 04 plot
#   left with one visit (no interval). Then: state changes in transitions that are not modelled (< MIN_EVENTS).
# Checks: the final counts must equal 04's visit-interval table (stops otherwise).
# Inputs: data/raw/pep/PEP.gpkg (placette_mes), data/build/{state_classification, disturbance_exclusions, disturbance_covariates,
#   soil_covariates, climate_covariates, combined_covariates, stand_age}.rds, TABLE_DIR/{visit_intervals, transition_set}.rds
#   (+ TABLE_DIR/old_regen_excluded.rds in the S8 run only: a table with id_pe).
# Output: TABLE_DIR/exclusion_flow.csv (the table) and TABLE_DIR/exclusions_by_visit.rds (id_pe, no_mes, step, rule).
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(sf); library(here) })
source(here::here("R", "config.R"))
options(width = 200)

rd <- function(f) { x <- as.data.table(readRDS(f)); if ("no_mes" %in% names(x)) x[, no_mes := as.integer(no_mes)]; x }
raw   <- as.data.table(st_read(here("data", "raw", "pep", "PEP.gpkg"), layer = "placette_mes", quiet = TRUE))[, .(id_pe, no_mes = as.integer(no_mes))]
sc    <- rd(file.path(BUILD_DIR, "state_classification.rds"))[, .(id_pe, no_mes, state)]
excl  <- rd(file.path(BUILD_DIR, "disturbance_exclusions.rds"))
dist  <- rd(file.path(BUILD_DIR, "disturbance_covariates.rds"))[, .(id_pe, no_mes)]
soil  <- rd(file.path(BUILD_DIR, "soil_covariates.rds"))[, .(id_pe, no_mes)]
clim  <- rd(file.path(BUILD_DIR, "climate_covariates.rds"))[, .(id_pe, no_mes)]
comb  <- rd(file.path(BUILD_DIR, "combined_covariates.rds"))[, .(id_pe, no_mes)]
sa    <- as.data.table(readRDS(file.path(BUILD_DIR, "stand_age.rds")))
iv    <- as.data.table(readRDS(file.path(TABLE_DIR, "visit_intervals.rds")))
tr    <- as.data.table(readRDS(file.path(TABLE_DIR, "transition_set.rds")))
OLD_REGEN_PATH <- file.path(TABLE_DIR, "old_regen_excluded.rds")
old_regen <- if (file.exists(OLD_REGEN_PATH)) unique(as.data.table(readRDS(OLD_REGEN_PATH))$id_pe) else NULL
key <- function(x) paste(x$id_pe, x$no_mes)

# ---- attribute every raw visit to the first step that removes it ------------------------------------------------------------------------------
v <- unique(raw)[, `:=`(step = NA_integer_, rule = NA_character_)]
setkey(v, id_pe, no_mes)
drop_at <- function(sel, st, rl) { i <- which(is.na(v$step) & sel); if (length(i)) set(v, i, c("step", "rule"), list(st, rl)); invisible(length(i)) }
drop_at(!key(v) %in% key(sc), 1L, "02 state classification (02's own exclusions)")
ek <- key(excl)
for (rl in sort(unique(excl$rule))) drop_at(key(v) %in% ek[excl$rule == rl], if (grepl("drop-plot", rl)) 2L else 3L, rl)
drop_at(!key(v) %in% key(dist), 4L, "03 other (visit not in the disturbance table)")
drop_at(!key(v) %in% key(soil), 5L, "03z no soil covariates")
drop_at(!key(v) %in% key(clim), 5L, "03z no climate covariates")
drop_at(!key(v) %in% key(comb), 5L, "03z other (not in the combined table)")
drop_at(v$id_pe %in% sa[include == FALSE, id_pe] | !v$id_pe %in% sa$id_pe, 6L, "04 no stand age (W_j)")
if (!is.null(old_regen)) drop_at(v$id_pe %in% old_regen, 7L, "04 old regen (T4)")
v[is.na(step), n_left := .N, by = id_pe]
drop_at(!is.na(v$n_left) & v$n_left == 1L, 8L, "04 plot left with one visit (no interval)")
v[, n_left := NULL]

# ---- counts still in the data after each step ------------------------------------------------------------------------------------------------
v[sc, state := i.state, on = c("id_pe", "no_mes")]
count_left <- function(after) {                       # visits not removed by a step <= after
  x <- v[is.na(step) | step > after][order(id_pe, no_mes)]
  x[, nxt := shift(state, -1L), by = id_pe]
  n_int <- x[, .(n = .N - 1L), by = id_pe][, sum(n)]
  data.table(plots = uniqueN(x$id_pe), visits = nrow(x), intervals = if (after >= 1L) n_int else NA_integer_,
             state_changes = if (after >= 1L) x[!is.na(nxt) & !is.na(state), sum(state != nxt)] else NA_integer_)
}
steps <- data.table(step = 0:8, name = c("0 raw PEP visits (placette_mes)", "1 02 state classification", "2 03 drop-plot codes",
  "3 03 CHT/DT truncation", "4 03 other", "5 03z covariate join (soil, climate)", "6 04 no stand age", "7 04 old regen (T4)", "8 04 single-visit plots"))
flow <- rbindlist(lapply(steps$step, function(s) cbind(steps[step == s], count_left(s))))
flow[, `:=`(plots_removed = shift(plots) - plots, visits_removed = shift(visits) - visits,
            intervals_removed = shift(intervals) - intervals, state_changes_removed = shift(state_changes) - state_changes)]
if (is.null(old_regen)) flow[step == 7L, name := "7 04 old regen: none in the main run (S8)"]

# ---- check against 04 ------------------------------------------------------------------------------------------------------------------------
fin <- flow[step == 8L]
stopifnot("flow table does not reproduce 04's plots"         = fin$plots == uniqueN(iv$id_pe),
          "flow table does not reproduce 04's intervals"     = fin$intervals == nrow(iv),
          "flow table does not reproduce 04's state changes" = fin$state_changes == sum(iv$is_change))

# ---- state changes in transitions that are not modelled -----------------------------------------------------------------------------------------
ch <- iv[is_change == TRUE, .N, by = .(from = state_start, to = state_end)]
ch[tr, kept := i.kept, on = c("from", "to")]
not_mod <- data.table(step = 9L, name = sprintf("9 transitions not modelled (< %d state changes)", if (exists("MIN_EVENTS")) MIN_EVENTS else 15L),
                      plots = NA_integer_, visits = NA_integer_, intervals = NA_integer_, state_changes = ch[kept == TRUE, sum(N)],
                      plots_removed = NA_integer_, visits_removed = NA_integer_, intervals_removed = NA_integer_,
                      state_changes_removed = ch[kept == FALSE, sum(N)])
flow <- rbind(flow, not_mod)

by_rule <- v[!is.na(step), .(plots = uniqueN(id_pe), visits = .N), by = .(step, rule)][order(step, rule)]
cat("Exclusion flow (", RUN_TAG, "): numbers still in the data after each step, and what the step removed\n", sep = "")
print(flow[, .(name, plots, visits, intervals, state_changes, plots_removed, visits_removed, intervals_removed, state_changes_removed)], row.names = FALSE)
cat("\nRemoved visits by rule (a visit counts at the first step that removes it; plots = plots with >= 1 visit removed by that rule):\n")
print(by_rule, row.names = FALSE)
cat(sprintf("\nCheck: final plots / intervals / state changes = 04's visit-interval table (%d / %d / %d): OK\n", fin$plots, fin$intervals, fin$state_changes))
fwrite(rbind(flow[, .(table = "flow", name, plots, visits, intervals, state_changes, plots_removed, visits_removed, intervals_removed, state_changes_removed)],
             by_rule[, .(table = "removed_by_rule", name = paste(step, rule), plots, visits, intervals = NA_integer_, state_changes = NA_integer_,
                         plots_removed = plots, visits_removed = visits, intervals_removed = NA_integer_, state_changes_removed = NA_integer_)]),
       file.path(TABLE_DIR, "exclusion_flow.csv"))
saveRDS(v[!is.na(step), .(id_pe, no_mes, step, rule)], file.path(TABLE_DIR, "exclusions_by_visit.rds"))
cat("Saved:", file.path(TABLE_DIR, "exclusion_flow.csv"), "and exclusions_by_visit.rds\n")
