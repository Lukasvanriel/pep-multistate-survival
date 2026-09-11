# ==============================================================================
# 06_combine_covariates.R
#
# Combines the four visit-level covariate tables produced by 02, 03, 03b and
# 03c into one analysis-ready data.frame for the INLA multistate survival
# model (04_fit_inla.R). Originally the raw-PEP-wrangling stub
# (01_wrangle_pep.R); renamed/repurposed as the combination step since it
# necessarily runs AFTER 02/03/03b/03c, not before them -- the "01" prefix no
# longer reflected its place in the pipeline, hence the move to 06.
#
# Inputs (all keyed on id_pe + no_mes = plot + visit number):
#   - data/processed/state_classification.rds        (from 02)
#   - data/covariates/pep/disturbance_covariates.rds  (from 03)
#   - data/covariates/pep/soil_covariates.rds         (from 03b)
#   - data/covariates/pep/climate_covariates.rds      (from 03c)
# NOT used: data/covariates/pep/climate_by_year_raw.rds -- that's an
# intermediate per-calendar-year checkpoint cache from 03c, a different grain
# entirely (id_pe x year, not id_pe x visit), out of scope here.
#
# Join strategy -- ROW-LEVEL (per plot-observation):
#   Plain inner join on (id_pe, no_mes) across all four tables. If a single
#   visit is missing its covariates in any one table, only THAT visit is
#   dropped -- the rest of the plot's visits are kept. This replaces an
#   earlier, more conservative whole-plot version of this script that
#   dropped an entire plot if even one visit had a gap; that turned out to
#   be dropping far more data than intended (538 extra plots on top of the
#   already-applied upstream filtering) for a gap-risk that upstream
#   filtering already handles differently (see below).
#
#   Plantation plots are the one case that SHOULD be fully excluded,
#   "before" observations included -- there's no reliable way to tell
#   whether an early, unlabelled visit on an eventual plantation plot was
#   genuinely natural forest or just pre-dates the record that later
#   labelled it. That exclusion is NOT re-implemented here: it's handled
#   upstream in 02_classify_states.R (section 3), which now drops such
#   plots entirely from state_classification.rds, before included. Because
#   state_classification.rds is one of the four tables in this join, a
#   plantation plot simply has zero rows there and so cannot appear in the
#   combined output either -- no extra logic needed on this end.
#
# Planned steps:
#   1. Load the four input tables
#   2. Diagnostic gate: total_ba (state vs disturbance) and year
#      (disturbance vs soil vs climate) must agree wherever both are
#      present -- drop the duplicate copies once confirmed
#   3. Row-level inner join on (id_pe, no_mes) across all four tables
#   4. Assemble final combined table, save
#   5. Column-by-column diagnostics summary (for the paper's methods
#      section), save as CSV alongside the combined data
# ==============================================================================

library(data.table)
library(here)

## ---- 0. config ------------------------------------------------------------

STATE_PATH <- here("data", "processed", "state_classification.rds")
DIST_PATH  <- here("data", "covariates", "pep", "disturbance_covariates.rds")
SOIL_PATH  <- here("data", "covariates", "pep", "soil_covariates.rds")
CLIM_PATH  <- here("data", "covariates", "pep", "climate_covariates.rds")

OUT_PATH      <- here("data", "processed", "combined_covariates.rds")
SUMMARY_PATH  <- here("data-output", "data", "combined_covariates_summary.csv")

## ==============================================================================
## STEP 1. LOAD
## ==============================================================================

cat("Loading input tables...\n")
state <- as.data.table(readRDS(STATE_PATH))
dist  <- as.data.table(readRDS(DIST_PATH))
soil  <- as.data.table(readRDS(SOIL_PATH))
clim  <- as.data.table(readRDS(CLIM_PATH))

# normalize no_mes to integer everywhere (stored as double upstream, but
# integer-valued -- avoid float-equality surprises in the join)
state[, no_mes := as.integer(no_mes)]
dist[,  no_mes := as.integer(no_mes)]
soil[,  no_mes := as.integer(no_mes)]
clim[,  no_mes := as.integer(no_mes)]

tables <- list(state = state, dist = dist, soil = soil, clim = clim)
for (nm in names(tables)) {
  cat(sprintf("  %-6s %6d rows, %5d distinct plots\n",
              nm, nrow(tables[[nm]]), uniqueN(tables[[nm]]$id_pe)))
}

## ==============================================================================
## STEP 2. DIAGNOSTIC GATE -- DUPLICATE COLUMN AGREEMENT (total_ba, year)
##
## Checked on the RAW tables, before any join filtering -- these columns are
## expected to agree on every (id_pe, no_mes) pair the tables have in common,
## whether or not that pair ends up in the final row-level join.
## ==============================================================================

cat("\nChecking duplicate columns agree before dropping...\n")
chk_ba <- merge(state[, .(id_pe, no_mes, total_ba_state = total_ba)],
                 dist[, .(id_pe, no_mes, total_ba_dist = total_ba)],
                 by = c("id_pe", "no_mes"))
n_ba_mismatch <- chk_ba[abs(total_ba_state - total_ba_dist) > 1e-6, .N]
cat(sprintf("  total_ba (state vs disturbance): %d / %d overlapping rows mismatch\n",
            n_ba_mismatch, nrow(chk_ba)))
if (n_ba_mismatch > 0) {
  stop("total_ba disagrees between state_classification and disturbance_covariates ",
       "for overlapping rows -- investigate before dropping either copy.")
}

chk_year <- merge(dist[, .(id_pe, no_mes, year_dist = year)],
                   soil[, .(id_pe, no_mes, year_soil = year)],
                   by = c("id_pe", "no_mes"))
chk_year <- merge(chk_year, clim[, .(id_pe, no_mes, year_clim = year)],
                   by = c("id_pe", "no_mes"))
n_year_mismatch <- chk_year[year_dist != year_soil | year_soil != year_clim, .N]
cat(sprintf("  year (disturbance vs soil vs climate): %d / %d overlapping rows mismatch\n",
            n_year_mismatch, nrow(chk_year)))
if (n_year_mismatch > 0) {
  stop("year disagrees across disturbance/soil/climate covariates for overlapping ",
       "rows -- investigate before dropping duplicate copies.")
}

## drop the now-confirmed-redundant columns before merging (explicit keep-list,
## matching the select() style used throughout 02/03/03b/03c)
dist_f <- dist[, .(id_pe, no_mes, year, fire, harv, pest, tsf, tsh, sfp,
                    tsf_imputed, tsh_imputed, sfp_imputed,
                    tsf_imputed_prev, tsh_imputed_prev, sfp_imputed_prev)]
soil_f <- soil[, .(id_pe, no_mes, ph_humus, drainage_ordinal, soil_pc1, soil_pc2,
                    soil_pc_source, soil_pc_donor_id_pe, soil_pc_donor_dist_km)]
clim_f <- clim[, .(id_pe, no_mes, win_start, win_end, start_gs, end_gs, length_gs,
                    pcp_1, pcp_2, pcp_3, pcp_4, gdd_1, gdd_2, gdd_3, gdd_4,
                    an_meanT, an_minT, an_maxT, an_meanT_3, an_range_T_3,
                    cmi_1, cmi_2, cmi_3, cmi_4, cmi_5, cmi_6, cmi_7, cmi_8,
                    cmi_9, cmi_10, cmi_11, cmi_12, cmi_sum, n_years_used)]

## ==============================================================================
## STEP 3. ROW-LEVEL INNER JOIN
##
## Keeps a visit as long as it has covariates in all four tables; a visit
## missing from just one table is dropped on its own, without taking the
## rest of that plot's visits with it. (Plantation plots are the exception
## meant to disappear entirely -- see header comment: that's handled
## upstream in 02, not here.)
## ==============================================================================

cat("\nJoining (row-level, inner join on id_pe + no_mes)...\n")

step1 <- merge(state, dist_f, by = c("id_pe", "no_mes"))
cat(sprintf("  state (%d rows) x disturbance (%d rows) -> %d rows (%d dropped)\n",
            nrow(state), nrow(dist_f), nrow(step1), nrow(state) - nrow(step1)))

step2 <- merge(step1, soil_f, by = c("id_pe", "no_mes"))
cat(sprintf("  ...x soil (%d rows) -> %d rows (%d further dropped)\n",
            nrow(soil_f), nrow(step2), nrow(step1) - nrow(step2)))

combined <- merge(step2, clim_f, by = c("id_pe", "no_mes"))
cat(sprintf("  ...x climate (%d rows) -> %d rows (%d further dropped)\n",
            nrow(clim_f), nrow(combined), nrow(step2) - nrow(combined)))
setorder(combined, id_pe, no_mes)

## ==============================================================================
## STEP 4. ASSEMBLE + SAVE
## ==============================================================================

cat(sprintf("\nFinal combined table: %d rows, %d columns, %d distinct plots\n",
            nrow(combined), ncol(combined), uniqueN(combined$id_pe)))

if (nrow(combined) != uniqueN(combined, by = c("id_pe", "no_mes"))) {
  stop("Combined table has duplicate (id_pe, no_mes) keys -- join produced a ",
       "cartesian blow-up somewhere. Investigate before saving.")
}

dir.create(dirname(OUT_PATH), showWarnings = FALSE, recursive = TRUE)
saveRDS(combined, OUT_PATH)
cat(sprintf("Saved: %s\n", OUT_PATH))

## ==============================================================================
## STEP 5. COLUMN-BY-COLUMN DIAGNOSTICS (for methods section)
## ==============================================================================

cat("\nBuilding column summary...\n")

summarize_column <- function(x, name) {
  n <- length(x)
  n_missing <- sum(is.na(x))
  pct_missing <- round(100 * n_missing / n, 2)
  if (is.numeric(x)) {
    xx <- x[!is.na(x)]
    data.table(
      variable = name, type = "numeric", n = n, n_missing = n_missing,
      pct_missing = pct_missing,
      mean = if (length(xx)) round(mean(xx), 4) else NA_real_,
      sd   = if (length(xx)) round(sd(xx), 4)   else NA_real_,
      min  = if (length(xx)) round(min(xx), 4)  else NA_real_,
      max  = if (length(xx)) round(max(xx), 4)  else NA_real_,
      n_unique = uniqueN(x[!is.na(x)]),
      top_values = NA_character_
    )
  } else {
    tab <- sort(table(x), decreasing = TRUE)
    top <- paste(sprintf("%s=%d", names(tab)[seq_len(min(5, length(tab)))],
                          as.integer(tab)[seq_len(min(5, length(tab)))]),
                 collapse = "; ")
    data.table(
      variable = name, type = "categorical", n = n, n_missing = n_missing,
      pct_missing = pct_missing, mean = NA_real_, sd = NA_real_,
      min = NA_real_, max = NA_real_, n_unique = uniqueN(x[!is.na(x)]),
      top_values = top
    )
  }
}

summary_dt <- rbindlist(lapply(names(combined), function(nm) summarize_column(combined[[nm]], nm)))

cat("\n--- COMBINED DATASET SUMMARY ---\n")
print(summary_dt)

dir.create(dirname(SUMMARY_PATH), showWarnings = FALSE, recursive = TRUE)
fwrite(summary_dt, SUMMARY_PATH)
cat(sprintf("\nSaved column summary: %s\n", SUMMARY_PATH))

## ---- notes for methods write-up -------------------------------------------
# - tsf/tsh/sfp (and their imputed variants) are NA by design when that
#   disturbance is not currently active -- high missingness there is not a
#   data-quality problem, don't report it as one without this caveat.
# - soil_pc_donor_id_pe / soil_pc_donor_dist_km are NA except for the small
#   subset of plots filled via nearest-donor imputation (soil_pc_source ==
#   "nearest_donor") -- see 03b for the imputation method.
# - start_gs/end_gs/length_gs/pcp_*/gdd_*/an_*/cmi_* share a small (<1%)
#   missingness footprint: plots whose 10-yr trailing climate window fell
#   partly or fully outside the 1960-2022 data availability (n_years_used
#   captures how many real years actually went into each average).
# - This combine step joins at the ROW (plot-visit) level: a visit missing
#   covariates in any one input table is dropped on its own, without taking
#   the rest of that plot's visits with it. Some retained plots therefore
#   have a gap in their visit sequence (e.g. visits 1, 2, 4 present but 3
#   missing) -- if the multistate model treats consecutive retained visits
#   as adjacent transitions, a gap like this silently spans more real time
#   than a normal interval. Check data/processed/combined_covariates.rds
#   per-plot no_mes sequences for gaps before fitting if this matters for
#   the model structure.
# - Plantation plots are the deliberate exception: they are fully absent
#   from this output (all visits, "before" included), because
#   02_classify_states.R drops them entirely upstream. See that script's
#   section 3 for the rationale and which plots/rows that affects.
