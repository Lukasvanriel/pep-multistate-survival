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
#   - data/build/state_classification.rds        (from 02)
#   - data/build/disturbance_covariates.rds      (from 03)
#   - data/build/soil_covariates.rds             (from 03b)
#   - data/build/climate_covariates.rds          (from 03c)
# NOT used: data/build/climate_by_year_raw.rds -- that's an
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
#   6. Transition diagnostics: counts of transitions between states,
#      disturbance occurrence per transition, and inter-visit interval
#      distributions (overall and per transition) -- saved as CSVs for the
#      methods section
#   7. Transition intensity matrix figure: from_state x to_state heatmap,
#      cells with n < 15 masked black -- saved as PNG for the methods section
# ==============================================================================

library(data.table)
library(here)
library(ggplot2)

## ---- 0. config ------------------------------------------------------------

source(here::here("R", "config.R"))   # project paths (data/build etc.)
dir.create(BUILD_SUM_DIR, recursive = TRUE, showWarnings = FALSE)
STATE_PATH <- file.path(BUILD_DIR, "state_classification.rds")
DIST_PATH  <- file.path(BUILD_DIR, "disturbance_covariates.rds")
SOIL_PATH  <- file.path(BUILD_DIR, "soil_covariates.rds")
CLIM_PATH  <- file.path(BUILD_DIR, "climate_covariates.rds")

OUT_PATH      <- file.path(BUILD_DIR, "combined_covariates.rds")
SUMMARY_PATH  <- file.path(BUILD_SUM_DIR, "combined_covariates_summary.csv")

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
                    tsf_imputed_prev, tsh_imputed_prev, sfp_imputed_prev,
                    tsf_imputed_mid, tsh_imputed_mid, sfp_imputed_mid,   # added 2026-10-26 (TSD_SOURCE "mid")
                    fire_dating, harv_dating, pest_dating)]
soil_f <- soil[, .(id_pe, no_mes, ph_humus, drainage_ordinal, soil_pc1, soil_pc2,
                    soil_pc_source, soil_pc_donor_id_pe, soil_pc_donor_dist_km)]
stopifnot("cmi_gs missing from climate_covariates.rds: rerun 03c (step 4b)" = "cmi_gs" %in% names(clim))
clim_f <- clim[, .(id_pe, no_mes, win_start, win_end, start_gs, end_gs, length_gs,
                    pcp_1, pcp_2, pcp_3, pcp_4, gdd_1, gdd_2, gdd_3, gdd_4,
                    an_meanT, an_minT, an_maxT, an_meanT_3, an_range_T_3,
                    cmi_1, cmi_2, cmi_3, cmi_4, cmi_5, cmi_6, cmi_7, cmi_8,
                    cmi_9, cmi_10, cmi_11, cmi_12, cmi_gs, n_years_used)]
## Growing-season CMI (cmi_gs): built in 03c step 4b since 2026-10-24 = monthly CMI weighted by the share of each month inside the
## plot-year's NRCan growing season (the window of gdd_3), averaged over the same 10-yr window. It replaces the fixed May-Oct sum of
## 2026-10-22. Annual CMI (cmi_sum) is NOT carried: 73% of its between-plot variance came from Oct-Apr (winter precipitation).

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

## ==============================================================================
## STEP 6. TRANSITION DIAGNOSTICS (for methods section)
##
## One row per observed transition: a plot's state at visit i vs. its state
## at the next RETAINED visit. "Next retained" is not necessarily "no_mes +
## 1" -- the row-level join (STEP 3) can leave gaps in a plot's visit
## sequence, so a transition can silently span more than one nominal
## measurement interval. no_mes_gap flags this below.
## ==============================================================================

cat("\nBuilding transition diagnostics...\n")

setorder(combined, id_pe, no_mes)
combined[, `:=`(
  .to_state  = shift(state, type = "lead"),
  .to_no_mes = shift(no_mes, type = "lead"),
  .to_year   = shift(year, type = "lead"),
  .to_fire   = shift(fire, type = "lead"),
  .to_harv   = shift(harv, type = "lead"),
  .to_pest   = shift(pest, type = "lead")
), by = id_pe]

trans <- combined[!is.na(.to_state), .(
  id_pe, no_mes, from_state = state, to_state = .to_state,
  no_mes_gap = .to_no_mes - no_mes,
  dt_years   = .to_year - year,
  fire = .to_fire, harv = .to_harv, pest = .to_pest
)]

combined[, c(".to_state", ".to_no_mes", ".to_year", ".to_fire", ".to_harv", ".to_pest") := NULL]

n_trans <- nrow(trans)
n_gap   <- trans[no_mes_gap > 1, .N]
cat(sprintf("  %d transitions, from %d plots\n", n_trans, uniqueN(trans$id_pe)))
cat(sprintf("  %d / %d transitions (%.1f%%) skip at least one no_mes (row-level join gap)\n",
            n_gap, n_trans, 100 * n_gap / n_trans))

## -- 6a. transition counts (from_state -> to_state) -------------------------
transition_counts <- trans[, .N, by = .(from_state, to_state)]
transition_counts[, pct := round(100 * N / n_trans, 2)]
setorder(transition_counts, -N)

cat("\n--- TRANSITION COUNTS (from_state -> to_state) ---\n")
print(transition_counts)

## -- 6b. disturbance occurrence per transition -------------------------------
cat(sprintf("\n  overall disturbance occurrence: fire %.1f%%, harv %.1f%%, pest %.1f%%\n",
            100 * mean(trans$fire > 0, na.rm = TRUE),
            100 * mean(trans$harv > 0, na.rm = TRUE),
            100 * mean(trans$pest > 0, na.rm = TRUE)))

disturbance_by_transition <- trans[, .(
  n        = .N,
  n_fire   = sum(fire > 0, na.rm = TRUE),
  pct_fire = round(100 * mean(fire > 0, na.rm = TRUE), 2),
  n_harv   = sum(harv > 0, na.rm = TRUE),
  pct_harv = round(100 * mean(harv > 0, na.rm = TRUE), 2),
  n_pest   = sum(pest > 0, na.rm = TRUE),
  pct_pest = round(100 * mean(pest > 0, na.rm = TRUE), 2)
), by = .(from_state, to_state)]
setorder(disturbance_by_transition, -n)

cat("\n--- DISTURBANCE OCCURRENCE BY TRANSITION ---\n")
cat("(fire/harv/pest read from the DESTINATION visit -- i.e. what the\n")
cat(" disturbance covariates record as having occurred since the previous visit)\n")
print(disturbance_by_transition)

## -- 6c. inter-visit interval (years): overall, then by transition ----------
interval_stats <- function(x) {
  x <- x[!is.na(x)]
  data.table(
    n = length(x), mean = round(mean(x), 2), sd = round(sd(x), 2),
    min = min(x), p10 = round(quantile(x, .10), 1), p25 = round(quantile(x, .25), 1),
    median = round(median(x), 1), p75 = round(quantile(x, .75), 1),
    p90 = round(quantile(x, .90), 1), max = max(x)
  )
}

interval_overall <- interval_stats(trans$dt_years)
cat("\n--- INTER-VISIT INTERVAL (years), ALL TRANSITIONS ---\n")
print(interval_overall)

interval_by_transition <- trans[, interval_stats(dt_years), by = .(from_state, to_state)]
setorder(interval_by_transition, from_state, to_state)
cat("\n--- INTER-VISIT INTERVAL (years), BY TRANSITION ---\n")
print(interval_by_transition)

## -- save for the methods section -------------------------------------------
TRANSITION_COUNTS_PATH <- file.path(BUILD_SUM_DIR, "transition_counts.csv")
DISTURBANCE_TRANS_PATH <- file.path(BUILD_SUM_DIR, "transition_disturbance_summary.csv")
INTERVAL_OVERALL_PATH  <- file.path(BUILD_SUM_DIR, "interval_overall_summary.csv")
INTERVAL_BY_TRANS_PATH <- file.path(BUILD_SUM_DIR, "interval_by_transition_summary.csv")

fwrite(transition_counts, TRANSITION_COUNTS_PATH)
fwrite(disturbance_by_transition, DISTURBANCE_TRANS_PATH)
fwrite(interval_overall, INTERVAL_OVERALL_PATH)
fwrite(interval_by_transition, INTERVAL_BY_TRANS_PATH)
cat(sprintf("\nSaved: %s\n", TRANSITION_COUNTS_PATH))
cat(sprintf("Saved: %s\n", DISTURBANCE_TRANS_PATH))
cat(sprintf("Saved: %s\n", INTERVAL_OVERALL_PATH))
cat(sprintf("Saved: %s\n", INTERVAL_BY_TRANS_PATH))

## ==============================================================================
## STEP 7. TRANSITION INTENSITY MATRIX (figure for methods section)
##
## Heatmap of transition_counts: from_state (rows/x) x to_state (cols/y),
## filled by log10(n) on a single-hue sequential ramp. Cells with n < 15 --
## including (from, to) pairs that never co-occur, n = 0 -- are masked solid
## black: too few observed transitions to trust a rate or coefficient
## estimated from that cell in the multistate model.
## ==============================================================================

cat("\nBuilding transition intensity matrix...\n")

state_levels <- sort(unique(c(trans$from_state, trans$to_state)))

mat_dt <- CJ(from_state = state_levels, to_state = state_levels)
mat_dt <- merge(mat_dt, transition_counts[, .(from_state, to_state, N)],
                by = c("from_state", "to_state"), all.x = TRUE)
mat_dt[is.na(N), N := 0L]
mat_dt[, masked := N < 15]
mat_dt[, label_color := fifelse(masked | N >= 50, "white", "#0b0b0b")]
mat_dt[, from_state := factor(from_state, levels = state_levels)]
mat_dt[, to_state := factor(to_state, levels = rev(state_levels))]

seq_ramp <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")

p_matrix <- ggplot(mat_dt, aes(x = from_state, y = to_state)) +
  geom_tile(data = mat_dt[masked == FALSE], aes(fill = log10(N)),
            color = "white", linewidth = 0.5) +
  geom_tile(data = mat_dt[masked == TRUE], fill = "#0b0b0b",
            color = "white", linewidth = 0.5) +
  geom_text(aes(label = N, color = label_color), size = 2.6) +
  scale_color_identity() +
  scale_fill_gradientn(colors = seq_ramp, name = "n transitions\n(log10 scale)",
                        breaks = log10(c(1, 10, 100, 1000)),
                        labels = c(1, 10, 100, 1000)) +
  labs(title = "Transition intensity matrix",
       subtitle = sprintf("n = %d transitions total; %d / %d cells masked (n < 15)",
                           n_trans, sum(mat_dt$masked), nrow(mat_dt)),
       x = "From state", y = "To state") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold"))

FIGURE_PATH <- file.path(BUILD_SUM_DIR, "transition_intensity_matrix.png")
dir.create(dirname(FIGURE_PATH), showWarnings = FALSE, recursive = TRUE)
ggsave(FIGURE_PATH, p_matrix, width = 9, height = 7.5, dpi = 200)
cat(sprintf("Saved: %s\n", FIGURE_PATH))

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
#   than a normal interval. Check data/build/combined_covariates.rds
#   per-plot no_mes sequences for gaps before fitting if this matters for
#   the model structure.
# - Plantation plots are the deliberate exception: they are fully absent
#   from this output (all visits, "before" included), because
#   02_classify_states.R drops them entirely upstream. See that script's
#   section 3 for the rationale and which plots/rows that affects.
# - STEP 6's no_mes_gap quantifies the row-level join gap directly: a
#   transition with no_mes_gap > 1 skips at least one nominal visit, so its
#   dt_years spans more real time than a normal interval. Check the console
#   output's gap rate before treating all transitions as equal-opportunity
#   in the multistate model.
