# ==============================================================================
# pep_soil_climate_covariates.R
#
# Builds the soil covariate for PEP permanent plots (climate: still TODO).
#
# DECISION HISTORY (see project discussion for full reasoning):
#   - Started from every soil-related variable in station_sol/station_pe;
#     narrowed to ph_humus and cl_drai (drainage class) because (a) they're
#     the best-covered soil variables, (b) Brice et al. 2020's own model
#     uses exactly these two (PH_HUMUS, DRAIN), and (c) drainage and
#     epmatorg are too collinear (r=0.69) to both be useful, while
#     ph_horizb is nearly orthogonal to drainage (r=0.02) and mostly
#     redundant with ph_humus's role -- see conversation for the full
#     variable-selection trail.
#   - Soil is only measured at a subset of visits, and checking WHERE
#     showed measurement is concentrated on the plot's LAST visit (96% of
#     multi-visit plots with any real measurement have it on their last
#     visit, only 24% on their first) -- so a plain forward-only fill
#     barely helps. Instead: take each plot's single most recent REAL
#     (non-NA) value and use it for every visit of that plot. A real
#     measurement is never overwritten by this rule; it just decides what
#     fills the gaps.
#   - ph_humus and drainage are correlated (r=0.22 raw, 0.22 after fill --
#     the fill doesn't change this) in the direction the minerotrophic/
#     ombrotrophic wetland literature predicts: poorer drainage tends to
#     go with LESS acidic (higher) humus pH, because poorly-drained sites
#     fed by lateral groundwater/seepage bring in base cations, while
#     well-drained sites in this climate get leached and acidify
#     (podzolization). Checked this isn't a coding error (drainage codes
#     and pH field definitions confirmed against the official PEP data
#     dictionary), isn't a Simpson's-paradox artifact of mixing organic
#     and mineral soils (both subgroups show the same-direction, weaker
#     correlation individually), but IS partly leverage-driven: Spearman's
#     rho (0.11) is about half of Pearson's r (0.22), and dropping the
#     most extreme drainage classes or winsorizing pH tails pulls r down
#     to 0.17-0.19. Treat the relationship as real but weaker/more
#     tail-driven than a raw Pearson r would suggest.
#   - Because pH and drainage are positively correlated with each other,
#     but the tree-species-suitability literature (Carteron et al. 2020;
#     sugar maple range-limit studies) shows high pH and good drainage are
#     BOTH favorable for temperate species -- i.e. the two variables pull
#     in OPPOSITE directions relative to species performance despite
#     moving together in the raw data -- the two orthogonal combinations
#     of these variables have different, both useful, interpretations.
#     Rather than pick one, we keep BOTH axes (soil_pc1, soil_pc2) and
#     decide what to do with them once we see how they behave in the
#     model.
#   - Of the 11,393 analysis-ready plots (checked carefully after an
#     earlier miscount -- see project discussion): 170 have no real
#     ph_humus ever, 85 have no real drainage ever, 78 are missing BOTH
#     entirely, and 99 are missing exactly ONE of the two (92 have
#     drainage but not pH, 7 have pH but not drainage) -- 177 plots
#     total missing at least one variable (170+85-78=177, cross-checked).
#     These are not a random subset: fewer visits on average, and
#     markedly less disturbed (fire ~4% vs ~14%, pest ~3% vs ~13%
#     plot-ever rate) than the rest of the dataset -- reason unknown
#     (possibly an older-cohort/protocol artifact: median first visit
#     1974 vs 1976), noted for later interpretation. Rejected an
#     environmental-covariate-similarity approach (the actual digital
#     soil mapping literature standard -- Mansuy et al. NRCan; PeerJ
#     Quebec ecoforest soil-texture paper, which also uses PEP data but
#     for texture, not pH/drainage, and doesn't publish reusable
#     predictions) as disproportionate effort for this scale of problem.
#     Instead: nearest-neighbor donor fallback, but applied to the FINAL
#     PC axes (soil_pc1/soil_pc2) rather than to the two raw variables
#     separately -- one mechanism that fixes both the 78-missing-both
#     and the 99-missing-one cases at once, since both simply show up as
#     "PC score not computed" after step 5. See step 5b below. Every
#     borrowed row is flagged (soil_pc_source / soil_pc_donor_id_pe /
#     soil_pc_donor_dist_km) for a sensitivity check once transition
#     models exist: refit with these plots' PC scores imputed vs. with
#     their rows dropped entirely, and compare the soil coefficient.
#
# Planned steps:
#   1. Load ph_humus (station_sol) and cl_drai (station_pe), recode
#      drainage to an ordinal 0 (excessif/best) - 6 (tres mauvais/worst)  [DONE]
#   2. Merge onto every plot visit (placette_mes)                        [DONE]
#   3. Last-REAL-value broadcast per plot, per variable                  [DONE]
#   4. Diagnostic: coverage before/after, and how often the adopted last
#      real value contradicts an earlier-agreeing real value             [DONE]
#   5. PCA on the two standardized variables -> soil_pc1, soil_pc2, with
#      an enforced sign convention (see interpretation notes below)      [DONE]
#   5b. Nearest-donor fallback ON THE PC AXES for any row still missing
#       soil_pc1/soil_pc2 (covers BOTH plots missing one raw variable
#       and plots missing both -- one mechanism instead of two)          [DONE]
#   6. Save id_pe / no_mes / id_pe_mes / year / ph_humus / drainage /
#      soil_pc1 / soil_pc2 / soil_pc_source / soil_pc_donor_id_pe /
#      soil_pc_donor_dist_km to disk                                     [DONE]
#   7. Climate covariates                                                [TODO]
#
# INTERPRETATION of the two output axes (soil_pc1, soil_pc2):
#   soil_pc1 ("agreement" axis, ~61% of the two variables' joint variance)
#     - loadings: ph_humus +0.707, drainage +0.707 (same sign, enforced)
#     - HIGH soil_pc1 = high pH AND poor drainage together -- minerotrophic,
#       richer wet sites (fens/swamps fed by groundwater seepage)
#     - LOW soil_pc1  = low pH AND good drainage together -- leached,
#       acidic, well-drained upland/podzolic sites
#     - This is the dominant statistical axis, but it reflects a
#       wetland-richness/paludification gradient, NOT a simple
#       temperate-vs-boreal or "good vs bad soil" story.
#   soil_pc2 ("contrast" axis, ~39% of the joint variance)
#     - loadings: ph_humus +0.707, drainage -0.707 (opposite signs, enforced)
#     - HIGH soil_pc2 = high pH AND good drainage together -- matches the
#       literature's description of conditions favorable to temperate
#       species (e.g. sugar maple) at their northern/boreal range limit
#     - LOW soil_pc2  = low pH AND poor drainage together -- matches the
#       documented soil stress (acidic + waterlogged) temperate species
#       face pushing into boreal soils (Carteron et al. 2020 and others)
#     - This axis explains less raw variance than PC1, but is the one
#       that actually lines up with the species-performance mechanism.
#
# Output: data/build/soil_covariates.rds
# ==============================================================================

library(sf)
library(data.table)
library(here)

source(here::here("R", "config.R"))   # project paths (data/build etc.)
PEP_PATH <- here("data", "raw", "pep", "PEP.gpkg")

# ==============================================================================
# STEP 1. LOAD AND RECODE
# ==============================================================================

cat("Loading data...\n")

placette <- as.data.table(st_read(PEP_PATH, layer = "placette_mes", quiet = TRUE))
sol      <- as.data.table(st_read(PEP_PATH, layer = "station_sol",  quiet = TRUE))
station  <- as.data.table(st_read(PEP_PATH, layer = "station_pe",   quiet = TRUE))

placette[, year := as.integer(format(date_sond, "%Y"))]

sol     <- sol[, .(id_pe, id_pe_mes, ph_humus)]
station <- station[, .(id_pe, id_pe_mes, cl_drai)]

# Drainage codes (confirmed against the official PEP data dictionary,
# CLASSE_DE_DRAINAGE): 00=Excessif ... 60-64=Tres mauvais, in steps of 10,
# with a units digit modifier (lateral drainage / frozen horizon / etc)
# that doesn't change the severity class. Code 16 ("Complexe") is NOT part
# of this ordinal scale -- it means drainage varies within the plot and
# can't be assigned a single class, so it's deliberately left NA here
# rather than naively integer-divided into class 1.
drai_ordinal <- function(code) {
  c <- suppressWarnings(as.integer(code))
  fcase(
    is.na(c),                 NA_integer_,
    c == 0,                   0L,
    c >= 10 & c <= 14,        1L,
    c >= 20 & c <= 24,        2L,
    c >= 30 & c <= 34,        3L,
    c >= 40 & c <= 44,        4L,
    c >= 50 & c <= 54,        5L,
    c >= 60 & c <= 64,        6L,
    default = NA_integer_
  )
}
station[, drainage_ordinal := drai_ordinal(cl_drai)]
station[, cl_drai := NULL]

cat(sprintf("  placette_mes: %d rows (%d plots)\n", nrow(placette), uniqueN(placette$id_pe)))
cat(sprintf("  station_sol:  %d rows\n", nrow(sol)))
cat(sprintf("  station_pe:   %d rows\n", nrow(station)))

# ==============================================================================
# STEP 2. MERGE ONTO EVERY PLOT VISIT
# ==============================================================================

dt <- placette[, .(id_pe, no_mes, id_pe_mes, year)]
dt <- merge(dt, sol,     by = c("id_pe", "id_pe_mes"), all.x = TRUE)
dt <- merge(dt, station, by = c("id_pe", "id_pe_mes"), all.x = TRUE)

setorder(dt, id_pe, no_mes)

soil_vars <- c("ph_humus", "drainage_ordinal")

dt_raw <- copy(dt[, c("id_pe", "id_pe_mes", "no_mes", soil_vars), with = FALSE])
setnames(dt_raw, soil_vars, paste0(soil_vars, "_raw"))

# ==============================================================================
# STEP 3. LAST-REAL-VALUE BROADCAST (per plot, per variable)
# ==============================================================================

last_real <- function(x) {
  observed <- x[!is.na(x)]
  if (length(observed) == 0L) return(x)
  rep(observed[length(observed)], length(x))
}

dt[, (soil_vars) := lapply(.SD, last_real), by = id_pe, .SDcols = soil_vars]

cat("\nCoverage BEFORE / AFTER last-real-value broadcast:\n")
for (v in soil_vars) {
  before <- mean(!is.na(dt_raw[[paste0(v, "_raw")]])) * 100
  after  <- mean(!is.na(dt[[v]])) * 100
  cat(sprintf("  %-18s before %5.1f%%  after %5.1f%%\n", v, before, after))
}

# ==============================================================================
# STEP 4. DIAGNOSTIC -- DOES THE ADOPTED (LAST REAL) VALUE CONTRADICT AN
#          EARLIER CONSENSUS?
# ==============================================================================

diag_results <- data.table(variable = character(), n_multi_real = integer(),
                            n_weird = integer(), pct_weird = numeric())

for (v in soil_vars) {
  raw_col <- paste0(v, "_raw")
  sub <- dt_raw[!is.na(get(raw_col)), .(id_pe, val = get(raw_col))]
  per_plot <- sub[, .(vals = list(val)), by = id_pe]
  per_plot <- per_plot[lengths(vals) >= 2]

  is_weird <- vapply(per_plot$vals, function(vals) {
    n        <- length(vals)
    earlier  <- vals[seq_len(n - 1)]
    last_val <- vals[n]
    length(unique(earlier)) == 1 && earlier[1] != last_val
  }, logical(1))

  diag_results <- rbind(diag_results, data.table(
    variable     = v,
    n_multi_real = nrow(per_plot),
    n_weird      = sum(is_weird),
    pct_weird    = if (nrow(per_plot) > 0) round(100 * sum(is_weird) / nrow(per_plot), 1) else NA_real_
  ))
}

cat("\nPlots where the adopted (last real) value contradicts an earlier-agreeing value:\n")
print(diag_results)

# ==============================================================================
# STEP 5. PCA -> soil_pc1, soil_pc2
#
# Fit PCA on the complete-case rows (both variables non-NA after the fill),
# then enforce a deterministic sign convention -- prcomp's sign is
# arbitrary (an artifact of the underlying SVD), so without this the
# reported loadings could flip from run to run with no substantive change.
# Convention: ph_humus always loads POSITIVELY on both PC1 and PC2. Given
# ph_humus and drainage are positively correlated, this automatically
# forces drainage to also load positively on PC1 (the "agreement" axis)
# and negatively on PC2 (the "contrast" axis) -- see header for what each
# axis means.
# ==============================================================================

complete <- dt[!is.na(ph_humus) & !is.na(drainage_ordinal)]
pca_fit  <- prcomp(complete[, .(ph_humus, drainage_ordinal)], scale. = TRUE)

if (pca_fit$rotation["ph_humus", "PC1"] < 0) pca_fit$rotation[, "PC1"] <- -pca_fit$rotation[, "PC1"]
if (pca_fit$rotation["ph_humus", "PC2"] < 0) pca_fit$rotation[, "PC2"] <- -pca_fit$rotation[, "PC2"]

cat("\nPCA loadings (sign-enforced: ph_humus positive on both axes):\n")
print(pca_fit$rotation)
cat("\nVariance explained:\n")
print(summary(pca_fit)$importance)

# project ALL complete-case rows using the (possibly sign-flipped) loadings
center <- pca_fit$center
scale_ <- pca_fit$scale
rot    <- pca_fit$rotation
z <- scale(complete[, .(ph_humus, drainage_ordinal)], center = center, scale = scale_)
scores <- z %*% rot

dt[, soil_pc1 := NA_real_]
dt[, soil_pc2 := NA_real_]
dt[!is.na(ph_humus) & !is.na(drainage_ordinal), `:=`(
  soil_pc1 = scores[, "PC1"],
  soil_pc2 = scores[, "PC2"]
)]

cat(sprintf("\nsoil_pc1: n=%d, range [%.2f, %.2f]\n", sum(!is.na(dt$soil_pc1)), min(dt$soil_pc1, na.rm=TRUE), max(dt$soil_pc1, na.rm=TRUE)))
cat(sprintf("soil_pc2: n=%d, range [%.2f, %.2f]\n", sum(!is.na(dt$soil_pc2)), min(dt$soil_pc2, na.rm=TRUE), max(dt$soil_pc2, na.rm=TRUE)))

# ==============================================================================
# STEP 5b. NEAREST-DONOR FALLBACK ON THE PC AXES
#
# Any row still missing soil_pc1/soil_pc2 at this point is missing it for
# ONE reason only: ph_humus and/or drainage_ordinal wasn't available for
# that plot even after the last-real-value broadcast (step 3), so PCA in
# step 5 had nothing to project. This covers two different plot
# populations at once -- plots missing BOTH raw variables entirely, and
# plots missing exactly ONE of them -- without needing separate fallback
# logic for each, because both simply show up here as "no PC score yet".
#
# Fallback: borrow the nearest plot's own real PC scores directly, rather
# than imputing the raw variable(s) and re-running PCA. These are
# equivalent when a donor plot's own (ph_humus, drainage_ordinal) pair is
# used, since the PCA projection is a fixed linear transform -- but
# working at the PC-score level means one fallback mechanism handles
# every remaining-NA case uniformly. Distance uses the "placette" layer's
# point geometry (one row per id_pe, already projected & metric --
# NAD83 / Quebec Lambert, EPSG:32198), so Euclidean distance is meters,
# no reprojection needed. Every borrowed row is flagged (soil_pc_source /
# soil_pc_donor_id_pe / soil_pc_donor_dist_km) for a sensitivity check
# once transition models exist: refit with these plots' PC scores
# imputed vs. with their rows dropped entirely, and compare the soil
# coefficient.
# ==============================================================================

dt[, soil_pc_source := "measured"]
dt[, soil_pc_donor_id_pe := NA_character_]
dt[, soil_pc_donor_dist_km := NA_real_]

pc_missing_plots <- dt[, .(all_na = all(is.na(soil_pc1))), by = id_pe][all_na == TRUE, id_pe]
cat(sprintf("\nPlots with NO PC score on any visit (still missing after step 5): %d\n",
            length(pc_missing_plots)))

if (length(pc_missing_plots) > 0) {
  coords <- st_read(PEP_PATH, layer = "placette", quiet = TRUE)[, c("id_pe", "geom")]

  donor_ids     <- dt[!is.na(soil_pc1), unique(id_pe)]
  donor_coords  <- coords[coords$id_pe %in% donor_ids, ]
  missing_coords <- coords[coords$id_pe %in% pc_missing_plots, ]

  dist_mat <- st_distance(missing_coords, donor_coords)  # meters (rows=missing, cols=donor)

  nearest_idx  <- apply(dist_mat, 1, which.min)
  nearest_dist <- apply(dist_mat, 1, min)

  donor_lookup <- unique(dt[!is.na(soil_pc1), .(id_pe, soil_pc1, soil_pc2)], by = "id_pe")

  fill_map <- data.table(
    id_pe         = missing_coords$id_pe,
    donor_id_pe   = donor_coords$id_pe[nearest_idx],
    donor_dist_km = as.numeric(nearest_dist) / 1000
  )
  fill_map <- merge(fill_map, donor_lookup, by.x = "donor_id_pe", by.y = "id_pe")

  cat(sprintf("Nearest-donor distance (km): median %.2f, max %.2f\n",
              median(fill_map$donor_dist_km), max(fill_map$donor_dist_km)))

  dt[fill_map, on = "id_pe", `:=`(
    soil_pc1              = i.soil_pc1,
    soil_pc2               = i.soil_pc2,
    soil_pc_source          = "nearest_donor",
    soil_pc_donor_id_pe     = i.donor_id_pe,
    soil_pc_donor_dist_km   = i.donor_dist_km
  )]
}

cat(sprintf("\nsoil_pc1 AFTER PC-axis fallback: n=%d / %d non-NA\n", sum(!is.na(dt$soil_pc1)), nrow(dt)))
cat(sprintf("soil_pc2 AFTER PC-axis fallback: n=%d / %d non-NA\n", sum(!is.na(dt$soil_pc2)), nrow(dt)))
cat(sprintf("Rows filled via PC-axis nearest donor: %d / %d\n",
            sum(dt$soil_pc_source == "nearest_donor"), nrow(dt)))


# ==============================================================================
# STEP 6. SAVE
# ==============================================================================

out <- dt[, .(id_pe, no_mes, id_pe_mes, year, ph_humus, drainage_ordinal,
              soil_pc1, soil_pc2, soil_pc_source, soil_pc_donor_id_pe, soil_pc_donor_dist_km)]

dir.create(BUILD_DIR, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(BUILD_DIR, "soil_covariates.rds")
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
saveRDS(out, out_path)
cat(sprintf("\nSaved %d rows to %s\n", nrow(out), out_path))
