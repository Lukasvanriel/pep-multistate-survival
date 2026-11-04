# ==============================================================================
# 03c_build_climate_covariates.R
#
# Builds climate covariates (gdd, cmi) for PEP permanent plots, as two
# separate variable families (NOT PCA-condensed like the soil covariates in
# 03b). Companion to 03b_build_soil_climate_covariates.R, kept in its own
# script per instruction.
#
# Source: NRCan / McKenney et al. (2011, BAMS 92(12):1611-1622) ANUSPLIN
# gridded North American monthly climate data, 60 arcsec (~2km) resolution,
# hosted at:
#   https://ftp.maps.canada.ca/pub/nrcan_rncan/Climate-archives_Archives-climatologiques/NAM_monthly/monthly_by_var/
#
# Two product families pulled from that mirror, confirmed live 1960-2022
# (404 for 2023 onward, checked directly against the server):
#   - sg{year}_60arcsec.zip  -> 16 vars: start_gs, end_gs, length_gs,
#       pcp_1-4, gdd_1-4, an_meanT, an_minT, an_maxT, an_meanT_3,
#       an_range_T_3. Documented at cfs.nrcan.gc.ca/projects/3/9.
#       gdd_1-4 = growing-degree-days over 4 periods relative to a defined
#       growing season (base temperature not published anywhere we could
#       find; inferred ~5C by consistency with the documented growing-season
#       onset threshold -- this is a caveat, not a confirmed fact).
#   - cmi{year}_60arcsec.zip -> cmi_1..cmi_12 (monthly) + cmi_sum (annual).
#       Bundled documentation (cmi60_en.txt / cmi60.xml inside the zip)
#       cites the same McKenney et al. (2011) + McKenney et al. (2013,
#       Forestry Chronicle 89(5):659-663) as the source. CMI = monthly
#       precip minus potential evapotranspiration (simplified
#       Penman-Monteith); positive = wetter than average, negative = drier.
#
# NOTE 2026-10-22: the model uses growing-season CMI (cmi_gs = sum of cmi_5..cmi_10, built in 03z), NOT cmi_sum (annual).
#   cmi_sum is still written here (raw product) but dropped in 03z. Brice et al. (2020) used May-Sept CMI (r = 0.997 with May-Oct).
#
# Methodology (matches Brice et al. 2020, Global Change Biology, methods
# text -- see numericalecology.com/Reprints/Brice_et_al_Global_Change_
# Biology_2020.pdf): each climate variable is averaged over the 10-year
# period prior to each plot's measurement year, to reduce the effect of
# inter-annual variability. NOTE: Brice et al. themselves only use mean
# growing-season temperature + May-Sept CMI, not gdd -- gdd here is this
# project's own addition, independent of that precedent.
#
# Known limitation, decided explicitly (not defaulted into): server data
# stops at 2022. Plot measurements in 2023-2025 (whose 10-yr trailing
# window would need 2023/2024) get a TRUNCATED window instead -- i.e. fewer
# than 10 real years averaged, never a repeated/duplicated year. This is
# recorded per plot-visit via n_years_used so it's fully auditable
# downstream, never silently hidden.
#
# Planned steps:
#   1. Determine the full set of plot-visits (id_pe, no_mes, year) that
#      matter -- i.e. those retained in state_classification.rds
#   2. Determine the set of calendar years actually needed (union of each
#      visit's truncated 10-yr trailing window)
#   3. For each needed year: download sg + cmi zips (skip if already done),
#      extract all variables at every plot location, append to a long-format
#      checkpoint table (id_pe x year x variable), delete the raw rasters
#      once extracted (they're large; only the extracted values are kept)
#   4. For each plot-visit: average each variable over its own truncated
#      window, from the long-format table
#   5. Save final wide covariate table
# ==============================================================================

library(sf)
library(terra)
library(curl)
library(data.table)
library(here)

## ---- 0. config ---------------------------------------------------------

BASE_URL <- "https://ftp.maps.canada.ca/pub/nrcan_rncan/Climate-archives_Archives-climatologiques/NAM_monthly/monthly_by_var"
MIN_YEAR_AVAILABLE <- 1960L   # confirmed via HEAD-check (1958-1962 all HTTP 200)
MAX_YEAR_AVAILABLE <- 2022L   # confirmed via HEAD-check (2023-2025 all HTTP 404)
TRAIL_WINDOW <- 10L           # years prior to measurement, per Brice et al. methods

RAW_DIR   <- here("data", "climate_raw")           # scratch, deleted per-year after extraction
CKPT_PATH <- here("data", "covariates", "pep", "climate_by_year_raw.rds")
OUT_PATH  <- here("data", "covariates", "pep", "climate_covariates.rds")

dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUT_PATH), recursive = TRUE, showWarnings = FALSE)

options(timeout = 600)

sg_names  <- c("start_gs", "end_gs", "length_gs",
               "pcp_1", "pcp_2", "pcp_3", "pcp_4",
               "gdd_1", "gdd_2", "gdd_3", "gdd_4",
               "an_meanT", "an_minT", "an_maxT", "an_meanT_3", "an_range_T_3")
cmi_month_names <- paste0("cmi_", 1:12)

## ---- 1. plot-visits that matter -----------------------------------------

sc <- readRDS(here("data", "processed", "state_classification.rds"))
dc <- readRDS(here("data", "covariates", "pep", "disturbance_covariates.rds"))

stopifnot("year" %in% colnames(dc))

visits <- merge(
  sc[, .(id_pe, no_mes)],
  dc[, .(id_pe, no_mes, year)],
  by = c("id_pe", "no_mes"),
  all.x = TRUE
)
visits <- as.data.table(visits)

n_missing_year <- sum(is.na(visits$year))
if (n_missing_year > 0) {
  cat(sprintf(
    "*** FLAG: %d of %d plot-visits in state_classification.rds have no matching year in disturbance_covariates.rds -- these will be dropped from the climate covariate output. Inspect before proceeding if this number is large. ***\n",
    n_missing_year, nrow(visits)
  ))
}
visits <- visits[!is.na(year)]

cat(sprintf("Plot-visits needing climate covariates: %d (%d distinct plots)\n",
            nrow(visits), uniqueN(visits$id_pe)))
cat("Measurement year range:", paste(range(visits$year), collapse = " - "), "\n")

## ---- 2. years actually needed -------------------------------------------

# truncated trailing window per visit (option (a): fewer real years near the
# data ceiling, never a duplicated/reused year)
visits[, win_start := pmax(MIN_YEAR_AVAILABLE, year - TRAIL_WINDOW)]
visits[, win_end   := pmin(MAX_YEAR_AVAILABLE, year - 1L)]

if (any(visits$win_end < visits$win_start)) {
  bad <- visits[win_end < win_start]
  cat(sprintf(
    "*** FLAG: %d plot-visits have an empty trailing window (measurement year too early/late relative to data availability 1960-2022). Their climate covariates will be NA. ***\n",
    nrow(bad)
  ))
  print(bad[, .(id_pe, no_mes, year)])
}

years_needed <- sort(unique(unlist(
  Map(function(a, b) if (a <= b) a:b else integer(0), visits$win_start, visits$win_end)
)))
cat(sprintf("Distinct calendar years needed: %d (%d to %d)\n",
            length(years_needed), min(years_needed), max(years_needed)))

## ---- 3. plot locations (fixed across years) ------------------------------

pep <- st_read(here("data", "raw", "pep", "PEP.gpkg"), layer = "placette")
pep <- pep[pep$id_pe %in% unique(visits$id_pe), ]
cat(sprintf("Plot locations loaded: %d\n", nrow(pep)))

## ---- 4. per-year download + extract, with checkpointing -----------------

download_zip <- function(info, year) {
  url <- sprintf("%s/%s%d_60arcsec.zip", BASE_URL, info, year)
  dest_dir <- file.path(RAW_DIR, info, year)
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  zip_path <- file.path(dest_dir, "data.zip")
  if (!file.exists(zip_path)) {
    ok <- tryCatch({
      download.file(url, zip_path, mode = "wb", method = "libcurl", quiet = TRUE)
      TRUE
    }, error = function(e) {
      cat(sprintf("  DOWNLOAD FAILED for %s %d: %s\n", info, year, conditionMessage(e)))
      FALSE
    })
    if (!ok || !file.exists(zip_path) || file.info(zip_path)$size == 0) return(NULL)
  }
  unzip(zip_path, exdir = dest_dir)
  dest_dir
}

load_sg <- function(dest_dir) {
  files <- list.files(dest_dir, pattern = "sg60_(0[1-9]|1[0-6])\\.tif$",
                       full.names = TRUE, recursive = TRUE)
  if (length(files) != 16) {
    cat(sprintf("  *** FLAG: expected 16 sg files, found %d in %s ***\n", length(files), dest_dir))
    return(NULL)
  }
  # order by the numeric suffix, not lexical sort (avoids 10 sorting before 2)
  ord <- order(as.integer(sub(".*sg60_(\\d+)\\.tif$", "\\1", files)))
  rast(files[ord])
}

load_cmi <- function(dest_dir) {
  m_files <- list.files(dest_dir, pattern = "cmi60_(0[1-9]|1[0-2])\\.tif$",
                         full.names = TRUE, recursive = TRUE)
  s_files <- list.files(dest_dir, pattern = "cmi60_sum\\.tif$",
                         full.names = TRUE, recursive = TRUE)
  if (length(m_files) != 12 || length(s_files) != 1) {
    cat(sprintf("  *** FLAG: expected 12 monthly + 1 sum cmi files, found %d + %d in %s ***\n",
                length(m_files), length(s_files), dest_dir))
    return(NULL)
  }
  ord <- order(as.integer(sub(".*cmi60_(\\d+)\\.tif$", "\\1", m_files)))
  list(monthly = rast(m_files[ord]), sum = rast(s_files))
}

# resume support: load whatever's already been extracted
if (file.exists(CKPT_PATH)) {
  long_dt <- readRDS(CKPT_PATH)
  done_years <- unique(long_dt$year)
  cat(sprintf("Resuming: %d years already in checkpoint.\n", length(done_years)))
} else {
  long_dt <- data.table()
  done_years <- integer(0)
}

pts_ll <- vect(pep[, "id_pe"])

for (yr in years_needed) {
  if (yr %in% done_years) next
  cat(sprintf("\n== year %d ==\n", yr))

  sg_dir  <- download_zip("sg", yr)
  cmi_dir <- download_zip("cmi", yr)

  if (is.null(sg_dir) || is.null(cmi_dir)) {
    cat(sprintf("  Skipping %d -- download failed for sg and/or cmi.\n", yr))
    next
  }

  r_sg  <- load_sg(sg_dir)
  r_cmi <- load_cmi(cmi_dir)

  if (is.null(r_sg) || is.null(r_cmi)) {
    cat(sprintf("  Skipping %d -- raster load failed.\n", yr))
    next
  }

  pts_sg  <- project(pts_ll, crs(r_sg))
  pts_cmi <- project(pts_ll, crs(r_cmi$monthly))

  vals_sg      <- extract(r_sg, pts_sg)
  vals_cmi_m   <- extract(r_cmi$monthly, pts_cmi)
  vals_cmi_sum <- extract(r_cmi$sum, pts_cmi)

  colnames(vals_sg)    <- c("ID", sg_names)
  colnames(vals_cmi_m) <- c("ID", cmi_month_names)
  colnames(vals_cmi_sum) <- c("ID", "cmi_sum")

  year_dt <- data.table(
    id_pe = pep$id_pe,
    year  = yr,
    vals_sg[, -1],
    vals_cmi_m[, -1],
    cmi_sum = vals_cmi_sum[, -1]
  )

  n_na <- sum(!complete.cases(year_dt))
  if (n_na > 0) {
    cat(sprintf("  *** FLAG: %d plots have at least one NA variable in %d ***\n", n_na, yr))
  }

  long_dt <- rbind(long_dt, year_dt)
  saveRDS(long_dt, CKPT_PATH)  # checkpoint after every year -- safe to interrupt/resume

  # done with the raw rasters for this year -- delete to save disk, the
  # extracted values are what's kept
  unlink(file.path(RAW_DIR, "sg", yr), recursive = TRUE)
  unlink(file.path(RAW_DIR, "cmi", yr), recursive = TRUE)

  cat(sprintf("  done: %d plots, %d NA-flagged, checkpoint saved (%d years total so far)\n",
              nrow(year_dt), n_na, uniqueN(long_dt$year)))
}

cat(sprintf("\nAll years processed. Checkpoint has %d years, %d rows.\n",
            uniqueN(long_dt$year), nrow(long_dt)))

## ---- 5. per-visit truncated trailing average -----------------------------

value_cols <- setdiff(colnames(long_dt), c("id_pe", "year"))

compute_visit_avg <- function(pe, y0, y1) {
  sub <- long_dt[id_pe == pe & year >= y0 & year <= y1]
  if (nrow(sub) == 0) {
    out <- as.list(setNames(rep(NA_real_, length(value_cols)), value_cols))
    out$n_years_used <- 0L
    return(out)
  }
  out <- as.list(sub[, lapply(.SD, mean, na.rm = TRUE), .SDcols = value_cols])
  out$n_years_used <- nrow(sub)
  out
}

cat("\nComputing per-visit trailing averages...\n")
visit_results <- vector("list", nrow(visits))
for (i in seq_len(nrow(visits))) {
  visit_results[[i]] <- compute_visit_avg(visits$id_pe[i], visits$win_start[i], visits$win_end[i])
}
visit_results_dt <- rbindlist(visit_results)

out <- cbind(visits[, .(id_pe, no_mes, year, win_start, win_end)], visit_results_dt)

cat(sprintf("\nFinal table: %d plot-visits\n", nrow(out)))
cat("n_years_used distribution:\n")
print(table(out$n_years_used))

n_incomplete <- sum(out$n_years_used < TRAIL_WINDOW)
cat(sprintf(
  "\n%d of %d plot-visits (%.1f%%) have a TRUNCATED window (< %d years) -- these are the ones close to either end of the 1960-2022 data availability. Expected and documented, not a bug.\n",
  n_incomplete, nrow(out), 100 * n_incomplete / nrow(out), TRAIL_WINDOW
))

n_all_na <- sum(out$n_years_used == 0)
if (n_all_na > 0) {
  cat(sprintf(
    "*** FLAG: %d plot-visits have n_years_used == 0 -- completely outside the 1960-2022 window. These have NA climate covariates. Inspect: ***\n",
    n_all_na
  ))
  print(out[n_years_used == 0, .(id_pe, no_mes, year)])
}

saveRDS(out, OUT_PATH)
cat(sprintf("\nSaved: %s\n", OUT_PATH))
