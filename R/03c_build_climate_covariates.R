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
# NOTE 2026-10-24: the model uses GROWING-SEASON-MATCHED CMI (cmi_gs, built here in step 4b): monthly CMI weighted by the share of each
#   month inside the plot-year's NRCan growing season [start_gs, end_gs] (= the window of gdd_3), then averaged like every column.
#   Replaces the fixed May-Oct sum of 2026-10-22 (plot-level r = 0.95 with it; correlation with gdd_3 -0.61 instead of -0.80).
#   cmi_sum (annual) is still written here (raw product) but dropped in 03z. Brice et al. (2020) used May-Sept CMI.
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

source(here::here("R", "config.R"))   # project paths (data/build etc.)
RAW_DIR   <- here("data", "raw", "climate_raw")    # scratch, deleted per-year after extraction
dir.create(BUILD_DIR, showWarnings = FALSE, recursive = TRUE)
CKPT_PATH <- file.path(BUILD_DIR, "climate_by_year_raw.rds")
OUT_PATH  <- file.path(BUILD_DIR, "climate_covariates.rds")

dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUT_PATH), recursive = TRUE, showWarnings = FALSE)

options(timeout = 600)

sg_names  <- c("start_gs", "end_gs", "length_gs",
               "pcp_1", "pcp_2", "pcp_3", "pcp_4",
               "gdd_1", "gdd_2", "gdd_3", "gdd_4",
               "an_meanT", "an_minT", "an_maxT", "an_meanT_3", "an_range_T_3")
cmi_month_names <- paste0("cmi_", 1:12)

## ---- 1. plot-visits that matter -----------------------------------------

sc <- readRDS(file.path(BUILD_DIR, "state_classification.rds"))
dc <- readRDS(file.path(BUILD_DIR, "disturbance_covariates.rds"))

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

# A1 (2026-11-04): the checkpoint is checked per PLOT-YEAR, not per year. Plots added to the analysis after the checkpoint was
# saved (5 plots on 2026-11-03) are extracted for the years their windows need; plots/years already in the checkpoint are not touched.
# Safety copies of the old files (made once; never overwritten): *.pre_A1.rds, used to verify that old visits are unchanged.
for (f in c(CKPT_PATH, OUT_PATH)) {
  bak <- sub("\\.rds$", ".pre_A1.rds", f)
  if (file.exists(f) && !file.exists(bak)) {
    stopifnot("could not copy the pre-A1 backup" = file.copy(f, bak, overwrite = FALSE))
    cat(sprintf("Backup written: %s\n", basename(bak)))
  }
}

if (file.exists(CKPT_PATH)) {
  long_dt <- readRDS(CKPT_PATH)
  cat(sprintf("Resuming: checkpoint has %d years, %d plots, %d plot-years.\n",
              uniqueN(long_dt$year), uniqueN(long_dt$id_pe), nrow(long_dt)))
} else {
  long_dt <- data.table()
}

stopifnot("duplicate id_pe in PEP placette layer" = !anyDuplicated(pep$id_pe),
          "visits with a plot that has no location in PEP.gpkg" = all(visits$id_pe %in% pep$id_pe))

# every plot-year the visit windows need (step 2), minus what the checkpoint holds
need <- unique(visits[win_start <= win_end, .(year = seq.int(win_start, win_end)), by = .(id_pe, no_mes)][, .(id_pe, year)])
miss <- if (nrow(long_dt) == 0L) need else need[!long_dt[, .(id_pe, year)], on = c("id_pe", "year")]
years_todo <- sort(unique(miss$year))
cat(sprintf("Plot-years needed: %d | missing from the checkpoint: %d (%d plots, %d years)\n",
            nrow(need), nrow(miss), uniqueN(miss$id_pe), length(years_todo)))
if (nrow(miss) > 0) print(miss[, .(n_years = .N, from = min(year), to = max(year)), by = id_pe][order(id_pe)])

for (yr in years_todo) {
  ids_yr <- miss[year == yr, id_pe]
  pep_yr <- pep[pep$id_pe %in% ids_yr, ]
  stopifnot(nrow(pep_yr) == length(ids_yr))
  pts_ll <- vect(pep_yr[, "id_pe"])           # ONLY the plots missing in this year
  cat(sprintf("\n== year %d (%d plots) ==\n", yr, length(ids_yr)))

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
    id_pe = pep_yr$id_pe,
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
stopifnot("duplicate (id_pe, year) in the checkpoint" = !anyDuplicated(long_dt, by = c("id_pe", "year")))

## ---- 4b. growing-season-matched CMI per plot-year (decided 2026-10-24) --
## cmi_gs = sum_m cmi_m x (days of month m inside [start_gs, end_gs]) / (days in month m).
## start_gs / end_gs = NRCan growing season (period 3 = the window of gdd_3), inclusive day-of-year, per plot-year.
## Month boundaries from each year's own calendar (in leap years every day-of-year after February shifts by one).
## Assumption: CMI is spread evenly within a partial month (only the season's first and last months are partial).
## Months with weight 0 do not count, so an NA in a month outside the season does not make cmi_gs NA.
## Computed BEFORE the 10-yr averaging: step 5 averages cmi_gs like every other column. The checkpoint file is not modified.
## Tested 2026-10-24 on all 717,835 plot-years (see PIPELINE_NOTES): weighted days == length_gs in every row; hand example reproduced.
cmi_gs_matched <- function(year, start_gs, end_gs, cmi_mat) {
  stopifnot(ncol(cmi_mat) == 12L, length(year) == nrow(cmi_mat), length(start_gs) == length(year), length(end_gs) == length(year))
  leap  <- (year %% 4L == 0L & year %% 100L != 0L) | year %% 400L == 0L
  md    <- rbind(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31),
                 c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))
  mdays  <- md[leap + 1L, , drop = FALSE]                                # n x 12
  mend   <- t(apply(md, 1, cumsum))[leap + 1L, , drop = FALSE]
  mstart <- mend - mdays + 1
  inside <- pmax(pmin(mend, end_gs) - pmax(mstart, start_gs) + 1, 0)    # matrix first: keeps dim; vectors recycle per row
  w <- inside / mdays                                                    # NA where start/end is NA
  cm <- as.matrix(cmi_mat); cm[!is.na(w) & w == 0] <- 0
  list(cmi_gs = rowSums(cm * w), w = w)
}
gsm <- cmi_gs_matched(long_dt$year, long_dt$start_gs, long_dt$end_gs, long_dt[, paste0("cmi_", 1:12), with = FALSE])
long_dt[, cmi_gs := gsm$cmi_gs]
cat(sprintf("\nGS-matched CMI: %d plot-years, %d NA (missing season dates or in-season CMI); %.1f%% of plot-years have season days outside May-Oct.\n",
            nrow(long_dt), sum(is.na(long_dt$cmi_gs)), 100 * mean(rowSums(gsm$w[, c(1:4, 11:12), drop = FALSE]) > 0, na.rm = TRUE)))
rm(gsm)

## A1 guard (2026-11-04): every needed plot-year must be present (a failed download is not skipped silently) and have no NA in the
## averaged columns (not only gdd_3 / cmi_gs). mean(na.rm = TRUE) in step 5 would otherwise average fewer years without a signal.
used <- long_dt[need, on = c("id_pe", "year"), nomatch = NULL]
if (nrow(used) != nrow(need)) {
  print(need[!long_dt[, .(id_pe, year)], on = c("id_pe", "year")][, .(n_years = .N, from = min(year), to = max(year)), by = id_pe])
  stop(sprintf("A1: %d needed plot-years are missing from the extraction (failed download?). Rerun 03c.", nrow(need) - nrow(used)))
}
vcols <- setdiff(names(long_dt), c("id_pe", "year"))     # every column that step 5 averages (16 sg, cmi_1-12, cmi_sum, cmi_gs)
na_by_col <- colSums(is.na(used[, ..vcols]))
if (any(na_by_col > 0)) {
  cat("NA per column among the used plot-years:\n"); print(na_by_col[na_by_col > 0])
  bad_py <- used[rowSums(is.na(used[, ..vcols])) > 0, .(id_pe, year)]
  print(head(bad_py, 30))
  stop(sprintf("A1: %d used plot-years have NA in the averaged columns (%d plots).", nrow(bad_py), uniqueN(bad_py$id_pe)))
}
cat(sprintf("A1 guard: all %d needed plot-years present, 0 NA in the %d averaged columns.\n", nrow(used), length(vcols)))
rm(used)

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

## A1 guard (2026-11-04): stop unless EVERY visit used its full (already 1960-2022-clipped) window.
bad <- out[n_years_used != (win_end - win_start + 1L) | n_years_used < 1L]
if (nrow(bad) > 0) {
  print(head(bad[, .(id_pe, no_mes, year, win_start, win_end, n_years_used)], 30))
  stop(sprintf("A1: %d plot-visits do not have n_years_used == win_end - win_start + 1 (>= 1). climate_covariates.rds NOT written.", nrow(bad)))
}
cat(sprintf("A1 guard: all %d plot-visits have n_years_used == window length.\n", nrow(out)))
na_out <- colSums(is.na(out[, ..value_cols]))
if (any(na_out > 0)) { print(na_out[na_out > 0]); stop("A1: NA in the averaged columns of the visit table. climate_covariates.rds NOT written.") }
cat(sprintf("A1 guard: 0 NA in the %d averaged columns of the %d visits.\n", length(value_cols), nrow(out)))

saveRDS(out, OUT_PATH)
cat(sprintf("\nSaved: %s\n", OUT_PATH))
