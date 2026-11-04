# ==============================================================================
# R/tests/03x_check_against_old.R -- check the new disturbance build against the old one on the old code's principles.
# P1 back-fill: an undated code takes the year of the first later dated code of its type that is possible (year <= visit).
# P2 a dated code shows its own year at its own visit.
# P3 a recorded code is flagged at its own visit.
# For every code where the NEW build does not follow the principle, the reason is classified by the agreed rule
# (from the event table of the new build); "UNEXPLAINED" must be inspected. Old-build counts are shown for comparison.
# Run from the project root after 03:  Rscript R/tests/03x_check_against_old.R [old_build.rds]
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here) })
source(here("R", "config.R")); source(here("R", "03a_disturbance_events.R"))
args <- commandArgs(TRUE)
OLD <- if (length(args)) args[1] else here("data", "build", "disturbance_covariates_build_9298c98.rds")
ex <- parse(here("R", "03_build_covariates.R"))
for (e in ex) { s <- paste(deparse(e), collapse = ""); if (grepl("^(fire_o|fire_p|harv_o|harv_p|pest_o|pest_p) <-", s)) eval(e) }
nw  <- as.data.table(readRDS(here("data", "build", "disturbance_covariates.rds"))); setkey(nw, id_pe, no_mes)
old <- as.data.table(readRDS(OLD)); setkey(old, id_pe, no_mes)
E <- readRDS(here("data", "build", "disturbance_events.rds")); ev <- E$events; rc <- E$records; vv <- E$visits; setkey(vv, id_pe, no_mes)
stopifnot(identical(nw[, .(id_pe, no_mes)], vv[, .(id_pe, no_mes)]))
# visits present in only one build (e.g. CHT/DT truncation keeping other plots) are reported and left out of the comparison
k_new <- nw[, paste(id_pe, no_mes)]; k_old <- old[, paste(id_pe, no_mes)]
cat(sprintf("visits only in the new build: %d (%d plots) | only in the old build: %d (%d plots)\n",
    sum(!k_new %in% k_old), uniqueN(nw$id_pe[!k_new %in% k_old]), sum(!k_old %in% k_new), uniqueN(old$id_pe[!k_old %in% k_new])))
keep <- k_new %in% k_old; nw <- nw[keep]; vv <- vv[keep]; old <- old[paste(id_pe, no_mes) %in% k_new]
stopifnot(identical(nw[, .(id_pe, no_mes)], old[, .(id_pe, no_mes)]))
nw <- cbind(nw, vv[, .(origine, an_origine, perturb, an_perturb, fire_event, harv_event, pest_event)])
nw[, vi := seq_len(.N), by = id_pe]

# raw codes -> records (type, level, year) exactly as step 0
lst <- list(fire_o = fire_o, fire_p = fire_p, harv_o = harv_o, harv_p = harv_p, pest_o = pest_o, pest_p = pest_p)
raw <- rbindlist(lapply(c("fire", "harv", "pest"), function(T) rbind(
  nw[!is.na(origine) & origine %in% lst[[paste0(T, "_o")]], .(id_pe, no_mes, vi, year, type = T, level = 2L, yr = as.integer(an_origine), code = origine)],
  nw[!is.na(perturb) & perturb %in% lst[[paste0(T, "_p")]], .(id_pe, no_mes, vi, year, type = T, level = 1L, yr = as.integer(an_perturb), code = perturb)])))
setorder(raw, id_pe, vi)
val <- function(b, T) { ts <- c(fire = "tsf", harv = "tsh", pest = "sfp")[[T]]
  b[, .(id_pe, no_mes, flag = get(T), real = as.numeric(get(ts)), ev_shown = if (paste0(T, "_event") %in% names(b)) get(paste0(T, "_event")) else NA_integer_)] }
look <- function(x, b) rbindlist(lapply(c("fire", "harv", "pest"), function(T) merge(x[type == T], val(b, T), by = c("id_pe", "no_mes"))))

# event of each record in the NEW build (records are after A1: one per visit and type)
rc2 <- rc[, .(id_pe, no_mes, type, rec_level = level, rec_yr = yr, yr_src, n_raw, event, action)]
ev2 <- ev[, .(id_pe, event, e_level = level, e_yr = yr, dropped, start_no_mes)]
classify <- function(z) {
  z[, zid := .I]
  z <- merge(z, rc2, by = c("id_pe", "no_mes", "type"), all.x = TRUE)
  # several records at one visit and type (different intervals): keep the record of this code (same year, else undated)
  z[, pri := fifelse(!is.na(yr) & !is.na(rec_yr) & rec_yr == yr, 1L, fifelse(is.na(yr) & action %like% "pair|continuation|undated|backfill", 2L, 3L))]
  setorder(z, zid, pri); z <- z[!duplicated(zid)][, c("zid", "pri") := NULL]
  z <- merge(z, ev2, by = c("id_pe", "event"), all.x = TRUE)
  z <- merge(z, nw[, .(id_pe, no_mes, ba_event)], by = c("id_pe", "no_mes"))
  z[, reason := fcase(
    !is.na(yr) & yr > year, "year after the visit (event after the measurement)",
    action == "pair_leftover", "pair: undated code is a leftover of a tracked event",
    n_raw > 1 & (is.na(yr) | rec_yr != yr | rec_level != level), "A1: same visit, same interval -> highest intensity",
    dropped %in% TRUE, "stand-origin selection: event dropped",
    !is.na(e_yr) & !is.na(rec_yr) & (e_yr != rec_yr | e_level != rec_level), "A2: later record in the same interval (highest / higher year)",
    !is.na(e_yr) & e_yr == year & !ba_event & flag == 0, "same-year rule: event after the measurement",
    flag == 0 & !is.na(start_no_mes) & start_no_mes <= no_mes, "reset: ended by a stand-origin event known to be later",
    !is.na(ev_shown) & ev_shown != event, "C(i): a more recent event of this type is shown",
    default = "UNEXPLAINED")]
  z
}

# ---- P1 back-fill -----------------------------------------------------------------------------------------------
u <- raw[is.na(yr)]; d <- raw[!is.na(yr), .(id_pe, type, dvi = vi, dyr = yr)]
m <- merge(u, d, by = c("id_pe", "type"), allow.cartesian = TRUE)[dvi > vi & dyr <= year][order(id_pe, type, vi, dvi)]
m <- unique(m, by = c("id_pe", "type", "vi", "level"))[, expect := year - dyr]
for (nm in c("old", "new")) { z <- look(m, if (nm == "old") old else nw)
  cat(sprintf("P1 %s: %d undated codes with a later possible dated code | shows exactly that year: %d | another year: %d | no year: %d | not flagged: %d\n",
      nm, nrow(z), z[flag > 0 & !is.na(real) & real == expect, .N], z[flag > 0 & !is.na(real) & real != expect, .N], z[flag > 0 & is.na(real), .N], z[flag == 0, .N])) }
z <- look(m, nw)[!(flag > 0 & !is.na(real) & real == expect)]
z <- classify(z); cat("  new build, codes not showing that year -- reasons:\n"); print(z[, .(codes = .N, plots = uniqueN(id_pe)), by = reason][order(-codes)], row.names = FALSE)
P1_unexpl <- z[reason == "UNEXPLAINED"]

# ---- P2 dated code shows its own year at its own visit -----------------------------------------------------------
dd <- raw[!is.na(yr)][, expect := year - yr]
for (nm in c("old", "new")) { z <- look(dd, if (nm == "old") old else nw)
  cat(sprintf("\nP2 %s: %d dated codes | own visit shows that year: %d | another year: %d | no year: %d | not flagged: %d\n",
      nm, nrow(z), z[flag > 0 & !is.na(real) & real == expect, .N], z[flag > 0 & !is.na(real) & real != expect, .N], z[flag > 0 & is.na(real), .N], z[flag == 0, .N])) }
z <- look(dd, nw)[!(flag > 0 & !is.na(real) & real == expect)]
z <- classify(z); cat("  new build, dated codes not showing their year -- reasons:\n"); print(z[, .(codes = .N, plots = uniqueN(id_pe)), by = reason][order(-codes)], row.names = FALSE)
P2_unexpl <- z[reason == "UNEXPLAINED"]

# ---- P3 recorded code flagged at its own visit ----------------------------------------------------------------------
for (nm in c("old", "new")) { z <- look(raw, if (nm == "old") old else nw)
  cat(sprintf("\nP3 %s: %d codes | flagged at their visit: %d | not flagged: %d\n", nm, nrow(z), z[flag > 0, .N], z[flag == 0, .N])) }
z <- look(raw, nw)[flag == 0]; z <- classify(z)
cat("  new build, codes not flagged at their visit -- reasons:\n"); print(z[, .(codes = .N, plots = uniqueN(id_pe)), by = reason][order(-codes)], row.names = FALSE)
P3_unexpl <- z[reason == "UNEXPLAINED"]

un <- unique(rbind(P1_unexpl[, .(id_pe, type, no_mes, check = "P1")], P2_unexpl[, .(id_pe, type, no_mes, check = "P2")], P3_unexpl[, .(id_pe, type, no_mes, check = "P3")]))
cat(sprintf("\nUNEXPLAINED in total: %d codes in %d plots\n", nrow(un), uniqueN(un$id_pe)))
if (nrow(un)) { cat("examples:\n"); print(head(un, 15)) }
