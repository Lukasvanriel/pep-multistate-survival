# ==============================================================================
# 02b_stand_age.R -- stand age at the first visit (W_j) per plot.
# Moved 2026-10-26 from 04_prepare_inla.R step 1 (user decision), logic UNCHANGED, so that 03 can use the stand age
# (TSD of stand-origin events undated at the first visit) and 04 reads the same table.
# Only difference with the old 04 step 1: the visits come from 02 (state_classification.rds + visit year from
# placette_mes) instead of combined_covariates.rds, because 02b runs before 03/03z. Checked against the W_j table of the
# inttest run (see docs/PIPELINE_NOTES.md).
# cl_age source: pee_ori_sond first, gaps filled from station_pe (per visit).
# Rules: 1 first visit | 2 later visit, back-calculated (never 0 unless first state is regen) | 3 time since
# fire/harvest at the first visit | 4 first state regen AND later exits regen -> 0 | otherwise EXCLUDED.
# Output: data/build/stand_age.rds (id_pe, year_1, state_1, W_j, rule, any_exit, include)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(sf); library(here) })
source(here::here("R", "config.R"))
PEP_PATH <- here("data", "raw", "pep", "PEP.gpkg")

sc <- as.data.table(readRDS(file.path(BUILD_DIR, "state_classification.rds")))
pl <- as.data.table(st_read(PEP_PATH, layer = "placette_mes", quiet = TRUE))[, .(id_pe, no_mes = as.integer(no_mes), year = as.integer(format(date_sond, "%Y")))]
dt <- merge(sc[, .(id_pe, no_mes = as.integer(no_mes), state)], pl, by = c("id_pe", "no_mes"))
setorder(dt, id_pe, no_mes)
cat(sprintf("02 visits with a year: %d visits, %d plots\n\n", nrow(dt), uniqueN(dt$id_pe)))

# ---- from here: 04_prepare_inla.R step 1, unchanged except the final save path ----------------------------------------
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

miss <- wj[is.na(W_j), id_pe]
r2 <- src[id_pe %in% miss & !is.na(cl_num), .SD[1L], by = id_pe][, .(id_pe, cl_later = cl_num, year_later = year)]
wj[, `:=`(W_r2 = NA_real_, cl_later = NA_real_)]
wj[r2, on = "id_pe", `:=`(W_r2 = i.cl_later - (i.year_later - year_1), cl_later = i.cl_later)]
n_nonpos <- wj[!is.na(W_r2) & W_r2 <= 0, .(regen = sum(state_1 == "regen"), non_regen = sum(state_1 != "regen"))]
wj[!is.na(W_r2) & W_r2 <= 0, W_r2 := fifelse(state_1 == "regen", 0, cl_later)]
cat(sprintf("Rule 2 back-calculation <= 0: %d plots with first state regen -> W_j 0; %d non-regen -> later class value\n",
            n_nonpos$regen, n_nonpos$non_regen))
wj[is.na(W_j) & !is.na(W_r2), `:=`(W_j = W_r2, rule = 2L)]

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

ex <- dt[, .(any_exit = any(state[-.N] == "regen" & state[-1L] != "regen")), by = id_pe]
wj[ex, on = "id_pe", any_exit := i.any_exit]
wj[is.na(W_j) & state_1 == "regen" & any_exit %in% TRUE, `:=`(W_j = 0, rule = 4L)]

wj[, include := !is.na(W_j)]
wj[, c("W_r2", "cl_later", "W_r3") := NULL]

cat("W_j rule breakdown:\n"); print(wj[, .N, by = .(rule, include)][order(rule)])
cat(sprintf("\nIncluded: %d / %d plots (%.2f%%)\n", sum(wj$include), nrow(wj), 100 * mean(wj$include)))
saveRDS(wj, file.path(BUILD_DIR, "stand_age.rds"))
cat("Saved", file.path(BUILD_DIR, "stand_age.rds"), "\n")
