# WP1 / A1 check (read-only). Run from the project root AFTER 03c (and optionally 03z):  Rscript R/tests/wp1_A1_check.R
# Compares data/build/climate_covariates.rds with the pre-A1 copy (climate_covariates.pre_A1.rds, made by 03c).
suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("R", "config.R"))
new <- as.data.table(readRDS(file.path(BUILD_DIR, "climate_covariates.rds")))
old <- as.data.table(readRDS(file.path(BUILD_DIR, "climate_covariates.pre_A1.rds")))
five <- c("7109604001", "7209002301", "7209002302", "7209605702", "7501600102")
cat(sprintf("rows old %d / new %d | plots old %d / new %d | same columns: %s\n", nrow(old), nrow(new), uniqueN(old$id_pe), uniqueN(new$id_pe),
            identical(names(old), names(new))))
cat("NA gdd_3 / cmi_gs / start_gs / end_gs  old:", sapply(old[, .(gdd_3, cmi_gs, start_gs, end_gs)], function(x) sum(is.na(x))),
    " new:", sapply(new[, .(gdd_3, cmi_gs, start_gs, end_gs)], function(x) sum(is.na(x))), "\n")
cat("n_years_used == window length in new:", new[, all(n_years_used == win_end - win_start + 1L)], "\n")
m <- merge(old, new, by = c("id_pe", "no_mes"), suffixes = c(".o", ".n"))
stopifnot(nrow(m) == nrow(old), nrow(m) == nrow(new))
vc <- setdiff(names(old), c("id_pe", "no_mes"))
differs <- m[, Reduce(`|`, lapply(vc, function(v) { a <- get(paste0(v, ".o")); b <- get(paste0(v, ".n")); xor(is.na(a), is.na(b)) | (!is.na(a) & !is.na(b) & a != b) }))]
cat(sprintf("visits that differ (any column): %d | of these in the 5 plots: %d | outside the 5 plots (MUST be 0): %d\n",
            sum(differs), sum(differs & m$id_pe %in% five), sum(differs & !(m$id_pe %in% five))))
cat("\nthe 5 plots, new values (with the percentile of gdd_3 / cmi_gs among all visits of the same measurement year):\n")
new[, `:=`(p_gdd = frank(gdd_3) / .N, p_cmi = frank(cmi_gs) / .N), by = year]
print(new[id_pe %in% five, .(id_pe, no_mes, year, win_start, win_end, n_years_used, gdd_3 = round(gdd_3, 1), cmi_gs = round(cmi_gs, 1),
                             soil_note = NA_character_, p_gdd = round(p_gdd, 2), p_cmi = round(p_cmi, 2))][order(id_pe, no_mes)][, soil_note := NULL])
cat("\nplot pair 7209002301 vs 7209002302 (should be near-identical per visit):\n")
print(merge(new[id_pe == "7209002301", .(no_mes, year, gdd_3_01 = gdd_3, cmi_gs_01 = cmi_gs)],
            new[id_pe == "7209002302", .(no_mes, year, gdd_3_02 = gdd_3, cmi_gs_02 = cmi_gs)], by = c("no_mes", "year"), all = TRUE))
