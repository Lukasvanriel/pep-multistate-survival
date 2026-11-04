# ==============================================================================
# R/tests/03w_twenty_plots_expected.R -- the 20 plots inspected by the user on 2026-10-26, with the values expected
# under the agreed rules, DERIVED BY HAND before running the code (main settings: TSD_SOURCE mid, first-visit ages
# harvest 4.5 / pest 32, stand ages from 02b). Compares data/build/disturbance_covariates.rds with them.
# Each cell: flag|TSD (mid), "0" = not flagged. Run from the project root after 03.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here) })
E <- list(   # plot = list(type = c(per visit))
  "7000400301" = list(harv = c("0","0","0","1|6","1|17","1|31")),
  "7000409902" = list(harv = c("0","0","1|3.5","1|16.5","1|28.5","1|34.5")),
  "7000700401" = list(harv = c("1|4.5","1|12.5","1|19.5","2|4.5","2|19.5","2|32.5")),
  "7009800801" = list(harv = c("0","0","2|6","1|7","1|20")),   # CT 1981 + CT 1986 (same interval, equal -> 1986); EPC 1996 other interval -> own event, most recent shown
  "7100600901" = list(harv = c("1|4.5","1|3.5","1|17.5","1|28.5","1|41.5")),
  "7109603102" = list(harv = c("0","0","1|7","1|18","1|30")),
  "7209404502" = list(harv = c("0","0","1|7","1|17","1|30")),
  "7209900402" = list(harv = c("2|10","2|17","2|32","2|41","2|53")),   # CP back-filled 1962, CT 1962 same interval -> heavy (A2 highest)
  "7309500101" = list(pest = c("1|32","1|39","0","0","0"), harv = c("0","0","2|10","2|22","2|34")),
  "7506900901" = list(harv = c("0","0","1|6","1|14","1|28")),
  "7601604601" = list(harv = c("2|30","2|38","2|43","2|57","2|69")),
  "7601700502" = list(harv = c("2|22","2|30","2|43","2|54")),
  "7609504402" = list(harv = c("2|70","2|80","2|91","2|97","2|110")),
  "7609507202" = list(harv = c("2|30","2|41","2|52","2|61","2|75")),
  "8901903402" = list(fire = c("1|30","1|42","1|52")),
  "8909202602" = list(pest = c("0","1|4","1|12","1|25")),
  "8909901202" = list(pest = c("1|32","1|5")),
  "9109609601" = list(harv = c("1|4.5","1|10.5","1|16.5","1|13")),
  "9500701601" = list(harv = c("2|10","2|18","2|10")),
  "9509900202" = list(harv = c("0","1|7","1|20")))
nw <- as.data.table(readRDS(here("data", "build", "disturbance_covariates.rds"))); setkey(nw, id_pe, no_mes)
# The expected values are for pee_ori_sond codes only. With station codes (USE_STATION_CODES = TRUE, main) these plots change by
# design (station-only codes; first-visit harvest age 4.5 -> 4.0 y because station drop-plot/CHT codes change the reference set):
# run this check on an S20 build (USE_STATION_CODES = FALSE) to test the pee logic; the station merge is checked separately.
if (isTRUE(readRDS(here("data", "build", "disturbance_events.rds"))$use_station_codes))
  cat("NOTE: this build includes station codes -- mismatches below are expected; the expectations are for pee codes only (S20).\n\n")
fmt <- function(f, t) ifelse(f == 0, "0", paste0(f, "|", ifelse(is.na(t), "NA", format(t, trim = TRUE, drop0trailing = TRUE))))
res <- rbindlist(lapply(names(E), function(p) rbindlist(lapply(names(E[[p]]), function(T) {
  ts <- c(fire = "tsf_imputed_mid", harv = "tsh_imputed_mid", pest = "sfp_imputed_mid")[[T]]
  x <- nw[id_pe == p]; got <- fmt(x[[T]], x[[ts]]); exp <- E[[p]][[T]]
  if (length(got) != length(exp)) stop(p, " ", T, ": ", length(got), " visits in the build, ", length(exp), " expected")
  data.table(plot = p, type = T, visit = x$no_mes, year = x$year, expected = exp, got = got, ok = got == exp) }))))
print(res[, .(plot, type, visit, year, expected, got, ok = ifelse(ok, "", "MISMATCH"))], nrows = 200, row.names = FALSE)
cat(sprintf("\n%d of %d cells as expected; plots fully as expected: %d of %d\n", sum(res$ok), nrow(res),
            res[, all(ok), by = plot][V1 == TRUE, .N], uniqueN(res$plot)))
