# ==============================================================================
# R/tests/03u_compare_events_real_plots.R -- compare the NEW disturbance event logic (R/03a_disturbance_events.R)
# with the CURRENT build (data/build/disturbance_covariates.rds) on real plots. Diagnostic only; writes nothing.
# Run from the project root:  Rscript R/tests/03u_compare_events_real_plots.R
#
# Inputs are rebuilt here for the chosen plots exactly as 03 steps 1-5a do (visit year from placette_mes, pee_ori_sond codes,
# ignore codes nulled, live basal area from dendro_arbres). Plots with drop-plot or CHT/DT codes are skipped (03 steps 4-5
# would change their rows). Code lists are read from R/03_build_covariates.R itself, so they cannot drift.
# Check built in: year and total_ba rebuilt here must equal those in the current build for the same plot.
#
# NOTE: the current build was made with the committed 03 of 2026-09-11 (git 9298c98), i.e. what noint/inttest used.
# Columns: F/H/P = fire / harvest (0/1/2) / pest flag; tsd = TSD lower bound in years (*_imputed; = real TSD when dated).
# A "*" marks a visit where flag or TSD changed.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(sf); library(here) })
source(here("R", "config.R")); source(here("R", "03a_disturbance_events.R"))
PEP_PATH <- if (nzchar(Sys.getenv("PEP_GPKG"))) Sys.getenv("PEP_GPKG") else here("data", "raw", "pep", "PEP.gpkg")

# plots to show: the examples discussed 2026-10-25/26 (broken in the current build), then plots where nothing should change
SHOW <- c(
  "7009501001",  # D: clear-cut seen 1992 erased by a back-filled fire 1947 (and the 14 Sep version)
  "8901804302",  # D: EL at v1, CT 1955 re-recorded later (EL is at the first visit -> not flagged, T3d)
  "7200500902",  # A: CP 1983 then CT 1984 in the same window
  "9509702501",  # A: CT 1974 then CT 1973 in the same window (old: a second clear-cut)
  "7000400301",  # F3: undated harvest first seen 1997 -> TSD stuck at 0
  "7509500602",  # CE erased by a back-filled ES 1981
  "8901804502")  # T3d + B: undated EL at the first visit (not flagged), EL again at v4 with a BA loss (new event)
N_CONTROL <- 4L   # + this many random plots whose codes are all dated within their own interval (simple cases);
                  # set it large (e.g. 1e5) to check ALL simple plots: 396 of 396 identical in the cloud check of 2026-10-26
SEED <- 11L

# ---- code lists from 03 ------------------------------------------------------------------------------------------------
ex <- parse(here("R", "03_build_covariates.R"))
for (e in ex) { s <- paste(deparse(e), collapse = "")
  if (grepl("^(ignore_origine|ignore_perturb|dropplot_origine|dropplot_perturb|droppart_origine|fire_o|fire_p|harv_o|harv_p|pest_o|pest_p) <-", s)) eval(e) }
CODES <- list(fire_o = fire_o, fire_p = fire_p, harv_o = harv_o, harv_p = harv_p, pest_o = pest_o, pest_p = pest_p)

# ---- raw inputs (as 03 steps 1-5a) ------------------------------------------------------------------------------------
cat("Reading", PEP_PATH, "...\n")
pl  <- as.data.table(st_read(PEP_PATH, layer = "placette_mes", quiet = TRUE))[, .(id_pe, no_mes, year = as.integer(format(date_sond, "%Y")))]
pee <- as.data.table(st_read(PEP_PATH, layer = "pee_ori_sond", quiet = TRUE))[, .(id_pe, no_mes, origine, an_origine = as.integer(an_origine),
                                                                                   perturb, an_perturb = as.integer(an_perturb))]
d <- merge(pee, pl, by = c("id_pe", "no_mes"), all.y = TRUE); setorder(d, id_pe, no_mes)
bad <- unique(d[origine %in% c(dropplot_origine, droppart_origine) | perturb %in% dropplot_perturb, id_pe])
d <- d[!id_pe %in% bad]
d[origine %in% ignore_origine, `:=`(origine = NA_character_, an_origine = NA_integer_)]
d[perturb %in% ignore_perturb, `:=`(perturb = NA_character_, an_perturb = NA_integer_)]
old <- as.data.table(readRDS(here("data", "build", "disturbance_covariates.rds")))
d <- d[id_pe %in% old$id_pe]

# random controls: every code dated, dated within its own interval (previous visit < year <= visit), one code per visit
d[, prev_year := shift(year), by = id_pe]
simple <- d[, .(ok = any(!is.na(origine) | !is.na(perturb)) &
                  all(is.na(origine) | (!is.na(an_origine) & !is.na(prev_year) & an_origine > prev_year & an_origine < year)) &
                  all(is.na(perturb) | (!is.na(an_perturb) & !is.na(prev_year) & an_perturb > prev_year & an_perturb < year)) &
                  !any(!is.na(origine) & !is.na(perturb))), by = id_pe][ok == TRUE, id_pe]
set.seed(SEED); ctrl <- sample(simple, min(N_CONTROL, length(simple)))
plots <- c(intersect(SHOW, d$id_pe), ctrl)
cat(sprintf("Showing %d plots (%d examples, %d random 'simple' controls out of %d simple plots)\n", length(plots), length(intersect(SHOW, d$id_pe)), length(ctrl), length(simple)))
if (length(setdiff(SHOW, d$id_pe))) cat("Not available (drop-plot / CHT-DT / not in build):", setdiff(SHOW, d$id_pe), "\n")

dend <- as.data.table(st_read(PEP_PATH, layer = "dendro_arbres", quiet = TRUE,
                              query = sprintf("SELECT id_pe, no_mes, etat, st_ha FROM dendro_arbres WHERE id_pe IN (%s)", paste0("'", plots, "'", collapse = ","))))
ba <- dend[etat %in% c("10", "12", "30", "32", "40", "42", "50", "52"), .(total_ba = sum(st_ha, na.rm = TRUE)), by = .(id_pe, no_mes)]
d <- merge(d[id_pe %in% plots], ba, by = c("id_pe", "no_mes"), all.x = TRUE); d[is.na(total_ba), total_ba := 0]; setorder(d, id_pe, no_mes)

# ---- run and print ----------------------------------------------------------------------------------------------------
fmt <- function(f, t) ifelse(f == 0, "0", paste0(f, ifelse(is.na(t), "|NA", paste0("|", t))))
for (p in plots) {
  x <- d[id_pe == p]
  r <- disturbance_plot(x$year, x$total_ba, x$origine, x$an_origine, x$perturb, x$an_perturb, CODES, BA_REGEN, BA_LOSS_THRESH)$visits
  o <- old[id_pe == p][order(no_mes)]
  stopifnot("rebuilt visits differ from the build" = identical(o$year, x$year), "rebuilt BA differs" = isTRUE(all.equal(o$total_ba, x$total_ba)))
  chg <- o$fire != r$fire | o$harv != r$harv | o$pest != r$pest |
         !mapply(identical, as.numeric(o$tsf_imputed), as.numeric(r$tsf_imputed)) |
         !mapply(identical, as.numeric(o$tsh_imputed), as.numeric(r$tsh_imputed)) |
         !mapply(identical, as.numeric(o$sfp_imputed), as.numeric(r$sfp_imputed))
  tag <- if (p %in% ctrl) "control" else "example"
  cat(sprintf("\n=== %s (%s): %s\n", p, tag, if (any(chg)) sprintf("%d of %d visits changed", sum(chg), nrow(x)) else "IDENTICAL"))
  print(data.table(v = x$no_mes, year = x$year, BA = round(x$total_ba, 1), ba_ev = as.integer(r$ba_event),
                   origine = ifelse(is.na(x$origine), "", paste(x$origine, ifelse(is.na(x$an_origine), "-", x$an_origine))),
                   perturb = ifelse(is.na(x$perturb), "", paste(x$perturb, ifelse(is.na(x$an_perturb), "-", x$an_perturb))),
                   `old F` = fmt(o$fire, o$tsf_imputed), `old H` = fmt(o$harv, o$tsh_imputed), `old P` = fmt(o$pest, o$sfp_imputed),
                   `new F` = fmt(r$fire, r$tsf_imputed), `new H` = fmt(r$harv, r$tsh_imputed), `new P` = fmt(r$pest, r$sfp_imputed),
                   ` ` = ifelse(chg, "*", "")), row.names = FALSE)
}
