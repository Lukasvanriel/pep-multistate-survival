# ==============================================================================
# classify_pep_states.R
#
# Classifies each plot-measurement in the PEP permanent plot network into
# a forest state based on:
#   1. Plantation exclusion (EPO/PIS presence, pee_ori_sond$origine, or
#      station_pe$origine carrying a plantation code at ANY visit -> whole
#      plot dropped, before included -- see section 3)
#   2. Data-quality exclusion (BA=0 + bad status/terrain code, or a plot
#      that is barren -- no trees/saplings/seedlings -- across its entire
#      recorded history)
#   3. Regeneration rule: total_ba < 5 m2/ha (Brice et al. 2020)
#   4. Dominant live basal area species (applied to non-regen)
#   5. Drop any observation whose dominant species has no state mapping
#      (residual "unclassified_species")
#   6. Drop any plot left with fewer than 2 observations
#
# States:
#   fir, black_spruce, sugar_maple, red_maple, yellow_birch, white_birch,
#   trembling_aspen, jack_pine, cedar, white_spruce,
#   other_temperate, other_pioneer, regen
#
# NOTE: "unclassified" no longer exists as an output state under this rule
# set -- every retained observation gets a real state. See section 6 for why.
#
# Output: data.table with columns id_pe, no_mes, state, total_ba
# ==============================================================================
library(sf)
library(data.table)
library(here)
PEP_PATH <- here("data", "raw", "pep", "PEP.gpkg")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
cat("Loading data...\n")
placette <- as.data.table(st_read(PEP_PATH, layer = "placette_mes",   quiet = TRUE))
station  <- as.data.table(st_read(PEP_PATH, layer = "station_pe",     quiet = TRUE))
dendro   <- as.data.table(st_read(PEP_PATH, layer = "dendro_arbres",  quiet = TRUE))
gaules   <- as.data.table(st_read(PEP_PATH, layer = "dendro_gaules",  quiet = TRUE))
semis    <- as.data.table(st_read(PEP_PATH, layer = "station_semis",  quiet = TRUE))
ori_sond <- as.data.table(st_read(PEP_PATH, layer = "pee_ori_sond",   quiet = TRUE))

placette[, year := as.integer(format(as.Date(date_sond), "%Y"))]

cat(sprintf("  placette_mes:  %d rows\n", nrow(placette)))
cat(sprintf("  station_pe:    %d rows\n", nrow(station)))
cat(sprintf("  dendro_arbres: %d rows\n", nrow(dendro)))
cat(sprintf("  pee_ori_sond:  %d rows\n", nrow(ori_sond)))

# ==============================================================================
# 2. DEFINITIONS
# ==============================================================================
LIVE_ETAT <- c("10", "12", "30", "32", "40", "42", "50", "52")

# Plantation-indicator species (Brice et al. 2020): presence at ANY visit
# disqualifies the whole plot. PICABI (Norway spruce) = EPO, PINSYL (Scots
# pine) = PIS in PEP essence codes.
PLANTATION_SPECIES <- c("EPO", "PIS")

# Plantation-origin codes: a code from this list appearing at ANY visit, in
# EITHER pee_ori_sond$origine or station_pe$origine, disqualifies the WHOLE
# plot (every visit, before and after -- see section 3 for rationale).
PLANTATION_ORIGINE <- c("P", "PL", "PLR", "PLN", "PRR", "ENS", "ENM")

# Bad plot status / non-forest terrain: only disqualifying when paired with
# total_ba == 0 (see section 6). AL, DH, DS were checked individually this
# round and confirmed to carry real forest when BA > 0 -- do not drop them
# outright.
BAD_STATUS <- c("AB", "SR", "NT", "DE")
BAD_COTER  <- c("RO", "AL", "DH", "DS", "EAU", "LTE")

# Species-to-state mapping. ERP (striped maple), NOC (butternut), ORT (rock
# elm) added to other_temperate -- Brice et al.'s own temperate group
# (ACEPEN, JUGCIN, ULMTHO) includes these. PRV, AUR, MAS, MEH are
# deliberately left unmapped: they fall through to "unclassified_species"
# and get dropped in section 9 only if one of them ends up dominant --
# checked against pre-filtering them out entirely, the difference was 3
# observations out of ~50,000, not worth the extra complexity.
SPECIES_TO_STATE <- list(
  fir             = "SAB",
  black_spruce    = "EPN",
  sugar_maple     = "ERS",
  red_maple       = "ERR",
  yellow_birch    = "BOJ",
  white_birch     = "BOP",
  trembling_aspen = "PET",
  jack_pine       = "PIG",
  cedar           = "THO",
  white_spruce    = "EPB",
  other_temperate = c("HEG", "TIL", "FRA", "FRN", "CHR", "OSV", "ORA",
                      "ERA", "PRU", "ERE", "ERN", "CET", "CHB", "CHG",
                      "CHE", "SAL", "BOG", "PIB", "PIR", "SOA", "SOD",
                      "EPR", "PID", "FRP", "ERP", "NOC", "ORT"),
  other_pioneer   = c("PEG", "PEB", "PED", "PEH", "PRP", "MEL")
)

# ==============================================================================
# 3. PLANTATION EXCLUSION
#
# Three independent signals feed one whole-plot drop:
#   (a) EPO/PIS live-tree presence at any visit (dendro_arbres)
#   (b) plantation-origin code in pee_ori_sond$origine at any visit
#   (c) plantation-origin code in station_pe$origine at any visit
#
# RESOLVED (was previously an open question -- see former section 9 note):
# station_pe$origine is a per-visit field, separate from pee_ori_sond, that
# was being loaded (for co_ter) but never checked for plantation codes. A
# plot whose ONLY plantation signal is a station_pe$origine code on one
# later visit still gets the WHOLE plot dropped here, "before" observations
# included -- there is no reliable way to tell whether an earlier,
# unlabelled visit was genuinely natural forest or simply pre-dates
# whichever record-keeping produced the later plantation label, so we don't
# try to split the difference. Verified against the raw geopackage: of the
# ~565 plots that carry a plantation-origin code in station_pe$origine at
# some visit, 468 were already caught by (a)/(b); this section adds the
# remaining ~97 plots (~390 raw plot-measurement rows) that (a)/(b) alone
# were missing.
# ==============================================================================
cat("\nApplying plantation exclusions (whole-plot, before included)...\n")
epo_pis_plots     <- unique(
  dendro[etat %in% LIVE_ETAT & essence %in% PLANTATION_SPECIES & st_ha > 0, id_pe]
)
origine_plots     <- unique(ori_sond[origine %in% PLANTATION_ORIGINE, id_pe])
station_origine_plots <- unique(station[origine %in% PLANTATION_ORIGINE, id_pe])
plantation_plots  <- Reduce(union, list(epo_pis_plots, origine_plots, station_origine_plots))
cat(sprintf("  Plots dropped for EPO/PIS presence: %d\n", length(epo_pis_plots)))
cat(sprintf("  Plots dropped for plantation origine (pee_ori_sond, any visit): %d\n", length(origine_plots)))
cat(sprintf("  Plots dropped for plantation origine (station_pe, any visit): %d\n", length(station_origine_plots)))
cat(sprintf("  Of which newly caught only by station_pe (not by EPO/PIS or pee_ori_sond): %d\n",
            length(setdiff(station_origine_plots, union(epo_pis_plots, origine_plots)))))
cat(sprintf("  Union (whole-plot drop, all 3 signals): %d\n", length(plantation_plots)))

placette <- placette[!id_pe %in% plantation_plots]
station  <- station[!id_pe %in% plantation_plots]
dendro   <- dendro[!id_pe %in% plantation_plots]
gaules   <- gaules[!id_pe %in% plantation_plots]
semis    <- semis[!id_pe %in% plantation_plots]

# ==============================================================================
# 4. BARREN-PLOT IDENTIFICATION (before any other data-quality drop)
#
# A plot is barren if, across its ENTIRE recorded history (not just what
# survives sections 3a/3b), it never has live trees, saplings, or seedlings.
# Computed on the full, unfiltered layers so a plot's whole history counts,
# even visits already removed by section 3.
# ==============================================================================
cat("\nIdentifying barren plots...\n")
dendro_full <- as.data.table(st_read(PEP_PATH, layer = "dendro_arbres", quiet = TRUE))
gaules_full <- as.data.table(st_read(PEP_PATH, layer = "dendro_gaules", quiet = TRUE))
semis_full  <- as.data.table(st_read(PEP_PATH, layer = "station_semis", quiet = TRUE))

ba_by_plot   <- dendro_full[etat %in% LIVE_ETAT, .(any_ba = sum(st_ha, na.rm = TRUE) > 0), by = id_pe]
sap_by_plot  <- gaules_full[, .(any_sap = sum(nb_tige, na.rm = TRUE) > 0), by = id_pe]
seed_by_plot <- semis_full[, .(any_seed = sum(nb_semis, na.rm = TRUE) > 0), by = id_pe]

all_plots_full <- unique(as.data.table(st_read(PEP_PATH, layer = "placette_mes", quiet = TRUE))[, .(id_pe)])
plot_ever <- merge(all_plots_full, ba_by_plot, by = "id_pe", all.x = TRUE)
plot_ever <- merge(plot_ever, sap_by_plot, by = "id_pe", all.x = TRUE)
plot_ever <- merge(plot_ever, seed_by_plot, by = "id_pe", all.x = TRUE)
plot_ever[is.na(any_ba), any_ba := FALSE]
plot_ever[is.na(any_sap), any_sap := FALSE]
plot_ever[is.na(any_seed), any_seed := FALSE]

barren_plots <- plot_ever[any_ba == FALSE & any_sap == FALSE & any_seed == FALSE, id_pe]
cat(sprintf("  Barren plots (no trees/saplings/seedlings, ever): %d\n", length(barren_plots)))

# ==============================================================================
# 5. COMPUTE LIVE BASAL AREA PER PLOT-MEASUREMENT
# ==============================================================================
cat("\nComputing live basal area...\n")
live_ba <- dendro[etat %in% LIVE_ETAT,
                  .(total_ba = sum(st_ha, na.rm = TRUE)),
                  by = .(id_pe, no_mes)]
all_meas    <- placette[, .(id_pe, no_mes)]
all_meas_ba <- merge(all_meas, live_ba, by = c("id_pe", "no_mes"), all.x = TRUE)
all_meas_ba[is.na(total_ba), total_ba := 0]

meas_full <- merge(all_meas_ba,
                   station[, .(id_pe, no_mes, co_ter)],
                   by = c("id_pe", "no_mes"))
meas_full <- merge(meas_full,
                   placette[, .(id_pe, no_mes, statut_mes, year)],
                   by = c("id_pe", "no_mes"))
cat(sprintf("  Total plot-measurements (post plantation exclusion): %d\n", nrow(meas_full)))

# ==============================================================================
# 6. DATA-QUALITY DROP
#
# Drop an observation if total_ba == 0 AND (bad statut_mes OR bad co_ter
# OR the plot is barren). BA=0 alone is not disqualifying -- it's a normal
# regen value (see section 8); it's only dropped in combination with a
# real data-quality signal.
# ==============================================================================
cat("\nApplying data-quality drop...\n")
meas_full[, problematic := (statut_mes %in% BAD_STATUS) |
            (co_ter %in% BAD_COTER) |
            (id_pe %in% barren_plots)]
meas_full[, drop_dq := (total_ba == 0) & problematic]
cat(sprintf("  Dropped (BA=0 + bad status/terrain/barren): %d\n", meas_full[drop_dq == TRUE, .N]))
meas_full <- meas_full[drop_dq == FALSE]

# ==============================================================================
# 7. ASSIGN DOMINANT STATE
# ==============================================================================
cat("\nAssigning dominant state...\n")
state_lookup <- rbindlist(lapply(names(SPECIES_TO_STATE), function(s) {
  data.table(essence = SPECIES_TO_STATE[[s]], state_group = s)
}))
dendro_live <- dendro[etat %in% LIVE_ETAT]
ba_by_species <- dendro_live[, .(species_ba = sum(st_ha, na.rm = TRUE)),
                             by = .(id_pe, no_mes, essence)]
dominant_species <- ba_by_species[, .SD[which.max(species_ba)], by = .(id_pe, no_mes)]
setnames(dominant_species, c("essence", "species_ba"),
         c("dominant_essence", "dominant_ba"))
dominant <- merge(dominant_species, state_lookup,
                  by.x = "dominant_essence", by.y = "essence",
                  all.x = TRUE)
dominant[is.na(state_group), state_group := "unclassified_species"]
setnames(dominant, "state_group", "dominant_state")

# ==============================================================================
# 8. COMBINE: REGEN THRESHOLD FIRST, THEN DOMINANCE
#
# Regen is now a pure basal-area threshold (Brice et al. 2020), with no
# dependency on cl_age/origine/perturb -- those fields were shown this
# session to be unreliable and circular.
# ==============================================================================
cat("\nCombining...\n")
classified <- merge(meas_full[, .(id_pe, no_mes, total_ba, year)],
                    dominant,
                    by = c("id_pe", "no_mes"),
                    all.x = TRUE)
classified[, state := fcase(
  total_ba < 5,            "regen",
  !is.na(dominant_state),  dominant_state,
  default =                "unclassified_species"
)]

# ==============================================================================
# 9. FINAL CLEANUP
#
# 9a. Drop residual unclassified_species (unmapped dominant species -- PRV/
#     AUR/MAS/MEH when dominant, or any other unmapped code).
# 9b. Drop plots left with fewer than 2 observations.
# ==============================================================================
cat("\nFinal cleanup...\n")
n_unclassified_species <- classified[state == "unclassified_species", .N]
cat(sprintf("  Dropping unclassified_species: %d\n", n_unclassified_species))
classified <- classified[state != "unclassified_species"]

obs_counts <- classified[, .N, by = id_pe]
single_obs_plots <- obs_counts[N < 2, id_pe]
cat(sprintf("  Plots dropped for <2 remaining observations: %d\n", length(single_obs_plots)))
classified <- classified[!id_pe %in% single_obs_plots]

# ==============================================================================
# 10. DIAGNOSTICS
# ==============================================================================
cat("\n--- STATE DISTRIBUTION ---\n")
print(classified[, .N, by = state][order(-N)])
cat(sprintf("\nTotal retained: %d plot-measurements, %d plots\n",
            nrow(classified), uniqueN(classified$id_pe)))

# ==============================================================================
# 11. SAVE
# ==============================================================================
out <- classified[, .(id_pe, no_mes, state, total_ba)]
dir.create("data/processed/", showWarnings = FALSE, recursive = TRUE)
saveRDS(out, "data/processed/state_classification.rds")
cat("\nSaved: data/processed/state_classification.rds\n")
cat(sprintf("  %d plot-measurements\n", nrow(out)))
