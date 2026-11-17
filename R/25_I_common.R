# ==============================================================================
# 25_I_common.R -- shared definitions for the final checks I1-I16 (2026-11-11). Sourced, not run.
# Builds on 13_fig_helpers.R (ORDER, PAL, save_fig, vi, ...). OUT = data-output/<tag>/figures/I/
# State thermal index STI: mean z-GDD of the visits ending in each state (as 12, 16, 22).
# Shade tolerance (Niinemets & Valladares 2006, Appendix A): the ten states as in 12; other_temperate 3.76 and other_pioneer 1.18 = plot-frequency
#   weighted means of the member species' scores (I16, user 2026-11-11; weight = plot visits in the state where the species is dominant; species
#   without a published score dropped: other_temperate 2.0% of visits (willow, Sorbus decora), other_pioneer 12.9% (pin cherry)). regen unscored (user).
# Groups (as 19 / 22): boreal, temperate, pioneer, regen. Calendar periods (as 22b): year_end <=1990, 1991-2000, 2001-2010, >2010.
# ==============================================================================
source(here::here("R", "13_fig_helpers.R"))
OUT <- file.path(FIG_DIR, "I"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
scz <- readRDS(file.path(FIT_DIR, "scaling.rds")); zf <- function(x, v) (x - scz[var == v, mean]) / scz[var == v, sd]
STI <- vi[, .(STI = mean(zf(gdd_3, "gdd_3"))), by = .(state = state_end)]; STI <- setNames(STI$STI, STI$state)[ORDER]
TOL <- c(fir = 5.01, black_spruce = 4.08, white_spruce = 4.15, jack_pine = 1.36, cedar = 3.45, sugar_maple = 4.76, red_maple = 3.44,
         yellow_birch = 3.17, white_birch = 1.54, trembling_aspen = 1.21, other_pioneer = 1.18, other_temperate = NA, regen = NA)   # other_temperate unscored (user 2026-11-16)
source(here::here("R", "00_groups.R"))   # GRP: thermal groups
GCOL <- c(all = "black", boreal = "#08306B", pioneer = "#E6B800", regen = "#777777", temperate = "#B2182B")
period_of <- function(y) factor(cut(y, c(-Inf, 1990, 2000, 2010, Inf), labels = c("≤ 1990", "1991–2000", "2001–2010", "> 2010")))
