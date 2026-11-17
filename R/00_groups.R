# 00_groups.R -- state groupings used for REPORTING results (fits are always at the 13-state level; groups only aggregate outputs).
# Single definition: every script sources this file instead of defining its own copy.
#
# Thermal groups (user-defined): boreal, pioneer, temperate, regen.
GRP <- c(black_spruce = "boreal", white_spruce = "boreal", fir = "boreal", jack_pine = "boreal", sugar_maple = "temperate", red_maple = "temperate",
         yellow_birch = "temperate", cedar = "temperate", other_temperate = "temperate", white_birch = "pioneer", trembling_aspen = "pioneer",
         other_pioneer = "pioneer", regen = "regen")
#
# Shade-tolerance groups (Niinemets & Valladares 2006 scores of the dominant species; user 2026-11-09): intolerant (1.2-1.5), intermediate (3.2-3.5),
# tolerant (4.1-5.0), regen its own group. other_temperate: NA everywhere (user 2026-11-16) -- its species span 1.9-4.8, so it has no single tolerance.
TOL_GROUP <- c(fir = "tolerant", black_spruce = "tolerant", white_spruce = "tolerant", jack_pine = "intolerant", cedar = "intermediate",
               sugar_maple = "tolerant", red_maple = "intermediate", yellow_birch = "intermediate", white_birch = "intolerant",
               trembling_aspen = "intolerant", other_pioneer = "intolerant", other_temperate = NA_character_, regen = "regen")
#
# Shade-tolerance scores (Niinemets & Valladares 2006, Appendix A; vector copied here unchanged from R/25_I_common.R (TOL), 2026-11-17).
# other_pioneer 1.18 = plot-frequency weighted mean of its member species (I16). other_temperate and regen: NA (no score).
TOL_SCORE <- c(fir = 5.01, black_spruce = 4.08, white_spruce = 4.15, jack_pine = 1.36, cedar = 3.45, sugar_maple = 4.76, red_maple = 3.44,
               yellow_birch = 3.17, white_birch = 1.54, trembling_aspen = 1.21, other_pioneer = 1.18, other_temperate = NA, regen = NA)
