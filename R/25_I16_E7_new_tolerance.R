# ==============================================================================
# 25_I16_E7_new_tolerance.R   (no INLA; <1 min; after 23_extended_quantities.R; run from the project root)
# I16: E7 (successional direction) with plot-frequency-weighted shade tolerance for the pooled classes (Niinemets & Valladares 2006 scores of the
#   member species, weighted by their frequency in the visits of that class; unscored species dropped): other temperate 3.76 (2.0% of visits unscored),
#   other pioneer 1.18 (12.9% unscored); regen unscored (user 2026-11-11). Tolerance groups: other temperate -> intermediate (3.76 lies between the
#   intermediate states, 3.17-3.45, and the lowest tolerant state, 4.08). Everything else as 23b's E7; original E7 read from figures/E for comparison.
# ==============================================================================
src <- readLines(here::here("R", "23b_extended_summaries.R")); eval(parse(text = src[6:grep("^base <- setNames", src)]))
TOL["other_pioneer"] <- 1.18   # other_temperate stays unscored and ungrouped (user 2026-11-16)
i0 <- grep("^E7 <- list\\(\\); E7g <- list\\(\\)", src); i1 <- grep("^E7 <- rbindlist\\(E7\\)", src) - 1; eval(parse(text = src[i0:i1]))
E7 <- rbindlist(E7); E7g <- rbindlist(E7g); OUTI <- file.path(FIG_DIR, "I")
O7 <- fread(file.path(OUTE, "E7_successional_direction.csv"))
C <- merge(O7[, .(scen, horizon, quantity, old = est, old_lo = lo, old_hi = hi)], E7[, .(scen, horizon, quantity, new = est, new_lo = lo, new_hi = hi)], by = c("scen", "horizon", "quantity"))
fwrite(C, file.path(OUTI, "I16_E7_old_vs_new_tolerance.csv")); fwrite(E7g, file.path(OUTI, "I16_E7_tolerance_group_flows.csv"))
options(width = 220)
f <- function(e, l, h) sprintf("%.3f [%.3f, %.3f]", e, l, h)
print(C[horizon == 50 & quantity %in% c("share of moves scored for tolerance", "mean tolerance shift of scored moves", "net (more - less tolerant)"),
        .(scen, quantity, old = f(old, old_lo, old_hi), new = f(new, new_lo, new_hi))][order(quantity, match(scen, SCN))])
