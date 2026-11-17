# ==============================================================================
# 25_I1_raw_thermophilization.R   (no INLA; seconds; run from the project root)
# I1: the observed, model-free record. Per observed visit interval: change in the state thermal index (STI) and in shade tolerance between the
#   start and end visit, per decade (change / interval length x 10); stays count as 0. Tolerance: intervals with regen at either end excluded (regen
#   unscored). Means by origin group x calendar period (end visit), plus all intervals; 95% CI from a plot bootstrap (200).
# Caveat carried: STI is computed from the same plots (as 12); calendar periods differ in which plots were measured (A1e).
# ==============================================================================
source(here::here("R", "25_I_common.R")); set.seed(25)
d <- vi[, .(id_pe, a = state_start, b = state_end, len = year_end - year_start, per = period_of(year_end))]
d[, `:=`(grp = GRP[a], dsti = (STI[b] - STI[a]) / len * 10, dtol = (TOL[b] - TOL[a]) / len * 10)]
cells <- rbind(d[, .(grp = "all", id_pe, per, dsti, dtol)], d[, .(grp, id_pe, per, dsti, dtol)])
cells <- rbind(cells, cells[, .(grp, id_pe, per = "all periods", dsti, dtol)])
est <- function(x) x[, .(n = .N, sti = mean(dsti), tol = mean(dtol, na.rm = TRUE), n_tol = sum(!is.na(dtol))), by = .(grp, per)]
E0 <- est(cells); ids <- unique(d$id_pe); cs <- split(cells, by = "id_pe")
B <- rbindlist(lapply(1:200, function(r) est(rbindlist(cs[sample(ids, replace = TRUE)]))[, r := r]))
CI <- B[, .(sti_lo = quantile(sti, .025), sti_hi = quantile(sti, .975), tol_lo = quantile(tol, .025, na.rm = TRUE), tol_hi = quantile(tol, .975, na.rm = TRUE)), by = .(grp, per)]
R1 <- merge(E0, CI, by = c("grp", "per")); R1[, per := factor(per, levels = c(levels(d$per), "all periods"))]
R1[, grp := factor(grp, levels = names(GCOL))]; setorder(R1, grp, per)
L <- melt(R1, id.vars = c("grp", "per", "n"), measure.vars = list(c("sti", "tol"), c("sti_lo", "tol_lo"), c("sti_hi", "tol_hi")),
          variable.name = "q", value.name = c("est", "lo", "hi"))
L[, q := factor(c("thermal index (ΔSTI per decade)", "shade tolerance (per decade)")[q], levels = c("thermal index (ΔSTI per decade)", "shade tolerance (per decade)"))]
p <- ggplot(L, aes(per, est, colour = grp)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(ymin = lo, ymax = hi, group = grp), width = 0, position = position_dodge(width = 0.6)) +
  geom_point(aes(group = grp), position = position_dodge(width = 0.6), size = 1.8) +
  scale_colour_manual(values = GCOL, name = "origin group") + facet_wrap(~ q, scales = "free_y") +
  labs(x = "period of the end visit", y = "observed change per decade", title = "I1  Observed change in state thermal index and shade tolerance",
       subtitle = "every observed visit interval (stays count as 0); 95% CI from a plot bootstrap; tolerance excludes intervals involving regen")
save_fig(p, "I1_raw_thermophilization", 12, 5, R1)
# I1b (user 2026-11-11): by latitude band instead of domain (no domain field in the data); tertiles of plot latitude, all periods.
con <- DBI::dbConnect(RSQLite::SQLite(), "data/raw/pep/PEP.gpkg", flags = RSQLite::SQLITE_RO)
pl <- unique(data.table::as.data.table(DBI::dbGetQuery(con, "SELECT id_pe, latitude FROM placette")), by = "id_pe"); DBI::dbDisconnect(con)
qs <- quantile(pl[id_pe %in% ids, latitude], c(1/3, 2/3), na.rm = TRUE)
pl[, band := cut(latitude, c(-Inf, qs, Inf), labels = sprintf(c("south (< %.1f°N)", "middle (%.1f–%.1f°N)", "north (> %.1f°N)"), c(qs[1], qs[1], qs[2]), c(NA, qs[2], NA))[1:3])]
pl[, band := factor(band, levels = levels(band))]
c2 <- merge(rbind(d[, .(grp = "all", id_pe, dsti, dtol)], d[, .(grp, id_pe, dsti, dtol)]), pl[, .(id_pe, band)], by = "id_pe")
est2 <- function(x) x[, .(n = .N, sti = mean(dsti), tol = mean(dtol, na.rm = TRUE)), by = .(grp, band)]
E2 <- est2(c2); cs2 <- split(c2, by = "id_pe")
B2 <- rbindlist(lapply(1:200, function(r) est2(rbindlist(cs2[sample(ids, replace = TRUE)]))))
R2 <- merge(E2, B2[, .(sti_lo = quantile(sti, .025), sti_hi = quantile(sti, .975), tol_lo = quantile(tol, .025, na.rm = TRUE), tol_hi = quantile(tol, .975, na.rm = TRUE)), by = .(grp, band)], by = c("grp", "band"))
R2[, grp := factor(grp, levels = names(GCOL))]; setorder(R2, grp, band)
L2 <- melt(R2, id.vars = c("grp", "band", "n"), measure.vars = list(c("sti", "tol"), c("sti_lo", "tol_lo"), c("sti_hi", "tol_hi")), variable.name = "q", value.name = c("est", "lo", "hi"))
L2[, q := factor(levels(L$q)[q], levels = levels(L$q))]
p2 <- ggplot(L2, aes(band, est, colour = grp)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(ymin = lo, ymax = hi, group = grp), width = 0, position = position_dodge(width = 0.6)) +
  geom_point(aes(group = grp), position = position_dodge(width = 0.6), size = 1.8) +
  scale_colour_manual(values = GCOL, name = "origin group") + facet_wrap(~ q, scales = "free_y") +
  labs(x = "latitude band (tertiles of plots)", y = "observed change per decade", title = "I1b  Observed change by latitude band",
       subtitle = "all periods; 95% CI from a plot bootstrap; tolerance excludes intervals involving regen")
save_fig(p2, "I1b_raw_thermophilization_latitude", 12, 5, R2)
cat("I1 done\n")
