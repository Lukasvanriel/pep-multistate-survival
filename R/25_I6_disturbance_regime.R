# ==============================================================================
# 25_I6_disturbance_regime.R   (no INLA; seconds; run from the project root)
# I6: the disturbance regime behind the disturbance x climate comparisons (SI figure + Methods sentence).
#  (a) co-occurrence: share of visit intervals with each combination of active flags (fire, pest, light harvest, heavy harvest; as in the fit).
#  (b) climate position of NEW events: GDD, CMI and soil (z, as in the fit) of the interval in which an event starts (disturbance_events.rds,
#      not dropped, start visit = end visit of the interval), by type and level, vs all intervals; mean and 95% CI (plot bootstrap, 200).
#  (c) return intervals: years between consecutive dated events in the same plot (any type; same type), retained plots only.
# ==============================================================================
source(here::here("R", "25_I_common.R")); set.seed(26)
ev <- as.data.table(readRDS(file.path(BUILD_DIR, "disturbance_events.rds"))$events)[dropped == FALSE & id_pe %in% vi$id_pe]
ev[, kind := fcase(type == "harv" & level == 1, "light harvest", type == "harv" & level == 2, "heavy harvest",
                   type == "fire", "fire", type == "pest", "pest")]   # pest intensities merged, as in the model (user 2026-11-11)
cat("events in retained plots:\n"); print(ev[, .N, by = kind][order(-N)])
# (a) co-occurrence of active flags per interval
f <- vi[, .(fire = fire > 0, pest = pest > 0, light = harv == 1, heavy = harv == 2)]
f[, combo := fifelse(!fire & !pest & !light & !heavy, "none", paste0(fifelse(fire, "fire+", ""), fifelse(pest, "pest+", ""), fifelse(light, "light harvest+", ""), fifelse(heavy, "heavy harvest+", "")))]
f[, combo := sub("\\+$", "", combo)]; CO <- f[, .N, by = combo][, share := round(N / sum(N), 3)][order(-N)]
cat("\nactive flag combinations per interval:\n"); print(CO); fwrite(CO, file.path(OUT, "I6_flag_combinations.csv"))
# (b) climate position of new events
iv <- vi[, .(id_pe, no_mes = visit_end, z_gdd = zf(gdd_3, "gdd_3"), z_cmi = zf(cmi_gs, "cmi_gs"), z_soil = zf(soil_pc1, "soil_pc1"))]
ne <- merge(ev[, .(id_pe, no_mes = start_no_mes, kind)], iv, by = c("id_pe", "no_mes"))
allv <- iv[, .(id_pe, kind = "all intervals", z_gdd, z_cmi, z_soil)]
X <- rbind(ne[, .(id_pe, kind, z_gdd, z_cmi, z_soil)], allv)
m <- function(x) x[, .(n = .N, gdd = mean(z_gdd), cmi = mean(z_cmi), soil = mean(z_soil, na.rm = TRUE)), by = kind]
M0 <- m(X); ids <- unique(X$id_pe); xs <- split(X, by = "id_pe")
B <- rbindlist(lapply(1:200, function(r) m(rbindlist(xs[sample(ids, replace = TRUE)]))))
CL <- melt(M0, id.vars = c("kind", "n"), variable.name = "cov", value.name = "est")
BL <- melt(B, id.vars = c("kind", "n"), variable.name = "cov")[, .(lo = quantile(value, .025), hi = quantile(value, .975)), by = .(kind, cov)]
CL <- merge(CL, BL, by = c("kind", "cov")); fwrite(CL, file.path(OUT, "I6_climate_position.csv"))
KL <- c("all intervals", "fire", "pest", "light harvest", "heavy harvest")
CL[, `:=`(kind = factor(kind, levels = rev(KL)), cov = factor(c(gdd = "GDD", cmi = "CMI", soil = "soil PC1")[as.character(cov)], levels = c("GDD", "CMI", "soil PC1")))]
KC <- c(`all intervals` = "black", fire = "#D95F02", pest = "#7570B3", `light harvest` = "#66A61E", `heavy harvest` = "#1B6E3C")
p <- ggplot(CL, aes(est, kind, colour = kind)) + geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y") + geom_point(size = 2) + scale_colour_manual(values = KC, guide = "none") +
  facet_wrap(~ cov, nrow = 1) + scale_y_discrete(labels = function(x) paste0(x, " (", M0$n[match(x, M0$kind)], ")")) +
  labs(x = "mean z-score (SD from the mean)", y = NULL, title = "I6  Climate and soil where disturbances happen",
       subtitle = "intervals in which a new event starts, vs all intervals; 95% CI from a plot bootstrap; in brackets: number of events")
save_fig(p, "I6_disturbance_climate_position", 11, 3.6, NULL)
# (c) return intervals (dated events only)
de <- ev[!is.na(yr)][order(id_pe, yr)]
ri_any <- de[, .(gap = as.numeric(diff(yr)), pair = paste(head(kind, -1), "→", tail(kind, -1))), by = id_pe][gap > 0]
ri_same <- de[, .(gap = as.numeric(diff(yr))), by = .(id_pe, kind)][gap > 0]
cat("\nreturn intervals (years), dated events:\n")
print(rbind(ri_any[, .(what = "any → any", n = .N, median = median(gap), p25 = quantile(gap, .25), p75 = quantile(gap, .75))],
            ri_same[, .(what = paste("same type:", kind[1]), n = .N, median = median(gap), p25 = quantile(gap, .25), p75 = quantile(gap, .75)), by = kind][, kind := NULL]))
print(ri_any[, .N, by = pair][order(-N)][1:8])
fwrite(ri_any, file.path(OUT, "I6_return_intervals.csv"))
cat("I6 done\n")
