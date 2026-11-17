# ==============================================================================
# 25_I13b_summary.R   (no INLA; after 25_I13_disturbance_regime.R; run from the project root)
# I13 summary: check (simulated, no events) vs I8's deterministic equilibrium; equilibrium under the disturbance regime (S, SC) vs observed; figure.
# ==============================================================================
source(here::here("R", "25_I_common.R"))
I8 <- fread(file.path(OUT, "I8_equilibrium_gradient.csv"))[, .(bin, state, I8 = model, obs, zg)]
SM <- rbindlist(lapply(list.files(file.path(OUT, "cache_I13"), full.names = TRUE), readRDS))
X <- merge(dcast(SM, bin + state ~ version, value.var = "share"), I8, by = c("bin", "state"))
fwrite(X, file.path(OUT, "I13_equilibrium_disturbance.csv"))
cat("check: max |simulated no-event - I8| per bin:\n"); print(X[, .(max_abs = round(max(abs(none - I8)), 3), TVD = round(0.5 * sum(abs(none - I8)), 3)), by = bin])
TV <- X[, .(TVD_undisturbed_I8 = 0.5 * sum(abs(I8 - obs)), TVD_regime_S = 0.5 * sum(abs(S - obs)), TVD_regime_SC = 0.5 * sum(abs(SC - obs))), by = .(bin, zg)]
cat("\nTVD vs observed:\n"); print(TV[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])
G <- X[, lapply(.SD, sum), by = .(bin, zg, grp = GRP[state]), .SDcols = c("obs", "I8", "S", "SC", "SC_noharv", "SC_lightonly", "SC_heavyonly", "SC_harvonly")]
fwrite(G, file.path(OUT, "I13_group_shares.csv"))   # incl. no-harvest version (fire and pest only)
cat("\ngroup shares obs | undisturbed | regime S | regime SC:\n")
print(dcast(G[, .(bin, grp, v = sprintf("%.2f | %.2f | %.2f | %.2f", obs, I8, S, SC))], bin ~ grp, value.var = "v"))
GL <- melt(G, id.vars = c("bin", "zg", "grp"), variable.name = "src")[!src %in% c("S", "SC_noharv", "SC_lightonly", "SC_heavyonly", "SC_harvonly")]
GL[, `:=`(src = factor(c(obs = "observed", I8 = "equilibrium, undisturbed (I8)", SC = "equilibrium, disturbance regime")[as.character(src)],
  levels = c("observed", "equilibrium, undisturbed (I8)", "equilibrium, disturbance regime")), xlab = sprintf("%s\n(GDD %+.1f SD)", bin, zg))]
p <- ggplot(GL, aes(xlab, value, colour = grp, group = interaction(grp, src), linetype = src)) + geom_line() + geom_point(size = 1.4) +
  scale_colour_manual(values = GCOL[c("boreal", "pioneer", "regen", "temperate")], name = NULL) + scale_linetype_manual(values = c("solid", "dotted", "dashed"), name = NULL) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(x = "plots binned by mean GDD (quintiles)", y = "share of stands", title = "I13  Observed composition vs equilibrium with and without the observed disturbance regime",
       subtitle = "disturbance rates by state and climate (Poisson); overlaps as in the data; point estimates; a direction, not a forecast")
save_fig(p, "I13_equilibrium_disturbance", 9, 5, NULL)
