# ==============================================================================
# 25_I8_equilibrium_gradient.R   (no INLA; ~1-2 min; run from the project root)
# I8: long-run (1000-y, undisturbed) composition along the climate gradient vs the observed composition. Plots in five bins by quintile of their mean
#   z-GDD (over their visit intervals); per bin: equilibrium of 17's engine (as C5; point estimates) with GDD, CMI and soil at the bin's medians (as I7),
#   vs the observed composition of the visits in the bin (state at the end of each interval, as C5). Distance = total variation distance (TVD).
# Not a forecast: no disturbance (B15 / I13 adds a disturbance regime), climate fixed, beyond the 55-y record.
# ==============================================================================
src <- readLines(here::here("R", "17_figures_C.R")); eval(parse(text = src[1:(grep("^# ---- propagation engine", src) - 1)]))
scl <- readRDS(file.path(FIT_DIR, "scaling.rds")); zz <- function(x, v) (x - scl[var == v, mean]) / scl[var == v, sd]
# 17's propagate(), generalised to CMI and soil (point estimates, all starts)
prop_eq <- function(zg, zc, zs, tmax = 1000, dt = 1) { nA <- tmax + 1L; ages <- 0:(nA - 1)
  TR <- lapply(S, function(o) lapply(EX[[o]]$lb, function(lb) { X <- PAR[[lb]][1L, , drop = FALSE]; a <- X[, "alpha"]
    list(eta = X[, "b0y"] + gcol(X, "z_gdd_3") * zg + gcol(X, "z_cmi_gs") * zc + gcol(X, "z_soil_pc1") * zs, inc = (ages + 1)^a - ages^a) })); names(TR) <- S
  M <- array(0, c(nS, nS, nA)); for (i in seq_len(nS)) M[i, i, 1] <- 1
  for (st in seq_len(tmax)) { New <- array(0, dim(M))
    for (s in seq_len(nS)) { ex <- EX[[s]]; Ms <- M[, s, , drop = FALSE]; dim(Ms) <- c(nS, nA)
      if (!length(ex$lb)) { New[, s, 2:nA] <- New[, s, 2:nA] + Ms[, 1:(nA - 1)]; next }
      dH <- lapply(TR[[s]], function(tr) exp(tr$eta) * tr$inc); tot <- Reduce(`+`, dH); pl <- 1 - exp(-tot)
      New[, s, 2:nA] <- New[, s, 2:nA] + sweep(Ms, 2, 1 - pl, "*")[, 1:(nA - 1)]
      for (k in seq_along(dH)) { pk <- ifelse(tot > 0, dH[[k]] / tot, 0) * pl; New[, ex$to[k], 1] <- New[, ex$to[k], 1] + Ms %*% pk } }
    M <- New }
  occ <- apply(M, 1:2, sum); list(eq = colMeans(occ), spread = apply(occ, 2, function(v) diff(range(v)))) }
source(here::here("R", "00_groups.R"))   # GRP: thermal groups
v <- vi[, .(id_pe, state = state_end, z_g = zz(gdd_3, "gdd_3"), z_c = zz(cmi_gs, "cmi_gs"), z_s = zz(soil_pc1, "soil_pc1"))]
pm <- v[, .(pg = mean(z_g)), by = id_pe]; pm[, bin := cut(pg, quantile(pg, 0:5 / 5), include.lowest = TRUE, labels = paste0("Q", 1:5))]
v <- merge(v, pm[, .(id_pe, bin)], by = "id_pe")
B <- v[, .(zg = median(z_g), zc = median(z_c), zs = median(z_s, na.rm = TRUE), n_plots = uniqueN(id_pe), n_visits = .N), by = bin][order(bin)]
res <- rbindlist(lapply(seq_len(nrow(B)), function(i) { r <- prop_eq(B$zg[i], B$zc[i], B$zs[i])
  data.table(bin = B$bin[i], state = S, model = r$eq, spread = r$spread) }))
obs <- v[, .N, by = .(bin, state)][, obs := N / sum(N), by = bin]
R <- merge(res, obs[, .(bin, state, obs)], by = c("bin", "state"), all.x = TRUE); R[is.na(obs), obs := 0]
R <- merge(R, B, by = "bin"); fwrite(R, file.path(FIG_DIR, "I", "I8_equilibrium_gradient.csv"))
TV <- R[, .(TVD = 0.5 * sum(abs(model - obs)), max_spread = max(spread)), by = .(bin, zg)]
G <- R[, .(model = sum(model), obs = sum(obs)), by = .(bin, zg, grp = GRP[state])]
cat("bins:\n"); print(B); cat("\nTVD equilibrium vs observed per bin:\n"); print(TV)
cat("\ngroup shares (model equilibrium vs observed):\n"); print(dcast(G, bin + zg ~ grp, value.var = c("obs", "model"), fun.aggregate = function(x) round(x, 2)))
# figure: group shares along the bins
source(here::here("R", "25_I_common.R"))
GL <- melt(G, id.vars = c("bin", "zg", "grp"), variable.name = "src")[, `:=`(src = factor(c(obs = "observed", model = "equilibrium (undisturbed)")[as.character(src)],
  levels = c("observed", "equilibrium (undisturbed)")), xlab = sprintf("%s\n(GDD %+.1f SD)", bin, zg))]
p <- ggplot(GL, aes(xlab, value, colour = grp, group = interaction(grp, src), linetype = src)) + geom_line() + geom_point(size = 1.6) +
  scale_colour_manual(values = GCOL[c("boreal", "pioneer", "regen", "temperate")], name = NULL) + scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(x = "plots binned by mean GDD (quintiles)", y = "share of stands", title = "I8  Observed composition vs the equilibrium of the rates along the gradient",
       subtitle = "equilibrium: 1000 y, undisturbed, climate and soil at the bin's medians; a direction, not a forecast")
save_fig(p, "I8_equilibrium_gradient", 9, 5, NULL)
cat("I8 done\n")
