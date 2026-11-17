# ==============================================================================
# 17b_figures_C_plots.R   (run after 17_figures_C.R has printed "all scenarios done"; from the project root)
# Figures C1-C5 from figures/C/cache/. C4 (pathways) and C5 (long-run composition) are computed here with POINT estimates only.
# C4: 20,000 simulated stands per start state, 1-y steps, 100 y, undisturbed, GDD z = 0; the first three moves are recorded.
# C5: mass propagation (as 17) to 1000 y, undisturbed, GDD z = -1/0/+1; compared with the share of all observed visits per state.
# ==============================================================================
src <- readLines(here::here("R", "17_figures_C.R")); eval(parse(text = src[1:(grep("^# ---- scenario list", src) - 1)]))
olab <- function(x) gsub("_", " ", x); fo <- function(x) factor(x, levels = ORDER)
OCC <- rbindlist(lapply(list.files(CACHE, "^occ_.*rds$", full.names = TRUE), readRDS))
FPS <- rbindlist(lapply(list.files(CACHE, "^fp_.*rds$", full.names = TRUE), readRDS))
share <- fread(file.path(FIG_DIR, "B", "B_flags_share_assumed.csv"))
scz <- readRDS(file.path(FIT_DIR, "scaling.rds")); zf <- function(x, v) (x - scz[var == v, mean]) / scz[var == v, sd]
STI <- vi[, .(STI = mean(zf(gdd_3, "gdd_3"))), by = .(state = state_end)]; STI <- setNames(STI$STI, STI$state)
DLAB <- c(none = "undisturbed", fire = "fire", pest = "pest", harv_light = "light harvest", harv_heavy = "heavy harvest")
OUT <- OUTC

# ---- C1 occupancy -----------------------------------------------------------------------------------------------------------------------------------------
occ_fig <- function(sc, ttl, labs_fun = olab) { x <- OCC[scen == sc][, `:=`(start = fo(start), state = fo(state))]
  ggplot(x, aes(t, est, fill = state)) +
    geom_area(colour = "white", linewidth = 0.1) + geom_vline(xintercept = 55, linetype = 2, linewidth = 0.4) +
    scale_fill_manual(values = PAL, labels = olab, name = "state at time t") + scale_y_continuous(labels = function(v) paste0(100 * v, "%")) +
    facet_wrap(~ start, ncol = 5, labeller = as_labeller(labs_fun)) +
    labs(x = "years since entering the start state", y = "probability", title = ttl,
         subtitle = "several moves, returns included; dashed = end of record (55 y)") }
for (sc in c("none_z-1", "none_z+0", "none_z+1")) save_fig(occ_fig(sc, sprintf("C1  State over time, by start state (undisturbed, %s)", gsd(as.integer(sub("none_z", "", sc))))),
  paste0("C1_occupancy_", sc), 16, 9.5, OCC[scen == sc, .(start, state, t, est, lo, hi)])
for (d in DIST) { sc <- paste0(d, "_z0"); fl <- share[dist == d]; lf <- function(x) { s <- fl[match(x, origin)]
    paste0(olab(x), "\n", sprintf("%.0f%% assumed", 100 * s$share_assumed), ifelse(s$flag, " *", "")) }
  save_fig(occ_fig(sc, sprintf("C1  State over time after %s at t = 0", DLAB[[d]]), lf) +
             labs(caption = "panel: share of the start state's exit rate assumed unaffected by the disturbance (* = > 25%)"),
           paste0("C1_occupancy_", sc), 16, 10, OCC[scen == sc, .(start, state, t, est, lo, hi)]) }

# C1d (user 2026-11-10): one figure per start state, all scenarios side by side (undisturbed GDD -1 / 0 / +1, then each disturbance at t = 0, GDD 0).
#   The disturbance carries over after state changes, so no panel is empty; when the START state has no term for that disturbance the panel title says so.
SCC <- c(`none_z-1` = "undisturbed\nGDD -1 SD", `none_z+0` = "undisturbed\nGDD 0", `none_z+1` = "undisturbed\nGDD +1 SD",
         fire_z0 = "fire", pest_z0 = "pest", harv_light_z0 = "light harvest", harv_heavy_z0 = "heavy harvest")
DKEY <- c(fire_z0 = "fire", pest_z0 = "pest", harv_light_z0 = "harv_light", harv_heavy_z0 = "harv_heavy")
for (o in ORDER) { x <- OCC[start == o & scen %in% names(SCC)]
  lab_o <- SCC; for (k in names(DKEY)) if (!share[origin == o & dist == DKEY[[k]], any_term]) lab_o[[k]] <- paste0(SCC[[k]], "\n(no term for this state)")
  x[, `:=`(panel = factor(lab_o[scen], levels = lab_o), state = fo(state))]
  p <- ggplot(x, aes(t, est, fill = state)) + geom_area(colour = "white", linewidth = 0.1) + geom_vline(xintercept = 55, linetype = 2, linewidth = 0.4) +
    scale_fill_manual(values = PAL, labels = olab, name = "state at time t") + scale_y_continuous(labels = function(v) paste0(100 * v, "%")) +
    facet_wrap(~ panel, nrow = 1) +
    labs(x = "years since t = 0", y = "probability", title = sprintf("C1  State over time, starting in %s", olab(o)),
         subtitle = "several moves, returns included; disturbance at t = 0 carries over after state changes (GDD 0); exits without a term for it keep their undisturbed rate; dashed = end of record (55 y)")
  save_fig(p, paste0("C1d_occupancy_", o), 17, 3.8, NULL) }

# ---- C2 thermal index of the expected composition -------------------------------------------------------------------------------------------------------------
draws_cti <- function(sc) { x <- OCC[scen == sc]; x[, sti := STI[state]] ; x[, .(cti = sum(est * sti), cti_lo = NA_real_), by = .(start, t, zg)] }
ct <- OCC[dist == "none"][, sti := STI[state]][, .(cti = sum(est * sti)), by = .(start, t, zg)]
ct[, `:=`(start = fo(start), dcti = cti - STI[as.character(start)])]
p <- ggplot(ct, aes(t, dcti, colour = zg, group = zg)) + geom_vline(xintercept = 55, linetype = 2, linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") + geom_line() + scale_colour_gradient2(low = "#2166AC", mid = "grey60", high = "#B2182B", name = "GDD (SD)") +
  facet_wrap(~ start, ncol = 5, labeller = as_labeller(olab)) +
  labs(x = "years since start", y = "thermal shift of the expected composition (ΔSTI)", title = "C2  Thermal shift of the expected composition over time (undisturbed)",
       subtitle = "> 0 = toward warmer-affinity states; dashed = end of record (55 y)")
save_fig(p, "C2_thermal_shift_over_time", 16, 9.5, ct)
# at 50 y with CI: from the draw-level occupancy is not cached; CI from the occupancy CIs is not additive, so the CI here uses the 2.5/97.5 pct of the
# occupancy-weighted index computed per draw is NOT available -> point estimates only, stated in the subtitle
c50 <- ct[t == 50]
p <- ggplot(c50, aes(zg, dcti)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey50") + geom_line(colour = "#B2182B") + geom_point(size = 1.2) +
  facet_wrap(~ start, ncol = 5, labeller = as_labeller(olab)) + labs(x = "GDD (SD from the mean)", y = "thermal shift after 50 y (ΔSTI)",
  title = "C2b  Thermal shift of the expected composition after 50 y, along GDD (undisturbed)")
save_fig(p, "C2b_thermal_shift_50y_vs_gdd", 16, 9.5, c50)

# ---- C3 first passage -----------------------------------------------------------------------------------------------------------------------------------------
fp50 <- FPS[t == 50 & start != target]
for (zz in c(-1, 0, 1)) { x <- fp50[zg == zz]; g <- CJ(start = ORDER, target = ORDER); g <- merge(g, x, by = c("start", "target"), all.x = TRUE)
  g[, `:=`(start = factor(start, levels = rev(ORDER)), target = factor(target, levels = ORDER))]
  g[, txt := fifelse(is.na(est), "", sprintf("%.0f\n[%.0f-%.0f]", 100 * est, 100 * lo, 100 * hi))]; g[, tc := fifelse(!is.na(est) & est > 0.55, "white", "black")]
  p <- ggplot(g, aes(target, start)) + geom_tile(aes(fill = est), colour = "white") + geom_text(aes(label = txt, colour = tc), size = 1.9, lineheight = 0.85) + scale_colour_identity() +
    scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 1), na.value = "grey95", labels = function(v) paste0(100 * v, "%"), name = "probability") +
    scale_x_discrete(labels = olab, position = "top") + scale_y_discrete(labels = olab) + coord_equal() + theme(axis.text.x = element_text(angle = 45, hjust = 0), panel.grid = element_blank()) +
    labs(x = "state reached", y = "start state", title = sprintf("C3  Probability of reaching each state within 50 y (undisturbed, %s)", gsd(zz)),
         subtitle = "% [95% CI]; any route")
  save_fig(p, sprintf("C3_first_passage_50y_gdd%+d", zz), 9, 8.5, x) }
med <- FPS[start != target, .(t_half = { k <- which(est >= 0.5)[1]; if (is.na(k)) NA_real_ else t[k] }, p100 = est[t == max(t)]), by = .(start, target, zg)]
fwrite(med, file.path(OUT, "C3_first_passage_median_time.csv"))

# ---- C4 pathways (simulation, point estimates) -------------------------------------------------------------------------------------------------------------
set.seed(1); NS <- 20000L
sim_paths <- function(o, zg = 0) { st <- rep(match(o, S), NS); age <- rep(0, NS); path <- rep(o, NS); nmov <- rep(0L, NS)
  for (yr in 1:100) { st0 <- st; for (s in unique(st0)) { ii <- which(st0 == s); ex <- EX[[s]]; if (!length(ex$lb)) next   # states at the start of the year only
      dH <- sapply(ex$lb, function(lb) { X <- PAR[[lb]][1, , drop = FALSE]; exp(X[, "b0y"] + gcol(X, "z_gdd_3") * zg) * ((age[ii] + 1)^X[, "alpha"] - age[ii]^X[, "alpha"]) })
      if (!is.matrix(dH)) dH <- matrix(dH, nrow = length(ii)); tot <- rowSums(dH); mv <- runif(length(ii)) < 1 - exp(-tot)
      if (any(mv)) { jj <- ii[mv]; pk <- dH[mv, , drop = FALSE] / tot[mv]; k <- apply(pk, 1, function(p) sample.int(length(p), 1, prob = p))
        nst <- ex$to[k]; keep3 <- nmov[jj] < 3; path[jj[keep3]] <- paste(path[jj[keep3]], S[nst[keep3]], sep = " > "); nmov[jj] <- nmov[jj] + 1L
        st[jj] <- nst; age[jj] <- -1 } }
    age <- age + 1 }
  data.table(start = o, path = path, n_moves = nmov) }
pw <- rbindlist(lapply(S, sim_paths)); fwrite(pw[, .N, by = .(start, path)][, share := N / NS][order(start, -share)], file.path(OUT, "C4_pathways.csv"))
top <- pw[, .N, by = .(start, path)][, share := N / NS][order(-share), head(.SD, 6), by = start][, start := fo(start)]
top[, path_lab := gsub("_", " ", path)]
p <- ggplot(top, aes(share, reorder(paste(start, path_lab, sep = "|"), share))) + geom_col(fill = "grey40") +
  geom_text(aes(label = sprintf("%.0f%%", 100 * share)), hjust = -0.1, size = 2.2) + scale_y_discrete(labels = function(x) sub("^[^|]*\\|", "", x)) +
  scale_x_continuous(labels = function(v) paste0(100 * v, "%"), expand = expansion(mult = c(0, 0.15))) + facet_wrap(~ start, scales = "free_y", ncol = 3, labeller = as_labeller(olab)) +
  labs(x = "share of simulated stands", y = NULL, title = "C4  Most frequent pathways within 100 y (undisturbed)",
       subtitle = "first three moves; 20,000 simulated stands per start state; a single state = no move") + theme(axis.text.y = element_text(size = 6))
save_fig(p, "C4_pathways", 15, 14, NULL)

# ---- C5 long-run composition ---------------------------------------------------------------------------------------------------------------------------------
lr <- rbindlist(lapply(c(-1, 0, 1), function(zz) { r <- propagate(zg = zz, tmax = 1000, rows = 1L); o <- r$occ[1, , , 1000]
  data.table(state = S, est = colMeans(o), spread_across_starts = apply(o, 2, function(v) diff(range(v))), zg = zz) }))
obs <- vi[, .N, by = .(state = state_end)][, est := N / sum(N)][, .(state, est, zg = NA_real_, spread_across_starts = NA_real_)]
lr2 <- rbind(lr[, src := paste0("model, ", gsd(zg))], obs[, src := "observed visits"], fill = TRUE)
lr2[, `:=`(state = fo(state), src = factor(src, levels = c("observed visits", paste0("model, ", gsd(c(-1, 0, 1))))))]
fwrite(lr2, file.path(OUT, "C5_long_run.csv")); cat(sprintf("C5: max spread across start states at 1000 y = %.4f (0 = converged)\n", max(lr$spread_across_starts)))
p <- ggplot(lr2, aes(est, state, fill = src)) + geom_col(position = position_dodge(width = 0.8), width = 0.75) + scale_y_discrete(labels = olab, limits = rev) +
  scale_fill_manual(values = c(`observed visits` = "grey30", `model, GDD -1 SD` = "#2166AC", `model, GDD 0` = "grey65", `model, GDD +1 SD` = "#B2182B"), name = NULL) +
  scale_x_continuous(labels = function(v) paste0(100 * v, "%")) +
  labs(x = "share", y = NULL, title = "C5  Long-run composition implied by the rates vs observed visits",
       subtitle = "undisturbed, after 1000 y; a property of the rates, not a forecast")
save_fig(p, "C5_long_run_composition", 10, 7, NULL)
cat("C figures done\n")
