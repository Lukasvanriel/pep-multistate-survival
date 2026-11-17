# ==============================================================================
# 25_I7_home_climate.R   (no INLA; resumable: re-run until "all I7 runs done"; run from the project root)
# I7 (user 2026-11-11): C1 (occupancy) and D1 (resistance / recovery / reorganisation) with each start state at its OWN climate instead of z = 0.
# Home climate, two definitions (user: compare both):
#   (a) "state"  = median z of GDD, CMI, soil PC1 over all visit intervals starting in the state (the climate its exits are modelled at);
#   (b) "exits"  = the same median over the intervals in which the stand left the state (where its transitions actually happen).
# Engine = 18's propagate_hist (history flags; 1-y steps; disturbance at t = 0 carried over, TSD = t; untermed exits at their undisturbed rate),
#   generalised to CMI and soil and run one start state at a time; after a move the new state is evaluated at the START state's home climate (the plot
#   does not move). Time in state at t = 0: 0 y (main, as C/D). Scenarios: undisturbed, fire, pest, light, heavy harvest. Draws as 18 (50 + point).
# ==============================================================================
src <- readLines(here::here("R", "18_figures_D1_recovery.R")); eval(parse(text = src[1:(grep("^SCD <- CJ", src) - 1)]))
OUTI <- file.path(FIG_DIR, "I"); CI7 <- file.path(OUTI, "cache_I7"); dir.create(CI7, recursive = TRUE, showWarnings = FALSE); t0 <- Sys.time()
scl <- readRDS(file.path(FIT_DIR, "scaling.rds")); zz <- function(x, v) (x - scl[var == v, mean]) / scl[var == v, sd]
H <- vi[, .(z_g = zz(gdd_3, "gdd_3"), z_c = zz(cmi_gs, "cmi_gs"), z_s = zz(soil_pc1, "soil_pc1"), chg = is_change), by = .(state = state_start)]
HOME <- rbind(H[, .(def = "state", zg = median(z_g), zc = median(z_c), zs = median(z_s, na.rm = TRUE), n = .N), by = state],
              H[chg == TRUE, .(def = "exits", zg = median(z_g), zc = median(z_c), zs = median(z_s, na.rm = TRUE), n = .N), by = state])
fwrite(HOME, file.path(OUTI, "I7_home_climate.csv"))
prop_home <- function(st0, zg, zc, zs, dist, tmax = TMAX) {
  D <- length(ROWS); nA <- tmax + 1L; ages <- 0:(nA - 1); X0 <- match(st0, S)
  TR <- lapply(S, function(o) lapply(EX[[o]]$lb, function(lb) { X <- PAR[[lb]][ROWS, , drop = FALSE]; a <- X[, "alpha"]
    eta <- X[, "b0y"] + gcol(X, "z_gdd_3") * zg + gcol(X, "z_cmi_gs") * zc + gcol(X, "z_soil_pc1") * zs; btsd <- rep(0, D)
    if (dist != "none" && TERMS[[dist]][1] %in% colnames(X)) { tm <- TERMS[[dist]]; eta <- eta + X[, tm[1]] + gcol(X, tm[2]) * zg; btsd <- gcol(X, tm[3]) }
    list(eta = eta, btsd = btsd, inc = outer(a, ages + 1, function(aa, x) x^aa) - outer(a, ages, function(aa, x) x^aa)) })); names(TR) <- S
  M <- array(0, c(D, 1, nS, 3, nA)); M[, 1, X0, 1, 1] <- 1; out <- array(NA_real_, c(D, 1, nS, 3, tmax))
  for (st in seq_len(tmax)) { ft <- f_tsd(st - 0.5); New <- array(0, dim(M))
    for (s in seq_len(nS)) { ex <- EX[[s]]; Ms <- M[, , s, , , drop = FALSE]; dim(Ms) <- c(D, 1, 3, nA)
      if (!length(ex$lb)) { New[, , s, , 2:nA] <- New[, , s, , 2:nA] + Ms[, , , 1:(nA - 1)]; next }
      dH <- lapply(TR[[s]], function(tr) exp(tr$eta + tr$btsd * ft) * tr$inc); tot <- Reduce(`+`, dH); pl <- 1 - exp(-tot)
      stayM <- sweep(Ms, c(1, 4), 1 - pl, "*"); New[, , s, , 2:nA] <- New[, , s, , 2:nA] + stayM[, , , 1:(nA - 1)]
      for (k in seq_along(dH)) { pk <- ifelse(tot > 0, dH[[k]] / tot, 0) * pl; fl <- apply(sweep(Ms, c(1, 4), pk, "*"), 1:3, sum); dim(fl) <- c(D, 1, 3)
        j <- ex$to[k]
        if (j == iR) New[, , j, 3, 1] <- New[, , j, 3, 1] + apply(fl, 1:2, sum)
        else { New[, , j, 2, 1] <- New[, , j, 2, 1] + fl[, , 1] + fl[, , 2]; New[, , j, 3, 1] <- New[, , j, 3, 1] + fl[, , 3] } } }
    M <- New; out[, , , , st] <- apply(M, 1:4, sum) }
  out }
RUNS <- CJ(def = c("state", "exits"), start = S, dist = c("none", DIST)); RUNS[, id := paste("i7", def, start, dist, sep = "_")]
for (i in seq_len(nrow(RUNS))) { f <- file.path(CI7, paste0(RUNS$id[i], ".rds")); if (file.exists(f)) next
  if (as.numeric(Sys.time() - t0, units = "secs") > TIME_BUDGET) { cat("time budget reached; re-run to continue\n"); quit(save = "no") }
  h <- HOME[def == RUNS$def[i] & state == RUNS$start[i]]; o <- prop_home(RUNS$start[i], h$zg, h$zc, h$zs, RUNS$dist[i])
  occ <- apply(o, c(1, 3, 5), sum)                                                        # [draw, state, t]
  O <- data.table(state = rep(S, times = TMAX), t = rep(seq_len(TMAX), each = nS), est = as.vector(occ[1, , ]),
                  lo = as.vector(apply(occ[-1, , , drop = FALSE], 2:3, quantile, 0.025)), hi = as.vector(apply(occ[-1, , , drop = FALSE], 2:3, quantile, 0.975)))
  saveRDS(list(occ = O[, `:=`(def = RUNS$def[i], start = RUNS$start[i], dist = RUNS$dist[i])], flags = o[, 1, , , c(50, 100), drop = FALSE]), f)
  cat("  ", RUNS$id[i], "\n") }
cat("all I7 runs done\n")
# ---- summary figure: still in the start state at 50 y, mean climate (z = 0, from C1) vs home climate (both definitions) ----------------------------
source(here::here("R", "25_I_common.R"))
O <- rbindlist(lapply(list.files(CI7, full.names = TRUE), function(f) readRDS(f)$occ))[t == 50]
c0 <- rbindlist(lapply(c(none = "none_z+0", fire = "fire_z0", pest = "pest_z0", harv_light = "harv_light_z0", harv_heavy = "harv_heavy_z0"), function(s)
  fread(file.path(FIG_DIR, "C", sprintf("C1_occupancy_%s.csv", s)))[t == 50]), idcol = "dist")[, def := "mean"]
A <- rbind(O[, .(def, dist, start, state, est, lo, hi)], c0[, .(def, dist, start, state, est, lo, hi)]); fwrite(A, file.path(OUT, "I7_occupancy_50y.csv"))
DF <- c(mean = "mean climate (z = 0)", state = "home: state median", exits = "home: median at exits")
DS <- c(none = "undisturbed", harv_light = "light harvest", harv_heavy = "heavy harvest", pest = "pest", fire = "fire")
P <- A[start == state & dist %in% c("none", "harv_light", "harv_heavy")][, `:=`(def = factor(DF[def], levels = DF), dist = factor(DS[dist], levels = DS),
                                                                             start = factor(start, levels = rev(ORDER)))]
p <- ggplot(P, aes(est, start, colour = def)) + geom_errorbar(aes(xmin = lo, xmax = hi, group = def), width = 0, orientation = "y", position = position_dodge(width = 0.6)) +
  geom_point(aes(group = def), position = position_dodge(width = 0.6), size = 1.6) + facet_wrap(~ dist, nrow = 1) + scale_y_discrete(labels = lab) +
  scale_colour_manual(values = c("grey60", "black", "#E08214"), name = NULL) + scale_x_continuous(labels = function(v) paste0(100 * v, "%")) +
  labs(x = "still in the start state at 50 y", y = NULL, title = "I7  Persistence at mean climate vs at each state's own climate",
       subtitle = "several moves (as C1); 95% CI")
save_fig(p, "I7_home_climate_persistence", 12, 5.5, NULL)
