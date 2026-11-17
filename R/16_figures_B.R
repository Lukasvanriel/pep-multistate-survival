# ==============================================================================
# 16_figures_B.R   (no INLA; ~1-2 min; run from the project root: Rscript R/16_figures_B.R)
# Exploratory figure set B: first move out of a state (semi-Markov, clock reset), derived from the fitted hazards. NOT manuscript figures.
# Needs: 08 (pooled_<tag>_final.csv, pooled_<tag>_final_cov.rds, fit_validity), 05/07 (mi_raw*), 04 (visit_intervals.rds), fits/scaling.rds.
# Writes: data-output/<tag>/figures/B/  (png/pdf + csv)
#
# Model (variant 0, years): exit k of origin i has hazard h_k(u) = alpha_k u^(alpha_k - 1) exp(eta_k(u)), eta_k = b0y_k + beta_k'x(u),
#   b0y = b0 (decade scale) - alpha log(10). All kept exits of an origin compete: S_i(t) = exp(-sum_k H_k(t)); CIF_k(t) = int_0^t h_k S_i.
#   Numerically on a 0.1-y grid: dH = exp(eta(u_mid)) (u_i^alpha - u_{i-1}^alpha); CIF increment = dH * S at the step midpoint.
#   Transitions dropped by the event rule (< 15 events) have hazard 0 (limitation).
# Uncertainty: N_DRAW joint draws per transition from N(pooled INLA mean, W + (1 + 1/M) B) (W: pooled covariance file; B: between-imputation
#   covariance of the INLA means of the 30 used fits, as 09/14); transitions independent (separate fits). Point = pooled estimates; CI = 2.5/97.5 pct.
# Scenarios (user, 2026-11-08/09: "make both, choose later"): CMI, soil z = 0 unless varied; GDD z = -1/0/+1 (and a gradient); disturbance
#   none / fire / pest / light / heavy harvest, entering at the start of the sojourn (TSD = time in state, main) or held at TSD = 10 y (check).
#   Exits WITHOUT that disturbance's term keep their undisturbed rate (assumption, option c, user 2026-11-09); each origin x disturbance carries the
#   share of its baseline exit rate (undisturbed, z = 0, at 10 y) that comes from such exits; share > 0.25 is FLAGGED in the figures.
#   Horizons: 50 y (inside the 55-y record; main), curves to 100 y (beyond 55 y = extrapolation, shaded).
# Quantities: P(stay), CIF to each destination, restricted mean residence time to 50 y (RMST50), median residence time, next-state probability
#   P(next = j | moved by 50 y), expected thermal shift of the first move = sum_j P(next = j) (STI_j - STI_origin)  (STI as in 12: same plots).
# B0 model check: undisturbed visit intervals with an observed entry into the current state (time in state at start from the entry midpoint);
#   observed vs predicted share changing to a modelled destination during the interval (point estimates), split by time since the entry.
# ==============================================================================
source(here::here("R", "13_fig_helpers.R"))
source(here::here("R", "covariate_spec.R")); f_tsd <- tsd_fun(TSD_FORM)
OUT <- file.path(FIG_DIR, "B"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
N_DRAW <- 500L; DT <- 0.1; TMAX <- 100; HOR <- 50; FLAG_SHARE <- 0.25; set.seed(20261109)
DIST <- c("fire", "pest", "harv_light", "harv_heavy"); DL <- c(none = "undisturbed", fire = "fire", pest = "pest", harv_light = "light harvest", harv_heavy = "heavy harvest")
TERMS <- list(fire = c("fire_f", "gdd_x_fire", "tsd_fire"), pest = c("pest_f", "gdd_x_pest", "tsd_pest"),
              harv_light = c("harv_light", "gdd_x_harv_light", "tsd_harv_light"), harv_heavy = c("harv_heavy", "gdd_x_harv_heavy", "tsd_harv_heavy"))
STATES_O <- ORDER

# ---- parameter draws per transition ------------------------------------------------------------------------------------------------------------------------
cvl <- readRDS(file.path(RESULTS_DIR, sprintf("pooled_%s_final_cov.rds", RUN_TAG)))
fv  <- fread(file.path(RESULTS_DIR, sprintf("fit_validity_%s.csv", RUN_TAG)))[used == TRUE, .(label, m)]
mir <- rbindlist(lapply(file.path(FIT_DIR, sprintf(c("mi_raw_%s.csv", "mi_raw_redraw_%s.csv"), RUN_TAG)), function(f) if (file.exists(f)) fread(f)[, .(label, m, parameter, est)]))
mir <- unique(mir[fv, on = c("label", "m"), nomatch = NULL], by = c("label", "m", "parameter"))
PAR <- list()
for (lb in modelled) {
  nm <- names(cvl[[lb]]$mean); pe <- setNames(pooled[label == lb, pooled_est], pooled[label == lb, parameter])[nm]; stopifnot(!anyNA(pe))
  Q <- as.matrix(dcast(mir[label == lb], m ~ parameter, value.var = "est")[, ..nm]); B <- stats::cov(Q); M <- nrow(Q)
  Tm <- cvl[[lb]]$W[nm, nm] + (1 + 1 / M) * B; Tm <- (Tm + t(Tm)) / 2
  ev <- eigen(Tm, symmetric = TRUE); L <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0))) 
  X <- rbind(pe, sweep(matrix(rnorm(N_DRAW * length(nm)), N_DRAW) %*% t(L), 2, pe, "+")); colnames(X) <- nm   # row 1 = point estimate
  stopifnot(all(X[, "alpha"] > 0))
  X <- cbind(X, b0y = X[, "(Intercept)"] - X[, "alpha"] * log(10)); PAR[[lb]] <- X
}
cat(sprintf("draws: %d per transition (+ point), %d transitions\n", N_DRAW, length(PAR)))
exits_of <- function(o) trow[from == o][order(match(to, ORDER)), label]
gcol <- function(X, nmv) if (nmv %in% colnames(X)) X[, nmv] else rep(0, nrow(X))   # one value per draw

# ---- state thermal index (as 12) and flags ----------------------------------------------------------------------------------------------------------------
sc <- readRDS(file.path(FIT_DIR, "scaling.rds")); zf <- function(x, v) (x - sc[var == v, mean]) / sc[var == v, sd]
STI <- vi[, .(STI = mean(zf(gdd_3, "gdd_3"))), by = .(state = state_end)]; STI <- setNames(STI$STI, STI$state)
h10 <- sapply(modelled, function(lb) { X <- PAR[[lb]]; X[1, "alpha"] * 10^(X[1, "alpha"] - 1) * exp(X[1, "b0y"]) })
has_term <- function(lb, d) TERMS[[d]][1] %in% colnames(PAR[[lb]])
share <- CJ(origin = STATES_O, dist = DIST)
share[, c("share_assumed", "any_term") := { ex <- exits_of(origin); ht <- vapply(ex, has_term, logical(1), d = dist)
  list(sum(h10[ex][!ht]) / sum(h10[ex]), any(ht)) }, by = .(origin, dist)]
share[, flag := share_assumed > FLAG_SHARE]; fwrite(share, file.path(OUT, "B_flags_share_assumed.csv"))

# ---- engine ------------------------------------------------------------------------------------------------------------------------------------------------
TG <- seq(0, TMAX, by = DT); UM <- (TG[-1] + TG[-length(TG)]) / 2; iH <- which.min(abs(TG - HOR)) - 1L
run_origin <- function(o, zg = 0, zc = 0, zs = 0, dist = "none", mode = "entry", keep_curves = FALSE) {
  ex <- exits_of(o); dH <- list()
  for (lb in ex) { X <- PAR[[lb]]; a <- X[, "alpha"]
    eta <- X[, "b0y"] + gcol(X, "z_gdd_3") * zg + gcol(X, "z_cmi_gs") * zc + gcol(X, "z_soil_pc1") * zs
    tsd_part <- 0
    if (dist != "none" && TERMS[[dist]][1] %in% colnames(X)) { tm <- TERMS[[dist]]
      eta <- eta + X[, tm[1]] + gcol(X, tm[2]) * zg
      ft <- if (mode == "entry") f_tsd(UM) else rep(f_tsd(10), length(UM)); tsd_part <- outer(gcol(X, tm[3]), ft) }
    pw <- outer(a, TG, function(aa, tt) tt^aa); inc <- pw[, -1, drop = FALSE] - pw[, -ncol(pw), drop = FALSE]
    dH[[lb]] <- exp(eta + tsd_part) * inc }                                    # rows = draws (1 = point), cols = steps
  Htot <- Reduce(`+`, dH); Hc <- t(apply(Htot, 1, cumsum)); Hprev <- cbind(0, Hc[, -ncol(Hc), drop = FALSE])
  Smid <- exp(-(Hprev + Hc) / 2); S <- exp(-Hc)
  cif <- lapply(dH, function(m) t(apply(m * Smid, 1, cumsum)))
  stay <- S[, iH]; c50 <- sapply(cif, function(m) m[, iH]); if (!is.matrix(c50)) c50 <- matrix(c50, nrow = 1)
  rmst <- rowSums(Smid[, seq_len(iH), drop = FALSE]) * DT
  med <- apply(S, 1, function(s) { k <- which(s <= 0.5)[1]; if (is.na(k)) NA_real_ else TG[k + 1] })
  nxt <- c50 / (1 - stay); dest <- sub(".*->", "", ex)
  shift <- drop(nxt %*% (STI[dest] - STI[o]))
  out <- list(dest = dest, stay = stay, cif50 = c50, rmst = rmst, med = med, nxt = nxt, shift = shift)
  if (keep_curves) { step <- seq(10, ncol(S), by = 10); out$curves <- data.table(t = TG[step + 1], stay = S[1, step],
                        do.call(cbind, setNames(lapply(cif, function(m) m[1, step]), dest))) }
  out
}
qs <- function(x) list(est = x[1], lo = unname(quantile(x[-1], 0.025, na.rm = TRUE)), hi = unname(quantile(x[-1], 0.975, na.rm = TRUE)))
tidy <- function(r, o, keys) {
  base <- c(list(origin = o), keys)
  rbind(as.data.table(c(base, list(quantity = "stay50", dest = NA_character_), qs(r$stay))),
        as.data.table(c(base, list(quantity = "rmst50", dest = NA_character_), qs(r$rmst))),
        as.data.table(c(base, list(quantity = "median", dest = NA_character_), qs(r$med))),
        as.data.table(c(base, list(quantity = "shift50", dest = NA_character_), qs(r$shift))),
        rbindlist(lapply(seq_along(r$dest), function(j) rbind(as.data.table(c(base, list(quantity = "cif50", dest = r$dest[j]), qs(r$cif50[, j]))),
                                                             as.data.table(c(base, list(quantity = "next50", dest = r$dest[j]), qs(r$nxt[, j])))))))
}

# ---- run all scenarios -------------------------------------------------------------------------------------------------------------------------------------
res <- list(); draws <- list(); curves <- list(); t0 <- Sys.time()
for (o in STATES_O) {
  for (zg in c(-1, 0, 1)) { r <- run_origin(o, zg = zg, keep_curves = TRUE); k <- list(   # curves for GDD -1/0/+1 (user 2026-11-10, B1 per-origin figure)
    scen = "none", dist = "none", mode = "none", zg = zg, zc = 0, zs = 0)
    res[[length(res) + 1]] <- tidy(r, o, k); draws[[paste(o, "none", "none", zg)]] <- r[c("stay", "rmst", "shift")]
    curves[[length(curves) + 1]] <- r$curves[, `:=`(origin = o, dist = "none", zg = zg)] }
  for (v in c("zc", "zs")) for (zz in c(-1, 1)) { a <- list(o = o); a[[v]] <- zz; r <- do.call(run_origin, a)
    k <- list(scen = v, dist = "none", mode = "none", zg = 0, zc = if (v == "zc") zz else 0, zs = if (v == "zs") zz else 0)
    res[[length(res) + 1]] <- tidy(r, o, k); draws[[paste(o, v, zz)]] <- r[c("stay", "rmst", "shift")] }
  for (zg in seq(-2, 2, by = 0.25)) { r <- run_origin(o, zg = zg)
    res[[length(res) + 1]] <- tidy(r, o, list(scen = "gdd_gradient", dist = "none", mode = "none", zg = zg, zc = 0, zs = 0)) }
  for (d in DIST) { if (!share[origin == o & dist == d, any_term]) next
    for (md in c("entry", "fixed10")) for (zg in c(0, 1)) { r <- run_origin(o, zg = zg, dist = d, mode = md, keep_curves = (zg == 0 && md == "entry"))
      res[[length(res) + 1]] <- tidy(r, o, list(scen = "dist", dist = d, mode = md, zg = zg, zc = 0, zs = 0))
      draws[[paste(o, d, md, zg)]] <- r[c("stay", "rmst", "shift")]
      if (zg == 0 && md == "entry") curves[[length(curves) + 1]] <- r$curves[, `:=`(origin = o, dist = d, zg = 0)] } }
}
res <- rbindlist(res, fill = TRUE); res <- merge(res, share[, .(origin, dist, share_assumed, flag)], by = c("origin", "dist"), all.x = TRUE)
fwrite(res, file.path(OUT, "B_quantities.csv")); cat(sprintf("scenarios done in %.0f s\n", as.numeric(Sys.time() - t0, units = "secs")))
curves <- rbindlist(curves, fill = TRUE); fwrite(curves, file.path(OUT, "B1_curves_point.csv"))
# check: P(stay) + sum CIF = 1 at 50 y (point estimates)
ck <- res[quantity %in% c("stay50", "cif50"), .(tot = sum(est)), by = .(origin, scen, dist, mode, zg, zc, zs)]
cat(sprintf("check stay + sum CIF at 50 y: range %.6f - %.6f\n", min(ck$tot), max(ck$tot))); stopifnot(all(abs(ck$tot - 1) < 1e-3))

# ---- differences (B6) and interactions on the probability scale (B7), per draw ---------------------------------------------------------------------------
dd <- function(a, b, q) { x <- draws[[a]][[q]] - draws[[b]][[q]]; qs(x) }
B6 <- rbindlist(lapply(STATES_O, function(o) { base <- paste(o, "none", "none", 0)
  rows <- list(c("GDD +1 SD", paste(o, "none", "none", 1)), c("CMI +1 SD", paste(o, "zc", 1)), c("soil PC1 +1 SD", paste(o, "zs", 1)))
  for (d in DIST) if (!is.null(draws[[paste(o, d, "entry", 0)]])) rows[[length(rows) + 1]] <- c(DL[[d]], paste(o, d, "entry", 0))
  rbindlist(lapply(rows, function(rr) rbindlist(lapply(c("rmst", "stay", "shift"), function(q) as.data.table(c(list(origin = o, contrast = rr[1], quantity = q), dd(rr[2], base, q)))))))
}))
B6 <- merge(B6, share[, .(origin, contrast = DL[dist], share_assumed, flag)], by = c("origin", "contrast"), all.x = TRUE); B6[is.na(flag), flag := FALSE]
fwrite(B6, file.path(OUT, "B6_climate_vs_disturbance.csv"))
B7 <- rbindlist(lapply(STATES_O, function(o) rbindlist(lapply(DIST, function(d) { k1 <- paste(o, d, "entry", 1); k0 <- paste(o, d, "entry", 0)
  if (is.null(draws[[k1]])) return(NULL); n1 <- paste(o, "none", "none", 1); n0 <- paste(o, "none", "none", 0)
  rbindlist(lapply(c("stay", "rmst", "shift"), function(q) { x <- (draws[[k1]][[q]] - draws[[k0]][[q]]) - (draws[[n1]][[q]] - draws[[n0]][[q]])
    as.data.table(c(list(origin = o, dist = d, quantity = q), qs(x))) })) }))))
B7 <- merge(B7, share[, .(origin, dist, share_assumed, flag)], by = c("origin", "dist")); fwrite(B7, file.path(OUT, "B7_interaction_probability_scale.csv"))

# ---- B0 model check --------------------------------------------------------------------------------------------------------------------------------------------
# Undisturbed visit intervals of plots whose entry into the current state was observed (entry = midpoint of the interval with the last observed
#   change); time in state at the interval start a = year_start - entry. Split: entered in the previous interval vs >= 2 intervals earlier
#   (2026-11-09 diagnosis: half of the changes right after a change are reversals A->B->A, which the model does not reproduce).
setorder(vi, id_pe, interval)
vi[, entry_mid := { e <- rep(NA_real_, .N); cur <- NA_real_; for (i in seq_len(.N)) { e[i] <- cur; if (is_change[i]) cur <- (year_start[i] + year_end[i]) / 2 }; e }, by = id_pe]
vi[, prev_change := shift(is_change, fill = FALSE), by = id_pe]
b0 <- vi[!is.na(entry_mid) & fire == 0 & pest == 0 & harv == 0 & state_start %in% STATES_O]
b0[, `:=`(a = year_start - entry_mid, L = year_end - year_start, zg = zf(gdd_3, "gdd_3"), zc = zf(cmi_gs, "cmi_gs"), zs = zf(soil_pc1, "soil_pc1"),
          grp = fifelse(prev_change, "entered in the previous interval", "entered >= 2 intervals earlier"))]
pred_int <- function(o, a, L, zg, zc, zs) { ex <- exits_of(o); u <- seq(a, a + L, length.out = 41)
  dH <- sapply(ex, function(lb) { X <- PAR[[lb]][1, , drop = FALSE]; e <- X[, "b0y"] + gcol(X, "z_gdd_3") * zg + gcol(X, "z_cmi_gs") * zc + gcol(X, "z_soil_pc1") * zs
    exp(e) * diff(u^X[, "alpha"]) }); if (!is.matrix(dH)) dH <- matrix(dH, nrow = 40)
  Ht <- cumsum(rowSums(dH)); Sm <- exp(-(c(0, Ht[-40]) + Ht) / 2); setNames(colSums(dH * Sm), sub(".*->", "", ex)) }
P <- rbindlist(lapply(seq_len(nrow(b0)), function(i) { r <- b0[i]; p <- pred_int(r$state_start, r$a, r$L, r$zg, r$zc, r$zs)
  data.table(i = i, origin = r$state_start, grp = r$grp, dest = names(p), p = p) }))
obs <- b0[, .(i = .I, origin = state_start, grp, prev_state = NA_character_, dest_obs = fifelse(is_change, state_end, NA_character_))]
obs[, modelled_dest := !is.na(dest_obs) & paste0(origin, "->", dest_obs) %in% modelled]
B0o <- merge(obs[, .(n = .N, obs_change = sum(modelled_dest)), by = .(origin, grp)], P[, .(pred_change = sum(p)), by = .(origin, grp)], by = c("origin", "grp"))
B0o[, `:=`(obs_p = obs_change / n, pred_p = pred_change / n)]
B0o[, c("obs_lo", "obs_hi") := { z <- 1.96; ph <- obs_p; den <- 1 + z^2 / n; ce <- (ph + z^2 / (2 * n)) / den; hw <- z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2)) / den; list(ce - hw, ce + hw) }]
B0t <- merge(obs[modelled_dest == TRUE, .(obs = .N), by = .(origin, grp, dest = dest_obs)], P[, .(pred = sum(p)), by = .(origin, grp, dest)], by = c("origin", "grp", "dest"), all = TRUE)
B0t[is.na(obs), obs := 0L]; fwrite(B0o, file.path(OUT, "B0_check_by_origin.csv")); fwrite(B0t, file.path(OUT, "B0_check_by_transition.csv"))
print(B0o[, .(n = sum(n), observed = sum(obs_change), predicted = round(sum(pred_change), 1)), by = grp][, ratio := round(observed / predicted, 2)][])
saveRDS(list(res = res, B6 = B6, B7 = B7, B0o = B0o, B0t = B0t, curves = curves, share = share, STI = STI), file.path(OUT, "B_results.rds"))
cat("computation done\n")
