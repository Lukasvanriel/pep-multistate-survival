# ==============================================================================
# 23_extended_quantities.R   (no INLA; resumable: re-run until it prints "all runs done"; run from the project root)
# Exploratory set E (user 2026-11-10: items 1-7). Same model, draws and conventions as 16/17 (variant 0, years; competing exits; dropped transitions
#   have hazard 0; CMI, soil z = 0; a disturbance enters at t = 0, carries over after state changes, TSD = t; exits without that disturbance's term keep
#   their undisturbed rate). Draws = 17's PAR (row 1 = point, 100 joint draws); CI = 2.5/97.5 pct.
# E1 P(next = j | leave at time t): single sojourn, shares of the exit hazards per 1-y step; undisturbed GDD -1/0/+1 and each disturbance at entry (z 0).
# E2 "how long" vs "which way": hazards of each state split into the total exit hazard (holding) and the destination shares (jump). Mixed runs use the
#   total from scenario H and the shares from scenario J. Shapley split (average of both orders): hold = [(P,B)-(B,B) + (P,P)-(B,P)]/2,
#   jump = total - hold. Perturbations: GDD +1 and -1 vs 0 (undisturbed); each disturbance vs undisturbed (z 0). Outcomes: first-move thermal shift
#   within 50 y (as B) and thermal index of the composition at 50 y (as C2), per start state.
# E3 elasticity: each transition's hazard x1.1 and x0.9 (intercept) and its alpha x1.1 and x0.9; point estimates; undisturbed, z 0; start = observed
#   first-visit composition (plot_info), time in state 0 at t = 0 (as D4 main); outcomes at 50 y: composition and its thermal index.
# E4 expected years in each state within 50 and 100 y (occupancy summed, returns included); scenarios as C1 (undisturbed -1/0/+1, each disturbance z 0).
# E5 expected number of state changes within 50/100 y; P(first return to the start state within t) after leaving (return flows absorbed), same scenarios.
# E6 convergence: point estimates, undisturbed, z -1/0/+1, from the observed first-visit composition, 1000 y; total variation distance to the 1000-y
#   composition; years until < 0.10.
# E7 successional direction: transition counts by (from, to) within 50/100 y (E5 runs) -> shade-tolerance shift of moves (Niinemets & Valladares 2006
#   scores, as 12; moves involving regen, other_pioneer, other_temperate are unscored and reported separately) and flows between tolerance groups (as 12;
#   other_temperate = "unclassified"). E1 also gives the tolerance shift of the next state as a function of exit time.
# ==============================================================================
src <- readLines(here::here("R", "17_figures_C.R")); eval(parse(text = src[1:(grep("^# ---- scenario list", src) - 1)]))
OUTE <- file.path(FIG_DIR, "E"); CE <- file.path(OUTE, "cache"); dir.create(CE, recursive = TRUE, showWarnings = FALSE); OUT <- OUTE
TIME_BUDGET <- 150; HOR <- c(50, 100)
sc <- readRDS(file.path(FIT_DIR, "scaling.rds")); zf <- function(x, v) (x - sc[var == v, mean]) / sc[var == v, sd]
STI <- vi[, .(STI = mean(zf(gdd_3, "gdd_3"))), by = .(state = state_end)]; STI <- setNames(STI$STI, STI$state)[S]
TOL <- c(fir = 5.01, black_spruce = 4.08, white_spruce = 4.15, jack_pine = 1.36, cedar = 3.45, sugar_maple = 4.76, red_maple = 3.44, yellow_birch = 3.17,
         white_birch = 1.54, trembling_aspen = 1.21, other_pioneer = 1.18, other_temperate = NA, regen = NA)[S]   # other_pioneer as 25_I_common (user 2026-11-16)
source(here::here("R", "00_groups.R")); TGRP <- TOL_GROUP[S]   # shade-tolerance groups, ordered as S
pinf <- readRDS(file.path(TABLE_DIR, "plot_info.rds")); W0 <- pinf[first_state %in% S, .N, by = first_state]; W0 <- setNames(W0$N / sum(W0$N), W0$first_state)[S]; W0[is.na(W0)] <- 0
saveRDS(list(STI = STI, TOL = TOL, TGRP = TGRP, W0 = W0, S = S), file.path(CE, "meta.rds"))
SCN <- list(none_zm1 = list(zg = -1, dist = "none"), none_z0 = list(zg = 0, dist = "none"), none_zp1 = list(zg = 1, dist = "none"),
            fire = list(zg = 0, dist = "fire"), pest = list(zg = 0, dist = "pest"), harv_light = list(zg = 0, dist = "harv_light"), harv_heavy = list(zg = 0, dist = "harv_heavy"))
PERT <- list(gdd_p1 = c("none_zp1", "none_z0"), gdd_m1 = c("none_zm1", "none_z0"), fire = c("fire", "none_z0"), pest = c("pest", "none_z0"),
             harv_light = c("harv_light", "none_z0"), harv_heavy = c("harv_heavy", "none_z0"))           # perturbed, base
over_time <- function() as.numeric(Sys.time() - t_start, units = "secs") > TIME_BUDGET
eta_of <- function(X, sc0, mod = NULL) { eta <- X[, "b0y"] + gcol(X, "z_gdd_3") * sc0$zg; btsd <- rep(0, nrow(X))
  if (sc0$dist != "none" && TERMS[[sc0$dist]][1] %in% colnames(X)) { tm <- TERMS[[sc0$dist]]; eta <- eta + X[, tm[1]] + gcol(X, tm[2]) * sc0$zg; btsd <- gcol(X, tm[3]) }
  a <- X[, "alpha"]; if (!is.null(mod)) { if (mod$type == "b0") eta <- eta + log(mod$f) else a <- a * mod$f }
  list(eta = eta, btsd = btsd, a = a) }

# ---- E1 / E2 first move: single sojourn on a 1-y grid ----------------------------------------------------------------------------------------------------------
sojourn <- function(o, scH, scJ = scH, rows = seq_len(ND), tmax = 100) { ex <- EX[[o]]; if (!length(ex$lb)) return(NULL)
  tg <- 0:tmax; um <- tg[-1] - 0.5
  mk <- function(sc0) lapply(ex$lb, function(lb) { X <- PAR[[lb]][rows, , drop = FALSE]; e <- eta_of(X, sc0)
    inc <- outer(e$a, tg[-1], function(aa, x) x^aa) - outer(e$a, tg[-length(tg)], function(aa, x) x^aa)
    exp(e$eta + outer(e$btsd, f_tsd(um))) * inc })
  dH <- mk(scH)
  if (!identical(scH, scJ)) { dJ <- mk(scJ); tH <- Reduce(`+`, dH); tJ <- Reduce(`+`, dJ); dH <- lapply(dJ, function(m) tH * ifelse(tJ > 0, m / tJ, 0)) }
  list(dH = dH, to = ex$to, t = tg[-1]) }
first_move <- function(sj, hor = 50) { tot <- Reduce(`+`, sj$dH); Hc <- t(apply(tot, 1, cumsum)); Sm <- exp(-(cbind(0, Hc[, -ncol(Hc), drop = FALSE]) + Hc) / 2)
  cif <- sapply(sj$dH, function(m) rowSums((m * Sm)[, seq_len(hor), drop = FALSE])); if (!is.matrix(cif)) cif <- matrix(cif, nrow = nrow(tot))
  pl <- 1 - exp(-Hc[, hor]); nxt <- cif / pmax(pl, 1e-12); list(pleave = pl, shift = drop(nxt %*% (STI[sj$to] - STI[S[1]] * 0)), nxt = nxt) }
f1 <- file.path(CE, "E1_E2_first_move.rds")
if (!file.exists(f1)) {
  E1 <- list(); FM <- list()
  for (o in S) { if (!length(EX[[o]]$lb)) next
    for (sn in names(SCN)) { sj <- sojourn(o, SCN[[sn]]); tot <- Reduce(`+`, sj$dH)
      for (k in seq_along(sj$to)) { sh <- ifelse(tot > 0, sj$dH[[k]] / tot, 0)
        E1[[length(E1) + 1]] <- data.table(origin = o, scen = sn, dest = S[sj$to[k]], t = sj$t, est = sh[1, ], lo = apply(sh[-1, ], 2, quantile, 0.025), hi = apply(sh[-1, ], 2, quantile, 0.975)) }
      dS <- Reduce(`+`, lapply(seq_along(sj$to), function(k) ifelse(tot > 0, sj$dH[[k]] / tot, 0) * (STI[sj$to[k]] - STI[o])))
      tolk <- TOL[sj$to] - TOL[o]; ok <- !is.na(tolk)
      dT <- if (any(ok)) Reduce(`+`, lapply(which(ok), function(k) sj$dH[[k]] * tolk[k])) / Reduce(`+`, sj$dH[which(ok)]) else NULL
      E1[[length(E1) + 1]] <- data.table(origin = o, scen = sn, dest = "<expected dSTI>", t = sj$t, est = dS[1, ], lo = apply(dS[-1, ], 2, quantile, 0.025), hi = apply(dS[-1, ], 2, quantile, 0.975))
      if (!is.null(dT)) E1[[length(E1) + 1]] <- data.table(origin = o, scen = sn, dest = "<expected dTOL, scored moves>", t = sj$t, est = dT[1, ], lo = apply(dT[-1, ], 2, quantile, 0.025, na.rm = TRUE), hi = apply(dT[-1, ], 2, quantile, 0.975, na.rm = TRUE))
      fm <- first_move(sj); FM[[paste(o, sn, sn)]] <- fm$shift - STI[o] * 0 }
    for (pn in names(PERT)) { P <- SCN[[PERT[[pn]][1]]]; B <- SCN[[PERT[[pn]][2]]]
      FM[[paste(o, PERT[[pn]][1], PERT[[pn]][2])]] <- first_move(sojourn(o, P, B))$shift
      FM[[paste(o, PERT[[pn]][2], PERT[[pn]][1])]] <- first_move(sojourn(o, B, P))$shift } }
  E1 <- rbindlist(E1); E1[, origin_sti := STI[origin]]
  saveRDS(list(E1 = E1, FM = FM), f1); cat("E1/E2 first move done\n") }

# ---- multi-step engine (mass over start x state x age; as 17 + mixing, counts, returns) ------------------------------------------------------------------------
propagate_x <- function(scH, scJ = scH, rows = seq_len(ND), tmax = 100, returns = FALSE, mods = NULL, light = FALSE, w0 = W0) {
  dt <- 1; nA <- tmax + 1L; nT <- tmax; D <- length(rows); ages <- 0:(nA - 1)
  prep <- function(sc0) lapply(S, function(o) lapply(EX[[o]]$lb, function(lb) { X <- PAR[[lb]][rows, , drop = FALSE]
    e <- eta_of(X, sc0, if (!is.null(mods) && mods$lb == lb) mods else NULL)
    list(eta = e$eta, btsd = e$btsd, inc = outer(e$a, ages + dt, function(aa, x) x^aa) - outer(e$a, ages, function(aa, x) x^aa)) }))
  TH <- prep(scH); TJ <- if (identical(scH, scJ)) NULL else prep(scJ)
  M <- array(0, c(D, nS, nS, nA)); for (i in seq_len(nS)) M[, i, i, 1] <- 1
  yrs <- array(0, c(D, nS, nS, length(HOR))); TI <- array(NA_real_, c(D, nS, nT)); MIX <- array(NA_real_, c(D, nS, nT))
  CNT <- if (light) NULL else array(0, c(D, nS, nS, nS, length(HOR)))         # [draw, start, from, to, horizon]
  RET <- if (returns) array(0, c(D, nS, nT)) else NULL
  for (st in seq_len(nT)) { ft <- f_tsd(st - 0.5); New <- array(0, c(D, nS, nS, nA))
    for (s in seq_len(nS)) { ex <- EX[[s]]; if (!length(ex$lb)) { New[, , s, 2:nA] <- New[, , s, 2:nA] + M[, , s, 1:(nA - 1)]; next }
      dH <- lapply(TH[[s]], function(tr) exp(tr$eta + tr$btsd * ft) * tr$inc)
      if (!is.null(TJ)) { dJ <- lapply(TJ[[s]], function(tr) exp(tr$eta + tr$btsd * ft) * tr$inc); tH <- Reduce(`+`, dH); tJ <- Reduce(`+`, dJ)
        dH <- lapply(dJ, function(m) tH * ifelse(tJ > 0, m / tJ, 0)) }
      tot <- Reduce(`+`, dH); pl <- 1 - exp(-tot)
      Ms <- M[, , s, , drop = FALSE]; dim(Ms) <- c(D, nS, nA)
      stayM <- sweep(Ms, c(1, 3), 1 - pl, "*"); New[, , s, 2:nA] <- New[, , s, 2:nA] + stayM[, , 1:(nA - 1)]
      for (k in seq_along(dH)) { pk <- ifelse(tot > 0, dH[[k]] / tot, 0) * pl; flow <- apply(sweep(Ms, c(1, 3), pk, "*"), c(1, 2), sum); j <- ex$to[k]
        if (!light) for (h in seq_along(HOR)) if (st <= HOR[h]) CNT[, , s, j, h] <- CNT[, , s, j, h] + flow
        if (returns) { RET[, j, st] <- RET[, j, st] + flow[, j]; flow[, j] <- 0 }
        New[, , j, 1] <- New[, , j, 1] + flow } }
    M <- New; occ <- apply(M, c(1, 2, 3), sum)                                   # [draw, start, state]
    for (h in seq_along(HOR)) if (st <= HOR[h]) yrs[, , , h] <- as.vector(yrs[, , , h]) + as.vector(occ) * dt   # as.vector: D = 1 drops a dim
    TI[, , st] <- apply(occ, c(1, 2), function(p) sum(p * STI)); MIX[, , st] <- apply(occ, c(1, 3), function(p) sum(p * w0)) }
  list(yrs = yrs, TI = TI, MIX = MIX, CNT = CNT, RET = RET) }

# ---- run list (resumable) ----------------------------------------------------------------------------------------------------------------------------------------
RUNS <- c(lapply(names(SCN), function(sn) list(id = paste0("base_", sn), H = sn, J = sn, returns = FALSE)),
          lapply(names(SCN), function(sn) list(id = paste0("ret_", sn), H = sn, J = sn, returns = TRUE)),
          unlist(lapply(names(PERT), function(pn) list(list(id = paste0("mixHP_", pn), H = PERT[[pn]][1], J = PERT[[pn]][2], returns = FALSE),
                                                       list(id = paste0("mixJP_", pn), H = PERT[[pn]][2], J = PERT[[pn]][1], returns = FALSE))), recursive = FALSE))
for (R in RUNS) { f <- file.path(CE, paste0(R$id, ".rds")); if (file.exists(f)) next; if (over_time()) { cat("time budget reached; re-run\n"); quit(save = "no") }
  t0 <- Sys.time(); r <- propagate_x(SCN[[R$H]], SCN[[R$J]], returns = R$returns, light = grepl("^mix", R$id))
  if (R$returns) r <- list(RET = r$RET)
  saveRDS(r, f); cat(sprintf("  %s (%.0f s)\n", R$id, as.numeric(Sys.time() - t0, units = "secs"))) }
# E3 elasticity: point estimates, 50 y, one file per transition (4 runs)
fe0 <- file.path(CE, "elas_base.rds"); if (!file.exists(fe0)) saveRDS(propagate_x(SCN$none_z0, rows = 1L, tmax = 50, light = TRUE)[c("MIX", "TI")], fe0)
for (lb in modelled) { f <- file.path(CE, paste0("elas_", gsub("->", "__", lb), ".rds")); if (file.exists(f)) next; if (over_time()) { cat("time budget reached; re-run\n"); quit(save = "no") }
  out <- list(); for (ty in c("b0", "alpha")) for (fa in c(1.1, 0.9)) { r <- propagate_x(SCN$none_z0, rows = 1L, tmax = 50, light = TRUE, mods = list(lb = lb, type = ty, f = fa))
    out[[paste(ty, fa)]] <- r$MIX[1, , 50] }
  saveRDS(out, f) }
# E6 convergence: point estimates, 1000 y
for (sn in c("none_zm1", "none_z0", "none_zp1")) { f <- file.path(CE, paste0("conv_", sn, ".rds")); if (file.exists(f)) next; if (over_time()) { cat("time budget reached; re-run\n"); quit(save = "no") }
  t0 <- Sys.time(); r <- propagate_x(SCN[[sn]], rows = 1L, tmax = 1000, light = TRUE); saveRDS(r$MIX[1, , ], f); cat(sprintf("  conv %s (%.0f s)\n", sn, as.numeric(Sys.time() - t0, units = "secs"))) }
cat("all runs done\n")
