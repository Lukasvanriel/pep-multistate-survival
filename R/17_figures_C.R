# ==============================================================================
# 17_figures_C.R   (no INLA; resumable: re-run until it prints "all scenarios done"; run from the project root)
# Exploratory set C: several moves. Probabilities implied by the fits under fixed scenarios (not a projection of real plots). NOT manuscript figures.
# Model as in 16_figures_B.R (variant 0, years; competing exits; transitions dropped by the event rule have hazard 0).
# Method: probability mass over (start state, current state, years in current state), propagated in DT-year steps from t = 0 (all mass in the start
#   state, time in state 0). Per step, exit k of state s at time in state a: dH_k = exp(eta_k(t)) ((a + DT)^alpha - a^alpha), eta at the step midpoint;
#   P(leave) = 1 - exp(-sum dH), P(-> k) = dH_k / sum dH * P(leave); at most one move per step (checked against DT = 0.25).
# Scenarios (user 2026-11-09): CMI, soil z = 0; GDD z in {-1, 0, 1} and a grid -2..2 (undisturbed); undisturbed or ONE disturbance at t = 0 that
#   CARRIES OVER after state changes, TSD = t updating (flags persist, as in the data build); exits without that disturbance term keep their
#   undisturbed rate (assumption; per-origin shares and flags in figures/B/B_flags_share_assumed.csv). Horizon 100 y (beyond 55 y = extrapolation).
# Uncertainty: N_DRAW joint parameter draws per transition (as 16), CI = 2.5/97.5 pct; point = pooled estimates.
# Outputs (figures/C/): C1 occupancy P(j, t | start i); C2 thermal index of the expected composition; C3 first passage P(reached j by t | i);
#   C4 pathways (simulation of individual stands, point estimates); C5 long-run composition (point, 1000 y, undisturbed) vs observed visits.
# ==============================================================================
source(here::here("R", "13_fig_helpers.R"))
source(here::here("R", "covariate_spec.R")); f_tsd <- tsd_fun(TSD_FORM)
OUTC <- file.path(FIG_DIR, "C"); CACHE <- file.path(OUTC, "cache"); dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
N_DRAW <- 100L; DT <- 1; TMAX <- 100; set.seed(20261109); TIME_BUDGET <- 150   # seconds per call (device limit); resumes from cache
DIST <- c("fire", "pest", "harv_light", "harv_heavy")
TERMS <- list(fire = c("fire_f", "gdd_x_fire", "tsd_fire"), pest = c("pest_f", "gdd_x_pest", "tsd_pest"),
              harv_light = c("harv_light", "gdd_x_harv_light", "tsd_harv_light"), harv_heavy = c("harv_heavy", "gdd_x_harv_heavy", "tsd_harv_heavy"))
S <- ORDER; nS <- length(S); t_start <- Sys.time()

# ---- parameter draws (as 16) ------------------------------------------------------------------------------------------------------------------------------
pf <- file.path(CACHE, "PAR.rds")
if (file.exists(pf)) PAR <- readRDS(pf) else {
  cvl <- readRDS(file.path(RESULTS_DIR, sprintf("pooled_%s_final_cov.rds", RUN_TAG)))
  fv  <- fread(file.path(RESULTS_DIR, sprintf("fit_validity_%s.csv", RUN_TAG)))[used == TRUE, .(label, m)]
  mir <- rbindlist(lapply(file.path(FIT_DIR, sprintf(c("mi_raw_%s.csv", "mi_raw_redraw_%s.csv"), RUN_TAG)), function(f) if (file.exists(f)) fread(f)[, .(label, m, parameter, est)]))
  mir <- unique(mir[fv, on = c("label", "m"), nomatch = NULL], by = c("label", "m", "parameter"))
  PAR <- list()
  for (lb in modelled) { nm <- names(cvl[[lb]]$mean); pe <- setNames(pooled[label == lb, pooled_est], pooled[label == lb, parameter])[nm]
    Q <- as.matrix(dcast(mir[label == lb], m ~ parameter, value.var = "est")[, ..nm]); Tm <- cvl[[lb]]$W[nm, nm] + (1 + 1 / nrow(Q)) * stats::cov(Q); Tm <- (Tm + t(Tm)) / 2
    ev <- eigen(Tm, symmetric = TRUE); L <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)))
    X <- rbind(pe, sweep(matrix(rnorm(N_DRAW * length(nm)), N_DRAW) %*% t(L), 2, pe, "+")); colnames(X) <- nm; stopifnot(all(X[, "alpha"] > 0))
    PAR[[lb]] <- cbind(X, b0y = X[, "(Intercept)"] - X[, "alpha"] * log(10)) }
  saveRDS(PAR, pf) }
ND <- nrow(PAR[[1]])                                                       # draws + point (row 1)
gcol <- function(X, v) if (v %in% colnames(X)) X[, v] else rep(0, nrow(X))
EX <- lapply(S, function(o) { lb <- trow[from == o][order(match(to, ORDER)), label]; list(lb = lb, to = match(sub(".*->", "", lb), S)) }); names(EX) <- S

# ---- propagation engine --------------------------------------------------------------------------------------------------------------------------------------
# returns occ[draw, start, state, step] (P in state at end of step) and, if absorbing = j, the mass absorbed in j (first passage)
propagate <- function(zg = 0, dist = "none", dt = DT, tmax = TMAX, absorb = NA_integer_, rows = seq_len(ND), keep = TRUE) {
  nA <- round(tmax / dt) + 1L; nT <- round(tmax / dt); D <- length(rows)
  ages <- (0:(nA - 1)) * dt
  # per transition: base eta (no TSD) and cumulative-hazard increments per age, [draw, age]
  TR <- lapply(S, function(o) lapply(EX[[o]]$lb, function(lb) { X <- PAR[[lb]][rows, , drop = FALSE]; a <- X[, "alpha"]
    eta <- X[, "b0y"] + gcol(X, "z_gdd_3") * zg; btsd <- rep(0, D)
    if (dist != "none" && TERMS[[dist]][1] %in% colnames(X)) { tm <- TERMS[[dist]]; eta <- eta + X[, tm[1]] + gcol(X, tm[2]) * zg; btsd <- gcol(X, tm[3]) }
    inc <- outer(a, ages + dt, function(aa, x) x^aa) - outer(a, ages, function(aa, x) x^aa)
    list(eta = eta, btsd = btsd, inc = inc) }))
  names(TR) <- S
  M <- array(0, c(D, nS, nS, nA)); for (i in seq_len(nS)) M[, i, i, 1] <- 1       # [draw, start, state, age]
  occ <- if (keep) array(NA_real_, c(D, nS, nS, nT)) else NULL; absorbed <- if (!is.na(absorb)) array(0, c(D, nS, nT)) else NULL
  for (st in seq_len(nT)) { tm <- (st - 0.5) * dt; ft <- f_tsd(tm)
    New <- array(0, c(D, nS, nS, nA))
    for (s in seq_len(nS)) { if (!is.na(absorb) && s == absorb) next
      ex <- EX[[s]]; if (!length(ex$lb)) { New[, , s, 2:nA] <- New[, , s, 2:nA] + M[, , s, 1:(nA - 1)]; next }
      dH <- lapply(TR[[s]], function(tr) exp(tr$eta + tr$btsd * ft) * tr$inc)  # [draw, age]
      tot <- Reduce(`+`, dH); pl <- 1 - exp(-tot); stay <- 1 - pl              # [draw, age]
      Ms <- M[, , s, , drop = FALSE]; dim(Ms) <- c(D, nS, nA)
      stayM <- sweep(Ms, c(1, 3), stay, "*"); New[, , s, 2:nA] <- New[, , s, 2:nA] + stayM[, , 1:(nA - 1)]
      for (k in seq_along(dH)) { pk <- ifelse(tot > 0, dH[[k]] / tot, 0) * pl; flow <- apply(sweep(Ms, c(1, 3), pk, "*"), c(1, 2), sum)   # [draw, start]
        j <- ex$to[k]; if (!is.na(absorb) && j == absorb) absorbed[, , st] <- absorbed[, , st] + flow else New[, , j, 1] <- New[, , j, 1] + flow } }
    if (!is.na(absorb)) New[, , absorb, ] <- 0
    M <- New; if (keep) occ[, , , st] <- apply(M, c(1, 2, 3), sum) }
  list(occ = occ, absorbed = absorbed, t = (1:nT) * dt)
}
qs3 <- function(a) list(est = a[1, , , , drop = FALSE], lo = apply(a[-1, , , , drop = FALSE], 2:4, quantile, 0.025), hi = apply(a[-1, , , , drop = FALSE], 2:4, quantile, 0.975))
to_dt <- function(occ, t, keys) { e <- occ[1, , , ]; lo <- apply(occ[-1, , , , drop = FALSE], 2:4, quantile, 0.025); hi <- apply(occ[-1, , , , drop = FALSE], 2:4, quantile, 0.975)
  g <- CJ(start = seq_len(nS), state = seq_len(nS), step = seq_along(t), sorted = FALSE)
  idx <- cbind(g$start, g$state, g$step); data.table(start = S[g$start], state = S[g$state], t = t[g$step], est = e[idx], lo = lo[idx], hi = hi[idx])[, names(keys) := keys][] }

# ---- scenario list (resumable) ----------------------------------------------------------------------------------------------------------------------------
SC <- rbind(data.table(id = sprintf("none_z%+g", c(-1, 0, 1)), zg = c(-1, 0, 1), dist = "none"),
            data.table(id = sprintf("%s_z0", DIST), zg = 0, dist = DIST),
            data.table(id = sprintf("grid_z%+.2f", setdiff(seq(-2, 2, 0.5), c(-1, 0, 1))), zg = setdiff(seq(-2, 2, 0.5), c(-1, 0, 1)), dist = "none"))
FP <- data.table(id = sprintf("fp_%s_z%+g", rep(S, each = 3), rep(c(-1, 0, 1), nS)), target = rep(seq_len(nS), each = 3), zg = rep(c(-1, 0, 1), nS))
for (i in seq_len(nrow(SC))) { f <- file.path(CACHE, paste0("occ_", SC$id[i], ".rds")); if (file.exists(f)) next
  if (as.numeric(Sys.time() - t_start, units = "secs") > TIME_BUDGET) { cat("time budget reached; re-run to continue\n"); quit(save = "no") }
  r <- propagate(zg = SC$zg[i], dist = SC$dist[i]); saveRDS(to_dt(r$occ, r$t, list(scen = SC$id[i], zg = SC$zg[i], dist = SC$dist[i])), f); cat("  occupancy", SC$id[i], "\n") }
for (i in seq_len(nrow(FP))) { f <- file.path(CACHE, paste0(FP$id[i], ".rds")); if (file.exists(f)) next
  if (as.numeric(Sys.time() - t_start, units = "secs") > TIME_BUDGET) { cat("time budget reached; re-run to continue\n"); quit(save = "no") }
  r <- propagate(zg = FP$zg[i], absorb = FP$target[i], keep = FALSE); cum <- aperm(apply(r$absorbed, c(1, 2), cumsum), c(2, 3, 1))   # [draw, start, step]
  e <- cum[1, , ]; lo <- apply(cum[-1, , , drop = FALSE], 2:3, quantile, 0.025); hi <- apply(cum[-1, , , drop = FALSE], 2:3, quantile, 0.975)
  g <- CJ(start = seq_len(nS), step = seq_along(r$t)); idx <- cbind(g$start, g$step)
  saveRDS(data.table(start = S[g$start], target = S[FP$target[i]], t = r$t[g$step], est = e[idx], lo = lo[idx], hi = hi[idx], zg = FP$zg[i]), f); cat("  first passage", FP$id[i], "\n") }
cat("all scenarios done\n")
