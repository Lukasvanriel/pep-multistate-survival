# ==============================================================================
# 18_figures_D1_recovery.R   (no INLA; resumable: re-run until "all D1 scenarios done"; then draws the figures; run from the project root)
# Exploratory set D1: what happens after a disturbance -- resistance, recovery, reorganization. NOT manuscript figures.
# Engine as 17_figures_C.R (variant 0, years, 1-y steps, competing exits, at most one move per step, disturbance at t = 0 carried over with TSD = t,
#   exits without that disturbance term at their undisturbed rate), extended with a history flag per stand:
#   flag 1 = never left the start state; 2 = left, never passed through regen; 3 = has passed through regen.
# Outcome at t for a stand that started in state X (X != regen):
#   resistance = in X, never left | recovery via regen = in X, passed through regen | return without regen = in X, left but no regen |
#   in regen | reorganization after regen = in another state, passed through regen (by group) | drift = in another state, no regen (by group).
# Groups (user 2026-11-09): boreal = black spruce, white spruce, fir, jack pine; temperate = sugar maple, red maple, yellow birch, cedar,
#   other temperate; pioneer = white birch, trembling aspen, other pioneer; regen.
# Scenarios: undisturbed / fire / pest / light / heavy harvest at t = 0; GDD z -1/0/+1 (CMI, soil z = 0); time in state at t = 0: 0 y (main,
#   consistent with C) and 40 y (established stand); both for every scenario (user 2026-11-09: show both). N_DRAW_D draws (first N_DRAW_D of 17's draws) + point.
# ==============================================================================
src <- readLines(here::here("R", "17_figures_C.R")); eval(parse(text = src[1:(grep("^# ---- propagation engine", src) - 1)]))
OUTD <- file.path(FIG_DIR, "D"); CD <- file.path(OUTD, "cache"); dir.create(CD, recursive = TRUE, showWarnings = FALSE)
N_DRAW_D <- 50L; ROWS <- 1:(N_DRAW_D + 1L); TIME_BUDGET <- 150; t0 <- Sys.time()
source(here::here("R", "00_groups.R"))   # GRP: thermal groups
iR <- match("regen", S); STARTS <- setdiff(seq_len(nS), iR); nI <- length(STARTS)

propagate_hist <- function(zg = 0, dist = "none", a0 = 0, tmax = TMAX) {
  D <- length(ROWS); nA <- tmax + a0 + 1L; ages <- 0:(nA - 1)
  TR <- lapply(S, function(o) lapply(EX[[o]]$lb, function(lb) { X <- PAR[[lb]][ROWS, , drop = FALSE]; a <- X[, "alpha"]
    eta <- X[, "b0y"] + gcol(X, "z_gdd_3") * zg; btsd <- rep(0, D)
    if (dist != "none" && TERMS[[dist]][1] %in% colnames(X)) { tm <- TERMS[[dist]]; eta <- eta + X[, tm[1]] + gcol(X, tm[2]) * zg; btsd <- gcol(X, tm[3]) }
    list(eta = eta, btsd = btsd, inc = outer(a, ages + 1, function(aa, x) x^aa) - outer(a, ages, function(aa, x) x^aa)) })); names(TR) <- S
  M <- array(0, c(D, nI, nS, 3, nA)); for (i in seq_len(nI)) M[, i, STARTS[i], 1, a0 + 1] <- 1
  out <- array(NA_real_, c(D, nI, nS, 3, tmax))                              # occupancy by state x flag (summed over age)
  for (st in seq_len(tmax)) { ft <- f_tsd(st - 0.5); New <- array(0, dim(M))
    for (s in seq_len(nS)) { ex <- EX[[s]]; Ms <- M[, , s, , , drop = FALSE]; dim(Ms) <- c(D, nI, 3, nA)
      if (!length(ex$lb)) { New[, , s, , 2:nA] <- New[, , s, , 2:nA] + Ms[, , , 1:(nA - 1)]; next }
      dH <- lapply(TR[[s]], function(tr) exp(tr$eta + tr$btsd * ft) * tr$inc); tot <- Reduce(`+`, dH); pl <- 1 - exp(-tot)
      stayM <- sweep(Ms, c(1, 4), 1 - pl, "*"); New[, , s, , 2:nA] <- New[, , s, , 2:nA] + stayM[, , , 1:(nA - 1)]
      for (k in seq_along(dH)) { pk <- ifelse(tot > 0, dH[[k]] / tot, 0) * pl; fl <- apply(sweep(Ms, c(1, 4), pk, "*"), 1:3, sum)   # [draw, start, flag]
        j <- ex$to[k]
        if (j == iR) New[, , j, 3, 1] <- New[, , j, 3, 1] + apply(fl, 1:2, sum)
        else { New[, , j, 2, 1] <- New[, , j, 2, 1] + fl[, , 1] + fl[, , 2]; New[, , j, 3, 1] <- New[, , j, 3, 1] + fl[, , 3] } } }
    M <- New; out[, , , , st] <- apply(M, 1:4, sum) }
  out
}
summarise_hist <- function(out, keys) { D <- dim(out)[1]; tmax <- dim(out)[5]; res <- list()
  for (i in seq_len(nI)) { X <- STARTS[i]; others <- setdiff(seq_len(nS), c(X, iR))
    cat_fun <- list(resistance = function(o) o[, i, X, 1, ], recovery_via_regen = function(o) o[, i, X, 3, ], return_without_regen = function(o) o[, i, X, 2, ],
                    in_regen = function(o) apply(o[, i, iR, , , drop = FALSE], c(1, 5), sum))
    for (g in c("boreal", "temperate", "pioneer")) { og <- others[GRP[S[others]] == g]; if (!length(og)) next
      cat_fun[[paste0("reorganized_", g)]] <- local({ og <- og; function(o) apply(o[, i, og, 3, , drop = FALSE], c(1, 5), sum) })
      cat_fun[[paste0("drift_", g)]] <- local({ og <- og; function(o) apply(o[, i, og, 2, , drop = FALSE], c(1, 5), sum) }) }
    for (cn in names(cat_fun)) { v <- matrix(cat_fun[[cn]](out), nrow = D)
      res[[length(res) + 1]] <- data.table(start = S[X], outcome = cn, t = seq_len(tmax), est = v[1, ], lo = apply(v[-1, , drop = FALSE], 2, quantile, 0.025),
                                           hi = apply(v[-1, , drop = FALSE], 2, quantile, 0.975))[, names(keys) := keys] } }
  rbindlist(res) }
SCD <- CJ(dist = c("none", DIST), zg = c(-1, 0, 1), a0 = c(0, 40))   # both times in state at the disturbance for all GDD (user 2026-11-09)
SCD[, id := sprintf("d1_%s_z%+d_a%d", dist, zg, a0)]
for (i in seq_len(nrow(SCD))) { f <- file.path(CD, paste0(SCD$id[i], ".rds")); if (file.exists(f)) next
  if (as.numeric(Sys.time() - t0, units = "secs") > TIME_BUDGET) { cat("time budget reached; re-run to continue\n"); quit(save = "no") }
  o <- propagate_hist(zg = SCD$zg[i], dist = SCD$dist[i], a0 = SCD$a0[i])
  saveRDS(summarise_hist(o, list(dist = SCD$dist[i], zg = SCD$zg[i], a0 = SCD$a0[i])), f); cat("  D1", SCD$id[i], sprintf("(%.0f s)\n", as.numeric(Sys.time() - t0, units = "secs"))) }
cat("all D1 scenarios done\n")
