# ==============================================================================
# R/tests/wp4_C1_check.R   (WP4: C1 piecewise P1, C5 weights; no INLA)
#   Rscript R/tests/wp4_C1_check.R          synthetic tables (runs anywhere): the real prepare_fit_inputs() on small made-up tables
#   Rscript R/tests/wp4_C1_check.R real     ALSO on the real tables (TABLE_DIR, after 04 has written the fv_* columns): base rows, onsets, pieces
# Checks (a failed check prints what it has and stops; bounds are never changed):
#   1  prepare_fit_inputs: base rows (X1) = first-visit values, X0 = X1 without disturbance, onset = W_j - first-visit TSD
#   2  pieces: grid from the first onset, EXTRA breakpoint at the exact onset of a later type (user decision 2026-11-04), onset <= 0 from time 0,
#      TSD at each piece midpoint, flags of the last piece = first-visit flags, contiguous cover of [0, end of the weighting bracket]
#   3  identity: plots without an onset before W_j give the same candidates as the WP3 p1_candidates fed the same first-visit lambda (same seed)
#   4  C1 unit test: 100,000 piecewise exit draws vs the exact piecewise survival, rebuilt independently of the production code (KS p > 0.01; exit type
#      chi-square p > 0.01); 4 cases, a separate seed each; every p-value is printed with its case
#   5  weights: p1_weights_pw (log space) vs an independent numerical integration of the piecewise hazard
#   6  p1_entries runs end to end (mean and draw), reproducible with the seed; fallback table by first state
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("R", "config.R")); source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
REAL <- "real" %in% commandArgs(TRUE)
fail <- function(msg, x = NULL) { cat("\n!! CHECK FAILED:", msg, "\n"); if (!is.null(x)) print(x, digits = 6)
  stop("check failed; bounds were NOT changed. Report this output to the user before changing anything.", call. = FALSE) }
ok <- function(cond, msg, x = NULL) if (!isTRUE(cond)) fail(msg, x)
set.seed(20261104)

# ---- synthetic tables ----------------------------------------------------------------------------------------------------------------------
st <- c("regen", "a", "b", "c")
tr <- data.table(from = c("regen", "regen", "a", "a", "a", "b", "b", "b", "c", "c"), to = c("a", "b", "b", "c", "regen", "a", "c", "regen", "a", "regen"))
tr[, `:=`(kept = TRUE, fire_term = c(T,F,T,F,T,F,T,F,T,F), pest_term = c(F,F,T,F,F,F,T,F,F,F), harv_term = ifelse(seq_len(.N) %% 2 == 1, "3level", "none"),
          n_harv1 = 5L, n_harv2 = 5L, n_undist = 20L)]
mk_plot <- function(id, W, first_state, fv, nint = 3, change_at = NA) {
  gap <- sample(5:10, nint, TRUE); ts <- W + c(0, cumsum(gap))[seq_len(nint)]; te <- ts + gap
  s_end <- rep(first_state, nint); s_start <- s_end
  if (!is.na(change_at)) { s_end[change_at:nint] <- if (first_state == "a") "b" else "a"; if (change_at < nint) s_start[(change_at + 1):nint] <- s_end[change_at] }
  iv <- data.table(id_pe = id, interval = seq_len(nint), t_start = ts, t_end = te, state_start = s_start, state_end = s_end,
                   gdd_3 = rnorm(nint, 1500, 150), cmi_gs = rnorm(nint, 5, 3), soil_pc1 = rnorm(nint),
                   fire = 0L, harv = 0L, pest = 0L, tsd_fire_y = 0, tsd_pest_y = 0, tsd_harv_y = 0)
  # the flags of the first visit continue in the observed intervals (TSD grows with the calendar)
  if (fv$fire) iv[, `:=`(fire = 1L, tsd_fire_y = fv$tf + (t_end - W))]
  if (fv$harv > 0) iv[, `:=`(harv = fv$harv, tsd_harv_y = fv$th + (t_end - W))]
  if (fv$pest) iv[, `:=`(pest = 1L, tsd_pest_y = fv$tp + (t_end - W))]
  pi <- data.table(id_pe = id, W_j = W, first_state = first_state, first_year = 1980L, n_visits = nint + 1L, n_changes = as.integer(!is.na(change_at)),
                   fv_gdd_3 = rnorm(1, 1500, 150), fv_cmi_gs = rnorm(1, 5, 3), fv_soil_pc1 = rnorm(1),
                   fv_fire = as.integer(fv$fire), fv_harv = as.integer(fv$harv), fv_pest = as.integer(fv$pest),
                   fv_tsd_fire_y = if (fv$fire) fv$tf else 0, fv_tsd_pest_y = if (fv$pest) fv$tp else 0, fv_tsd_harv_y = if (fv$harv > 0) fv$th else 0)
  iv[, is_change := state_start != state_end]; list(iv = iv, pi = pi)
}
F0 <- list(fire = 0, harv = 0, pest = 0)
specs <- list(
  list(id = 1, W = 40, s = "a", fv = list(fire = 1, harv = 1, pest = 0, tf = 12.5, th = 14.6)),           # onsets 27.5 (fire), 25.4 (light harvest)
  list(id = 2, W = 30, s = "b", fv = list(fire = 1, harv = 0, pest = 0, tf = 35, th = 0, tp = 0)),          # onset -5 (before stand origin)
  list(id = 3, W = 45, s = "a", fv = list(fire = 0, harv = 2, pest = 1, th = 3.3, tp = 0.4), change = 2),   # heavy harvest + pest, change observed
  list(id = 4, W = 25, s = "a", fv = list(fire = 1, harv = 0, pest = 0, tf = 0, th = 0)),                   # TSD 0: onset = W_j, never on before the first visit
  list(id = 5, W = 50, s = "c", fv = F0), list(id = 6, W = 33.5, s = "a", fv = F0, change = 3), list(id = 7, W = 60, s = "b", fv = F0),
  list(id = 8, W = 20, s = "regen", fv = F0))
for (i in 9:40) specs[[i]] <- list(id = i, W = runif(1, 15, 90), s = sample(st, 1), fv = if (i %% 3 == 0) list(fire = 0, harv = sample(1:2, 1), pest = 0, th = runif(1, 0, 40)) else F0)
pl <- lapply(specs, function(s) mk_plot(s$id, s$W, s$s, s$fv, nint = sample(2:4, 1), change_at = if (is.null(s$change)) NA else s$change))
d <- tempfile("wp4"); dir.create(d)
saveRDS(rbindlist(lapply(pl, `[[`, "iv")), file.path(d, "visit_intervals.rds")); saveRDS(rbindlist(lapply(pl, `[[`, "pi")), file.path(d, "plot_info.rds"))
saveRDS(tr, file.path(d, "transition_set.rds"))
inp <- prepare_fit_inputs(d)
K <- nrow(inp$td); P <- length(inp$terms)
mk_theta <- function(ok_all = TRUE) { th <- list(alpha = runif(K, 1.1, 2.1), intercept = log(1 / runif(K, 25, 70)), beta = matrix(rnorm(K * P, 0, 0.25), K, P, dimnames = list(inp$td$label, inp$terms)),
                                               ok = rep(TRUE, K)); th$beta[!inp$mask] <- 0; th }
theta <- mk_theta()
pid <- inp$p1$id_pe; jj <- function(id) match(id, pid)
f_tsd <- tsd_fun()

# ---- 1. base rows and onsets ------------------------------------------------------------------------------------------------------------------
pi_all <- readRDS(file.path(d, "plot_info.rds"))
for (j in seq_along(pid)) {
  p <- pi_all[id_pe == pid[j]]; sc <- inp$scaling
  z <- function(v, x) (x - sc[var == v, mean]) / sc[var == v, sd]
  ok(abs(inp$X1[j, "z_gdd_3"] - z("gdd_3", p$fv_gdd_3)) < 1e-12 && abs(inp$X1[j, "z_cmi_gs"] - z("cmi_gs", p$fv_cmi_gs)) < 1e-12 &&
       abs(inp$X1[j, "z_soil_pc1"] - z("soil_pc1", p$fv_soil_pc1)) < 1e-12, sprintf("plot %s: X1 climate/soil != scaled first-visit values", pid[j]))
  ok(inp$X1[j, "fire_f"] == p$fv_fire && inp$X1[j, "pest_f"] == p$fv_pest && inp$X1[j, "harv_light"] == (p$fv_harv == 1) && inp$X1[j, "harv_heavy"] == (p$fv_harv == 2),
     sprintf("plot %s: X1 flags != first-visit flags", pid[j]))
  ok(abs(inp$X1[j, "tsd_fire"] - f_tsd(p$fv_tsd_fire_y) * p$fv_fire) < 1e-12 && abs(inp$X1[j, "tsd_pest"] - f_tsd(p$fv_tsd_pest_y) * p$fv_pest) < 1e-12,
     sprintf("plot %s: X1 TSD columns", pid[j]))
  ok(all(inp$X0[j, c("fire_f", "pest_f", "harv_light", "harv_heavy", "gdd_x_fire", "gdd_x_pest", "gdd_x_harv_light", "gdd_x_harv_heavy", "tsd_fire", "tsd_pest", "tsd_harv_light", "tsd_harv_heavy")] == 0) &&
       all(inp$X0[j, c("z_gdd_3", "z_cmi_gs", "z_soil_pc1")] == inp$X1[j, c("z_gdd_3", "z_cmi_gs", "z_soil_pc1")]), sprintf("plot %s: X0", pid[j]))
  ok(identical(is.na(inp$on$fire[j]), p$fv_fire == 0L) && (p$fv_fire == 0L || abs(inp$on$fire[j] - (p$W_j - p$fv_tsd_fire_y)) < 1e-12), sprintf("plot %s: fire onset", pid[j]))
  ok(identical(is.na(inp$on$harv[j]), p$fv_harv == 0L) && (p$fv_harv == 0L || abs(inp$on$harv[j] - (p$W_j - p$fv_tsd_harv_y)) < 1e-12), sprintf("plot %s: harvest onset", pid[j]))
  ok(identical(is.na(inp$on$pest[j]), p$fv_pest == 0L) && (p$fv_pest == 0L || abs(inp$on$pest[j] - (p$W_j - p$fv_tsd_pest_y)) < 1e-12), sprintf("plot %s: pest onset", pid[j]))
}
cat(sprintf("1. base rows / onsets: %d synthetic plots, all equal the first-visit values; onset = W_j - TSD\n", length(pid)))

# ---- 2. pieces ---------------------------------------------------------------------------------------------------------------------------------
L0 <- lambda_matrix(theta, inp$X0); Liv <- lambda_matrix(theta, inp$Xiv)
pcs <- lapply(seq_along(pid), function(j) p1_pieces(j, theta, inp, L0, Liv, f_tsd))
for (j in seq_along(pid)) {          # contiguous cover of [0, obs_upper], pre-pieces end at W_j, post-pieces = observed intervals
  pc <- pcs[[j]]; W <- inp$p1$W_j[j]
  ok(abs(pc$p0[1] - 0) < 1e-12 && all(abs(pc$p0[-1] - pc$p1[-length(pc$p1)]) < 1e-9) && all(pc$p1 > pc$p0), sprintf("plot %s: pieces not contiguous", pid[j]))
  ok(abs(pc$p1[pc$npre] - W) < 1e-9 && abs(tail(pc$p1, 1) - inp$p1$obs_upper[j]) < 1e-9, sprintf("plot %s: pieces do not end at W_j / obs_upper", pid[j]))
}
j <- jj(1); pc <- pcs[[j]]; W <- 40
exp_brk <- sort(unique(round(c(0, seq(25.4, 40, 1), 27.5, 40), 9)))
ok(isTRUE(all.equal(c(pc$p0[1:pc$npre], W), exp_brk)), "plot 1: breakpoints", list(got = c(pc$p0[1:pc$npre], W), want = exp_brk))
mid <- (pc$p0 + pc$p1)[1:pc$npre] / 2
# TSD at piece midpoints: re-derive from the rules and compare through lambda (lambda_matrix of an independently built row)
build_row <- function(j, m, on_f, on_h, lvl) {
  x <- inp$X0[j, ]; zg <- x[["z_gdd_3"]]
  if (!is.na(on_f) && m >= on_f) { x[["fire_f"]] <- 1; x[["gdd_x_fire"]] <- zg; x[["tsd_fire"]] <- f_tsd(m - on_f) }
  if (!is.na(on_h) && m >= on_h) { x[[lvl]] <- 1; x[[paste0("gdd_x_", lvl)]] <- zg; x[[paste0("tsd_", lvl)]] <- f_tsd(m - on_h) }
  x }
for (q in seq_len(pc$npre)) {
  x <- build_row(j, mid[q], 27.5, 25.4, "harv_light")
  ok(max(abs(lambda_matrix(theta, matrix(x, 1, dimnames = list(NULL, inp$terms))) - pc$L[q, ])) < 1e-12, sprintf("plot 1 piece %d: lambda != independent build", q))
}
nf <- sum(mid < 27.5 & mid >= 25.4); ok(nf >= 1, "plot 1: no piece with harvest on and fire off")        # the split at 27.5 gives a harvest-only piece
ok(any(abs(pc$p0 - 27.5) < 1e-9) && any(abs(pc$p1 - 27.5) < 1e-9), "plot 1: no breakpoint at the exact fire onset 27.5")
ok(abs(mid[pc$npre] - 39.7) < 1e-9, "plot 1: last piece [39.4, 40] must have midpoint 39.7 (cut at W_j, own midpoint)")
# last piece of every flagged plot: flags equal the first-visit flags
for (jx in seq_along(pid)) { pc <- pcs[[jx]]; if (pc$npre == 0) next
  fl <- inp$X1[jx, c("fire_f", "pest_f", "harv_light", "harv_heavy")]; if (all(fl == 0)) next
  on <- c(inp$on$fire[jx], inp$on$pest[jx], inp$on$harv[jx]); if (all(on[!is.na(on)] >= inp$p1$W_j[jx])) next        # onset at W_j: never on before the first visit
  ok(!identical(pc$L[pc$npre, ], L0[jx, ]), sprintf("plot %s: flagged at the first visit but the last pre-piece equals the base row", pid[jx])) }
j <- jj(2); pc <- pcs[[j]]                   # onset -5: on from time 0, grid from 0
ok(isTRUE(all.equal(c(pc$p0[1:pc$npre], 30), 0:30)), "plot 2: grid from 0 expected for an onset <= 0")
x <- build_row(j, 0.5, 30 - 35, NA, "harv_light"); ok(abs(x[["tsd_fire"]] - log1p(5.5)) < 1e-12, "plot 2: TSD at the first midpoint must be 0.5 - (-5) = 5.5")
ok(max(abs(lambda_matrix(theta, matrix(x, 1, dimnames = list(NULL, inp$terms))) - pcs[[j]]$L[1, ])) < 1e-12, "plot 2: first piece lambda")
j <- jj(4); ok(pcs[[j]]$npre == 1L && identical(pcs[[j]]$L[1, ], L0[j, ]), "plot 4 (TSD 0, onset = W_j): one base piece expected")
j <- jj(6); pc <- pcs[[j]]; ok(abs(tail(pc$p1, 1) - inp$p1$obs_upper[j]) < 1e-9 && tail(pc$p1, 1) < inp$iv[id_pe == 6, max(t_end)] + 1e-9, "plot 6: pieces stop at the end of the change interval")
cat("2. pieces: grid, extra breakpoint at the later onset, onset <= 0, TSD at midpoints, cover and ends: OK\n")

# ---- 3. identity with the WP3 candidates for plots without an onset before W_j ------------------------------------------------------------------
# reference lambda = the base row X0 (first-visit climate/soil, no disturbance). For an unflagged plot X0 == X1; a plot flagged with TSD 0 (onset = W_j)
# has its flag off during the whole of P1's history by design (C1.3), so X0, not X1, is its P1 row.
Lx1 <- L0; worst <- 0; n3 <- 0L
for (jx in seq_along(pid)) {
  on <- c(inp$on$fire[jx], inp$on$pest[jx], inp$on$harv[jx]); on <- on[!is.na(on)]
  if (length(on) && any(on < inp$p1$W_j[jx])) next
  set.seed(100 + jx); old <- p1_candidates(inp$start_idx, inp$p1$s_idx[jx], inp$p1$W_j[jx], theta$alpha, Lx1[jx, ], inp$out_k, inp$td$to_idx, 200)
  set.seed(100 + jx); new <- p1_candidates_pw(pcs[[jx]], inp$start_idx, inp$p1$s_idx[jx], inp$p1$W_j[jx], theta$alpha, inp$out_k, inp$td$to_idx, theta$ok, 200)
  ok(length(old) == length(new), sprintf("plot %s: number of candidates differs (%d vs %d)", pid[jx], length(old), length(new)))
  if (length(old)) worst <- max(worst, max(abs(old - new))); n3 <- n3 + 1L
}
ok(worst < 1e-9, sprintf("candidates differ from p1_candidates: max %.3e", worst))
cat(sprintf("3. identity: %d plots without an onset before W_j, same seed: identical candidates (max |difference| %.1e)\n", n3, worst))

# ---- 4. KS test of piecewise exit draws ---------------------------------------------------------------------------------------------------------
# The DRAWS come from the production code (p1_pieces -> p1_next_exit, pieces of pcs). The EXACT side is written here without any production function:
# the piece boundaries and the hazard of every piece are rebuilt from the raw first-visit values of plot_info, the scaling and theta (not from inp$X0,
# not from pc), and the survival is the closed-form piecewise Weibull S(t) = exp(-sum_k sum_pieces lambda_kq^alpha_k [(t_q1 - e)^alpha_k - (t_q0 - e)^alpha_k]),
# t_q0 = max(a_q, e), t_q1 = min(b_q, t) (pieces beyond t or before e contribute 0). Each case has its own seed.
NDRAW <- 100000L
stopifnot("check 4 builds log1p TSD columns by hand" = TSD_FORM == "log1p")
ind_pieces <- function(id) {          # list(a, b, lam = pieces x K) for plot `id`, from raw values
  p <- pi_all[id_pe == id]; W <- p$W_j; sc <- inp$scaling; trm <- inp$terms
  onset <- c(fire = if (p$fv_fire == 1L) W - p$fv_tsd_fire_y, pest = if (p$fv_pest == 1L) W - p$fv_tsd_pest_y, harv = if (p$fv_harv > 0L) W - p$fv_tsd_harv_y)
  lvl <- if (p$fv_harv == 2L) "harv_heavy" else "harv_light"
  g0 <- if (length(onset)) max(0, min(onset)) else Inf
  bk <- sort(unique(c(0, if (g0 < W) seq(g0, W, by = 1), onset[onset > 0 & onset < W]))); bk <- bk[bk < W - 1e-9]; bk <- c(bk, W)
  a <- bk[-length(bk)]; b <- bk[-1]; mid <- (a + b) / 2
  zg <- (p$fv_gdd_3 - sc[var == "gdd_3", mean]) / sc[var == "gdd_3", sd]
  lam <- t(sapply(seq_along(a), function(q) {
    x <- setNames(numeric(length(trm)), trm)
    x[["z_gdd_3"]] <- zg; x[["z_cmi_gs"]] <- (p$fv_cmi_gs - sc[var == "cmi_gs", mean]) / sc[var == "cmi_gs", sd]
    x[["z_soil_pc1"]] <- (p$fv_soil_pc1 - sc[var == "soil_pc1", mean]) / sc[var == "soil_pc1", sd]
    for (ty in names(onset)) if (mid[q] >= onset[[ty]]) {
      nm <- switch(ty, fire = "fire", pest = "pest", harv = sub("harv_", "", lvl))
      fl <- switch(ty, fire = "fire_f", pest = "pest_f", harv = lvl); ts <- switch(ty, fire = "tsd_fire", pest = "tsd_pest", harv = paste0("tsd_", lvl))
      gx <- switch(ty, fire = "gdd_x_fire", pest = "gdd_x_pest", harv = paste0("gdd_x_", lvl))
      x[[fl]] <- 1; x[[gx]] <- zg; x[[ts]] <- log(mid[q] - onset[[ty]] + 1) }
    exp(theta$intercept + drop(theta$beta %*% x)) }))
  list(a = a, b = b, lam = lam, W = W)
}
S_exact <- function(ip, ks, e, t) {   # vectorised in t; total survival of the exits ks of the start state, sojourn entered at e
  H <- numeric(length(t))
  for (q in seq_along(ip$a)) { t0 <- max(ip$a[q], e); t1 <- pmax(pmin(ip$b[q], t), t0)
    for (k in ks) H <- H + ip$lam[q, k]^theta$alpha[k] * ((t1 - e)^theta$alpha[k] - (t0 - e)^theta$alpha[k]) }
  exp(-H)
}
cases <- list(list(id = 1, e = 0), list(id = 1, e = 5.3), list(id = 2, e = 0), list(id = 3, e = 11.7)); KS_LOG <- list()
for (ci in seq_along(cases)) { case <- cases[[ci]]
  jx <- jj(case$id); pc <- pcs[[jx]]; W <- inp$p1$W_j[jx]; e <- case$e; ks <- inp$out_k[[inp$start_idx]]; ip <- ind_pieces(case$id)
  ok(max(abs(ip$lam - pc$L[1:pc$npre, ])) < 1e-10 * max(ip$lam) && isTRUE(all.equal(c(ip$a, ip$W), c(pc$p0[1:pc$npre], W))),
     sprintf("plot %d: independent pieces differ from p1_pieces", case$id))
  p0 <- pc$p0[1:pc$npre]; p1e <- pc$p1[1:pc$npre]; LA <- pc$LA[1:pc$npre, , drop = FALSE]
  set.seed(20261104L + ci)
  dr <- t(vapply(seq_len(NDRAW), function(i) p1_next_exit(p0, p1e, LA, theta$alpha, ks, e, -log(runif(length(ks)))), c(0, 0)))
  hit <- dr[, 1] < W; FW <- 1 - S_exact(ip, ks, e, W)
  u <- (1 - S_exact(ip, ks, e, dr[hit, 1])) / FW                       # conditional CDF given an exit before W_j
  pks <- suppressWarnings(ks.test(u, "punif"))$p.value                  # ties: runif has 32-bit granularity; harmless
  # exit type: P(k, exit before W) = integral of h_k(s) S_tot(s) ds, piece by piece with integrate()
  pk <- sapply(ks, function(k) sum(sapply(seq_along(ip$a), function(q) { t0 <- max(ip$a[q], e); if (ip$b[q] <= t0) return(0)
    integrate(function(s) theta$alpha[k] * ip$lam[q, k]^theta$alpha[k] * (s - e)^(theta$alpha[k] - 1) * S_exact(ip, ks, e, s), t0, ip$b[q], rel.tol = 1e-10)$value })))
  pk <- c(pk, 1 - sum(pk)); obs <- c(table(factor(dr[hit, 2], levels = ks)), sum(!hit)); pchi <- chisq.test(obs, p = pk)$p.value
  KS_LOG[[ci]] <- data.table(case = ci, plot = case$id, entry = e, seed = 20261104L + ci, n_draws = NDRAW, n_exit = sum(hit), ks_p = pks, chisq_p = pchi)
  cat(sprintf("4. case %d: plot %d, entry %.1f, seed %d: %d draws, %d exit before W_j=%.1f | KS p = %.3f | exit type chi-square p = %.3f\n", ci, case$id, e, 20261104L + ci, NDRAW, sum(hit), W, pks, pchi))
  ok(pks > 0.01, sprintf("KS p = %.4f <= 0.01 (case %d: plot %d, entry %.1f)", pks, ci, case$id, e))
  ok(pchi > 0.01, sprintf("exit-type chi-square p = %.4f <= 0.01 (case %d: plot %d, entry %.1f)", pchi, ci, case$id, e), rbind(observed = obs, expected = round(pk * sum(obs), 1)))
}

# ---- 5. weights vs independent numerical integration --------------------------------------------------------------------------------------------
for (id in c(1, 3, 6, 7, 2)) {
  jx <- jj(id); pc <- pcs[[jx]]; sj <- inp$p1$s_idx[jx]; ks <- inp$out_k[[sj]]; W <- inp$p1$W_j[jx]
  lo <- inp$p1$obs_lower[jx]; hi <- inp$p1$obs_upper[jx]
  cand <- sort(runif(12, 0, W * 0.98))
  Hnum <- function(e, a, b) {    # integral of the total hazard of the first-state exits, piece by piece with integrate()
    s <- 0; for (q in seq_along(pc$p0)) { lq <- max(pc$p0[q], a, e); uq <- min(pc$p1[q], b); if (uq <= lq) next
      for (k in ks) s <- s + integrate(function(t) theta$alpha[k] * pc$L[q, k]^theta$alpha[k] * (t - e)^(theta$alpha[k] - 1), lq, uq, rel.tol = 1e-10)$value }; s }
  lw <- sapply(cand, function(e) if (!is.na(lo)) -Hnum(e, W, lo) + log(-expm1(-Hnum(e, lo, hi))) else -Hnum(e, W, hi))
  w_ref <- exp(lw - max(lw)); w_ref <- w_ref / sum(w_ref)
  w_new <- p1_weights_pw(cand, sj, lo, hi, W, pc, theta$alpha, inp$out_k, theta$ok)
  ok(!is.null(w_new) && max(abs(w_new - w_ref)) < 1e-7, sprintf("plot %d: weights differ from the numerical integral (max %.2e)", id, if (is.null(w_new)) NA else max(abs(w_new - w_ref))))
  cat(sprintf("5. plot %d (%s): %d pieces, weights equal the numerical integral (max |diff| %.1e)\n", id, if (is.na(lo)) "no change observed" else "change observed", length(pc$p0), max(abs(w_new - w_ref))))
}

# ---- 6. end to end ---------------------------------------------------------------------------------------------------------------------------------
set.seed(5); e1 <- p1_entries(theta, inp, M = 100, mode = "draw"); set.seed(5); e2 <- p1_entries(theta, inp, M = 100, mode = "draw")
ok(identical(e1, e2), "p1_entries(draw) is not reproducible with the seed")
em <- p1_entries(theta, inp, M = 100, mode = "mean")
ok(all(e1$entry >= 0 & e1$entry <= inp$p1$W_j) && all(em$entry >= 0 & em$entry <= inp$p1$W_j), "entries outside [0, W_j]")
fb <- p1_fallback_by_state(e1, inp, m = 3L); print(fb)
ok(sum(fb$n_plots) == nrow(inp$p1) && sum(fb$n_fallback) == sum(e1$fallback), "p1_fallback_by_state totals")
th_bad <- theta; th_bad$ok[1:2] <- FALSE; eb <- p1_entries(th_bad, inp, M = 50, mode = "draw"); ok(nrow(eb) == nrow(inp$p1), "p1_entries with failed transitions")
# no ok transition at all: nothing can leave regen, so a plot that starts in regen keeps entry 0 (no change observed) or falls back (change observed:
# weight 0); every plot that starts in another state has no candidate and falls back (same behaviour as the WP3 code)
th_none <- theta; th_none$ok[] <- FALSE; en <- p1_entries(th_none, inp, M = 20, mode = "draw")
fs <- inp$p1$first_state; chg <- !is.na(inp$p1$obs_lower)
ok(all(en$fallback[fs != "regen"]) && all(en$fallback[fs == "regen" & chg]) && all(!en$fallback[fs == "regen" & !chg] & en$entry[fs == "regen" & !chg] == 0),
   "no ok transition: unexpected fallback pattern", data.table(fs, chg, fb = en$fallback, entry = en$entry)[, .N, by = .(fs, chg, fb, entry)])
cat(sprintf("6. p1_entries: draw reproducible, %d of %d plots without candidate (random theta), failed-transition and no-transition cases run\n", sum(e1$fallback), nrow(e1)))

# ---- 7. the C5 shift script's function on the synthetic tables (the real run is in WP5) ----------------------------------------------------------
WP4_C5_SHIFT_NO_MAIN <- TRUE; source(here("R", "tests", "wp4_C5_shift.R"))
rs <- c5_shift(inp, theta, which(inp$p1$first_state != "regen"), M_PATHS = 300L)
ok(nrow(rs) > 0 && all(is.finite(rs$entry_old)) && all(is.finite(rs$entry_new)), "c5_shift returned no / non-finite entries")
cat(sprintf("7. c5_shift on %d synthetic plots: mean old - corrected = %.2f y (random theta: information only here)\n", nrow(rs), mean(rs$entry_old - rs$entry_new)))
invisible(capture.output(c5_report(copy(rs))))

# ---- 8. lambda_matrix: clean X gives the matrix of before; an NA in X stops and names the column(s) and the number of rows --------------------------------
lam_before <- function(theta, X) { eta <- X %*% t(theta$beta); lam <- exp(sweep(eta, 2, theta$intercept, "+")); lam[, !theta$ok] <- 0; lam[is.na(lam)] <- 0; lam }   # the WP3 / WP4 version
th_f <- theta; th_f$ok[2] <- FALSE; th_f$intercept[2] <- NA_real_; th_f$beta[2, ] <- 0                                                                       # a failed fit, as extract_theta returns it
for (X in list(inp$Xiv, inp$X0, inp$X1)) ok(identical(lambda_matrix(th_f, X), lam_before(th_f, X)) && identical(lambda_matrix(theta, X), lam_before(theta, X)), "lambda_matrix differs from the version before on a clean X")
Xna <- inp$Xiv; Xna[c(3, 9), "z_gdd_3"] <- NA; Xna[9, "tsd_fire"] <- NA
msg <- tryCatch({ lambda_matrix(theta, Xna); NA_character_ }, error = function(e) conditionMessage(e))
ok(!is.na(msg) && grepl("z_gdd_3 (2 rows)", msg, fixed = TRUE) && grepl("tsd_fire (1 rows)", msg, fixed = TRUE) && grepl(sprintf("2 of %d rows have an NA", nrow(Xna)), msg, fixed = TRUE),
   "lambda_matrix did not stop with the expected message on an NA in X", msg)
cat(sprintf("8. lambda_matrix: identical to the version before on clean X (also with a failed fit); an NA in X stops: \"%s\"\n", msg))

# ---- information: first-visit flags that are off at the interval-1 end visit (per type; no change to P1) -------------------------------------------
off_at_end1 <- function(ir) {
  i1 <- ir$iv[interval == 1L]; fvq <- ir$p1[match(i1$id_pe, id_pe)]; jf <- match(i1$id_pe, ir$p1$id_pe)
  f <- list(fire = list(fv = ir$X1[jf, "fire_f"] == 1, off = i1$fire == 0L), pest = list(fv = ir$X1[jf, "pest_f"] == 1, off = i1$pest == 0L),
            harvest = list(fv = ir$on$harv_lvl[jf] > 0L, off = i1$harv == 0L))
  r <- rbindlist(lapply(names(f), function(t) data.table(type = t, flagged_at_first_visit = sum(f[[t]]$fv), off_at_interval1_end = sum(f[[t]]$fv & f[[t]]$off))))
  r <- rbind(r, data.table(type = "any type", flagged_at_first_visit = sum(Reduce(`|`, lapply(f, `[[`, "fv"))),
                           off_at_interval1_end = sum(Reduce(`|`, lapply(f, function(x) x$fv & x$off)))))
  r[, share_off := off_at_interval1_end / pmax(flagged_at_first_visit, 1L)]; r[]
}
cat("\nInformation (no change to P1): plots flagged at the first visit whose flag is OFF at the interval-1 end visit:\n"); print(off_at_end1(inp))

# ---- real tables (Mac) ------------------------------------------------------------------------------------------------------------------------------
if (REAL) {
  cat("\n== real tables:", TABLE_DIR, "==\n")
  ir <- prepare_fit_inputs(TABLE_DIR); pir <- readRDS(file.path(TABLE_DIR, "plot_info.rds")); setDT(pir)
  ok(!anyNA(pir[, grep("^fv_", names(pir)), with = FALSE]), "NA in fv_ columns of the real plot_info")
  fl <- ir$p1[, id_pe][ir$X1[, "fire_f"] + ir$X1[, "pest_f"] + ir$X1[, "harv_light"] + ir$X1[, "harv_heavy"] > 0]
  cat(sprintf("plots %d | flagged at the first visit %d | unflagged %d (plan, main build: 6,330 unflagged)\n", nrow(ir$p1), length(fl), nrow(ir$p1) - length(fl)))
  pw <- pir[match(ir$p1$id_pe, id_pe)]
  ok(max(abs(ir$on$fire - (ir$p1$W_j - pw$fv_tsd_fire_y)), na.rm = TRUE) < 1e-9 && max(abs(ir$on$harv - (ir$p1$W_j - pw$fv_tsd_harv_y)), na.rm = TRUE) < 1e-9 &&
       max(abs(ir$on$pest - (ir$p1$W_j - pw$fv_tsd_pest_y)), na.rm = TRUE) < 1e-9, "real: onsets != W_j - first-visit TSD")
  sc <- ir$scaling; samp <- sample(nrow(ir$p1), 300)
  for (v in c("gdd_3", "cmi_gs", "soil_pc1")) ok(max(abs(ir$X1[samp, paste0("z_", v)] - (pw[[paste0("fv_", v)]][samp] - sc[var == v, mean]) / sc[var == v, sd])) < 1e-9, paste("real: X1", v))
  cat(sprintf("onset summary (flagged types): fire %d, pest %d, harvest %d | onset <= 0: %d, onset in (W_j-1, W_j]: %d\n", sum(!is.na(ir$on$fire)), sum(!is.na(ir$on$pest)), sum(!is.na(ir$on$harv)),
              sum(unlist(ir$on[, .(fire, pest, harv)]) <= 0, na.rm = TRUE), sum(unlist(ir$on[, .(fire, pest, harv)]) > ir$p1$W_j - 1, na.rm = TRUE)))
  Kr <- nrow(ir$td); Pr <- length(ir$terms); thr <- list(alpha = runif(Kr, 1.1, 2), intercept = log(1 / runif(Kr, 25, 70)), beta = matrix(rnorm(Kr * Pr, 0, .2), Kr, Pr, dimnames = list(ir$td$label, ir$terms)), ok = rep(TRUE, Kr))
  thr$beta[!ir$mask] <- 0; L0r <- lambda_matrix(thr, ir$X0); Livr <- lambda_matrix(thr, ir$Xiv)
  jf <- match(sample(fl, min(300, length(fl))), ir$p1$id_pe); np <- integer(length(jf))
  for (i in seq_along(jf)) { pc <- p1_pieces(jf[i], thr, ir, L0r, Livr); np[i] <- pc$npre
    ok(all(abs(pc$p0[-1] - pc$p1[-length(pc$p1)]) < 1e-9) && abs(pc$p1[pc$npre] - ir$p1$W_j[jf[i]]) < 1e-9 && abs(tail(pc$p1, 1) - ir$p1$obs_upper[jf[i]]) < 1e-9, paste("real: pieces of plot", ir$p1$id_pe[jf[i]])) }
  cat("\nReal tables, information (no change to P1): plots flagged at the first visit whose flag is OFF at the interval-1 end visit:\n"); print(off_at_end1(ir))
  cat(sprintf("pieces before the first visit in %d flagged plots: median %d, max %d; all contiguous from 0 to the end of the weighting bracket\n", length(jf), as.integer(median(np)), max(np)))
}
cat("\nALL CHECKS PASSED\n")
