# ==============================================================================
# R/tests/wp4_C5_toy.R   (WP4, C5 check)   toy: regen -> A -> B, exact posterior mean of the entry time into A by integration
# Toy supplied by the user 2026-11-04 (the plan's 27.3 / 37.3 / 27.4 came from a version that cannot be found): entry into A ~ Weibull(a1, l1),
# sojourn in A ~ Weibull(a2, l2), first visit W = 60 (the plot is in A), the exit to B is observed in (lo, hi] = (70, 80].
# PASS: the corrected weights give the exact posterior mean within 3 Monte Carlo SE, and the current (double-counting) weights are clearly outside
# (> 5 SE). Run twice: (1) the toy's own weight formulas, (2) the production functions p1_candidates_pw / p1_weights_pw (05a) on the same toy;
# the old p1_weights is run for the "current" value. The Monte Carlo SE of a self-normalised mean is sqrt(sum(w^2 (x - xbar)^2)) / sum(w).
#   Rscript R/tests/wp4_C5_toy.R          (from the project root; no INLA)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("R", "config.R")); source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
fail <- function(msg) stop("CHECK FAILED: ", msg, " (bounds were NOT changed; report this output to the user)", call. = FALSE)
set.seed(1)
a1 <- 1.5; l1 <- 1/30; a2 <- 1.8; l2 <- 1/40; W <- 60; lo <- 70; hi <- 80
SA <- function(d) exp(-(l2 * pmax(d, 0))^a2)
f1 <- function(e) a1 * l1^a1 * e^(a1 - 1) * exp(-(l1 * e)^a1)
post <- function(e) f1(e) * (SA(lo - e) - SA(hi - e))
Z <- integrate(post, 0, W)$value; exact <- integrate(function(e) e * post(e), 0, W)$value / Z
wm <- function(x, w) { xb <- sum(x * w) / sum(w); c(mean = xb, se = sqrt(sum(w^2 * (x - xb)^2)) / sum(w)) }
judge <- function(tag, cur, fix, ncand) {
  zc <- (cur[["mean"]] - exact) / cur[["se"]]; zf <- (fix[["mean"]] - exact) / fix[["se"]]
  cat(sprintf("%-26s candidates %6d | exact %.2f | current %.2f (SE %.2f, %.1f SE away) | corrected %.2f (SE %.2f, %.1f SE away)\n",
              tag, ncand, exact, cur[["mean"]], cur[["se"]], zc, fix[["mean"]], fix[["se"]], zf))
  if (abs(zf) >= 3) fail(sprintf("%s: corrected mean is %.1f SE from the exact mean (limit 3)", tag, abs(zf)))
  if (abs(zc) <= 5) fail(sprintf("%s: current weights are only %.1f SE from the exact mean (need > 5: the toy must show the bug)", tag, abs(zc)))
}

# (1) the toy's own formulas
M <- 2e5; e <- (-log(runif(M)))^(1 / a1) / l1; e <- e[e < W]
tA <- e + (-log(runif(length(e))))^(1 / a2) / l2; cand <- e[tA >= W]
w_now <- SA(lo - cand) - SA(hi - cand); w_fix <- w_now / SA(W - cand)
judge("toy formulas", wm(cand, w_now), wm(cand, w_fix), length(cand))

# (2) production functions: pieces = [0, W] (regen->A and A->B constant hazard) and the observed interval [W, hi]
set.seed(2)
L <- rbind(c(l1, l2), c(l1, l2)); al <- c(a1, a2)
pc <- list(p0 = c(0, W), p1 = c(W, hi), L = L, LA = sweep(L, 2, al, `^`), npre = 1L)
out_k <- list(1L, 2L, integer(0)); to_idx <- c(2L, 3L)       # states: 1 regen, 2 A, 3 B ; k = 1 regen->A, k = 2 A->B
cand2 <- p1_candidates_pw(pc, start_idx = 1L, s_j = 2L, W = W, alpha = al, out_k = out_k, to_idx = to_idx, ok = c(TRUE, TRUE), M = 1e5)
w_fix2 <- p1_weights_pw(cand2, 2L, lo, hi, W, pc, al, out_k, c(TRUE, TRUE))
w_cur2 <- p1_weights(cand2, 2L, lo, hi, al, c(l1, l2), out_k)       # the old function (survival from entry, no division by S(W - e))
judge("production functions", wm(cand2, w_cur2), wm(cand2, w_fix2), length(cand2))
cat("ALL CHECKS PASSED\n")
