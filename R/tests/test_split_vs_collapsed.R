# ==============================================================================
# test_split_vs_collapsed.R
#
# Does INLA (weibullsurv, exact times, left truncation) give the same result
# when a sojourn is ONE row group (collapsed) vs split into one row group per
# visit interval (segments with delayed entry)? Covariate is STATIC here, so
# the two must agree. Exit time is exact (post-imputation form), which isolates
# row construction from P1/P2. Two competing destinations test double counting.
#
# READ-ONLY use of the simulation folder (source only, nothing is written there).
# Times are in decades (as after to_time_since_entry with TIME_SCALE = 10).
# ==============================================================================
suppressPackageStartupMessages(library(INLA))
SIM_DIR <- "/Users/lukas/Desktop/RProjects/Quebec-forest-succession/R-scripts/Sim4stateExact/p1p2_v2_test"
source(file.path(SIM_DIR, "model_fitting.R"))   # split_by_transition(), fit_inla_models()

set.seed(2026)
N <- 4000
dests <- c("B", "C")
alpha <- c(B = 1.3, C = 0.9); b0 <- c(B = -2.2, C = -2.6); b1 <- c(B = 0.5, C = -0.4)
W <- sample(c(1, 3, 5), N, TRUE)          # stand age at first visit (decades)
x <- rnorm(N)                             # static covariate
Tdest <- sapply(dests, function(d) (rexp(N) / exp(b0[d] + b1[d] * x))^(1 / alpha[d]))
Tex   <- apply(Tdest, 1, min)             # exact exit time (entry = 0)
dest  <- dests[apply(Tdest, 1, which.min)]
keep  <- which(Tex > W)                   # left truncation: only survivors to W enter the sample
cat(sprintf("N simulated %d, in sample after truncation %d\n", N, length(keep)))

row_ <- function(i, to, status, tstop, trunc)
  data.frame(id = i, from = "A", to = to, status = status, Tstop = tstop, truncation = trunc, x = x[i])

collapsed <- list(); split_ <- list()
for (i in keep) {
  v <- W[i] + c(0, cumsum(runif(8, 0.8, 1.4))); v <- v[v <= W[i] + 4]   # visit times
  observed_exit <- Tex[i] <= max(v)
  if (observed_exit) {
    ks <- which(v >= Tex[i])[1]                     # first visit at/after the exit
    other <- setdiff(dests, dest[i])
    # collapsed: one event row group, truncation = W
    collapsed[[length(collapsed) + 1]] <- rbind(
      row_(i, dest[i], 1L, Tex[i], W[i]), row_(i, other, 0L, Tex[i], W[i]))
    # split: full censored segments before the event segment, then the event segment
    segs <- list()
    for (j in seq_len(ks - 2))
      segs[[length(segs) + 1]] <- rbind(row_(i, dests, 0L, v[j + 1], v[j]))
    segs[[length(segs) + 1]] <- rbind(
      row_(i, dest[i], 1L, Tex[i], v[ks - 1]), row_(i, other, 0L, Tex[i], v[ks - 1]))
    split_[[length(split_) + 1]] <- do.call(rbind, segs)
  } else {
    m <- length(v)
    collapsed[[length(collapsed) + 1]] <- row_(i, dests, 0L, v[m], W[i])
    segs <- lapply(seq_len(m - 1), function(j) row_(i, dests, 0L, v[j + 1], v[j]))
    split_[[length(split_) + 1]] <- do.call(rbind, segs)
  }
}
collapsed <- do.call(rbind, collapsed); split_ <- do.call(rbind, split_)

# ---- CHECK 1: row bookkeeping (no double counting of events) ----------------
cat("\nCHECK 1: rows / events\n")
print(rbind(collapsed = c(rows = nrow(collapsed), events = sum(collapsed$status == 1L),
                          ids_with_event = length(unique(collapsed$id[collapsed$status == 1L])),
                          max_events_per_id = max(table(collapsed$id[collapsed$status == 1L]))),
            split     = c(rows = nrow(split_), events = sum(split_$status == 1L),
                          ids_with_event = length(unique(split_$id[split_$status == 1L])),
                          max_events_per_id = max(table(split_$id[split_$status == 1L])))))
print(table(collapsed = collapsed$to[collapsed$status == 1L]))
print(table(split = split_$to[split_$status == 1L]))

# ---- CHECK 2: exact row-level log-likelihood at the TRUE parameters ---------
# Weibull (INLA variant 1): h = alpha*lam*t^(alpha-1), S = exp(-lam*t^alpha)
ll <- function(d) {
  a <- alpha[d$to]; lam <- exp(b0[d$to] + b1[d$to] * d$x)
  tapply(d$status * (log(a) + log(lam) + (a - 1) * log(d$Tstop)) - lam * d$Tstop^a + lam * d$truncation^a,
         d$to, sum)
}
cat("\nCHECK 2: total log-likelihood at true parameters, per destination\n")
print(rbind(collapsed = ll(collapsed), split = ll(split_), diff = ll(collapsed) - ll(split_)), digits = 12)

# ---- CHECK 3: INLA fits, same call as the pipeline --------------------------
fit_both <- function(d) fit_inla_models(split_by_transition(d, rep("A", 2), dests), "x", exact = TRUE)
mc <- fit_both(collapsed); ms <- fit_both(split_)
tab <- do.call(rbind, lapply(seq_along(dests), function(k) {
  g <- function(m) c(intercept = m$summary.fixed["(Intercept)", "mean"], x = m$summary.fixed["x", "mean"],
                     int_sd = m$summary.fixed["(Intercept)", "sd"], x_sd = m$summary.fixed["x", "sd"],
                     alpha = m$summary.hyperpar[1, "mean"], alpha_sd = m$summary.hyperpar[1, "sd"],
                     mlik = m$mlik[1, 1])
  a <- g(mc[[k]]); b <- g(ms[[k]])
  data.frame(dest = dests[k], quantity = names(a), collapsed = a, split = b, diff = a - b, row.names = NULL)
}))
cat("\nCHECK 3: INLA (true values: B intercept -2.2, x 0.5, alpha 1.3; C intercept -2.6, x -0.4, alpha 0.9)\n")
print(tab, digits = 8)
cat("\nMax |diff| over means/SDs (excluding mlik):", max(abs(tab$diff[tab$quantity != "mlik"])), "\n")
