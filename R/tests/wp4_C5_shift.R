# ==============================================================================
# R/tests/wp4_C5_shift.R   (WP4, C5 real-data check; RUN IN WP5 after the first full imputation of the new run; needs no INLA)
# How much do the OLD weights (survival from entry, S(lo - e) - S(hi - e): double counting of the survival to the first visit) shift the entry
# of the plots, compared with the CORRECTED weights (divided by S(W_j - e), C5)? Same candidates, same piecewise hazards (C1) and the SAME theta
# (theta_m of imputation m of the new run, read from extras_<tag>_m<mm>.rds); only the weights differ, so the difference is the C5 effect alone.
# Plots: first state not regen (the plan's 1,466 plots on the main build), M_PATHS paths each (plan: 1,000).
# PASS: the corrected entries are earlier on average (mean of old - corrected > 0). The size is reported for information only (plan, main build:
# 3.8 y when a change is observed later, 2.9 y when not; time in the first state at the first visit 18.5 -> 22.3 y and 24.0 -> 26.9 y).
#   Rscript R/tests/wp4_C5_shift.R [m = 1] [M_PATHS = 1000]        (project root; writes <FIT_DIR>/wp4_C5_shift_<tag>.csv)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("R", "config.R")); source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))

# plot-level weighted mean entry under old and corrected weights, for plots j (rows of inp$p1)
c5_shift <- function(inp, theta, plots, M_PATHS = 1000L, seed = 1L) {
  L0 <- lambda_matrix(theta, inp$X0); Liv <- lambda_matrix(theta, inp$Xiv); f_tsd <- tsd_fun()
  rbindlist(lapply(plots, function(j) {
    set.seed(seed * 100003L + j)
    pc <- p1_pieces(j, theta, inp, L0, Liv, f_tsd); W <- inp$p1$W_j[j]; sj <- inp$p1$s_idx[j]; lo <- inp$p1$obs_lower[j]; hi <- inp$p1$obs_upper[j]
    cand <- p1_candidates_pw(pc, inp$start_idx, sj, W, theta$alpha, inp$out_k, inp$td$to_idx, theta$ok, M_PATHS)
    if (!length(cand)) return(NULL)
    ks <- inp$out_k[[sj]]; ks <- ks[theta$ok[ks]]; H <- function(a, b) p1_H_between(pc, ks, theta$alpha, cand, a, b)
    lw_old <- if (!is.na(lo)) -H(0, lo) + log1mexp(H(lo, hi)) else -H(0, hi)       # S(lo - e) - S(hi - e)  /  S(hi - e): no division by S(W_j - e)
    wn <- function(lw) { if (anyNA(lw) || all(lw == -Inf)) return(NULL); w <- exp(lw - max(lw)); w / sum(w) }
    w_old <- wn(lw_old); w_new <- p1_weights_pw(cand, sj, lo, hi, W, pc, theta$alpha, inp$out_k, theta$ok)
    if (is.null(w_old) || is.null(w_new)) return(NULL)
    data.table(id_pe = inp$p1$id_pe[j], first_state = inp$p1$first_state[j], change_observed = !is.na(lo), W_j = W, n_cand = length(cand),
               entry_old = sum(cand * w_old), entry_new = sum(cand * w_new))
  }))
}
c5_report <- function(r) {
  r[, shift := entry_old - entry_new]
  s <- function(d) d[, .(n_plots = .N, mean_shift = mean(shift), median_shift = median(shift), p90_shift = quantile(shift, 0.9), max_shift = max(shift),
                         time_in_state_old = mean(W_j - entry_old), time_in_state_new = mean(W_j - entry_new))]
  cat("\nShift of the weighted mean entry, old - corrected (years; positive = old weights put the entry later):\n")
  print(rbind(cbind(group = "all", s(r)), cbind(group = "change observed", s(r[change_observed == TRUE])), cbind(group = "no change", s(r[change_observed == FALSE]))), digits = 3)
  cat("\nBy first state:\n"); print(r[, s(.SD), by = first_state][order(first_state)], digits = 3)
  mean(r$shift)
}

if (!exists("WP4_C5_SHIFT_NO_MAIN")) {
  a <- commandArgs(TRUE); m <- if (length(a) >= 1) as.integer(a[1]) else 1L; M_PATHS <- if (length(a) >= 2) as.integer(a[2]) else 1000L
  inp <- prepare_fit_inputs(TABLE_DIR)
  ex <- readRDS(file.path(FIT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, m))); tm <- ex$theta_m
  stopifnot("extras theta_m does not match the transitions of this table" = identical(tm$label, inp$td$label), "theta_m lacks intercept_dec (WP3 format)" = "intercept_dec" %in% names(tm))
  b0 <- tm$intercept_dec; bv <- as.matrix(tm[, inp$terms, with = FALSE]); bv[is.na(bv)] <- 0; okk <- tm$ok & !is.na(tm$alpha) & tm$alpha > 0
  theta <- list(alpha = tm$alpha, intercept = rep(NA_real_, nrow(tm)), beta = bv, ok = okk)
  cv <- v0_to_p1p2(b0[okk], bv[okk, , drop = FALSE], tm$alpha[okk]); theta$intercept[okk] <- cv$intercept; theta$beta[okk, ] <- cv$beta; theta$beta[!okk, ] <- 0
  plots <- which(inp$p1$first_state != "regen")
  cat(sprintf("imputation m=%d | %d plots with a first state other than regen | %d paths each\n", m, length(plots), M_PATHS))
  t0 <- proc.time()[["elapsed"]]; r <- c5_shift(inp, theta, plots, M_PATHS); cat(sprintf("done in %.0f s; %d plots with usable candidates\n", proc.time()[["elapsed"]] - t0, nrow(r)))
  ms <- c5_report(r); fwrite(r, file.path(FIT_DIR, sprintf("wp4_C5_shift_%s.csv", RUN_TAG)))
  if (!(ms > 0)) stop(sprintf("CHECK FAILED: corrected entries are NOT earlier on average (mean old - corrected = %.3f y). Report to the user; nothing was changed.", ms))
  cat(sprintf("\nPASS: corrected entries are earlier on average (mean shift %.2f y; plan, main build: 2.9-3.8 y, information only)\n", ms))
}
