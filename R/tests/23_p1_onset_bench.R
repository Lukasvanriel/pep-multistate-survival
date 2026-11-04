# ==============================================================================
# R/tests/23_p1_onset_bench.R   (2026-11-04)   TIMING TEST ONLY -- changes nothing in the pipeline; no INLA needed
# Benchmark: current P1 (constant covariates, one-formula draw) vs P1 with the disturbance flags switching on at their
# onset (onset on P1's clock = W_j - TSD at the first visit). Three versions of the onset P1 (MODE):
#   onset  = TSD fixed at its first-visit value (2-3 pieces per plot)
#   log    = pieces at onset + 1, 3, 7, 15, 31, 63, 127 y, TSD at the piece midpoint
#   yearly = yearly pieces after the first onset, TSD at the piece midpoint
# Uses the main-run tables and the theta draw of imputation m01 (variant-1 form, years; current model terms, TSD in decades).
# The onset P1 is a PROTOTYPE (plain R loops), for timing only.
# Run from the project root:   Rscript R/tests/23_p1_onset_bench.R        (env NB = plots per sample, default 300)
# ==============================================================================
suppressMessages({ library(data.table); library(here) })
source(here::here("R", "config.R"))   # TABLE_DIR, FIT_DIR, BUILD_DIR, RUN_TAG
SIM_M <- 200L; M_MI <- 30L; SD_MIN <- 1e-3; N_COV_SAMPLES <- 0L
source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
inp <- prepare_fit_inputs(TABLE_DIR)
ex  <- readRDS(file.path(FIT_DIR, sprintf("extras_%s_m01.rds", RUN_TAG)))
tm  <- ex$theta_m
theta <- list(alpha = tm$alpha, intercept = tm$intercept_years, ok = tm$ok,
              beta = as.matrix(tm[, inp$terms, with = FALSE]))
theta$beta[is.na(theta$beta)] <- 0
stopifnot(identical(tm$label, inp$td$label))
cpu_inla <- sum(vapply(ex$fits, function(f) if (is.null(f$cpu)) NA_real_ else f$cpu[["Total"]], 0), na.rm = TRUE)

# ---- first-visit disturbance info (onset on P1's clock = W_j - TSD at first visit) ----
cc <- as.data.table(readRDS(file.path(BUILD_DIR, "combined_covariates.rds")))
cc[, no_mes := as.integer(no_mes)]; setorder(cc, id_pe, no_mes)
fv <- cc[, .SD[1L], by = id_pe][, .(id_pe, fire, harv, pest, tf = tsf_imputed_mid, th = tsh_imputed_mid, tp = sfp_imputed_mid)]
p1 <- copy(inp$p1)[fv, on = "id_pe", `:=`(fire = i.fire, harv = i.harv, pest = i.pest, tf = i.tf, th = i.th, tp = i.tp)]
p1[, flagged := fire == 1L | harv > 0L | pest == 1L]
cat(sprintf("plots %d | flagged at first visit %d\n", nrow(p1), sum(p1$flagged)))

# ---- current P1 for one plot (as p1_entries' inner function, mode "draw") ----
L0 <- lambda_matrix(theta, inp$X)
one_current <- function(j) {
  cand <- p1_candidates(inp$start_idx, p1$s_idx[j], p1$W_j[j], theta$alpha, L0[j, ], inp$out_k, inp$td$to_idx, SIM_M)
  if (!length(cand)) return(NA_real_)
  w <- p1_weights(cand, p1$s_idx[j], p1$obs_lower[j], p1$obs_upper[j], theta$alpha, L0[j, ], inp$out_k)
  if (is.null(w)) NA_real_ else cand[sample.int(length(cand), 1L, prob = w)]
}

# ---- stepped P1 for one flagged plot ----
# pieces: [0, first onset) constant; then yearly pieces up to the end of the weighting bracket; covariates per piece from the
# base row (flags zeroed) plus every type that has started, with TSD at the piece midpoint (current model: TSD in decades).
flag_cols <- c("fire_f", "pest_f", "harv_light", "harv_heavy", "harv_any", "gdd_x_fire", "gdd_x_pest", "gdd_x_harv_light",
               "gdd_x_harv_heavy", "gdd_x_harv_any", "tsd_fire", "tsd_pest", "tsd_harv")
pieces_for <- function(j) {
  W <- p1$W_j[j]; tend <- max(W, p1$obs_upper[j], na.rm = TRUE)
  x0 <- inp$X[j, ]; x0[flag_cols] <- 0; zg <- x0[["z_gdd_3"]]
  on <- c(fire = if (p1$fire[j] == 1L) W - p1$tf[j] else NA, harv = if (p1$harv[j] > 0L) W - p1$th[j] else NA,
          pest = if (p1$pest[j] == 1L) W - p1$tp[j] else NA)
  first_on <- max(0, min(on, na.rm = TRUE))
  brk <- switch(MODE,
    yearly = unique(c(0, if (first_on > 0) first_on, seq(first_on, tend, by = 1)[-1])),
    log    = sort(unique(c(0, unlist(lapply(on[!is.na(on)], function(o) o + c(0, 1, 3, 7, 15, 31, 63, 127)))))),
    onset  = sort(unique(c(0, on[!is.na(on)]))))
  brk <- brk[brk >= 0 & brk < tend]; p0 <- brk; p1e <- c(brk[-1], tend)
  mid <- (p0 + p1e) / 2
  tsd_of <- function(o) if (MODE == "onset") rep(W - o, length(mid)) else mid - o   # onset: TSD fixed at the first-visit value
  Xp <- matrix(x0, length(p0), length(x0), byrow = TRUE, dimnames = list(NULL, names(x0)))
  if (!is.na(on["fire"])) { a <- mid >= on["fire"]; Xp[a, "fire_f"] <- 1; Xp[a, "gdd_x_fire"] <- zg; Xp[a, "tsd_fire"] <- tsd_of(on["fire"])[a] / 10 }
  if (!is.na(on["pest"])) { a <- mid >= on["pest"]; Xp[a, "pest_f"] <- 1; Xp[a, "gdd_x_pest"] <- zg; Xp[a, "tsd_pest"] <- tsd_of(on["pest"])[a] / 10 }
  if (!is.na(on["harv"])) { a <- mid >= on["harv"]; lv <- if (p1$harv[j] == 2L) "harv_heavy" else "harv_light"
    Xp[a, lv] <- 1; Xp[a, paste0("gdd_x_", lv)] <- zg; Xp[a, "harv_any"] <- 1; Xp[a, "gdd_x_harv_any"] <- zg; Xp[a, "tsd_harv"] <- tsd_of(on["harv"])[a] / 10 }
  list(p0 = p0, p1 = p1e, L = lambda_matrix(theta, Xp))   # pieces x K
}
# exact draw of the next exit for a sojourn entered at e, current time t_now, under piecewise-constant lambda (variant 1 form)
next_exit <- function(pc, ks, alpha, e, t_now, W) {
  jj <- which(pc$p1 > t_now & pc$p0 < W); if (!length(jj)) return(c(Inf, NA))
  a0 <- pmax(pc$p0[jj], t_now) - e; a1 <- pmin(pc$p1[jj], W) - e
  Ek <- -log(runif(length(ks)))
  best <- Inf; bk <- NA
  for (i in seq_along(ks)) { k <- ks[i]; al <- alpha[k]; la <- pc$L[jj, k]^al
    dH <- la * (a1^al - a0^al); cs <- cumsum(dH); q <- which(cs >= Ek[i])[1L]
    if (is.na(q)) next
    prev <- if (q > 1L) cs[q - 1L] else 0
    tt <- e + (a0[q]^al + (Ek[i] - prev) / la[q])^(1 / al)
    if (tt < best) { best <- tt; bk <- k } }
  c(best, bk)
}
p1_candidates_tv <- function(j, pc, M) {
  W <- p1$W_j[j]; s_j <- p1$s_idx[j]
  if (W <= 0) return(if (s_j == inp$start_idx) rep(0, M) else numeric(0))
  cand <- numeric(M); n <- 0L
  for (i in seq_len(M)) {
    state <- inp$start_idx; entry <- 0; t_now <- 0
    repeat {
      ks <- inp$out_k[[state]]; ks <- ks[theta$ok[ks]]; if (!length(ks)) break
      nx <- next_exit(pc, ks, theta$alpha, entry, t_now, W)
      if (!is.finite(nx[1]) || nx[1] >= W) break
      state <- inp$td$to_idx[nx[2]]; entry <- nx[1]; t_now <- entry }
    if (state == s_j) { n <- n + 1L; cand[n] <- entry } }
  cand[seq_len(n)]
}
# piecewise log-survival from each candidate entry to time t (all exits of s_j)
logS_tv <- function(pc, ks, alpha, cand, t) {
  o <- numeric(length(cand))
  for (k in ks) { al <- alpha[k]; la <- pc$L[, k]^al
    for (c in seq_along(cand)) { a0 <- pmax(pc$p0, cand[c]) - cand[c]; a1 <- pmin(pc$p1, t) - cand[c]; ok <- a1 > a0
      o[c] <- o[c] - sum(la[ok] * (a1[ok]^al - a0[ok]^al)) } }
  o
}
one_tv <- function(j) {
  pc <- pieces_for(j)
  cand <- p1_candidates_tv(j, pc, SIM_M); if (!length(cand)) return(NA_real_)
  ks <- inp$out_k[[p1$s_idx[j]]]; ks <- ks[theta$ok[ks]]; lo <- p1$obs_lower[j]; hi <- p1$obs_upper[j]
  w <- if (!is.na(lo)) { ll <- logS_tv(pc, ks, theta$alpha, cand, lo); lh <- logS_tv(pc, ks, theta$alpha, cand, hi); ifelse(ll > lh, exp(ll) - exp(lh), 0) } else exp(logS_tv(pc, ks, theta$alpha, cand, hi))
  if (!all(is.finite(w)) || sum(w) <= 0) NA_real_ else cand[sample.int(length(cand), 1L, prob = w)]
}

# ---- timing on the same sample of flagged plots ----
set.seed(1); N <- as.integer(Sys.getenv("NB", "300"))
fl <- sample(which(p1$flagged), N); un <- sample(which(!p1$flagged), N)
tm_cur_fl <- system.time(r_cur <- vapply(fl, one_current, 0))[["elapsed"]]
tm_cur_un <- system.time(vapply(un, one_current, 0))[["elapsed"]]
n_fl <- sum(p1$flagged); n_un <- sum(!p1$flagged)
cur_total <- tm_cur_fl / N * n_fl + tm_cur_un / N * n_un
cat(sprintf("\nsample %d flagged + %d unflagged plots, single core, M = %d paths per plot\n", N, N, SIM_M))
cat(sprintf("current: %.4f s per flagged plot, %.4f per unflagged | all plots %.2f min | fallbacks in sample %d\n",
            tm_cur_fl / N, tm_cur_un / N, cur_total / 60, sum(is.na(r_cur))))
for (MODE in c("onset", "log", "yearly")) {
  set.seed(2)
  tm_tv_fl <- system.time(r_tv <- vapply(fl, one_tv, 0))[["elapsed"]]
  npc <- vapply(fl, function(j) length(pieces_for(j)$p0), 0L)
  new_total <- tm_tv_fl / N * n_fl + tm_cur_un / N * n_un
  cat(sprintf("%-7s: %.4f s per flagged plot | pieces median %d, max %d | all plots %.2f min (x%.1f) | fallbacks in sample %d\n",
              MODE, tm_tv_fl / N, as.integer(median(npc)), max(npc), new_total / 60, new_total / cur_total, sum(is.na(r_tv))))
}
cat(sprintf("INLA time of imputation m01 (sum of Total over its 90 fits, your Mac): %.1f min\n", cpu_inla / 60))
