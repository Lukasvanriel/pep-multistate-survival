# ==============================================================================
# 05a_fit_funs.R   (functions for 05_fit_inla.R; no INLA needed to source this)
#
# Real-data adaptation of the simulation's P1 / P2 / theta-draw / pooling code
# (fit_p1_rungs_2a_2b.R, fit_p1_rung5_mi.R, p2_imputation_strategies.R,
# fit_p1p2_combined.R). The simulation folder is only READ; what is reused is
# copied here and adapted. Differences from the simulation:
#   * transitions = the KEPT transitions of transition_set (91), each with its
#     own term set (zeros in beta for terms a transition does not have)
#   * P1 weighting bracket from the real visits (no hidden-hop flag, no oracle)
#   * P1 finds no candidate  -> entry W_j/2, flagged (the simulation dropped the plot)
#   * transition with a failed fit (theta NA) -> hazard 0 in P1; P2: midpoint (Stage 1/2) or uniform draw (imputations), counted (WP5)
# Units (WP3, B3/B4): the models are fitted as Weibull VARIANT 0 on the decade scale, S(t) = exp(-e^eta (t/10)^alpha). The INLA
# intercept is kept on that scale everywhere (tidy_fit, mi_raw, screen, pooling). P1 and P2 work in the variant-1 form in YEARS,
# S(t) = exp(-(lambda1 t)^alpha); theta (extract_theta, draw_theta_m) holds the CONVERTED values (v0_to_p1p2) plus the raw variant-0 values.
# ==============================================================================

library(data.table)

# defaults only if the driver script has not set them (05_fit_inla.R sets them first)
if (!exists("SIM_M"))     SIM_M     <- 200L   # P1 candidate paths per plot (simulation value)
if (!exists("M_MI"))      M_MI      <- 30L    # imputations (simulation value)
if (!exists("N_WORKERS")) N_WORKERS <- 6L     # simulation value; adjust to the machine
if (!exists("SD_MIN")) SD_MIN <- 1e-3   # alpha posterior SD below this = near-degenerate fit (screening threshold)
if (!exists("N_COV_SAMPLES")) N_COV_SAMPLES <- 1000L   # posterior samples per fit for the joint covariance (0 = skip)
LOG_TIME_SCALE <- log(10)   # TIME_SCALE = 10 (04b_build_rows.R)
ALPHA_NAME     <- "alpha parameter for weibullsurv"
CAP            <- -log(0.001 * 2^-23)   # 22.85: INLA clamps F(truncation) at 1 - 0.001 * FLT_EPSILON (B6)

# ------------------------------------------------------------------------------
# fit_inla_models_pc() (2026-10-26, T5/T6): pipeline copy of the simulation's fit_inla_models() (model_fitting.R, read-only; both branches
# copied unchanged) with these changes: an explicit PC prior on the Weibull shape, pc.alphaw(ALPHA_PRIOR) (R/config.R), one thread per fit,
# and (WP3, B1/B2) the Weibull variant from WEIBULL_VARIANT and control.fixed = N(0, COEF_PRIOR_SD^2) on the coefficients, flat intercept.
# Checked 2026-10-26: with ALPHA_PRIOR = 5 (variant 1, default fixed priors) the fits were identical to the original (max difference 3e-7 SD).
# ------------------------------------------------------------------------------
fit_inla_models_pc <- function(data_list, rhs_formula, exact = FALSE, int_strategy = "grid") {
  stopifnot("ALPHA_PRIOR must be set in R/config.R" = exists("ALPHA_PRIOR"),
            "WEIBULL_VARIANT must be set in R/config.R" = exists("WEIBULL_VARIANT"),
            "COEF_PRIOR_SD must be set in R/config.R" = exists("COEF_PRIOR_SD"))
  K <- length(data_list)
  if (length(rhs_formula) == 1) rhs_formula <- rep(rhs_formula, K)
  stopifnot(length(rhs_formula) == K)
  lapply(seq_len(K), function(k) {
    d <- data_list[[k]]
    if (nrow(d) == 0) return(list(error = "no observed rows for this transition"))
    if (exact) {
      f_str <- sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", rhs_formula[k])
    } else {
      d$time_lower <- ifelse(d$status == 3L, d$Tstart, d$Tstop)
      f_str <- sprintf("inla.surv(time = time_lower, time2 = Tstop, event = status, truncation = truncation) ~ %s", rhs_formula[k])
    }
    tryCatch(INLA::inla(as.formula(f_str), family = "weibullsurv", data = d,
                        control.family = list(variant = WEIBULL_VARIANT, hyper = list(alpha = list(prior = "pc.alphaw", param = ALPHA_PRIOR))),
                        control.fixed = list(mean = 0, prec = 1 / COEF_PRIOR_SD^2, mean.intercept = 0, prec.intercept = 0),
                        control.compute = list(config = TRUE),
                        control.inla = list(strategy = "adaptive", int.strategy = int_strategy, control.vb = list(enable = FALSE)),
                        num.threads = "1:1"),
             error = function(e) list(error = conditionMessage(e)))
  })
}

# ------------------------------------------------------------------------------
# prepare_fit_inputs(): everything the stages need, built once.
# ------------------------------------------------------------------------------
prepare_fit_inputs <- function(proc_dir, subsample_n = NULL, seed = 1L) {
  iv <- readRDS(file.path(proc_dir, "visit_intervals.rds")); setDT(iv)
  pinfo <- readRDS(file.path(proc_dir, "plot_info.rds"));   setDT(pinfo)
  tr <- readRDS(file.path(proc_dir, "transition_set.rds")); setDT(tr)
  setorder(iv, id_pe, interval)
  scaling <- make_scaling(iv)                          # scaling from ALL intervals, also when subsampling
  iv <- add_model_covariates(iv, scaling)
  if (!is.null(subsample_n)) {
    set.seed(seed); keep <- sample(unique(iv$id_pe), subsample_n)
    iv <- iv[id_pe %in% keep]; pinfo <- pinfo[id_pe %in% keep]
  }
  states <- sort(unique(c(tr$from, tr$to)))
  td <- tr[kept == TRUE]
  td[, `:=`(k = .I, from_idx = match(from, states), to_idx = match(to, states), label = paste(from, to, sep = "->"))]
  mask <- term_mask(tr); stopifnot(identical(rownames(mask), td$label))
  terms <- all_terms()
  kmat <- matrix(NA_integer_, length(states), length(states), dimnames = list(states, states))
  kmat[cbind(td$from, td$to)] <- td$k
  out_k <- lapply(seq_along(states), function(s) td$k[td$from_idx == s])

  # P1 inputs, one line per plot in iv order: covariates of the FIRST interval, W_j, first state, bracket
  f1 <- iv[interval == 1L]
  chg <- iv[is_change == TRUE, .(lo = t_start[1L], hi = t_end[1L]), by = id_pe]
  lastv <- iv[, .(last_end = t_end[.N]), by = id_pe]
  p1 <- pinfo[match(f1$id_pe, id_pe), .(id_pe, W_j, first_state)]
  p1[, `:=`(last_end = lastv$last_end[match(id_pe, lastv$id_pe)], lo = chg$lo[match(id_pe, chg$id_pe)], hi = chg$hi[match(id_pe, chg$id_pe)])]
  p1[, `:=`(obs_lower = lo, obs_upper = fifelse(is.na(lo), last_end, hi), s_idx = match(first_state, states))]
  stopifnot(identical(p1$id_pe, f1$id_pe), !anyNA(p1$s_idx))
  X <- as.matrix(f1[, ..terms])    # LEGACY (end-visit values of interval 1): P1 no longer reads it since WP4 (C1); kept for the older diagnostics

  # WP4 / C1: first-visit covariates of 04 (plot_info fv_* columns), scaled with the interval scaling.
  #   X1 = first-visit values as observed (flags and TSD at the first visit); X0 = the same with every disturbance flag and TSD set to 0
  #   (climate, soil and GDD unchanged); on = onset of each flagged type on P1's clock = W_j - first-visit TSD (C1.3). harv_lvl = 1 light / 2 heavy.
  fvn <- c("fv_gdd_3", "fv_cmi_gs", "fv_soil_pc1", "fv_fire", "fv_harv", "fv_pest", "fv_tsd_fire_y", "fv_tsd_pest_y", "fv_tsd_harv_y")
  if (!all(fvn %in% names(pinfo))) stop("plot_info.rds has no first-visit columns (", paste(setdiff(fvn, names(pinfo)), collapse = ", "), "): rerun R/04_prepare_inla.R (WP4)")
  fv <- pinfo[match(p1$id_pe, id_pe)]
  stopifnot("A1.5: NA in a first-visit column of plot_info" = !anyNA(fv[, ..fvn]))
  d1 <- data.table(gdd_3 = fv$fv_gdd_3, cmi_gs = fv$fv_cmi_gs, soil_pc1 = fv$fv_soil_pc1, fire = fv$fv_fire, harv = fv$fv_harv, pest = fv$fv_pest,
                   tsd_fire_y = fv$fv_tsd_fire_y, tsd_pest_y = fv$fv_tsd_pest_y, tsd_harv_y = fv$fv_tsd_harv_y)
  X1 <- as.matrix(add_model_covariates(d1, scaling)[, ..terms])
  d0 <- copy(d1); d0[, `:=`(fire = 0L, harv = 0L, pest = 0L, tsd_fire_y = 0, tsd_pest_y = 0, tsd_harv_y = 0)]
  X0 <- as.matrix(add_model_covariates(d0, scaling)[, ..terms])
  on <- data.table(fire = fifelse(fv$fv_fire == 1L, p1$W_j - fv$fv_tsd_fire_y, NA_real_), pest = fifelse(fv$fv_pest == 1L, p1$W_j - fv$fv_tsd_pest_y, NA_real_),
                   harv = fifelse(fv$fv_harv > 0L, p1$W_j - fv$fv_tsd_harv_y, NA_real_), harv_lvl = fv$fv_harv)
  # observed intervals of each plot are contiguous rows of iv, the first starting at W_j (the weights run piecewise over them, C1.6)
  iv_first <- match(p1$id_pe, iv$id_pe); iv_n <- iv[, .N, by = id_pe]$N[match(p1$id_pe, unique(iv$id_pe))]
  same <- iv$id_pe[-1L] == iv$id_pe[-nrow(iv)]
  stopifnot("P1: intervals of a plot are not contiguous rows of iv" = identical(iv$id_pe[iv_first + iv_n - 1L], p1$id_pe),
            "P1: first interval does not start at W_j" = all(abs(iv$t_start[iv_first] - p1$W_j) < 1e-9),
            "P1: intervals of a plot are not contiguous in time" = all(abs(iv$t_start[-1L] - iv$t_end[-nrow(iv)])[same] < 1e-9))

  list(iv = iv, tr = tr, td = td, subsampled = !is.null(subsample_n), states = states, mask = mask, terms = terms, rhs = rhs_formulas(mask),
       kmat = kmat, out_k = out_k, p1 = p1, X = X, X1 = X1, X0 = X0, on = on, iv_first = iv_first, iv_n = iv_n, Xiv = as.matrix(iv[, ..terms]),
       start_idx = match("regen", states), scaling = scaling,
       keep_cols = c("fire", "harv", "pest", "year_end", terms))
}

# ------------------------------------------------------------------------------
# v0_to_p1p2() (B4): variant-0 fit (intercept b0 on the decade scale, coefficients beta, shape alpha) -> variant-1 form in YEARS used by
# P1 / P2. S = exp(-e^eta (t/10)^alpha) = exp(-(lambda1 t)^alpha) with log lambda1 = eta/alpha - log(10): intercept = b0/alpha - log 10,
# beta = beta/alpha. beta may be a vector (one fit) or a matrix with one row per fit (alpha then has one value per row).
# ------------------------------------------------------------------------------
v0_to_p1p2 <- function(b0, beta, alpha) {
  stopifnot("v0_to_p1p2: only Weibull variant 0 is supported (WEIBULL_VARIANT in R/config.R)" = identical(as.integer(WEIBULL_VARIANT), 0L))
  list(intercept = b0 / alpha - LOG_TIME_SCALE, beta = beta / alpha)
}

# ------------------------------------------------------------------------------
# theta: list(alpha, intercept [years, converted], beta [K x P, converted, 0 = term absent], ok [K],
#             intercept_dec [raw variant-0 intercept, decade scale], beta_v0 [K x P, raw variant-0 coefficients])
# extract_theta(): posterior means (conversion evaluated at the means; used for Stage-2 seeds only).
# ------------------------------------------------------------------------------
extract_theta <- function(models, inp, stat = "mean") {
  K <- nrow(inp$td); P <- length(inp$terms)
  alpha <- rep(NA_real_, K); b0 <- rep(NA_real_, K)
  beta0 <- matrix(0, K, P, dimnames = list(inp$td$label, inp$terms))
  for (k in seq_len(K)) {
    m <- models[[k]]; if (!is.null(m$error)) next
    a <- tryCatch(m$summary.hyperpar[ALPHA_NAME, stat], error = function(e) NA_real_)
    ic <- tryCatch(m$summary.fixed["(Intercept)", stat], error = function(e) NA_real_)
    alpha[k] <- if (length(a) == 1) a else NA_real_
    b0[k]    <- if (length(ic) == 1) ic else NA_real_
    for (tm in colnames(inp$mask)[inp$mask[k, ]]) {
      b <- tryCatch(m$summary.fixed[tm, stat], error = function(e) NA_real_)
      beta0[k, tm] <- if (length(b) == 1) b else NA_real_
    }
  }
  ok <- !is.na(alpha) & alpha > 0 & !is.na(b0) & !apply(is.na(beta0), 1, any)
  beta0[!ok, ] <- 0
  icpt <- rep(NA_real_, K); beta <- beta0
  if (any(ok)) { cv <- v0_to_p1p2(b0[ok], beta0[ok, , drop = FALSE], alpha[ok]); icpt[ok] <- cv$intercept; beta[ok, ] <- cv$beta }
  list(alpha = alpha, intercept = icpt, beta = beta, ok = ok, intercept_dec = b0, beta_v0 = beta0)
}

# seed > 0 and a single thread make INLA's sampler reproducible (seed = 0 means a clock seed)
default_sampler <- function(n, m, seed = 0L) INLA::inla.posterior.sample(n, m, seed = seed, num.threads = "1:1")

# One posterior draw per transition (variant-0 intercept, coefficients and alpha of the SAME draw, converted with that draw's alpha).
# A transition whose fit failed, or whose names are not found, falls back to the seed point estimate (theta_seed); if that
# is not ok either it stays ok = FALSE (hazard 0 / midpoint in the imputation). th$intercept_dec / th$beta_v0 hold the raw variant-0 values.
draw_theta_m <- function(models, inp, theta_seed, sampler = default_sampler) {
  K <- nrow(inp$td); P <- length(inp$terms)
  th <- list(alpha = theta_seed$alpha, intercept = theta_seed$intercept, beta = theta_seed$beta, ok = theta_seed$ok,
             intercept_dec = theta_seed$intercept_dec, beta_v0 = theta_seed$beta_v0)
  n_fallback <- 0L
  for (k in seq_len(K)) {
    m <- models[[k]]; if (!is.null(m$error)) { n_fallback <- n_fallback + 1L; next }
    samp <- sampler(1, m, seed = sample.int(1e6, 1L))[[1]]   # seed from R's RNG: reproducible after set.seed()
    ln <- rownames(samp$latent)
    ii <- grepl("^\\(Intercept\\)", ln)
    tms <- colnames(inp$mask)[inp$mask[k, ]]
    idx <- lapply(tms, function(tm) grepl(paste0("^", tm, "($|:)"), ln))
    a <- if (ALPHA_NAME %in% names(samp$hyperpar)) samp$hyperpar[[ALPHA_NAME]] else NA_real_
    if (!any(ii) || !all(vapply(idx, any, logical(1))) || !is.finite(a) || a <= 0) {
      n_fallback <- n_fallback + 1L; warning(sprintf("draw_theta_m(): names not found (or alpha not > 0) for %s, using seed estimate", inp$td$label[k])); next
    }
    b0 <- samp$latent[ii, 1][1]
    th$beta_v0[k, ] <- 0
    for (j in seq_along(tms)) th$beta_v0[k, tms[j]] <- samp$latent[idx[[j]], 1][1]
    cv <- v0_to_p1p2(b0, th$beta_v0[k, ], a)
    th$intercept_dec[k] <- b0; th$alpha[k] <- a
    th$intercept[k] <- cv$intercept; th$beta[k, ] <- cv$beta
    th$ok[k] <- TRUE
  }
  th$n_fallback <- n_fallback
  th
}

# lambda (years) per plot x transition from a covariate matrix; ok = FALSE -> 0
# WP4 review: the former `lam[is.na(lam)] <- 0` is gone. It was not needed for missing coefficients (a failed fit has ok = FALSE and is zeroed by the line below;
# an ok transition has no NA coefficient), and for an NA in X it silently gave ALL hazards of that row 0 (0 * NA = NA in X %*% t(beta)). An NA in X now stops.
lambda_matrix <- function(theta, X) {
  if (anyNA(X)) {
    nc <- colSums(is.na(X)); nc <- nc[nc > 0]
    stop(sprintf("lambda_matrix: NA in the covariate matrix X: column(s) %s; %d of %d rows have an NA (an NA would silently give every hazard of that row the value 0; fix the covariates, cf. A1)",
                 paste0(if (is.null(names(nc))) which(colSums(is.na(X)) > 0) else names(nc), " (", nc, " rows)", collapse = ", "), sum(rowSums(is.na(X)) > 0), nrow(X)), call. = FALSE)
  }
  eta <- X %*% t(theta$beta)
  lam <- exp(sweep(eta, 2, theta$intercept, "+"))
  lam[, !theta$ok] <- 0
  lam
}

# ------------------------------------------------------------------------------
# P1: candidate entries (ported from sim_entry_candidates_2b) and bracket weights
# (compute_weights_2b). Forward paths from regen at time 0 under the plot's lambda;
# a path contributes ONE candidate = the entry time of the state it occupies at W_j,
# only if that state is the plot's first observed state.
# ------------------------------------------------------------------------------
p1_candidates <- function(start_idx, s_j, W_j, alpha, lam, out_k, to_idx, M) {
  if (W_j <= 0) return(if (s_j == start_idx) rep(0, M) else numeric(0))
  cand <- numeric(M); n <- 0L
  for (i in seq_len(M)) {
    state <- start_idx; entry <- 0; t_now <- 0
    repeat {
      ks <- out_k[[state]]; ks <- ks[lam[ks] > 0]
      if (!length(ks)) break
      tc <- t_now + (-log(runif(length(ks))))^(1 / alpha[ks]) / lam[ks]
      w <- which.min(tc)
      if (tc[w] >= W_j) break
      state <- to_idx[ks[w]]; entry <- tc[w]; t_now <- entry
    }
    if (state == s_j) { n <- n + 1L; cand[n] <- entry }
  }
  cand[seq_len(n)]
}

p1_weights <- function(cand, s_j, lo, hi, alpha, lam, out_k) {
  ks <- out_k[[s_j]]; ks <- ks[lam[ks] > 0]
  logS <- function(dur) { o <- numeric(length(dur)); for (k in ks) o <- o - pmax(lam[k] * dur, 0)^alpha[k]; o }
  w <- if (!is.na(lo)) { ll <- logS(lo - cand); lh <- logS(hi - cand); ifelse(ll > lh, exp(ll) - exp(lh), 0) } else exp(logS(hi - cand))
  if (!all(is.finite(w)) || sum(w) <= 0) return(NULL)
  w / sum(w)
}

# ------------------------------------------------------------------------------
# WP4 (C1, C5): P1 with covariates that change over P1's history, and the corrected weights.
# The functions above (p1_candidates, p1_weights) are the WP3 versions with ONE lambda per plot; they are kept as the reference for the identity
# test (R/tests/wp4_C1_check.R) and are no longer used by p1_entries().
#
# Pieces (absolute P1 clock, years; time 0 = stand origin, W_j = first visit):
#   * before the first onset: one piece with the base covariates (X0: first-visit climate/soil, no disturbance) = the one-formula draw of before;
#   * from the first onset: yearly pieces (grid from the first onset; onset <= 0: grid from 0) up to W_j, plus an EXTRA breakpoint at the exact
#     onset of every other type (decided 2026-11-04, deviation from "yearly pieces": the piece splits in two, each with TSD at its own midpoint);
#   * a type is on in a piece when the piece midpoint >= its onset; its TSD is the piece midpoint minus the onset, through tsd_fun(TSD_FORM);
#   * after W_j: the observed intervals up to the end of the weighting bracket, each with its own end-visit row (as in the fitted rows).
# Hazard of exit k in a piece: lambda_k,piece^alpha_k * alpha_k * u^(alpha_k - 1) with u = time since entry (proportional hazards, entry clock),
# so H(e; a, b) = sum over pieces lambda^alpha * (u1^alpha - u0^alpha).
# ------------------------------------------------------------------------------
P1_FLAG_COLS <- list(fire = c(flag = "fire_f", gdd = "gdd_x_fire", tsd = "tsd_fire"),
                     pest = c(flag = "pest_f", gdd = "gdd_x_pest", tsd = "tsd_pest"),
                     harv_light = c(flag = "harv_light", gdd = "gdd_x_harv_light", tsd = "tsd_harv_light"),
                     harv_heavy = c(flag = "harv_heavy", gdd = "gdd_x_harv_heavy", tsd = "tsd_harv_heavy"))

# Pieces of plot j. L0, Liv = lambda matrices of the base rows (X0) and of the observed intervals (Xiv). Returns p0, p1, L (pieces x K), LA = L^alpha
# (0 for transitions that are not ok), npre = number of pieces before W_j (pieces 1..npre end at W_j; the rest are observed intervals).
p1_pieces <- function(j, theta, inp, L0, Liv, f_tsd = tsd_fun()) {
  W <- inp$p1$W_j[j]
  on <- c(fire = inp$on$fire[j], pest = inp$on$pest[j], harv = inp$on$harv[j])
  lv <- c(fire = "fire", pest = "pest", harv = if (isTRUE(inp$on$harv_lvl[j] == 2L)) "harv_heavy" else "harv_light")
  on <- on[!is.na(on)]
  K <- ncol(L0)
  if (W > 0) {
    first_on <- if (length(on)) max(0, min(on)) else NA_real_
    grid <- if (is.na(first_on) || first_on >= W) numeric(0) else seq(first_on, W, by = 1)
    brk <- sort(c(0, grid, unname(on[on > 0 & on < W]))); brk <- brk[brk < W - 1e-9]
    brk <- c(brk[c(TRUE, diff(brk) > 1e-9)], W)         # drop breakpoints closer than 1e-9 y; W itself stays exact
    p0 <- brk[-length(brk)]; p1e <- brk[-1L]
  } else { p0 <- numeric(0); p1e <- numeric(0) }
  npre <- length(p0)
  Lp <- if (npre) matrix(L0[j, ], npre, K, byrow = TRUE) else matrix(0, 0, K)   # W_j = 0: no piece before the first visit
  if (npre && length(on)) {
    mid <- (p0 + p1e) / 2
    Xp <- matrix(inp$X0[j, ], npre, length(inp$terms), byrow = TRUE, dimnames = list(NULL, inp$terms))
    zg <- inp$X0[j, "z_gdd_3"]; any_on <- rep(FALSE, npre)
    for (ty in names(on)) {
      a <- mid >= on[[ty]] - 1e-9; if (!any(a)) next
      cc <- P1_FLAG_COLS[[lv[[ty]]]]
      Xp[a, cc[["flag"]]] <- 1; Xp[a, cc[["gdd"]]] <- zg; Xp[a, cc[["tsd"]]] <- f_tsd(mid[a] - on[[ty]]); any_on <- any_on | a
    }
    if (any(any_on)) Lp[any_on, ] <- lambda_matrix(theta, Xp[any_on, , drop = FALSE])
  }
  r <- inp$iv_first[j] + seq_len(inp$iv_n[j]) - 1L
  r <- r[inp$iv$t_end[r] <= inp$p1$obs_upper[j] + 1e-9]
  L <- rbind(Lp, Liv[r, , drop = FALSE])
  LA <- matrix(0, nrow(L), K); okk <- which(theta$ok)
  LA[, okk] <- L[, okk, drop = FALSE]^rep(theta$alpha[okk], each = nrow(L))
  list(p0 = c(p0, inp$iv$t_start[r]), p1 = c(p1e, inp$iv$t_end[r]), L = L, LA = LA, npre = npre)
}

# Next exit of a sojourn entered at `entry`: for each competing exit k (columns ks of LA = lambda^alpha per piece) the time at which its cumulative
# hazard reaches E[k] (E = standard exponential draws, one per exit), found piece by piece; the earliest wins. Pieces p0/p1e must end at the
# horizon W; returns c(time, k), or c(Inf, 0) if no exit falls before the end of the last piece.
p1_next_exit <- function(p0, p1e, LA, alpha, ks, entry, E) {
  q <- which(p1e > entry); u0 <- pmax(p0[q] - entry, 0); u1 <- p1e[q] - entry
  best <- Inf; bk <- 0L
  for (ii in seq_along(ks)) {
    k <- ks[ii]; al <- alpha[k]; la <- LA[q, k]
    cs <- cumsum(la * (u1^al - u0^al)); h <- which(cs >= E[ii])[1L]
    if (is.na(h)) next
    t_exit <- entry + (u0[h]^al + (E[ii] - (if (h > 1L) cs[h - 1L] else 0)) / la[h])^(1 / al)
    if (t_exit < best) { best <- t_exit; bk <- k }
  }
  c(best, bk)
}

# Candidate entries: paths from regen at time 0 under the piecewise hazards, one candidate per path = entry time of the state occupied at W_j
# (kept only if that state is the plot's first state). Same random numbers as p1_candidates (one runif per competing exit per step); each exit
# time is the exact inverse of its own cumulative hazard, the earliest wins.
# WP5 jump limit (P1_MAX_JUMPS in R/config.R): a path with more jumps stops with an error (a real path has well under 50 over 120 y; more means
# non-finite or absurd hazards that cycle at time 0). ctx = list(id_pe, first_state, labels) only names the plot and the transitions in the message.
p1_jump_stop <- function(ctx, pc, np, lk, entry, state, out_k, ok, W) {
  lab <- if (is.null(ctx$labels)) as.character(lk) else ctx$labels[lk[lk > 0L]]
  ks <- out_k[[state]]; ks <- ks[ok[ks]]; q <- which(pc$p1[seq_len(np)] > entry)[1L]; if (is.na(q)) q <- np
  stop(sprintf("P1 jump limit: a candidate path of plot %s (first state %s, W_j %.2f) exceeded %d jumps (P1_MAX_JUMPS); entry now %.6g y. Last transitions: %s. Hazards lambda_k (per year) of the exits from the current state in the piece [%.4g, %.4g]: %s. Non-finite or absurd hazards (check the theta draw); no redraw, no capping, no fallback.",
               if (is.null(ctx$id_pe)) "?" else ctx$id_pe, if (is.null(ctx$first_state)) "?" else ctx$first_state, W, P1_MAX_JUMPS, entry,
               paste(lab, collapse = " | "), pc$p0[q], pc$p1[q],
               paste(sprintf("%s = %.4g", if (is.null(ctx$labels)) ks else ctx$labels[ks], pc$L[q, ks]), collapse = ", ")), call. = FALSE)
}
p1_candidates_pw <- function(pc, start_idx, s_j, W, alpha, out_k, to_idx, ok, M, ctx = NULL) {
  if (W <= 0) return(structure(if (s_j == start_idx) rep(0, M) else numeric(0), max_jumps = 0L))
  np <- pc$npre; p0 <- pc$p0[seq_len(np)]; p1e <- pc$p1[seq_len(np)]; LA <- pc$LA[seq_len(np), , drop = FALSE]
  cand <- numeric(M); n <- 0L; mx <- 0L      # mx = the largest number of jumps of any path simulated for this plot (attribute max_jumps of the result)
  for (i in seq_len(M)) {
    state <- start_idx; entry <- 0; nj <- 0L; lk <- integer(6L)
    repeat {
      ks <- out_k[[state]]; ks <- ks[ok[ks]]
      if (!length(ks)) break
      nx <- p1_next_exit(p0, p1e, LA, alpha, ks, entry, -log(runif(length(ks))))
      best <- nx[1L]; bk <- nx[2L]
      if (bk == 0L || best >= W) break
      state <- to_idx[bk]; entry <- best
      nj <- nj + 1L; lk <- c(lk[-1L], bk)
      if (nj > P1_MAX_JUMPS) p1_jump_stop(ctx, pc, np, lk, entry, state, out_k, ok, W)
    }
    mx <- max(mx, nj)
    if (state == s_j) { n <- n + 1L; cand[n] <- entry }
  }
  structure(cand[seq_len(n)], max_jumps = mx)
}

log1mexp <- function(h) ifelse(h > 0, ifelse(h < log(2), log(-expm1(-h)), log1p(-exp(-h))), -Inf)   # log(1 - exp(-h)), h >= 0

# Cumulative hazard of all exits ks of the first state between the absolute times a <= b, for each candidate entry e (a, b scalars).
p1_H_between <- function(pc, ks, alpha, cand, a, b) {
  nc <- length(cand); np <- length(pc$p0)
  U0 <- pmax(matrix(pmax(pc$p0, a), nc, np, byrow = TRUE) - cand, 0)
  U1 <- pmax(matrix(pmin(pc$p1, b), nc, np, byrow = TRUE) - cand, U0)
  H <- numeric(nc)
  for (k in ks) H <- H + drop((U1^alpha[k] - U0^alpha[k]) %*% pc$LA[, k])
  H
}

# C5: weight = [S(lo - e) - S(hi - e)] / S(W_j - e) with a change observed in (lo, hi]; S(hi - e) / S(W_j - e) without. In log space:
# S(lo)/S(W) = exp(-H(W, lo)); [S(lo) - S(hi)]/S(W) = exp(-H(W, lo)) * (1 - exp(-H(lo, hi))); S(hi)/S(W) = exp(-H(W, hi)).
p1_weights_pw <- function(cand, s_j, lo, hi, W, pc, alpha, out_k, ok) {
  ks <- out_k[[s_j]]; ks <- ks[ok[ks]]
  lw <- if (!is.na(lo)) -p1_H_between(pc, ks, alpha, cand, W, lo) + log1mexp(p1_H_between(pc, ks, alpha, cand, lo, hi)) else -p1_H_between(pc, ks, alpha, cand, W, hi)
  if (anyNA(lw) || any(lw == Inf) || all(lw == -Inf)) return(NULL)
  w <- exp(lw - max(lw)); w / sum(w)
}

# mode "mean": weighted mean (Stage 2); "draw": weighted draw (MI loop).
# Plots without a usable candidate get W_j/2 and fallback = TRUE.
# WP5 guard A: stops when a hazard of an ok transition is not finite (Inf or NaN) in the base rows or the interval rows. No redraw, no capping, no fallback.
# lambda_matrix() already sets the hazards of not-ok transitions to 0, so every non-finite value here belongs to an ok transition.
p1_check_hazards <- function(theta, inp, L0, Liv, m = NA_integer_) {
  b0 <- !is.finite(L0); bv <- !is.finite(Liv)
  if (!any(b0) && !any(bv)) return(invisible(TRUE))
  ks <- which(colSums(b0) > 0 | colSums(bv) > 0)
  raw <- if (!is.null(theta$beta_v0)) theta$beta_v0 else theta$beta
  det <- vapply(ks, function(k) sprintf("%s (alpha %.3g, max |beta| drawn %.4g, converted %.4g)", inp$td$label[k], theta$alpha[k],
                                        max(abs(raw[k, ])), max(abs(theta$beta[k, ]))), "")
  stop(sprintf("P1 guard A (%s): non-finite hazard in %d transition(s) with an ok fit: %s. Affected: %d of %d base rows (X0, one per plot), %d of %d interval rows (Xiv); %d Inf and %d NaN entries. Largest |beta| of these transitions: drawn %.4g, converted %.4g. No redraw, no capping, no fallback.",
               if (is.na(m)) "imputation not given (Stage 2 or direct call)" else sprintf("imputation m = %d", m), length(ks),
               paste(head(det, 10), collapse = "; ") , sum(rowSums(b0) > 0), nrow(L0), sum(rowSums(bv) > 0), nrow(Liv),
               sum(is.infinite(L0)) + sum(is.infinite(Liv)), sum(is.nan(L0)) + sum(is.nan(Liv)),
               max(abs(raw[ks, ])), max(abs(theta$beta[ks, ]))), call. = FALSE)
}

p1_entries <- function(theta, inp, M = SIM_M, mode = c("mean", "draw"), n_cores = 1L, m = NA_integer_) {
  mode <- match.arg(mode)
  L0 <- lambda_matrix(theta, inp$X0); Liv <- lambda_matrix(theta, inp$Xiv); f_tsd <- tsd_fun()
  p1_check_hazards(theta, inp, L0, Liv, m)
  W <- inp$p1$W_j; sj <- inp$p1$s_idx; lo <- inp$p1$obs_lower; hi <- inp$p1$obs_upper
  to_idx <- inp$td$to_idx
  one <- function(j) {
    pc <- p1_pieces(j, theta, inp, L0, Liv, f_tsd)
    cand <- p1_candidates_pw(pc, inp$start_idx, sj[j], W[j], theta$alpha, inp$out_k, to_idx, theta$ok, M,
                             ctx = list(id_pe = inp$p1$id_pe[j], first_state = inp$p1$first_state[j], labels = inp$td$label))
    mj <- attr(cand, "max_jumps"); if (is.null(mj)) mj <- NA_integer_
    if (!length(cand)) return(c(NA_real_, 0, mj))
    w <- p1_weights_pw(cand, sj[j], lo[j], hi[j], W[j], pc, theta$alpha, inp$out_k, theta$ok)
    if (is.null(w)) return(c(NA_real_, length(cand), mj))
    c(if (mode == "mean") sum(cand * w) else cand[sample.int(length(cand), 1L, prob = w)], length(cand), mj)
  }
  res <- parallel::mclapply(seq_along(W), one, mc.cores = n_cores)
  bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r), logical(1))
  if (any(bad)) stop("p1_entries: ", sum(bad), " plots failed in the workers: ", as.character(res[bad][[1]]))
  res <- do.call(rbind, res)
  fb <- is.na(res[, 1])
  entry <- ifelse(fb, W / 2, res[, 1])
  data.table(id_pe = inp$p1$id_pe, entry = entry, n_cand = res[, 2], fallback = fb, max_jumps = as.integer(res[, 3]))
}

# largest number of jumps of any P1 path and its plot (printed per imputation and for Stage 2; mi_diag columns p1_max_jumps, p1_max_jumps_plot)
p1_jump_summary <- function(p1) { i <- which.max(p1$max_jumps); list(max = p1$max_jumps[i], id_pe = p1$id_pe[i]) }

# P1 fallbacks (no usable candidate) by first state: n_plots, n_fallback, share (printed per imputation; compare with the main run's 124-196 in total)
p1_fallback_by_state <- function(p1, inp, m = NA_integer_) {
  d <- data.table(first_state = inp$p1$first_state, fallback = p1$fallback)
  r <- d[, .(n_plots = .N, n_fallback = sum(fallback)), by = first_state][order(-n_fallback, first_state)]
  r[, share := n_fallback / n_plots]
  if (!is.na(m)) r[, m := m]
  r[]
}

# ------------------------------------------------------------------------------
# P2 (moment_match_weibull_naive copied unchanged from p2_imputation_strategies.R; WP5: rule C2, draw C3)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# WP5: C2 (exit at the heavy disturbance for X->regen) and C3 (drawn exit). Plan: docs/PIPELINE_FIX_PLAN_2026-11-03.md, C2 / C3.
# All P2 functions below have the build_rows() signature p2_fun(a_dur, b_dur, from, to, r); r = row of inp$iv (the event interval).
# exit_source (one value per event): "rule" (C2), "draw" (closed-form C3 draw), "uniform" (no usable fit, uniform draw on (a, b)),
# "uniform_nonfinite" (usable fit but the closed form is not finite: uniform draw on (a, b), user decision 2026-11-04), "moment" (Stage 2 point), "moment_fallback" (Stage 2: the moment match fell back to the midpoint),
# "midpoint" (Stage 1; Stage 2 where the fit is unusable). `diag` is an environment: n_event, n_midpoint (kept names), n_rule_at_end, src[r].
# ------------------------------------------------------------------------------
P2_SOURCES <- c("rule", "draw", "uniform", "uniform_nonfinite", "moment", "moment_fallback", "midpoint")

p2_diag_init <- function(inp, diag) {
  diag$n_event <- 0L; diag$n_midpoint <- 0L; diag$n_rule_at_end <- 0L
  diag$src <- rep(NA_character_, nrow(inp$iv)); diag
}
p2_check_rule_inputs <- function(inp) {
  if (isTRUE(XREGEN_EXIT_RULE) && !all(c("D_heavy", "year_start", "t_start", "t_end") %in% names(inp$iv)))
    stop("C2: inp$iv has no D_heavy / year_start column (rerun R/04_prepare_inla.R, WP2 A9)")
}
# C2: exit time (duration since entry) of an X->regen event whose interval has a heavy-disturbance date, else NA.
#   T = t_start[r] + (D_heavy - year_start[r]) - E, with E = t_end[r] - b_dur (b_dur = t_end - E); clamped to [a_dur, b_dur].
#   A clamp that moves T by more than 1e-9 stops (it would mean A9 or the entry logic is broken). Applies to every theta (also failed fits, non-kept transitions).
p2_rule_T <- function(inp, diag, a_dur, b_dur, to, r) {
  if (!isTRUE(XREGEN_EXIT_RULE) || to != inp$states[inp$start_idx]) return(NA_real_)
  D <- inp$iv$D_heavy[r]; if (is.na(D)) return(NA_real_)
  E <- inp$iv$t_end[r] - b_dur
  Tr <- inp$iv$t_start[r] + (D - inp$iv$year_start[r]) - E
  Tc <- min(max(Tr, a_dur), b_dur)
  if (abs(Tc - Tr) > 1e-9) stop(sprintf("C2: the clamp moves the rule exit of plot %s interval %s by %.3g y (rule T %.6f, bracket [%.6f, %.6f])",
                                        inp$iv$id_pe[r], inp$iv$interval[r], Tc - Tr, Tr, a_dur, b_dur))
  if (abs(Tc - b_dur) < 1e-9) diag$n_rule_at_end <- diag$n_rule_at_end + 1L
  Tc
}
p2_note <- function(diag, r, src) { diag$n_event <- diag$n_event + 1L; diag$src[r] <- src; if (src == "midpoint") diag$n_midpoint <- diag$n_midpoint + 1L }

# Stage 1: rule date where the interval has D_heavy, midpoint otherwise (the rule needs no fitted parameters and Stage 1 seeds Stage 2).
make_p2_stage1 <- function(inp, diag = new.env()) {
  p2_diag_init(inp, diag); p2_check_rule_inputs(inp)
  function(a_dur, b_dur, from, to, r) {
    Tr <- p2_rule_T(inp, diag, a_dur, b_dur, to, r)
    if (!is.na(Tr)) { p2_note(diag, r, "rule"); return(Tr) }
    p2_note(diag, r, "midpoint"); (a_dur + b_dur) / 2
  }
}

# moment_match_weibull_info(): the moment match (copied from p2_imputation_strategies.R) returning list(T, fallback); fallback = TRUE where the
# midpoint was used (interval probability < 1e-10 or not finite, or the integral failed). moment_match_weibull_naive() = its T (values unchanged).
moment_match_weibull_info <- function(a, b, alpha, lambda) {
  S_a <- exp(-(lambda * a)^alpha)
  S_b <- exp(-(lambda * b)^alpha)
  prob_interval <- S_a - S_b
  if (!is.finite(prob_interval) || prob_interval < 1e-10) return(list(T = (a + b) / 2, fallback = TRUE))
  integrand <- function(t) (t^alpha) * alpha * lambda^alpha * t^(alpha - 1) * exp(-(lambda * t)^alpha) / prob_interval
  e_t_alpha <- tryCatch(integrate(integrand, lower = a, upper = b, rel.tol = 1e-6)$value, error = function(e) NA_real_)
  if (!is.finite(e_t_alpha) || e_t_alpha <= 0) return(list(T = (a + b) / 2, fallback = TRUE))
  list(T = max(a, min(b, e_t_alpha^(1 / alpha))), fallback = FALSE)
}
moment_match_weibull_naive <- function(a, b, alpha, lambda) moment_match_weibull_info(a, b, alpha, lambda)$T

# Stage 2 (the point for the seed models, C3 point 1): rule first (C2), then the moment-matched point for kept transitions with an ok fit,
# midpoint otherwise. A moment match that falls back to the midpoint inside the match is labelled "moment_fallback" (WP5 follow-up).
make_p2_moment <- function(theta, inp, diag = new.env()) {
  p2_diag_init(inp, diag); p2_check_rule_inputs(inp)
  function(a_dur, b_dur, from, to, r) {
    Tr <- p2_rule_T(inp, diag, a_dur, b_dur, to, r)
    if (!is.na(Tr)) { p2_note(diag, r, "rule"); return(Tr) }
    k <- inp$kmat[from, to]
    if (is.na(k) || !theta$ok[k]) { p2_note(diag, r, "midpoint"); return((a_dur + b_dur) / 2) }
    lam <- exp(theta$intercept[k] + sum(theta$beta[k, ] * inp$Xiv[r, ]))
    mm <- moment_match_weibull_info(a_dur, b_dur, theta$alpha[k], lam)
    p2_note(diag, r, if (mm$fallback) "moment_fallback" else "moment")
    mm$T
  }
}

# C3: closed-form draw from the Weibull conditional on (a, b], H(t) = (lam t)^alpha:  T = (H(a) - log1p(-U (-expm1(-D))))^(1/alpha) / lam, D = H(b) - H(a).
# Cumulative-hazard form (exact when S(a) underflows). Returns NA when not finite (H(a), H(b), D <= 0, T); otherwise T clamped to [a, b].
p2_draw_closed <- function(U, a, b, lam, alpha) {
  if (!is.finite(lam) || lam <= 0 || !is.finite(alpha) || alpha <= 0) return(NA_real_)
  Ha <- (lam * a)^alpha; Hb <- (lam * b)^alpha; D <- Hb - Ha
  if (!is.finite(Ha) || !is.finite(Hb) || !is.finite(D) || !(D > 0)) return(NA_real_)
  x <- Ha - log1p(-U * (-expm1(-D)))
  Tt <- x^(1 / alpha) / lam
  if (!is.finite(Tt)) return(NA_real_)
  min(max(Tt, a), b)
}

# Imputations (run_one_imputation): rule first (C2); else the closed-form draw for the realised transition; uniform draw on (a, b) where there is
# no usable fit ("uniform") or the closed form is not finite ("uniform_nonfinite"). Exactly ONE runif(1) per event, drawn first, also for rule events
# (the stream of the other events is then identical with XREGEN_EXIT_RULE TRUE and FALSE, until the first exit differs).
make_p2_draw <- function(theta, inp, diag = new.env()) {
  p2_diag_init(inp, diag); p2_check_rule_inputs(inp)
  function(a_dur, b_dur, from, to, r) {
    U <- runif(1)
    Tr <- p2_rule_T(inp, diag, a_dur, b_dur, to, r)
    if (!is.na(Tr)) { p2_note(diag, r, "rule"); return(Tr) }
    k <- inp$kmat[from, to]
    if (is.na(k) || !theta$ok[k]) { p2_note(diag, r, "uniform"); return(a_dur + U * (b_dur - a_dur)) }
    lam <- exp(theta$intercept[k] + sum(theta$beta[k, ] * inp$Xiv[r, ]))
    Td <- p2_draw_closed(U, a_dur, b_dur, lam, theta$alpha[k])
    if (is.na(Td)) { p2_note(diag, r, "uniform_nonfinite"); return(a_dur + U * (b_dur - a_dur)) }
    p2_note(diag, r, "draw"); Td
  }
}

# per-event table (id_pe, interval, from, to, exit_abs, floored, exit_source) of a build_rows() result and its diag environment
p2_event_table <- function(b, diag) {
  t <- b$segs[type == "event", .(id_pe, interval, iv_row, from, to = event_to, exit_abs = end_abs, floored)]
  t[, exit_source := diag$src[iv_row]]; t[, iv_row := NULL]
  stopifnot("exit_source missing for an event" = !anyNA(t$exit_source), all(t$exit_source %in% P2_SOURCES))
  t[]
}
# counts per transition and source (long); m = imputation number (NA for the stages)
p2_source_counts <- function(p2, m = NA_integer_) {
  r <- p2[, .(n = .N), by = .(transition = paste(from, to, sep = "->"), exit_source)][order(transition, exit_source)]
  r[, m := m]; setcolorder(r, c("m", "transition", "exit_source", "n")); r[]
}
# one printable line: exits per source; uniform_nonfinite is listed per transition
p2_source_line <- function(src) {
  tot <- src[, .(n = sum(n)), by = exit_source]; cnt <- function(s) { x <- tot$n[tot$exit_source == s]; if (length(x)) x else 0L }
  s <- paste(sprintf("%s %d", P2_SOURCES, vapply(P2_SOURCES, cnt, 1L)), collapse = ", ")
  un <- src[exit_source == "uniform_nonfinite"]
  if (nrow(un)) s <- paste0(s, " | uniform_nonfinite by transition: ", paste(sprintf("%s %d", un$transition, un$n), collapse = "; "))
  s
}
# diag columns of one imputation
p2_diag_cols <- function(src, diag) {
  cnt <- function(s) { x <- sum(src$n[src$exit_source == s]); as.integer(x) }
  list(p2_rule = cnt("rule"), p2_draw = cnt("draw"), p2_uniform = cnt("uniform"), p2_uniform_nonfinite = cnt("uniform_nonfinite"), p2_rule_at_end = diag$n_rule_at_end)
}

# ------------------------------------------------------------------------------
# fitting
# ------------------------------------------------------------------------------
# only_k (2026-10-22): fit only these transitions (indices into inp$td). Transitions are fitted independently (fit_inla_models is an lapply),
# so a transition's fit is the same call whether or not the others are fitted. Default = all transitions (unchanged behaviour).
fit_rows <- function(rows, inp, fit_fun, only_k = NULL) {
  ks <- if (is.null(only_k)) seq_len(nrow(inp$td)) else sort(unique(as.integer(only_k)))
  d <- rows[kept == TRUE]
  idx <- match(paste(d$from, d$to, sep = "->"), inp$td$label)
  stopifnot(!anyNA(idx))
  keep <- idx %in% ks; d <- d[keep]; idx <- idx[keep]
  cols <- c("from", "to", "status", "Tstop", "truncation", inp$terms)
  lst <- split(as.data.frame(d[, ..cols]), factor(idx, levels = ks))
  n_events <- vapply(lst, function(x) sum(x$status == 1L), integer(1))
  models <- fit_fun(unname(lst), rhs_formula = unname(inp$rhs[ks]), exact = TRUE)
  list(models = models, n_events = n_events, data = lst, ks = ks)
}

# Per-fit diagnostics kept in the output: error text, INLA mode status, log-marginal-likelihood,
# and any warning / failure lines of INLA's own log (collapsed). Small, saved for every fit.
fit_diagnostics <- function(m) {
  if (!is.null(m$error)) return(list(error_msg = as.character(m$error)[1], mode_status = NA_integer_, mlik = NA_real_, log_warn = NA_character_))
  ms <- tryCatch(as.integer(m$mode$mode.status), error = function(e) NA_integer_)
  ml <- tryCatch(as.numeric(m$mlik[1, 1]), error = function(e) NA_real_)
  lw <- tryCatch({ l <- m$logfile; l <- l[grepl("warn|fail|converge|not positive|singular|\\b(nan|inf)\\b", l, ignore.case = TRUE) & !grepl("(Relative|Absolute|Monte Carlo) error", l, ignore.case = TRUE)]
                   if (length(l)) paste(substr(trimws(utils::head(unique(l), 3)), 1, 120), collapse = " | ") else NA_character_ }, error = function(e) NA_character_)
  list(error_msg = NA_character_, mode_status = if (length(ms)) ms[1] else NA_integer_, mlik = ml, log_warn = lw)
}

# Compact extras per fit (kept for every fit of every imputation): full posterior summaries,
# marginals, evidence, optimiser trace, warnings and the JOINT posterior covariance of
# (Intercept, betas, alpha) from n_cov posterior samples. Intercept is on INLA's decade scale (variant 0);
# the years-scale intercept b0 - alpha*log(10) is computed in 08 from this covariance (B3).
fit_extras <- function(m, tms, sampler, n_cov = N_COV_SAMPLES, seed = 1L) {
  if (!is.null(m$error)) return(list(error = as.character(m$error)[1]))
  g <- function(expr) tryCatch(expr, error = function(e) NULL)
  ex <- list(summary_fixed = m$summary.fixed, summary_hyper = m$summary.hyperpar,
             internal_hyper = g(m$internal.summary.hyperpar), marg_fixed = g(m$marginals.fixed),
             marg_alpha = g(m$marginals.hyperpar[[ALPHA_NAME]]), mlik = g(m$mlik),
             mode_theta = g(m$mode$theta), mode_status = g(m$mode$mode.status), log_post_mode = g(m$mode$log.posterior.mode),
             opt_trace = g(m$misc$opt.trace), nfunc = g(m$misc$nfunc), warnings = g(m$misc$warnings),
             cpu = g(m$cpu.used), n_rows = g(nrow(m$summary.linear.predictor)), logfile_tail = g(utils::tail(m$logfile, 40)))
  if (n_cov > 0) ex$cov <- tryCatch({
    smp <- sampler(n_cov, m, seed = seed)
    S <- t(vapply(smp, function(x) { ln <- rownames(x$latent)
      c(x$latent[grepl("^\\(Intercept\\)", ln), 1][1], vapply(tms, function(tm) x$latent[grepl(paste0("^", tm, "($|:)"), ln), 1][1], numeric(1)),
        x$hyperpar[[ALPHA_NAME]]) }, numeric(length(tms) + 2L)))
    colnames(S) <- c("(Intercept)", tms, "alpha")
    list(mean = colMeans(S), cov = stats::cov(S), n = nrow(S))
  }, error = function(e) list(error = conditionMessage(e)))
  ex
}

tidy_fit <- function(models, inp, ks = seq_len(nrow(inp$td))) {   # models[[j]] belongs to transition ks[j]
  rbindlist(lapply(seq_along(ks), function(j) {
    k <- ks[j]; m <- models[[j]]; tms <- colnames(inp$mask)[inp$mask[k, ]]; pars <- c("alpha", "(Intercept)", tms)
    dgn <- fit_diagnostics(m)
    if (!is.null(m$error)) return(data.table(k = k, parameter = pars, est = NA_real_, sd = NA_real_, fit_error = TRUE,
                                             error_msg = dgn$error_msg, mode_status = dgn$mode_status, mlik = dgn$mlik, log_warn = dgn$log_warn))
    fx <- m$summary.fixed; hy <- m$summary.hyperpar
    g <- function(nm, col) if (nm %in% rownames(fx)) fx[nm, col] else NA_real_
    data.table(k = k, parameter = pars,
               est = c(hy[ALPHA_NAME, "mean"], g("(Intercept)", "mean"), vapply(tms, g, numeric(1), col = "mean")),   # intercept: INLA decade scale, variant 0 (B3)
               sd  = c(hy[ALPHA_NAME, "sd"],   g("(Intercept)", "sd"),                       vapply(tms, g, numeric(1), col = "sd")),
               fit_error = FALSE, error_msg = dgn$error_msg, mode_status = dgn$mode_status, mlik = dgn$mlik, log_warn = dgn$log_warn)
  }))[, label := inp$td$label[k]][]
}

# ------------------------------------------------------------------------------
# B6: cap diagnostic. fit_map_v0() = penalised maximum likelihood with the SAME variant-0 likelihood as the INLA fit (left-truncated
# Weibull, N(0, COEF_PRIOR_SD^2) penalty on the coefficients, flat intercept, NO truncation cap); H_at() = cumulative hazard at the
# truncation time of every row. Moved unchanged from tests/22_cap_test_new_rows.R (SIGMA_B -> COEF_PRIOR_SD).
# Diagnostics only, not a validity gate (user decision 2026-10-22).
# ------------------------------------------------------------------------------
fit_map_v0 <- function(d, terms) {      # h = alpha exp(eta) t^(alpha-1), H = exp(eta) t^alpha
  X <- cbind(1, as.matrix(d[, ..terms])); X[is.na(X)] <- 0; Tt <- d$Tstop; tr <- d$truncation; s <- d$status; lT <- log(Tt)
  nll <- function(th) { a <- exp(th[1]); eta <- drop(X %*% th[-1]); -(sum(s * (th[1] + eta + (a - 1) * lT)) - sum(exp(eta) * (Tt^a - tr^a))) + sum(th[-(1:2)]^2) / (2 * COEF_PRIOR_SD^2) }
  o <- optim(c(0, log(max(sum(s), 0.5) / sum(Tt - tr)), rep(0, length(terms))), nll, method = "BFGS", control = list(maxit = 10000, reltol = 1e-12))
  V <- tryCatch(solve(optimHess(o$par, nll)), error = function(e) NULL); se <- if (is.null(V)) rep(NA_real_, length(o$par)) else sqrt(pmax(diag(V), 0))
  a <- exp(o$par[1])
  list(tab = data.table(parameter = c("alpha", "(Intercept)", terms), est = c(a, o$par[-1]), sd = c(a * se[1], se[-1])), conv = o$convergence)
}
H_at <- function(d, terms, tab) {       # cumulative hazard at the truncation time of every row, at the given estimates
  est <- setNames(tab$est, tab$parameter); if (anyNA(est)) return(c(maxH = NA_real_, n_above = NA_real_, n_above_ev = NA_real_))
  X <- as.matrix(d[, ..terms]); X[is.na(X)] <- 0
  H <- exp(est[["(Intercept)"]] + drop(X %*% est[terms])) * d$truncation^est[["alpha"]]
  c(maxH = max(H), n_above = sum(H > CAP), n_above_ev = sum(H > CAP & d$status == 1L))
}
CAP_COLS <- c("maxH_inla", "n_above_inla", "maxH_map", "n_above_map", "z_map_max", "z_map_alpha", "z_map_intercept")
# One fit: d = its fitting data (data.frame), tms = its terms, tab = tidy_fit rows of the fit (parameter, est, sd; INLA means). Returns the
# CAP_COLS values. z_map_max = largest |INLA - MAP| / MAP SE over the term coefficients (not alpha, not the intercept); z_map_alpha and z_map_intercept the same for alpha and for the intercept.
# NA for a failed INLA fit or a failed MAP.
cap_diagnostics <- function(d, tms, tab) {
  na <- setNames(rep(NA_real_, length(CAP_COLS)), CAP_COLS)
  if (isTRUE(any(tab$fit_error)) || anyNA(tab$est)) return(na)
  d <- as.data.table(d)
  fm <- tryCatch(fit_map_v0(d, tms), error = function(e) NULL); if (is.null(fm)) return(na)
  hi <- H_at(d, tms, tab); hm <- H_at(d, tms, fm$tab)
  z <- merge(tab[, .(parameter, inla = est)], fm$tab[, .(parameter, map = est, se = sd)], by = "parameter")[, z := abs(inla - map) / se]
  zmax <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  c(maxH_inla = hi[["maxH"]], n_above_inla = hi[["n_above"]], maxH_map = hm[["maxH"]], n_above_map = hm[["n_above"]],
    z_map_max = zmax(z[!parameter %in% c("alpha", "(Intercept)"), z]), z_map_alpha = zmax(z[parameter == "alpha", z]),
    z_map_intercept = zmax(z[parameter == "(Intercept)", z]))
}

# Adds the CAP_COLS columns to tidy_fit() output: one value per fit (k), repeated over its parameter rows. data_list[[j]] = fitting data of transition ks[j].
add_cap_columns <- function(res, data_list, ks, inp) {
  cap <- rbindlist(lapply(seq_along(ks), function(j) { kk <- ks[j]
    as.data.table(as.list(c(k = kk, cap_diagnostics(data_list[[j]], colnames(inp$mask)[inp$mask[kk, ]], res[res$k == kk])))) }))
  merge(res, cap, by = "k")
}

# ------------------------------------------------------------------------------
# stages
# ------------------------------------------------------------------------------
# WP5: the Stage 1 rows, built in ONE place (run_stage1 and R/tests/wp5_recount_stage1.R both call it): entry W_j/2, C2 rule date for X->regen
# intervals with D_heavy, midpoint otherwise.
stage1_rows <- function(inp) {
  fe <- setNames(inp$p1$W_j / 2, inp$p1$id_pe)
  dg <- new.env()
  b <- build_rows(inp$iv, fe, inp$states, inp$tr, p2_fun = make_p2_stage1(inp, dg), keep_cols = inp$keep_cols)
  list(b = b, dg = dg, p2 = p2_event_table(b, dg))
}

run_stage1 <- function(inp, fit_fun) {
  s <- stage1_rows(inp); b <- s$b
  # B7: recount events on the final rows before the first fit (stops on a mismatch). Skipped for a subsample: inp$tr holds full-data counts.
  if (isTRUE(inp$subsampled)) cat("recount_events(): skipped (SUBSAMPLE_N is set; transition_set counts are full-data counts)\n") else {
    rc <- recount_events(b$rows, inp$tr); cat(sprintf("recount_events(): OK, %d kept transitions, counts equal transition_set, all >= %d events\n", nrow(rc), MIN_EVENTS)) }
  f <- fit_rows(b$rows, inp, fit_fun)
  list(theta = extract_theta(f$models, inp), tidy = tidy_fit(f$models, inp), n_events = f$n_events, n_floored = b$n_floored,
       p2 = s$p2, p2_src = p2_source_counts(s$p2), p2_rule_at_end = s$dg$n_rule_at_end)
}

run_stage2 <- function(inp, theta1, fit_fun, M = SIM_M, n_cores = 1L) {
  p1 <- p1_entries(theta1, inp, M, "mean", n_cores)
  dg <- new.env()
  b <- build_rows(inp$iv, setNames(p1$entry, p1$id_pe), inp$states, inp$tr,
                  p2_fun = make_p2_moment(theta1, inp, dg), keep_cols = inp$keep_cols)
  f <- fit_rows(b$rows, inp, fit_fun)
  p2 <- p2_event_table(b, dg)
  list(models = f$models, theta = extract_theta(f$models, inp), tidy = tidy_fit(f$models, inp), n_events = f$n_events,
       p1 = p1, p2_event = dg$n_event, p2_midpoint = dg$n_midpoint, n_floored = b$n_floored,
       p2 = p2, p2_src = p2_source_counts(p2), p2_rule_at_end = dg$n_rule_at_end)
}

run_one_imputation <- function(m, inp, models_seed, theta_seed, fit_fun, M = SIM_M, sampler = default_sampler, keep_models = FALSE, only_k = NULL) {
  # only_k: fit only these transitions. The imputation itself (theta draw for ALL transitions, P1, P2) is unchanged, so imputation m is
  # identical whatever only_k is; only the set of fitted transitions differs.
  th <- draw_theta_m(models_seed, inp, theta_seed, sampler)
  p1 <- p1_entries(th, inp, M, "draw", 1L, m = m)
  dg <- new.env()
  b <- build_rows(inp$iv, setNames(p1$entry, p1$id_pe), inp$states, inp$tr,
                  p2_fun = make_p2_draw(th, inp, dg), keep_cols = inp$keep_cols)
  f <- fit_rows(b$rows, inp, fit_fun, only_k)
  res <- tidy_fit(f$models, inp, f$ks); res[, m := m]
  ev <- data.table(k = f$ks, n_events = as.integer(f$n_events))
  res <- merge(res, ev, by = "k")
  # B6: per-fit cap diagnostic (one value per fit, repeated over its parameter rows like mlik); fixed column order for fwrite(append)
  res <- add_cap_columns(res, f$data, f$ks, inp)
  p2 <- p2_event_table(b, dg); p2_src <- p2_source_counts(p2, m)
  out <- list(p1_state = p1_fallback_by_state(p1, inp, m), p2_src = p2_src, res = res, diag = data.table(m = m, p1_fallback = sum(p1$fallback), p1_no_candidate = sum(p1$n_cand == 0), p1_max_jumps = p1_jump_summary(p1)$max, p1_max_jumps_plot = p1_jump_summary(p1)$id_pe,
                                    p2_events = dg$n_event, p2_midpoint = dg$n_midpoint, as.data.table(p2_diag_cols(p2_src, dg)),
                                    n_theta_fallback = th$n_fallback, n_fit_error = sum(vapply(f$models, function(x) !is.null(x$error), logical(1))),
                                    n_floored = b$n_floored, n_fitted = length(f$ks)))
  tms_k <- lapply(seq_len(nrow(inp$td)), function(k) colnames(inp$mask)[inp$mask[k, ]])
  out$extras <- list(
    theta_m = data.table(k = seq_len(nrow(inp$td)), label = inp$td$label, ok = th$ok, alpha = th$alpha, intercept_dec = th$intercept_dec,
                         as.data.table(th$beta_v0)),   # B4: the variant-0 draws (decade intercept, raw beta, alpha)
    p1 = p1[, .(id_pe, entry, n_cand, fallback)],
    p2 = p2,      # id_pe, interval, from, to, exit_abs, floored, exit_source (WP5, C2/C3)
    fits = setNames(lapply(seq_along(f$ks), function(j) fit_extras(f$models[[j]], tms_k[[f$ks[j]]], sampler, N_COV_SAMPLES, sample.int(1e6, 1L))), inp$td$label[f$ks]))
  # fitting data of failed or near-degenerate fits (alpha SD < SD_MIN), to rerun them with verbose = TRUE
  al <- res[parameter == "alpha"]
  bad_k <- al[fit_error == TRUE | (!is.na(sd) & sd < SD_MIN), k]
  if (length(bad_k)) out$bad_data <- setNames(lapply(bad_k, function(k) list(data = f$data[[match(k, f$ks)]], rhs = inp$rhs[[k]], label = inp$td$label[k])), inp$td$label[bad_k])
  if (keep_models) out$models <- f$models
  out
}

# ------------------------------------------------------------------------------
# Rubin pooling (pool_rubin copied from fit_p1_rung4_mi.R)
# ------------------------------------------------------------------------------
pool_rubin <- function(est, var) {
  ok <- !is.na(est) & !is.na(var); m_eff <- sum(ok)
  if (m_eff < 2) return(data.frame(pooled_est = if (m_eff == 1) est[ok] else NA_real_, W_bar = if (m_eff == 1) var[ok] else NA_real_,
                                   B = NA_real_, total_var = NA_real_, m_eff = m_eff))
  W_bar <- mean(var[ok]); B <- var(est[ok])
  data.frame(pooled_est = mean(est[ok]), W_bar = W_bar, B = B, total_var = W_bar + (1 + 1 / m_eff) * B, m_eff = m_eff)
}

# Fits flagged as computational failures: INLA error, or alpha posterior SD below sd_min (near-degenerate).
# A flagged fit is flagged for ALL its parameters in that imputation.
flag_fits <- function(raw, sd_min = SD_MIN) {
  fl <- raw[parameter == "alpha", .(m, k, label, flagged = fit_error == TRUE | is.na(sd) | sd < sd_min, alpha_est = est, alpha_sd = sd, error_msg)]
  fl[]
}

# Regular Rubin pooling (columns pooled_est, W_bar, B, total_var, m_eff, se, n_events_avg as before) plus
#  r (relative increase in variance), fmi, df_br (Barnard-Rubin df; complete-data df = events - parameters) and t-based 95% interval.
#  For alpha also pooling on the log scale (delta method, sd/est): log_est, log_se, est_from_log, ci_lo_log, ci_hi_log (back-transformed).
# exclude_flagged = TRUE drops flagged fits first (sensitivity); the default pools every fit that ran.
pool_all <- function(raw, exclude_flagged = FALSE, sd_min = SD_MIN) {
  d <- raw[!is.na(est) & !is.na(sd) & fit_error == FALSE]
  if (exclude_flagged) { fl <- flag_fits(raw, sd_min); d <- merge(d, fl[flagged == FALSE, .(m, k)], by = c("m", "k")) }
  npar <- raw[, .(npar = uniqueN(parameter)), by = k]
  out <- d[, { p <- pool_rubin(est, sd^2); c(as.list(p), list(se = sqrt(p$total_var), n_events_avg = mean(n_events))) }, by = .(k, label, parameter)][order(k)]
  out <- merge(out, npar, by = "k", sort = FALSE)
  out[, r := ifelse(W_bar > 0, (1 + 1 / m_eff) * B / W_bar, NA_real_)]
  out[, lam := r / (1 + r)]
  nu_old <- (out$m_eff - 1) / out$lam^2
  nu_com <- pmax(out$n_events_avg - out$npar, 1)
  nu_obs <- (nu_com + 1) / (nu_com + 3) * nu_com * (1 - out$lam)
  out[, df_br := ifelse(is.finite(nu_old) & !is.na(nu_old), 1 / (1 / nu_old + 1 / nu_obs), nu_obs)]
  out[, fmi := (r + 2 / (df_br + 3)) / (r + 1)]
  out[, `:=`(ci_lo = pooled_est - qt(0.975, df_br) * se, ci_hi = pooled_est + qt(0.975, df_br) * se)]
  al <- d[parameter == "alpha" & est > 0, { p <- pool_rubin(log(est), (sd / est)^2); list(log_est = p$pooled_est, log_se = sqrt(p$total_var), log_B = p$B, log_W = p$W_bar, log_m = p$m_eff) }, by = k]
  out <- merge(out, al, by = "k", all.x = TRUE, sort = FALSE)
  na <- out$parameter != "alpha"; for (v in c("log_est", "log_se", "log_B", "log_W", "log_m")) out[na, (v) := NA_real_]
  out[, `:=`(est_from_log = exp(log_est), ci_lo_log = exp(log_est - qt(0.975, df_br) * log_se), ci_hi_log = exp(log_est + qt(0.975, df_br) * log_se))]
  out[, c("npar", "lam") := NULL]
  setorder(out, k)
  out[]
}
