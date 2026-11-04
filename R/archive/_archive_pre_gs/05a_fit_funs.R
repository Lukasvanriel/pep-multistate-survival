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
#   * transition with a failed fit (theta NA) -> hazard 0 in P1, midpoint in P2, counted
# Units: P1 and P2 work in YEARS; the INLA intercept is on the decade scale and is
# corrected by -LOG_TIME_SCALE exactly as in the simulation. alpha and betas need no correction.
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
  X <- as.matrix(f1[, ..terms])

  list(iv = iv, tr = tr, td = td, states = states, mask = mask, terms = terms, rhs = rhs_formulas(mask),
       kmat = kmat, out_k = out_k, p1 = p1, X = X, Xiv = as.matrix(iv[, ..terms]),
       start_idx = match("regen", states), scaling = scaling,
       keep_cols = c("fire", "harv", "pest", "year_end", terms))
}

# ------------------------------------------------------------------------------
# theta: list(alpha, intercept [years], beta [K x P, 0 = term absent], ok [K])
# ------------------------------------------------------------------------------
extract_theta <- function(models, inp, stat = "mean") {
  K <- nrow(inp$td); P <- length(inp$terms)
  alpha <- rep(NA_real_, K); icpt <- rep(NA_real_, K)
  beta <- matrix(0, K, P, dimnames = list(inp$td$label, inp$terms))
  for (k in seq_len(K)) {
    m <- models[[k]]; if (!is.null(m$error)) next
    a <- tryCatch(m$summary.hyperpar[ALPHA_NAME, stat], error = function(e) NA_real_)
    ic <- tryCatch(m$summary.fixed["(Intercept)", stat], error = function(e) NA_real_)
    alpha[k] <- if (length(a) == 1) a else NA_real_
    icpt[k]  <- if (length(ic) == 1) ic - LOG_TIME_SCALE else NA_real_
    for (tm in colnames(inp$mask)[inp$mask[k, ]]) {
      b <- tryCatch(m$summary.fixed[tm, stat], error = function(e) NA_real_)
      beta[k, tm] <- if (length(b) == 1) b else NA_real_
    }
  }
  ok <- !is.na(alpha) & !is.na(icpt) & !apply(is.na(beta), 1, any)
  beta[!ok, ] <- 0
  list(alpha = alpha, intercept = icpt, beta = beta, ok = ok)
}

# seed > 0 and a single thread make INLA's sampler reproducible (seed = 0 means a clock seed)
default_sampler <- function(n, m, seed = 0L) INLA::inla.posterior.sample(n, m, seed = seed, num.threads = "1:1")

# One posterior draw per transition. A transition whose fit failed, or whose
# names are not found, falls back to the seed point estimate (theta_seed); if that
# is not ok either it stays ok = FALSE (hazard 0 / midpoint in the imputation).
draw_theta_m <- function(models, inp, theta_seed, sampler = default_sampler) {
  K <- nrow(inp$td); P <- length(inp$terms)
  th <- list(alpha = theta_seed$alpha, intercept = theta_seed$intercept, beta = theta_seed$beta, ok = theta_seed$ok)
  n_fallback <- 0L
  for (k in seq_len(K)) {
    m <- models[[k]]; if (!is.null(m$error)) { n_fallback <- n_fallback + 1L; next }
    samp <- sampler(1, m, seed = sample.int(1e6, 1L))[[1]]   # seed from R's RNG: reproducible after set.seed()
    ln <- rownames(samp$latent)
    ii <- grepl("^\\(Intercept\\)", ln)
    tms <- colnames(inp$mask)[inp$mask[k, ]]
    idx <- lapply(tms, function(tm) grepl(paste0("^", tm, "($|:)"), ln))
    if (!any(ii) || !all(vapply(idx, any, logical(1))) || !(ALPHA_NAME %in% names(samp$hyperpar))) {
      n_fallback <- n_fallback + 1L; warning(sprintf("draw_theta_m(): names not found for %s, using seed estimate", inp$td$label[k])); next
    }
    th$intercept[k] <- samp$latent[ii, 1][1] - LOG_TIME_SCALE
    th$beta[k, ] <- 0
    for (j in seq_along(tms)) th$beta[k, tms[j]] <- samp$latent[idx[[j]], 1][1]
    th$alpha[k] <- samp$hyperpar[[ALPHA_NAME]]
    th$ok[k] <- TRUE
  }
  th$n_fallback <- n_fallback
  th
}

# lambda (years) per plot x transition from a covariate matrix; ok = FALSE -> 0
lambda_matrix <- function(theta, X) {
  eta <- X %*% t(theta$beta)
  lam <- exp(sweep(eta, 2, theta$intercept, "+"))
  lam[, !theta$ok] <- 0
  lam[is.na(lam)] <- 0
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

# mode "mean": weighted mean (Stage 2); "draw": weighted draw (MI loop).
# Plots without a usable candidate get W_j/2 and fallback = TRUE.
p1_entries <- function(theta, inp, M = SIM_M, mode = c("mean", "draw"), n_cores = 1L) {
  mode <- match.arg(mode)
  L <- lambda_matrix(theta, inp$X)
  W <- inp$p1$W_j; sj <- inp$p1$s_idx; lo <- inp$p1$obs_lower; hi <- inp$p1$obs_upper
  to_idx <- inp$td$to_idx
  one <- function(j) {
    cand <- p1_candidates(inp$start_idx, sj[j], W[j], theta$alpha, L[j, ], inp$out_k, to_idx, M)
    if (!length(cand)) return(c(NA_real_, 0))
    w <- p1_weights(cand, sj[j], lo[j], hi[j], theta$alpha, L[j, ], inp$out_k)
    if (is.null(w)) return(c(NA_real_, length(cand)))
    c(if (mode == "mean") sum(cand * w) else cand[sample.int(length(cand), 1L, prob = w)], length(cand))
  }
  res <- parallel::mclapply(seq_along(W), one, mc.cores = n_cores)
  bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r), logical(1))
  if (any(bad)) stop("p1_entries: ", sum(bad), " plots failed in the workers: ", as.character(res[bad][[1]]))
  res <- do.call(rbind, res)
  fb <- is.na(res[, 1])
  entry <- ifelse(fb, W / 2, res[, 1])
  data.table(id_pe = inp$p1$id_pe, entry = entry, n_cand = res[, 2], fallback = fb)
}

# ------------------------------------------------------------------------------
# P2: moment-matched exit (copied unchanged from p2_imputation_strategies.R)
# ------------------------------------------------------------------------------
moment_match_weibull_naive <- function(a, b, alpha, lambda) {
  S_a <- exp(-(lambda * a)^alpha)
  S_b <- exp(-(lambda * b)^alpha)
  prob_interval <- S_a - S_b
  if (!is.finite(prob_interval) || prob_interval < 1e-10) return((a + b) / 2)
  integrand <- function(t) (t^alpha) * alpha * lambda^alpha * t^(alpha - 1) * exp(-(lambda * t)^alpha) / prob_interval
  e_t_alpha <- tryCatch(integrate(integrand, lower = a, upper = b, rel.tol = 1e-6)$value, error = function(e) NA_real_)
  if (!is.finite(e_t_alpha) || e_t_alpha <= 0) return((a + b) / 2)
  max(a, min(b, e_t_alpha^(1 / alpha)))
}

# p2_fun for build_rows(): moment match for kept transitions with an ok fit,
# midpoint otherwise. `diag` (an environment) counts both cases.
make_p2_moment <- function(theta, inp, diag = new.env()) {
  diag$n_event <- 0L; diag$n_midpoint <- 0L
  function(a_dur, b_dur, from, to, r) {
    diag$n_event <- diag$n_event + 1L
    k <- inp$kmat[from, to]
    if (is.na(k) || !theta$ok[k]) { diag$n_midpoint <- diag$n_midpoint + 1L; return((a_dur + b_dur) / 2) }
    lam <- exp(theta$intercept[k] + sum(theta$beta[k, ] * inp$Xiv[r, ]))
    moment_match_weibull_naive(a_dur, b_dur, theta$alpha[k], lam)
  }
}

# ------------------------------------------------------------------------------
# fitting
# ------------------------------------------------------------------------------
fit_rows <- function(rows, inp, fit_fun) {
  d <- rows[kept == TRUE]
  idx <- match(paste(d$from, d$to, sep = "->"), inp$td$label)
  stopifnot(!anyNA(idx))
  cols <- c("from", "to", "status", "Tstop", "truncation", inp$terms)
  lst <- split(as.data.frame(d[, ..cols]), factor(idx, levels = seq_len(nrow(inp$td))))
  n_events <- vapply(lst, function(x) sum(x$status == 1L), integer(1))
  models <- fit_fun(unname(lst), rhs_formula = unname(inp$rhs), exact = TRUE)
  list(models = models, n_events = n_events, data = lst)
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
# (Intercept, betas, alpha) from n_cov posterior samples. Intercept is on INLA's decade scale
# (subtract LOG_TIME_SCALE for years); covariances are unaffected by that shift.
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

tidy_fit <- function(models, inp) {
  rbindlist(lapply(seq_len(nrow(inp$td)), function(k) {
    m <- models[[k]]; tms <- colnames(inp$mask)[inp$mask[k, ]]; pars <- c("alpha", "(Intercept)", tms)
    dgn <- fit_diagnostics(m)
    if (!is.null(m$error)) return(data.table(k = k, parameter = pars, est = NA_real_, sd = NA_real_, fit_error = TRUE,
                                             error_msg = dgn$error_msg, mode_status = dgn$mode_status, mlik = dgn$mlik, log_warn = dgn$log_warn))
    fx <- m$summary.fixed; hy <- m$summary.hyperpar
    g <- function(nm, col) if (nm %in% rownames(fx)) fx[nm, col] else NA_real_
    data.table(k = k, parameter = pars,
               est = c(hy[ALPHA_NAME, "mean"], g("(Intercept)", "mean") - LOG_TIME_SCALE, vapply(tms, g, numeric(1), col = "mean")),
               sd  = c(hy[ALPHA_NAME, "sd"],   g("(Intercept)", "sd"),                       vapply(tms, g, numeric(1), col = "sd")),
               fit_error = FALSE, error_msg = dgn$error_msg, mode_status = dgn$mode_status, mlik = dgn$mlik, log_warn = dgn$log_warn)
  }))[, label := inp$td$label[k]][]
}

# ------------------------------------------------------------------------------
# stages
# ------------------------------------------------------------------------------
run_stage1 <- function(inp, fit_fun) {
  fe <- setNames(inp$p1$W_j / 2, inp$p1$id_pe)
  b <- build_rows(inp$iv, fe, inp$states, inp$tr, keep_cols = inp$keep_cols)
  f <- fit_rows(b$rows, inp, fit_fun)
  list(theta = extract_theta(f$models, inp), tidy = tidy_fit(f$models, inp), n_events = f$n_events, n_floored = b$n_floored)
}

run_stage2 <- function(inp, theta1, fit_fun, M = SIM_M, n_cores = 1L) {
  p1 <- p1_entries(theta1, inp, M, "mean", n_cores)
  dg <- new.env()
  b <- build_rows(inp$iv, setNames(p1$entry, p1$id_pe), inp$states, inp$tr,
                  p2_fun = make_p2_moment(theta1, inp, dg), keep_cols = inp$keep_cols)
  f <- fit_rows(b$rows, inp, fit_fun)
  list(models = f$models, theta = extract_theta(f$models, inp), tidy = tidy_fit(f$models, inp), n_events = f$n_events,
       p1 = p1, p2_event = dg$n_event, p2_midpoint = dg$n_midpoint, n_floored = b$n_floored)
}

run_one_imputation <- function(m, inp, models_seed, theta_seed, fit_fun, M = SIM_M, sampler = default_sampler, keep_models = FALSE) {
  th <- draw_theta_m(models_seed, inp, theta_seed, sampler)
  p1 <- p1_entries(th, inp, M, "draw", 1L)
  dg <- new.env()
  b <- build_rows(inp$iv, setNames(p1$entry, p1$id_pe), inp$states, inp$tr,
                  p2_fun = make_p2_moment(th, inp, dg), keep_cols = inp$keep_cols)
  f <- fit_rows(b$rows, inp, fit_fun)
  res <- tidy_fit(f$models, inp); res[, m := m]
  ev <- data.table(k = seq_along(f$n_events), n_events = as.integer(f$n_events))
  res <- merge(res, ev, by = "k")
  out <- list(res = res, diag = data.table(m = m, p1_fallback = sum(p1$fallback), p1_no_candidate = sum(p1$n_cand == 0),
                                    p2_events = dg$n_event, p2_midpoint = dg$n_midpoint,
                                    n_theta_fallback = th$n_fallback, n_fit_error = sum(vapply(f$models, function(x) !is.null(x$error), logical(1))),
                                    n_floored = b$n_floored))
  tms_k <- lapply(seq_len(nrow(inp$td)), function(k) colnames(inp$mask)[inp$mask[k, ]])
  out$extras <- list(
    theta_m = data.table(k = seq_len(nrow(inp$td)), label = inp$td$label, ok = th$ok, alpha = th$alpha, intercept_years = th$intercept,
                         as.data.table(th$beta)),
    p1 = p1[, .(id_pe, entry, n_cand, fallback)],
    p2 = b$segs[type == "event", .(id_pe, interval, from, to = event_to, exit_abs = end_abs, floored)],
    fits = setNames(lapply(seq_len(nrow(inp$td)), function(k) fit_extras(f$models[[k]], tms_k[[k]], sampler, N_COV_SAMPLES, sample.int(1e6, 1L))), inp$td$label))
  # fitting data of failed or near-degenerate fits (alpha SD < SD_MIN), to rerun them with verbose = TRUE
  al <- res[parameter == "alpha"]
  bad_k <- al[fit_error == TRUE | (!is.na(sd) & sd < SD_MIN), k]
  if (length(bad_k)) out$bad_data <- setNames(lapply(bad_k, function(k) list(data = f$data[[k]], rhs = inp$rhs[[k]], label = inp$td$label[k])), inp$td$label[bad_k])
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
