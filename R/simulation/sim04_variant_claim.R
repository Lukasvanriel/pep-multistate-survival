# ==============================================================================
# R/simulation/sim04_variant_claim.R   (2026-11-01)   RUN ON THE MACHINE WITH INLA (independent of sim01-03)
#
# QUESTION: does D. Rustand's observation (r-inla group, 15 Nov 2025: weibullsurv variant 1 gives "more stable results"
#  than variant 0) hold (a) in his own setting when replicated, and (b) when the setting moves toward ours?
#  His setting (as far as the thread shows; his code is not in this project): one interval-censored Weibull transition,
#  true shape 0.8, follow-up 50, N = 3000, fixed intervals of width 1-5, INTERCEPT ONLY, INLA's default alpha prior
#  (pc.alphaw 5) and an intercept prior N(0, precision 0.01). One run per setting.
#  Ours: imputed EXACT times, 4-13 coefficients, 15-100 events, pc.alphaw(0.1), coefficient prior N(0, 5^2).
# DESIGN (fully crossed, 20 replicates; all data simulated here, no project data):
#   N            3000 (his) / 100 / 40 (~75 / ~30 events: our events-per-coefficient range)
#   K            0 (his: intercept only) / 3 (true log HR 0.4, -0.3, 0.2) / 10 (the same 3 + 7 noise covariates, true 0)
#                K0 data are simulated WITHOUT covariate effects (first run, 2026-11-01 11:41-12:41, simulated them WITH the
#                3 effects and fitted intercept only = omitted-covariate misspecification; those K0 results are in
#                fits/_superseded_K0_with_effects/ and datasets/_superseded_K0_with_effects/, not used). K3/K10 unchanged.
#   observation  exact / interval-censored on a fixed grid of width 1, 3, 5 (event in (L, U]; T > 50 right-censored at 50)
#   variant      0 / 1
#   priors       "denis" (pc.alphaw 5; intercept N(0, prec 0.01); coefficients INLA default prec 0.001)
#                "pipeline" (pc.alphaw 0.1; flat intercept; coefficients N(0, 5^2))
#  Truth (proportional hazards): H(t | x) = exp(b0 + x'beta) t^0.8, b0 = log(log 2) - 0.8 log 20 (median 20 y at x = 0).
#  Variant 1 truth: alpha 0.8, coefficient beta / 0.8; log HR compared on the same scale (variant 1: alpha * beta).
#  All observation types and all fits of one replicate use the SAME simulated event times.
# MEASURES per cell: alpha relative bias (median, mean), coverage of the INLA 95% interval, median posterior SD, INLA
#  errors, frozen alpha-SD screen failures; log HR of the 3 true covariates (median error, coverage); share of noise
#  covariates whose interval excludes 0.
# OUTPUT data-output/simulation/sim04_variant_claim/: datasets/ (seed + times + covariates per replicate), fits/,
#  sim04_estimates.csv, sim04_alpha_summary.csv, sim04_coef_summary.csv, run_log.txt, progress.log
# Run from the project root:  Rscript R/simulation/sim04_variant_claim.R   (re-running skips finished replicates)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(tools) })
options(width = 200)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

FIT_ENGINE  <- Sys.getenv("SIM04_ENGINE", "inla")              # "mle" = local plumbing check only
N_WORKERS   <- as.integer(Sys.getenv("SIM04_WORKERS", "6"))
N_REP       <- as.integer(Sys.getenv("SIM04_NREP", "20"))
MASTER_SEED <- 20261102L
ALPHA_TRUE  <- 0.8; FOLLOW_UP <- 50; MEDIAN_T <- 20
B0_TRUE     <- log(log(2)) - ALPHA_TRUE * log(MEDIAN_T)
BETA_TRUE   <- c(0.4, -0.3, 0.2)                                # log HR of x1..x3; x4..x10 are noise (true 0)
NS <- c(3000L, 100L, 40L); KS <- c(0L, 3L, 10L); OBS <- c("exact", "ic1", "ic3", "ic5")
PRIORS <- list(denis = list(a_param = 5, prec = 0.001, prec_icpt = 0.01), pipeline = list(a_param = 0.1, prec = 1 / 25, prec_icpt = 0))
ALPHA_NAME <- "alpha parameter for weibullsurv"

OUT <- here("data-output", "simulation", paste0("sim04_variant_claim", if (FIT_ENGINE != "inla") paste0("_", FIT_ENGINE, "_check")))
for (s in c("datasets", "fits")) dir.create(file.path(OUT, s), recursive = TRUE, showWarnings = FALSE)
if (FIT_ENGINE == "inla") suppressPackageStartupMessages(library(INLA))
logf <- file.path(OUT, "run_log.txt"); progf <- file.path(OUT, "progress.log")
logl <- function(...) { s <- sprintf(...); cat(s, "\n"); cat(format(Sys.time()), s, "\n", file = logf, append = TRUE) }
prog <- function(console, ...) { s <- sprintf(...); if (console) cat(s, "\n"); cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), s, "\n", file = progf, append = TRUE) }
eta_text <- function(t_start, n_done_run, n_left) {
  if (n_done_run < 1) return("ETA: after the first one finishes")
  el <- as.numeric(difftime(Sys.time(), t_start, units = "mins")); rem <- el / n_done_run * n_left
  sprintf("elapsed %.1f min, ~%.0f min left (finish ~%s)", el, rem, format(Sys.time() + rem * 60, "%H:%M"))
}
self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); if (!length(self)) self <- here("R", "simulation", "sim04_variant_claim.R")
logl("==== start | engine %s | R %s | INLA %s", FIT_ENGINE, getRversion(), if (FIT_ENGINE == "inla") as.character(packageVersion("INLA")) else "not used")
logl("script %s md5 %s", self, md5sum(self))
logl("truth: alpha %.2f, b0 %.4f (median %g y), beta %s, follow-up %g | N %s | K %s | obs %s | reps %d", ALPHA_TRUE, B0_TRUE, MEDIAN_T,
     paste(BETA_TRUE, collapse = "/"), FOLLOW_UP, paste(NS, collapse = "/"), paste(KS, collapse = "/"), paste(OBS, collapse = "/"), N_REP)
capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))

# ---- data ---------------------------------------------------------------------------------------------------------
JOBS <- CJ(n = NS, k = KS, rep = seq_len(N_REP))
JOBS[, `:=`(seed = MASTER_SEED + 100000L * match(n, NS) + 10000L * match(k, KS) + rep, job_id = sprintf("N%d_K%d_rep%02d", n, k, rep))]
simulate <- function(n, k, seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * 10), n, 10); colnames(X) <- paste0("x", 1:10)     # always 10 columns drawn, so x1..x3 are identical across K
  eta <- B0_TRUE + if (k == 0) 0 else drop(X[, 1:3] %*% BETA_TRUE)  # K0: no covariate effects (homogeneous, as in Denis's setting; fixed 2026-11-01)
  T <- (-log(runif(n)) / exp(eta))^(1 / ALPHA_TRUE)                      # H(T) = exp(eta) T^alpha ~ Exp(1)
  data.table(id = seq_len(n), T = T, X[, seq_len(max(k, 1)), drop = FALSE])[, if (k == 0) .(id, T) else .SD]
}
observe <- function(d, obs) {
  d <- copy(d)
  if (obs == "exact") { d[, `:=`(status = as.integer(T <= FOLLOW_UP), time = pmin(T, FOLLOW_UP), time2 = pmin(T, FOLLOW_UP))]; return(d[]) }
  w <- as.numeric(sub("ic", "", obs)); grid <- sort(unique(c(seq(0, FOLLOW_UP, by = w), FOLLOW_UP)))
  ev <- d$T <= FOLLOW_UP; U <- grid[findInterval(d$T, grid, left.open = TRUE) + 1L]; L <- grid[findInterval(d$T, grid, left.open = TRUE)]
  d[, `:=`(status = fifelse(ev, 3L, 0L), time = fifelse(ev, L, FOLLOW_UP), time2 = fifelse(ev, U, FOLLOW_UP))]
  d[]
}
terms_k <- function(k) if (k == 0) character(0) else paste0("x", seq_len(k))

# ---- fits ---------------------------------------------------------------------------------------------------------
fit_inla_c <- function(d, terms, variant, pr) {
  rhs <- paste(c("1", terms), collapse = " + ")
  f <- as.formula(if (any(d$status == 3L)) sprintf("inla.surv(time = time, time2 = time2, event = status) ~ %s", rhs)
                  else sprintf("inla.surv(time = time, event = status) ~ %s", rhs))
  m <- tryCatch(INLA::inla(f, family = "weibullsurv", data = d,
                           control.family = list(variant = variant, hyper = list(alpha = list(prior = "pc.alphaw", param = pr$a_param))),
                           control.fixed = list(prec = pr$prec, mean.intercept = 0, prec.intercept = pr$prec_icpt), control.compute = list(config = TRUE),
                           control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)), num.threads = "1:1"),
                error = function(e) list(error = conditionMessage(e)))
  pars <- c("alpha", "(Intercept)", terms)
  if (!is.null(m$error)) return(data.table(parameter = pars, mean = NA_real_, sd = NA_real_, q025 = NA_real_, q975 = NA_real_, error = substr(m$error, 1, 150)))
  cols <- c("mean", "sd", "0.025quant", "0.975quant"); fx <- m$summary.fixed; hy <- m$summary.hyperpar
  g <- function(nm) if (nm == "alpha") unlist(hy[ALPHA_NAME, cols]) else if (nm %in% rownames(fx)) unlist(fx[nm, cols]) else rep(NA_real_, 4)
  v <- t(vapply(pars, g, numeric(4)))
  data.table(parameter = pars, mean = v[, 1], sd = v[, 2], q025 = v[, 3], q975 = v[, 4], error = NA_character_)
}
fit_map_c <- function(d, terms, variant, pr) {                       # plumbing-only stand-in (MAP, Wald), exact and interval-censored
  X <- if (length(terms)) cbind(1, as.matrix(d[, ..terms])) else matrix(1, nrow(d), 1); s <- d$status; t1 <- d$time; t2 <- d$time2
  Hf <- function(t, a, eta) if (variant == 1) exp(a * eta) * t^a else exp(eta) * t^a
  nll <- function(th) { a <- exp(th[1]); eta <- drop(X %*% th[-1]); H1 <- Hf(t1, a, eta); H2 <- Hf(t2, a, eta)
    lh <- if (variant == 1) th[1] + a * eta + (a - 1) * log(pmax(t1, 1e-12)) else th[1] + eta + (a - 1) * log(pmax(t1, 1e-12))
    ll <- sum(ifelse(s == 1, lh - H1, ifelse(s == 0, -H1, log(pmax(exp(-H1) - exp(-H2), 1e-300)))))
    -(ll - pr$prec / 2 * sum(th[-(1:2)]^2) - pr$prec_icpt / 2 * th[2]^2) }
  o <- tryCatch(optim(c(0, -2, rep(0, length(terms))), nll, method = "BFGS", control = list(maxit = 10000, reltol = 1e-12)), error = function(e) NULL)
  pars <- c("alpha", "(Intercept)", terms)
  if (is.null(o)) return(data.table(parameter = pars, mean = NA_real_, sd = NA_real_, q025 = NA_real_, q975 = NA_real_, error = "optim failed"))
  V <- tryCatch(solve(optimHess(o$par, nll)), error = function(e) NULL); se <- if (is.null(V)) rep(NA_real_, length(o$par)) else sqrt(pmax(diag(V), 0))
  data.table(parameter = pars, mean = c(exp(o$par[1]), o$par[-1]), sd = c(exp(o$par[1]) * se[1], se[-1]),
             q025 = c(exp(o$par[1] - 1.96 * se[1]), o$par[-1] - 1.96 * se[-1]), q975 = c(exp(o$par[1] + 1.96 * se[1]), o$par[-1] + 1.96 * se[-1]), error = NA_character_)
}
fit_c <- function(...) if (FIT_ENGINE == "inla") fit_inla_c(...) else fit_map_c(...)

run_job <- function(i) {
  jb <- JOBS[i]; f <- file.path(OUT, "fits", paste0(jb$job_id, ".rds")); if (file.exists(f)) return(invisible(NULL))
  d0 <- simulate(jb$n, jb$k, jb$seed); saveRDS(list(seed = jb$seed, data = d0), file.path(OUT, "datasets", paste0(jb$job_id, ".rds")))
  tm <- terms_k(jb$k); res <- list()
  for (ob in OBS) { d <- observe(d0, ob)
    for (v in 0:1) for (p in names(PRIORS)) {
      t1 <- proc.time()[["elapsed"]]; tb <- fit_c(d, tm, v, PRIORS[[p]])
      res[[length(res) + 1L]] <- tb[, `:=`(obs = ob, variant = v, prior = p, n_events = sum(d$status != 0), sec = proc.time()[["elapsed"]] - t1)] } }
  R <- rbindlist(res)[, `:=`(job_id = jb$job_id, n = jb$n, k = jb$k, rep = jb$rep, seed = jb$seed)]
  saveRDS(R, f)
  n_done <- sum(file.exists(file.path(OUT, "fits", paste0(JOBS$job_id, ".rds"))))
  prog(TRUE, "[sim04] %d/%d done: %s | events %d | INLA errors %d of %d | %s", n_done, nrow(JOBS), jb$job_id, R$n_events[1], R[parameter == "alpha" & !is.na(error), .N],
       length(OBS) * 2 * length(PRIORS), eta_text(T0, n_done - N0, nrow(JOBS) - n_done))
  invisible(NULL)
}
todo <- which(!file.exists(file.path(OUT, "fits", paste0(JOBS$job_id, ".rds"))))
N0 <- nrow(JOBS) - length(todo); T0 <- Sys.time()
prog(TRUE, "[sim04] start: %d of %d replicates to run, %d fits each (%d workers)", length(todo), nrow(JOBS), length(OBS) * 2 * length(PRIORS), N_WORKERS)
err <- unlist(mclapply(todo, function(i) tryCatch({ run_job(i); NULL }, error = function(e) paste(JOBS$job_id[i], conditionMessage(e))), mc.cores = N_WORKERS, mc.preschedule = FALSE))
if (length(err)) { logl("ERRORS:\n%s", paste(err, collapse = "\n")); stop("job errors, see run_log.txt") }

# ---- summaries -------------------------------------------------------------------------------------------------------
R <- rbindlist(lapply(file.path(OUT, "fits", paste0(JOBS$job_id, ".rds")), readRDS))
al <- R[parameter == "alpha", .(job_id, obs, variant, prior, a = mean, a_sd = sd)]
R <- merge(R, al, by = c("job_id", "obs", "variant", "prior"))
R[, failed := !is.na(error) | is.na(mean)]
R[, screen_fail := !failed & a_sd / (0.78 * a / sqrt(pmax(n_events, 1))) < 0.1]
R[, xi := suppressWarnings(as.integer(sub("^x", "", parameter)))]
R[, truth_loghr := fifelse(!is.na(xi), fifelse(xi <= 3, BETA_TRUE[pmin(xi, 3)], 0), NA_real_)]
R[, truth := fcase(parameter == "alpha", ALPHA_TRUE, !is.na(xi) & variant == 0L, truth_loghr, !is.na(xi) & variant == 1L, truth_loghr / ALPHA_TRUE, default = NA_real_)]
R[!is.na(xi), loghr := fifelse(variant == 1L, a * mean, mean)]
R[, covered := !failed & q025 <= truth & truth <= q975]
fwrite(R, file.path(OUT, "sim04_estimates.csv"))
A <- R[parameter == "alpha", .(fits = .N, failed = sum(failed), screen_fail = sum(screen_fail), events = round(mean(n_events)),
                               median_rel_bias_pct = round(100 * median((mean - ALPHA_TRUE) / ALPHA_TRUE, na.rm = TRUE), 1),
                               mean_rel_bias_pct = round(100 * mean((mean - ALPHA_TRUE) / ALPHA_TRUE, na.rm = TRUE), 1),
                               coverage = round(mean(covered[!failed]), 2), median_sd = signif(median(sd, na.rm = TRUE), 3)), keyby = .(n, k, obs, prior, variant)]
fwrite(A, file.path(OUT, "sim04_alpha_summary.csv"))
C <- R[!is.na(xi) & !failed, .(n_fits = .N, true_median_err_loghr = round(median((loghr - truth_loghr)[xi <= 3], na.rm = TRUE), 3),
                              true_coverage = round(mean(covered[xi <= 3]), 2),
                              noise_false_pos = if (any(xi > 3)) round(mean((q025 > 0 | q975 < 0)[xi > 3]), 3) else NA_real_), keyby = .(n, k, obs, prior, variant)]
fwrite(C, file.path(OUT, "sim04_coef_summary.csv"))
cat("\nalpha: relative bias % | coverage, variant 0 vs 1\n")
print(dcast(A, n + k + obs + prior ~ variant, value.var = c("median_rel_bias_pct", "coverage", "failed")), nrows = 200)
logl("outputs: sim04_estimates.csv (%d rows), sim04_alpha_summary.csv, sim04_coef_summary.csv", nrow(R))
logl("==== done")
