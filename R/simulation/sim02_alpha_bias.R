# ==============================================================================
# R/simulation/sim02_alpha_bias.R   (2026-11-01)   RUN ON THE MACHINE WITH INLA, AFTER sim01
#
# QUESTION: sim01 found that INLA's posterior for the Weibull shape alpha is biased LOW when there are few events per
#  coefficient (about -30% at 2-4 events per coefficient), while maximum likelihood on the same datasets is not. Which
#  of the candidate causes is it? A deliberately ROUGH, wide check (few datasets per variant) to zoom in afterwards.
#
# DATASETS: 3 sim01 cells where the bias is large and the fitted model is the true one (small harvest effect, 2-4 events
#  per coefficient in fit B): 10 datasets each = the first sim01 replicates with >= 15 events and no row above INLA's cap.
#  They are READ from sim01's saved outcomes (same data as sim01's fit B). Two variants need new data, generated from the
#  same truth with their own seeds: added noise covariates, and the same sojourns without left truncation.
#
# VARIANTS (all INLA weibullsurv; "ref" = sim01 fit B exactly: B terms, pc.alphaw(0.1), N(0, 5^2) on coefficients,
#  strategy adaptive, int.strategy grid, variant 1, time in decades):
#   prior        P_a0.01, P_a0.001 (alpha PC rate), P_aflat (flat prior on INLA's internal alpha scale = flat on log alpha),
#                P_s1, P_sdef (coefficient prior SD 1 / INLA default 31.6)
#   approximation A_laplace, A_gaussian (strategy), A_ccd, A_eb (int.strategy), A_finegrid (grid dz 0.25, diff.logdens 10)
#   coupling     C_years, C_centuries (exact time rescaling), C_variant0 (S = exp(-lambda t^alpha): coefficients no longer
#                multiplied by alpha), C_fixbeta (all coefficients fixed at the truth by an offset: only intercept + alpha),
#                C_fixicpt (intercept fixed at the truth by an offset: alpha + coefficients)
#   coefficients N_C (fewer terms: no harvest terms; ~true because the harvest effect is small), N_C3, N_C6 (C + 3 / 6 pure
#                noise covariates, true coefficient 0), N_Breduced (B without the per-level GDD / TSD terms of levels with
#                < 5 events in that dataset)
#   data         D_notrunc (same sojourns shifted to start at time 0, i.e. no left truncation; new outcomes, ref fit)
#  REFERENCES on every dataset: MAP (maximum likelihood with the same N(0, 5^2) coefficient penalty, no alpha prior; = the
#  "MLE" of the sim01 comparison) and, on 5 datasets, an exact posterior by random-walk Metropolis MCMC (same likelihood,
#  coefficient prior N(0, 5^2), flat intercept; alpha prior pc.alphaw(0.1) if INLA exposes its density, and flat on log alpha),
#  in BOTH parameterisations: variant 1 (pipeline: coefficients multiplied by alpha) and variant 0 (coefficients = log HR).
#  Local check 2026-11-01 (MCMC, 5 datasets, flat log-alpha prior): variant-1 posterior mean of alpha 0.2-0.5 below the
#  maximum-likelihood value, variant 0 within 0.02 of it -> the bias may come from the parameterisation x coefficient prior.
#
# READ-ONLY inputs: data-output/simulation/sim01_harvest_reference/ (designs, calibration, outcomes, datasets_index) and
#  data-output/main/results/pooled_main_final.csv. Nothing in sim01, the pipeline or the simulation folder is changed.
# REPRODUCIBILITY: new data from fixed seeds (sim01 seed + 5e6 / 7e6), outcomes saved; MCMC seeded; checkpoint per dataset;
#  run_log.txt (versions, md5s), progress.log (one line per dataset; tail -f it).
# OUTPUT data-output/simulation/sim02_alpha_bias/: fits/, outcomes/, alpha_results.csv (every dataset x variant),
#  alpha_summary.csv (per cell x variant), alpha_summary_overall.csv, mcmc.csv, run_log.txt, progress.log
# Run from the project root:  Rscript R/simulation/sim02_alpha_bias.R      (re-running skips finished datasets)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(tools) })
source(here("R", "config.R"))
options(width = 200)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# ---- settings -------------------------------------------------------------------------------------------------
FIT_ENGINE  <- Sys.getenv("SIM02_ENGINE", "inla")              # "mle" = local plumbing check only (no INLA)
N_WORKERS   <- as.integer(Sys.getenv("SIM02_WORKERS", "6"))
N_DS        <- as.integer(Sys.getenv("SIM02_NDS", "10"))         # datasets per cell
SIM01_DIR   <- here("data-output", "simulation", Sys.getenv("SIM02_SIM01", "sim01_harvest_reference"))
CELLS       <- c("black_spruce_to_fir__small__d0_30", "regen_to_yellow_birch__small__d0_5", "sugar_maple_to_regen__small__d0_10")
EFFECT_SMALL <- c(light = 0.2, heavy = 0.3)                      # sim01 "small" (log HR)
HARV_TERMS  <- c("harv_light", "harv_heavy", "gdd_x_harv_light", "gdd_x_harv_heavy", "tsd_harv_light", "tsd_harv_heavy")
SIGMA_B     <- 5
N_MCMC_DS   <- 5L; MCMC_ITER <- 40000L; MCMC_BURN <- 10000L; MCMC_SEED <- 20261101L
SEED_NOTRUNC <- 5000000L; SEED_NOISE <- 7000000L
TIME_SCALE <- 10; LOG_TIME_SCALE <- log(10); EPSILON <- 1e-6
ALPHA_NAME <- "alpha parameter for weibullsurv"
H_CAP <- -log(0.001 * 2^-23)

V <- rbindlist(list(
  list("ref",        "reference",     "B",  "pc", 0.1,   1 / SIGMA_B^2, "adaptive", "grid", NA, NA, 1, 1,   "none", "sim01"),
  list("P_a0.01",    "prior",         "B",  "pc", 0.01,  1 / SIGMA_B^2, "adaptive", "grid", NA, NA, 1, 1,   "none", "sim01"),
  list("P_a0.001",   "prior",         "B",  "pc", 0.001, 1 / SIGMA_B^2, "adaptive", "grid", NA, NA, 1, 1,   "none", "sim01"),
  list("P_aflat",    "prior",         "B",  "flat", NA,  1 / SIGMA_B^2, "adaptive", "grid", NA, NA, 1, 1,   "none", "sim01"),
  list("P_s1",       "prior",         "B",  "pc", 0.1,   1,             "adaptive", "grid", NA, NA, 1, 1,   "none", "sim01"),
  list("P_sdef",     "prior",         "B",  "pc", 0.1,   0.001,         "adaptive", "grid", NA, NA, 1, 1,   "none", "sim01"),
  list("A_laplace",  "approximation", "B",  "pc", 0.1,   1 / SIGMA_B^2, "laplace",  "grid", NA, NA, 1, 1,   "none", "sim01"),
  list("A_gaussian", "approximation", "B",  "pc", 0.1,   1 / SIGMA_B^2, "gaussian", "grid", NA, NA, 1, 1,   "none", "sim01"),
  list("A_ccd",      "approximation", "B",  "pc", 0.1,   1 / SIGMA_B^2, "adaptive", "ccd",  NA, NA, 1, 1,   "none", "sim01"),
  list("A_eb",       "approximation", "B",  "pc", 0.1,   1 / SIGMA_B^2, "adaptive", "eb",   NA, NA, 1, 1,   "none", "sim01"),
  list("A_finegrid", "approximation", "B",  "pc", 0.1,   1 / SIGMA_B^2, "adaptive", "grid", 0.25, 10, 1, 1, "none", "sim01"),
  list("C_years",    "coupling",      "B",  "pc", 0.1,   1 / SIGMA_B^2, "adaptive", "grid", NA, NA, 1, 10,  "none", "sim01"),
  list("C_centuries","coupling",      "B",  "pc", 0.1,   1 / SIGMA_B^2, "adaptive", "grid", NA, NA, 1, 0.1, "none", "sim01"),
  list("C_variant0", "coupling",      "B",  "pc", 0.1,   1 / SIGMA_B^2, "adaptive", "grid", NA, NA, 0, 1,   "none", "sim01"),
  list("C_fixbeta",  "coupling",      "B",  "pc", 0.1,   1 / SIGMA_B^2, "adaptive", "grid", NA, NA, 1, 1,   "betas", "sim01"),
  list("C_fixicpt",  "coupling",      "B",  "pc", 0.1,   1 / SIGMA_B^2, "adaptive", "grid", NA, NA, 1, 1,   "intercept", "sim01"),
  list("N_C",        "coefficients",  "C",  "pc", 0.1,   1 / SIGMA_B^2, "adaptive", "grid", NA, NA, 1, 1,   "none", "sim01"),
  list("N_C3",       "coefficients",  "C3", "pc", 0.1,   1 / SIGMA_B^2, "adaptive", "grid", NA, NA, 1, 1,   "none", "sim01"),
  list("N_C6",       "coefficients",  "C6", "pc", 0.1,   1 / SIGMA_B^2, "adaptive", "grid", NA, NA, 1, 1,   "none", "sim01"),
  list("N_Breduced", "coefficients",  "Bred","pc", 0.1,  1 / SIGMA_B^2, "adaptive", "grid", NA, NA, 1, 1,   "none", "sim01"),
  list("D_notrunc",  "data",          "B",  "pc", 0.1,   1 / SIGMA_B^2, "adaptive", "grid", NA, NA, 1, 1,   "none", "notrunc")))
setnames(V, c("variant", "group", "terms", "a_prior", "a_param", "beta_prec", "strategy", "int_strategy", "dz", "diff_logdens",
              "family_variant", "time_mult", "offset", "data"))

OUT <- here("data-output", "simulation", paste0("sim02_alpha_bias", if (FIT_ENGINE != "inla") paste0("_", FIT_ENGINE, "_check")))
for (s in c("fits", "outcomes")) dir.create(file.path(OUT, s), recursive = TRUE, showWarnings = FALSE)
if (FIT_ENGINE == "inla") suppressPackageStartupMessages(library(INLA))
logf <- file.path(OUT, "run_log.txt"); progf <- file.path(OUT, "progress.log")
logl <- function(...) { s <- sprintf(...); cat(s, "\n"); cat(format(Sys.time()), s, "\n", file = logf, append = TRUE) }
prog <- function(console, ...) { s <- sprintf(...); if (console) cat(s, "\n"); cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), s, "\n", file = progf, append = TRUE) }

self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); if (!length(self)) self <- here("R", "simulation", "sim02_alpha_bias.R")
logl("==== start | engine %s | R %s | data.table %s | INLA %s", FIT_ENGINE, getRversion(), packageVersion("data.table"),
     if (FIT_ENGINE == "inla") as.character(packageVersion("INLA")) else "not used")
logl("script %s md5 %s", self, md5sum(self))
for (f in c(file.path(SIM01_DIR, c("calibration.csv", "datasets_index.csv")), file.path(RESULTS_DIR, sprintf("pooled_%s_final.csv", RUN_TAG))))
  logl("input %s md5 %s", sub(here(), "", f, fixed = TRUE), md5sum(f))
logl("cells %s | %d datasets per cell | %d variants | MCMC on %d datasets (%d iterations, burn-in %d)", paste(CELLS, collapse = ", "), N_DS, nrow(V),
     N_MCMC_DS, MCMC_ITER, MCMC_BURN)
capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))

# ---- inputs from sim01 -------------------------------------------------------------------------------------------
cal <- fread(file.path(SIM01_DIR, "calibration.csv")); idx <- fread(file.path(SIM01_DIR, "datasets_index.csv"))
pooled <- fread(file.path(RESULTS_DIR, sprintf("pooled_%s_final.csv", RUN_TAG)))
stopifnot(all(CELLS %in% cal$cell_id))
DS <- idx[cell_id %in% CELLS & n_events >= 15 & n_rows_above_cap == 0][order(cell_id, rep)][, head(.SD, N_DS), by = cell_id]
DS <- merge(DS, cal[, .(cell_id, template, effect, d0_target, c_yr)], by = "cell_id")
DS[, ds_id := sprintf("%s_rep%02d", cell_id, rep)]
logl("datasets: %s", paste(DS[, .N, by = cell_id][, sprintf("%s %d", cell_id, N)], collapse = " | "))
designs <- lapply(setNames(unique(DS$template), unique(DS$template)), function(lb) readRDS(file.path(SIM01_DIR, "design", sprintf("design_%s.rds", gsub("->", "_to_", lb)))))
truth_for <- function(des) {                                  # sim01 truth_for() for the "small" effect, from the saved design
  p <- pooled[label == des$label]; a <- p[parameter == "alpha", pooled_est]
  b <- setNames(p[match(des$base, parameter), pooled_est], des$base); stopifnot(!anyNA(b))
  h <- c(harv_light = EFFECT_SMALL[["light"]], harv_heavy = EFFECT_SMALL[["heavy"]], gdd_x_harv_light = 0, gdd_x_harv_heavy = 0,
         tsd_harv_light = 0, tsd_harv_heavy = 0) / a
  list(alpha = a, beta = c(b, h))
}
seg_parts <- function(D, tru) { X <- as.matrix(D[, names(tru$beta), with = FALSE]); eta0 <- drop(X %*% tru$beta) + LOG_TIME_SCALE
  list(eta0 = eta0, b0 = exp(tru$alpha * eta0) * (D$Tstop^tru$alpha - D$truncation^tru$alpha)) }
expected_d0 <- function(c_yr, D, sp, alpha) { dH <- exp(alpha * c_yr) * sp$b0; cs <- ave(dH, D$soj, FUN = cumsum); p <- exp(-(cs - dH)) * (-expm1(-dH)); sum(p[D$harv == 0]) }
simulate_outcome <- function(D, sp, alpha, c_yr, seed) {      # identical to sim01
  set.seed(seed); n <- max(D$soj); u <- rexp(n)
  dH <- exp(alpha * c_yr) * sp$b0; cs <- ave(dH, D$soj, FUN = cumsum); prev <- cs - dH
  hit <- which(cs >= u[D$soj]); hit <- hit[!duplicated(D$soj[hit])]
  if (!length(hit)) return(data.table(soj = integer(0), row = integer(0), t = numeric(0)))
  rem <- u[D$soj[hit]] - prev[hit]; lam_a <- exp(alpha * (c_yr + sp$eta0[hit]))
  t <- (D$truncation[hit]^alpha + rem / lam_a)^(1 / alpha); t <- pmin(pmax(t, D$truncation[hit] + EPSILON), D$Tstop[hit])
  data.table(soj = D$soj[hit], row = hit, t = t)
}
regenerate_dataset <- function(D, outcome) {                  # identical to sim01
  cut <- rep(.Machine$integer.max, max(D$soj)); cut[outcome$soj] <- outcome$row
  keep <- seq_len(nrow(D)) <= cut[D$soj]; d <- D[keep]; d[, status := 0L]; ri <- match(outcome$row, which(keep))
  d[ri, `:=`(status = 1L, Tstop = outcome$t)]; d[]
}
notrunc_design <- function(D) { D <- copy(D); D[, s0 := truncation[1], by = soj]; D[, `:=`(truncation = truncation - s0, Tstop = Tstop - s0)][, s0 := NULL]; D[] }

# ---- fitting -------------------------------------------------------------------------------------------------------
terms_for <- function(set, des, info) {
  B <- c(des$base, HARV_TERMS)
  switch(set, B = B, C = des$base, C3 = c(des$base, paste0("noise", 1:3)), C6 = c(des$base, paste0("noise", 1:6)),
         Bred = setdiff(B, c(if (info$d_light < 5) c("gdd_x_harv_light", "tsd_harv_light"), if (info$d_heavy < 5) c("gdd_x_harv_heavy", "tsd_harv_heavy"))))
}
fit_inla_v <- function(d, terms, v, tru, c_yr) {
  d <- copy(d)
  if (v$time_mult != 1) d[, `:=`(Tstop = Tstop * v$time_mult, truncation = truncation * v$time_mult)]
  rhs <- c("1", terms)
  if (v$offset == "betas") { X <- as.matrix(d[, names(tru$beta), with = FALSE]); d[, off := drop(X %*% tru$beta)]; rhs <- c("1", "offset(off)") }
  if (v$offset == "intercept") { d[, off := c_yr + LOG_TIME_SCALE]; rhs <- c("-1", terms, "offset(off)") }
  f <- as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", paste(rhs, collapse = " + ")))
  hyp <- if (v$a_prior == "pc") list(alpha = list(prior = "pc.alphaw", param = v$a_param)) else list(alpha = list(prior = "flat", param = numeric(0)))
  ci <- list(strategy = v$strategy, int.strategy = v$int_strategy, control.vb = list(enable = FALSE))
  if (!is.na(v$dz)) ci$dz <- v$dz
  if (!is.na(v$diff_logdens)) ci$diff.logdens <- v$diff_logdens
  m <- tryCatch(INLA::inla(f, family = "weibullsurv", data = d, control.family = list(variant = v$family_variant, hyper = hyp),
                           control.fixed = list(prec = v$beta_prec, prec.intercept = 0), control.compute = list(config = TRUE),
                           control.inla = ci, num.threads = "1:1"), error = function(e) list(error = conditionMessage(e)))
  if (!is.null(m$error)) return(data.table(a_mean = NA_real_, a_sd = NA_real_, a_q025 = NA_real_, a_q500 = NA_real_, a_q975 = NA_real_,
                                           icpt_mean = NA_real_, mode_status = NA_integer_, error = substr(m$error, 1, 150)))
  h <- m$summary.hyperpar[ALPHA_NAME, ]; ic <- if ("(Intercept)" %in% rownames(m$summary.fixed)) m$summary.fixed["(Intercept)", "mean"] else NA_real_
  data.table(a_mean = h$mean, a_sd = h$sd, a_q025 = h$`0.025quant`, a_q500 = h$`0.5quant`, a_q975 = h$`0.975quant`, icpt_mean = ic,
             mode_status = tryCatch(as.integer(m$mode$mode.status)[1], error = function(e) NA_integer_), error = NA_character_)
}
# negative log posterior / MAP (no alpha prior unless lp_alpha given); variant 1, time in decades, optional offset
make_nlp <- function(d, terms, beta_prec, off = NULL, intercept = TRUE, lp_alpha = NULL, variant = 1) {   # variant 0: S = exp(-exp(eta) t^alpha)
  X <- as.matrix(d[, terms, with = FALSE]); if (intercept) X <- cbind(1, X)
  T <- d$Tstop; tr <- d$truncation; s <- d$status; lT <- log(T); o <- if (is.null(off)) 0 else off
  nb <- ncol(X) - intercept
  function(th) { la <- th[1]; a <- exp(la); b <- th[-1]; eta <- drop(X %*% b) + o
    ll <- if (variant == 1) sum(s * (la + a * eta + (a - 1) * lT)) - sum(exp(a * eta) * (T^a - tr^a)) else sum(s * (la + eta + (a - 1) * lT)) - sum(exp(eta) * (T^a - tr^a))
    pen <- if (nb > 0) beta_prec / 2 * sum(b[(intercept + 1):length(b)]^2) else 0
    -(ll - pen + (if (is.null(lp_alpha)) 0 else lp_alpha(la))) }
}
fit_map <- function(d, terms, beta_prec = 1 / SIGMA_B^2, variant = 1) {
  nlp <- make_nlp(d, terms, beta_prec, variant = variant)
  st <- c(0, log(max(sum(d$status), 0.5) / sum(d$Tstop - d$truncation)), rep(0, length(terms)))
  o <- tryCatch(optim(st, nlp, method = "BFGS", control = list(maxit = 10000, reltol = 1e-12)), error = function(e) NULL)
  if (is.null(o)) return(list(ok = FALSE))
  H <- tryCatch(optimHess(o$par, nlp), error = function(e) NULL); V <- tryCatch(solve(H), error = function(e) NULL)
  se <- if (is.null(V)) rep(NA_real_, length(o$par)) else sqrt(pmax(diag(V), 0))
  list(ok = TRUE, par = o$par, V = V, conv = o$convergence,
       row = data.table(a_mean = exp(o$par[1]), a_sd = exp(o$par[1]) * se[1], a_q025 = exp(o$par[1] - 1.96 * se[1]), a_q500 = exp(o$par[1]),
                        a_q975 = exp(o$par[1] + 1.96 * se[1]), icpt_mean = o$par[2], mode_status = o$convergence, error = NA_character_))
}
fit_variant <- function(d, terms, v, tru, c_yr) {
  if (FIT_ENGINE == "inla") return(fit_inla_v(d, terms, v, tru, c_yr))
  r <- fit_map(d, terms, v$beta_prec); if (!r$ok) return(data.table(a_mean = NA_real_, error = "map failed")); r$row   # plumbing check only
}

# ---- MCMC (exact posterior, random-walk Metropolis) ---------------------------------------------------------------
pc_logdens <- NULL
if (FIT_ENGINE == "inla" && exists("inla.pc.dalphaw", envir = asNamespace("INLA"))) {
  f_pc <- get("inla.pc.dalphaw", envir = asNamespace("INLA"))
  pc_logdens <- function(la) { a <- exp(la); v <- tryCatch(f_pc(a, lambda = 0.1, log = TRUE), error = function(e) NA_real_); v + la }   # + log Jacobian (prior on log alpha)
  if (!is.finite(pc_logdens(log(1.3)))) pc_logdens <- NULL
}
logl("MCMC alpha priors: %s", if (is.null(pc_logdens)) "flat on log alpha only (INLA's pc.alphaw density not available)" else "pc.alphaw(0.1) and flat on log alpha")
run_mcmc <- function(d, terms, lp_alpha, seed, variant = 1) {
  mp <- fit_map(d, terms, variant = variant); if (!mp$ok || is.null(mp$V)) return(NULL)
  nlp <- make_nlp(d, terms, 1 / SIGMA_B^2, lp_alpha = lp_alpha, variant = variant)
  set.seed(seed); k <- length(mp$par); Sig <- mp$V * 2.38^2 / k; L <- chol(Sig + diag(1e-10, k))
  th <- mp$par; cur <- -nlp(th); out <- numeric(MCMC_ITER - MCMC_BURN); acc <- 0L; keep <- matrix(NA_real_, MCMC_ITER, k)
  for (i in seq_len(MCMC_ITER)) {
    prop <- th + drop(rnorm(k) %*% L); lp <- -nlp(prop)
    if (is.finite(lp) && log(runif(1)) < lp - cur) { th <- prop; cur <- lp; if (i > MCMC_BURN) acc <- acc + 1L }
    keep[i, ] <- th
    if (i %in% c(2500L, 5000L, 7500L) && i <= MCMC_BURN) { S2 <- cov(keep[(i - 2499L):i, , drop = FALSE]); L <- tryCatch(chol(S2 * 2.38^2 / k + diag(1e-10, k)), error = function(e) L) }
    if (i > MCMC_BURN) out[i - MCMC_BURN] <- exp(th[1])
  }
  nb <- 40L; bm <- colMeans(matrix(out[seq_len(nb * floor(length(out) / nb))], ncol = nb)); ess <- length(out) * var(out) / (var(bm) * length(out) / nb)
  data.table(a_mean = mean(out), a_q025 = quantile(out, 0.025), a_q500 = median(out), a_q975 = quantile(out, 0.975), accept = acc / (MCMC_ITER - MCMC_BURN), ess = ess,
             map_alpha = exp(mp$par[1]))
}

# ---- one dataset: all variants -----------------------------------------------------------------------------------
noise_cols <- function(d, seed) { set.seed(seed); z <- matrix(rnorm(max(d$soj) * 6), ncol = 6); for (j in 1:6) set(d, j = paste0("noise", j), value = z[d$soj, j]); d }
run_dataset <- function(i) {
  ds <- DS[i]; f <- file.path(OUT, "fits", paste0(ds$ds_id, ".rds")); if (file.exists(f)) return(invisible(NULL))
  des <- designs[[ds$template]]; tru <- truth_for(des); D <- des$D
  oc <- readRDS(file.path(SIM01_DIR, "outcomes", paste0(ds$ds_id, ".rds")))
  d <- noise_cols(regenerate_dataset(D, oc), ds$seed + SEED_NOISE)
  ev <- d[status == 1L]; stopifnot(nrow(ev) == ds$n_events, ev[harv == 0, .N] == ds$d0)        # same data as sim01
  info <- list(d_light = ev[harv == 1, .N], d_heavy = ev[harv == 2, .N])
  # no-truncation data: same sojourns from time 0, intercept recalibrated to the same expected d0
  Dn <- notrunc_design(D); spn <- seg_parts(Dn, tru)
  c_n <- tryCatch(uniroot(function(c) log(expected_d0(c, Dn, spn, tru$alpha)) - log(ds$d0_target), c(-25, 0), extendInt = "upX", tol = 1e-10)$root, error = function(e) NA_real_)
  ocn <- simulate_outcome(Dn, spn, tru$alpha, c_n, ds$seed + SEED_NOTRUNC); saveRDS(ocn, file.path(OUT, "outcomes", paste0(ds$ds_id, "_notrunc.rds")))
  dn <- regenerate_dataset(Dn, ocn); evn <- dn[status == 1L]
  res <- list()
  for (j in seq_len(nrow(V))) {
    v <- V[j]; dd <- if (v$data == "notrunc") dn else d; cy <- if (v$data == "notrunc") c_n else ds$c_yr
    tm <- terms_for(v$terms, des, info)
    t0 <- proc.time()[["elapsed"]]; r <- fit_variant(dd, tm, v, tru, cy)
    res[[j]] <- cbind(data.table(variant = v$variant, group = v$group, n_coef = if (v$offset == "betas") 0L else length(tm), n_events = sum(dd$status),
                                 sec = proc.time()[["elapsed"]] - t0), r)
  }
  mp <- fit_map(d, terms_for("B", des, info)); mpn <- fit_map(dn, terms_for("B", des, info))
  res[[length(res) + 1L]] <- cbind(data.table(variant = "MAP_ref", group = "reference", n_coef = length(terms_for("B", des, info)), n_events = nrow(ev), sec = NA_real_), mp$row)
  res[[length(res) + 1L]] <- cbind(data.table(variant = "MAP_notrunc", group = "reference", n_coef = length(terms_for("B", des, info)), n_events = nrow(evn), sec = NA_real_), mpn$row)
  R <- rbindlist(res, fill = TRUE)
  R[, `:=`(ds_id = ds$ds_id, cell_id = ds$cell_id, template = ds$template, rep = ds$rep, alpha_true = tru$alpha, d0 = ds$d0, d_light = info$d_light, d_heavy = info$d_heavy)]
  saveRDS(R, f)
  n_done <- sum(file.exists(file.path(OUT, "fits", paste0(DS$ds_id, ".rds"))))
  prog(TRUE, "[sim02] dataset %d/%d done: %s | events %d | INLA errors %d of %d | ref alpha %.3f (true %.3f), MAP %.3f", n_done, nrow(DS), ds$ds_id, nrow(ev),
       R[!variant %like% "^MAP" & !is.na(error), .N], nrow(V), R[variant == "ref", a_mean], tru$alpha, R[variant == "MAP_ref", a_mean])
  invisible(NULL)
}

todo <- which(!file.exists(file.path(OUT, "fits", paste0(DS$ds_id, ".rds"))))
prog(TRUE, "[sim02] start: %d of %d datasets to run, %d variants each (%d workers)", length(todo), nrow(DS), nrow(V), N_WORKERS)
t0 <- Sys.time()
err <- unlist(mclapply(todo, function(i) tryCatch({ run_dataset(i); NULL }, error = function(e) paste(DS$ds_id[i], conditionMessage(e))), mc.cores = N_WORKERS, mc.preschedule = FALSE))
if (length(err)) { logl("ERRORS:\n%s", paste(err, collapse = "\n")); stop("dataset errors, see run_log.txt") }
logl("fits done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

# ---- MCMC on the first datasets ---------------------------------------------------------------------------------------
mc_f <- file.path(OUT, "mcmc.csv")
if (!file.exists(mc_f)) {
  mc_ds <- DS[, head(.SD, ceiling(N_MCMC_DS / length(CELLS))), by = cell_id][seq_len(N_MCMC_DS)]
  priors <- c(list(flat_log_alpha = function(la) 0), if (!is.null(pc_logdens)) list(pc_0.1 = pc_logdens))
  jobs <- CJ(i = seq_len(nrow(mc_ds)), p = names(priors), param = c(1L, 0L))    # parameterisation: variant 1 (pipeline) and variant 0
  MC <- rbindlist(mclapply(seq_len(nrow(jobs)), function(j) {
    ds <- mc_ds[jobs$i[j]]; des <- designs[[ds$template]]; tru <- truth_for(des)
    d <- regenerate_dataset(des$D, readRDS(file.path(SIM01_DIR, "outcomes", paste0(ds$ds_id, ".rds")))); ev <- d[status == 1L]
    info <- list(d_light = ev[harv == 1, .N], d_heavy = ev[harv == 2, .N])
    t0 <- proc.time()[["elapsed"]]; r <- run_mcmc(d, terms_for("B", des, info), priors[[jobs$p[j]]], MCMC_SEED + jobs$i[j], variant = jobs$param[j])
    prog(FALSE, "  MCMC %s prior %s variant %d done (%.0f s)", ds$ds_id, jobs$p[j], jobs$param[j], proc.time()[["elapsed"]] - t0)
    if (is.null(r)) return(NULL); cbind(data.table(ds_id = ds$ds_id, cell_id = ds$cell_id, alpha_prior = jobs$p[j], param_variant = jobs$param[j], alpha_true = tru$alpha), r)
  }, mc.cores = N_WORKERS, mc.preschedule = FALSE), fill = TRUE)
  fwrite(MC, mc_f); prog(TRUE, "[sim02] MCMC done: %d chains", nrow(MC))
}

# ---- collect and summarise -------------------------------------------------------------------------------------------
R <- rbindlist(lapply(file.path(OUT, "fits", paste0(DS$ds_id, ".rds")), readRDS), fill = TRUE)
R[, `:=`(rel_bias = (a_mean - alpha_true) / alpha_true, rel_bias_med = (a_q500 - alpha_true) / alpha_true,
         covered = a_q025 <= alpha_true & alpha_true <= a_q975, failed = !is.na(error) | is.na(a_mean))]
R[, alpha_sd_ratio := a_sd / (0.78 * a_mean / sqrt(n_events))]
fwrite(R, file.path(OUT, "alpha_results.csv"))
sm <- function(x) x[, .(n = .N, failed = sum(failed), screen_fail = sum(!failed & alpha_sd_ratio < 0.1, na.rm = TRUE), n_coef = round(mean(n_coef), 1),
                        events = round(mean(n_events)), median_rel_bias_pct = round(100 * median(rel_bias, na.rm = TRUE), 1),
                        mean_rel_bias_pct = round(100 * mean(rel_bias, na.rm = TRUE), 1), median_rel_bias_postmedian_pct = round(100 * median(rel_bias_med, na.rm = TRUE), 1),
                        coverage = round(mean(covered[!failed]), 2), median_sd = signif(median(a_sd, na.rm = TRUE), 3))]
S <- R[, sm(.SD), keyby = .(cell_id, group, variant)]; fwrite(S, file.path(OUT, "alpha_summary.csv"))
So <- R[, sm(.SD), keyby = .(group, variant)]; fwrite(So, file.path(OUT, "alpha_summary_overall.csv"))
cat("\nalpha: relative bias and coverage per variant (all cells pooled)\n"); print(So, nrows = 100)
if (file.exists(mc_f)) { MC <- fread(mc_f); if (nrow(MC)) { cat("\nMCMC (exact posterior) vs INLA ref / P_aflat and MAP on the same datasets\n")
  cmp <- merge(MC, dcast(R[variant %in% c("ref", "P_aflat", "C_variant0", "MAP_ref"), .(ds_id, variant, a_mean)], ds_id ~ variant, value.var = "a_mean"), by = "ds_id")
  print(cmp[, .(ds_id, param_variant, alpha_prior, inla_variant0 = round(C_variant0, 3), alpha_true = round(alpha_true, 3), mcmc_mean = round(a_mean, 3), mcmc_ci = sprintf("%.2f-%.2f", a_q025, a_q975), accept = round(accept, 2), ess = round(ess),
                inla_ref = round(ref, 3), inla_flat = round(P_aflat, 3), map = round(MAP_ref, 3))])
  fwrite(cmp, file.path(OUT, "mcmc_vs_inla.csv")) } }
logl("outputs: alpha_results.csv (%d rows), alpha_summary.csv, alpha_summary_overall.csv, mcmc.csv", nrow(R))
logl("==== done")
