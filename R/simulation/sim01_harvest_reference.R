# ==============================================================================
# R/simulation/sim01_harvest_reference.R   (2026-10-30)   RUN ON THE MACHINE WITH INLA
#
# ONE script, two tests. Questions (agreed 2026-10-30):
#  Q1  When few transitions happen WITHOUT harvest (sparse reference level), are the harvest estimates
#      (light, heavy; with their GDD interaction and TSD term per level) and their INLA 95% intervals trustworthy?
#  Q2  In the same situation, does LEAVING HARVEST OUT bias GDD / CMI / soil / alpha?
#  Q3  Does a N(0, SIGMA_B^2) prior on the coefficients change the real-data estimates?
#
# TEST 1 (simulation, Q1-Q2). Real design, exact times:
#  - Design = the real fitting rows of one stored imputation (DESIGN_M) of the main run, for the template transition
#    (from-state = template origin): every sojourn in that state with its segments (truncation, Tstop, covariates).
#    The real sojourn ends (exits to other states, last visit) are kept as censoring, so the amount and pattern of
#    status-0 rows (time at risk without the tested transition) is the real one. With exact times the likelihood of
#    the tested transition sees competing exits only as censoring, so they need not be simulated.
#  - A new exit time for the tested transition is drawn per sojourn from the TRUE piecewise Weibull hazard over its
#    segments (variant 1: S(t) = exp(-(lambda t)^alpha), t in decades, lambda = exp(eta)); a sojourn whose draw falls
#    beyond its real end is censored there.
#  - Truth: alpha and the non-harvest terms = the template's pooled main-run estimates; harvest is the 3-level coding
#    (light, heavy), each level with its own GDD interaction and its own TSD term (tsd_harv x level), true values from
#    the grid below (given on the log-hazard-ratio scale, beta = logHR / alpha). The intercept is solved so that the
#    EXPECTED number of events without harvest equals the d0 target (exact formula, no simulation needed).
#  - Fits per dataset: A harvest in, default coefficient prior (INLA: precision 0.001); B harvest in, N(0, SIGMA_B^2);
#    C harvest out, N(0, SIGMA_B^2). All: pc.alphaw(ALPHA_PRIOR) on alpha, one thread, same inla() call as the pipeline
#    (05a fit_inla_models_pc) plus control.fixed. The simulation's model_fitting.R cannot set either prior, so its
#    call is reproduced here; nothing in the simulation folder or the pipeline is changed.
# TEST 2 (real data, Q3): every kept transition, stored imputations REAL_M, default vs N(0, SIGMA_B^2). The default
#  refit is also compared with the main run's stored estimates (reproduction check).
#
# REPRODUCIBILITY: every dataset has its own seed (dataset_seed()); its simulated outcome (event row + time per sojourn)
#  is saved, and regenerate_dataset() rebuilds the exact fitting data from design + outcome; the script checks that
#  regenerating from the seed gives the identical outcome. Designs saved with md5. Checkpoint per cell (resumable).
#  Versions, settings and input md5s in run_log.txt.
# OUTPUT: data-output/simulation/sim01_harvest_reference/  (FIT_ENGINE "mle" = local check without INLA -> folder suffix)
#   design/ calibration.csv datasets_index.csv outcomes/ cells/ sim_estimates.csv sim_summary.csv sim_summary_by_d0.csv
#   sim_trends.pdf real/ real_prior_shift.csv real_prior_shift_summary.csv run_log.txt
# Run from the project root:  Rscript R/simulation/sim01_harvest_reference.R    (re-running skips finished cells)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(tools) })
source(here("R", "config.R"))
options(width = 200)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# ---- settings -------------------------------------------------------------------------------------------------
FIT_ENGINE  <- Sys.getenv("SIM01_ENGINE", "inla")     # "inla" (the real run); "mle" = local check of the machinery only
RUN_PARTS   <- strsplit(Sys.getenv("SIM01_PARTS", "sim,real"), ",")[[1]]
MASTER_SEED <- 20261030L
N_REP       <- as.integer(Sys.getenv("SIM01_NREP", "30"))      # env overrides exist only for quick local checks
N_WORKERS   <- as.integer(Sys.getenv("SIM01_WORKERS", "6"))
DESIGN_M    <- 1L                                      # stored main-run imputation used as the design
TEMPLATES   <- c("sugar_maple->regen", "regen->yellow_birch", "black_spruce->fir")
# True log hazard ratio vs no harvest (at mean GDD, TSD 0). "real" = the template's own crude log rate ratio per level on
# the design data (events / time at risk, real outcomes), so the real combination of few events without harvest and
# many with harvest is reproduced. Add  large = c(light = 1.5, heavy = 3.0)  for a fixed large effect (+1/3 run time).
EFFECTS     <- list(real     = "real",
                    small    = c(light =  0.2, heavy =  0.3),
                    negative = c(light = -0.5, heavy = -1.0))
D0_TARGETS  <- c(0.5, 2, 5, 10, 30)                    # EXPECTED events without harvest
D0_EXTRA    <- list("black_spruce->fir" = 100)         # extra target per template: the control's real reference (96 events)
TRUE_GDD_X_HARV <- 0                                   # true log-HR of the GDD x harvest-level interactions (per SD GDD)
TRUE_TSD_HARV   <- 0                                   # true log-HR change per decade since harvest, per level
SIGMA_B     <- 5                                       # prior SD of the coefficients in fits B and C
PREC_DEFAULT <- 0.001                                  # INLA's default coefficient precision (fit A)
REAL_M      <- 1:3                                     # stored imputations for test 2
HARV_TERMS  <- c("harv_light", "harv_heavy", "gdd_x_harv_light", "gdd_x_harv_heavy", "tsd_harv_light", "tsd_harv_heavy")
H_CAP       <- -log(0.001 * 2^-23)                     # 22.85: INLA clamps 1 - F at 0.001*FLT_EPSILON (inla-likelihood.c)
ARMS <- data.table(arm = c("A", "B", "C"), harvest = c(TRUE, TRUE, FALSE), beta_prec = c(PREC_DEFAULT, 1 / SIGMA_B^2, 1 / SIGMA_B^2))

OUT <- here("data-output", "simulation", paste0("sim01_harvest_reference", if (FIT_ENGINE != "inla") paste0("_", FIT_ENGINE, "_check")))
for (s in c("design", "outcomes", "cells", "real")) dir.create(file.path(OUT, s), recursive = TRUE, showWarnings = FALSE)
SIM_M <- 200L; M_MI <- 30L; N_COV_SAMPLES <- 0L
source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
if (FIT_ENGINE == "inla") suppressPackageStartupMessages(library(INLA))
logf <- file.path(OUT, "run_log.txt")
logl <- function(...) { s <- sprintf(...); cat(s, "\n"); cat(format(Sys.time()), s, "\n", file = logf, append = TRUE) }
# progress tracker: progress.log gets one line per finished dataset (3 fits) and per finished cell / transition;
# the console gets the per-cell / per-imputation lines with an estimated time remaining. Follow it with:  tail -f <OUT>/progress.log
progf <- file.path(OUT, "progress.log")
prog <- function(console, ...) { s <- sprintf(...); if (console) cat(s, "\n"); cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), s, "\n", file = progf, append = TRUE) }
eta_text <- function(t_start, n_done_run, n_left) {
  if (n_done_run < 1) return("ETA: after the first one finishes")
  el <- as.numeric(difftime(Sys.time(), t_start, units = "mins")); rem <- el / n_done_run * n_left
  sprintf("elapsed %.1f min, ~%.0f min left (finish ~%s)", el, rem, format(Sys.time() + rem * 60, "%H:%M"))
}

# ---- run log: versions, settings, input md5 ---------------------------------------------------------------------
self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (!length(self)) self <- here("R", "simulation", "sim01_harvest_reference.R")
inputs <- c(file.path(TABLE_DIR, c("visit_intervals.rds", "plot_info.rds", "transition_set.rds")),
            file.path(FIT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, sort(unique(c(DESIGN_M, REAL_M))))),
            file.path(FIT_DIR, sprintf("mi_raw_%s.csv", RUN_TAG)), file.path(RESULTS_DIR, sprintf("pooled_%s_final.csv", RUN_TAG)))
logl("==== start | engine %s | R %s | data.table %s | INLA %s", FIT_ENGINE, getRversion(), packageVersion("data.table"),
     if (FIT_ENGINE == "inla") as.character(packageVersion("INLA")) else "not used")
logl("script %s md5 %s", self, md5sum(self))
for (f in inputs) logl("input %s md5 %s", sub(here(), "", f, fixed = TRUE), md5sum(f))
logl("settings: MASTER_SEED %d N_REP %d DESIGN_M %d SIGMA_B %g ALPHA_PRIOR %g TRUE_GDD_X_HARV %g TRUE_TSD_HARV %g H_CAP %.3f",
     MASTER_SEED, N_REP, DESIGN_M, SIGMA_B, ALPHA_PRIOR, TRUE_GDD_X_HARV, TRUE_TSD_HARV, H_CAP)
logl("templates %s | d0 targets %s | effects %s", paste(TEMPLATES, collapse = ", "), paste(D0_TARGETS, collapse = ", "),
     paste(vapply(names(EFFECTS), function(n) if (identical(EFFECTS[[n]], "real")) "real(template-specific, see calibration.csv)" else
       sprintf("%s(light %g, heavy %g)", n, EFFECTS[[n]][["light"]], EFFECTS[[n]][["heavy"]]), character(1)), collapse = ", "))
capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))

# ---- shared inputs ------------------------------------------------------------------------------------------------
inp <- prepare_fit_inputs(TABLE_DIR)
pooled <- fread(file.path(RESULTS_DIR, sprintf("pooled_%s_final.csv", RUN_TAG)))
key_iv <- paste(inp$iv$id_pe, inp$iv$interval); b_abs <- inp$iv$t_end
rebuild <- function(m) {                                   # identical to R/tests/08f_failed_vs_valid.R
  ex <- readRDS(file.path(FIT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, m))); p1 <- as.data.table(ex$p1); p2 <- as.data.table(ex$p2); rm(ex)
  ex_row <- rep(NA_real_, nrow(inp$iv)); ii <- match(paste(p2$id_pe, p2$interval), key_iv); stopifnot(!anyNA(ii)); ex_row[ii] <- p2$exit_abs
  p2f <- function(a_dur, b_dur, from, to, r) { E <- b_abs[r] - b_dur; min(max(ex_row[r] - E, a_dur), b_dur) }
  build_rows(inp$iv, setNames(p1$entry, as.character(p1$id_pe)), inp$states, inp$tr, p2_fun = p2f, keep_cols = inp$keep_cols)$rows[kept == TRUE]
}

# ---- fitting (engine-specific) -------------------------------------------------------------------------------------
fit_inla <- function(d, terms, beta_prec) {
  f <- as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", paste(c("1", terms), collapse = " + ")))
  tryCatch(INLA::inla(f, family = "weibullsurv", data = d,
                      control.family = list(variant = 1, hyper = list(alpha = list(prior = "pc.alphaw", param = ALPHA_PRIOR))),
                      control.fixed = list(prec = beta_prec, prec.intercept = 0),
                      control.compute = list(config = TRUE),
                      control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)),
                      num.threads = "1:1"),
           error = function(e) list(error = conditionMessage(e)))
}
# MAP + Wald intervals with the same coefficient prior (no alpha prior): ONLY to check the simulation machinery without INLA.
fit_mle <- function(d, terms, beta_prec) {
  d <- d[complete.cases(d[, c("status", "Tstop", "truncation", terms), with = FALSE])]
  X <- cbind(1, as.matrix(d[, ..terms])); T <- d$Tstop; tr <- d$truncation; s <- d$status; lT <- log(T); ltr <- ifelse(tr > 0, log(pmax(tr, 1e-300)), 0)
  pen <- function(th) if (length(terms)) beta_prec / 2 * sum(th[-(1:2)]^2) else 0
  nll <- function(th) { a <- exp(th[1]); eta <- drop(X %*% th[-1]); -(sum(s * (th[1] + a * eta + (a - 1) * lT)) - sum(exp(a * eta) * (T^a - tr^a))) + pen(th) }
  gr <- function(th) { a <- exp(th[1]); eta <- drop(X %*% th[-1]); e <- exp(a * eta); Ta <- T^a; tra <- tr^a
    gb <- -(colSums(X * (s * a)) - colSums(X * (a * e * (Ta - tra)))); if (length(terms)) gb[-1] <- gb[-1] + beta_prec * th[-(1:2)]
    ga <- -(sum(s * (1 / a + eta + lT)) - sum(e * (eta * (Ta - tra) + Ta * lT - tra * ltr))); c(ga * a, gb) }
  o <- tryCatch(optim(c(0, log(max(sum(s), 0.5) / sum(T - tr)), rep(0, length(terms))), nll, gr, method = "BFGS", control = list(maxit = 5000, reltol = 1e-12)),
                error = function(e) NULL)
  if (is.null(o)) return(list(error = "optim failed"))
  V <- tryCatch(solve(optimHess(o$par, nll, gr)), error = function(e) NULL)
  se <- if (is.null(V)) rep(NA_real_, length(o$par)) else sqrt(pmax(diag(V), 0))
  est <- c(exp(o$par[1]), o$par[-1]); sdv <- c(exp(o$par[1]) * se[1], se[-1])
  lo <- c(exp(o$par[1] - 1.96 * se[1]), o$par[-1] - 1.96 * se[-1]); hi <- c(exp(o$par[1] + 1.96 * se[1]), o$par[-1] + 1.96 * se[-1])
  list(tab = data.table(parameter = c("alpha", "(Intercept)", terms), mean = est, sd = sdv, q025 = lo, q500 = est, q975 = hi), conv = o$convergence)
}
fit_one <- function(d, terms, beta_prec) if (FIT_ENGINE == "inla") fit_inla(d, terms, beta_prec) else fit_mle(d, terms, beta_prec)

# Tidy result: parameters on the fitted scale; intercept shifted to YEARS (as the pipeline's tidy_fit).
tidy_one <- function(m, terms) {
  pars <- c("alpha", "(Intercept)", terms)
  na <- data.table(parameter = pars, mean = NA_real_, sd = NA_real_, q025 = NA_real_, q500 = NA_real_, q975 = NA_real_)
  if (!is.null(m$error)) return(na[, `:=`(error = substr(as.character(m$error)[1], 1, 150), mode_status = NA_integer_)][])
  if (FIT_ENGINE == "inla") {
    cols <- c("mean", "sd", "0.025quant", "0.5quant", "0.975quant"); fx <- m$summary.fixed; hy <- m$summary.hyperpar
    g <- function(nm) if (nm == "alpha") unlist(hy[ALPHA_NAME, cols]) else if (nm %in% rownames(fx)) unlist(fx[nm, cols]) else rep(NA_real_, 5)
    v <- t(vapply(pars, g, numeric(5)))
    tab <- data.table(parameter = pars, mean = v[, 1], sd = v[, 2], q025 = v[, 3], q500 = v[, 4], q975 = v[, 5])
    ms <- tryCatch(as.integer(m$mode$mode.status)[1], error = function(e) NA_integer_)
  } else { tab <- copy(m$tab); ms <- as.integer(m$conv) }
  tab[parameter == "(Intercept)", c("mean", "q025", "q500", "q975") := lapply(.SD, function(x) x - LOG_TIME_SCALE), .SDcols = c("mean", "q025", "q500", "q975")]
  tab[, `:=`(error = NA_character_, mode_status = ms)][]
}

# ==================================================================================================================
# TEST 1: simulation
# ==================================================================================================================
base_terms <- function(lb) {
  k <- match(lb, inp$td$label); stopifnot(!is.na(k))
  c(paste0("z_", COV_CONT), if (inp$td$fire_term[k]) c("fire_f", "gdd_x_fire", "tsd_fire"), if (inp$td$pest_term[k]) c("pest_f", "gdd_x_pest", "tsd_pest"))
}
make_design <- function(rows, lb) {
  ft <- strsplit(lb, "->", fixed = TRUE)[[1]]
  d <- rows[from == ft[1] & to == ft[2]]; setorder(d, id_pe, sojourn, seg_id)
  d[, `:=`(tsd_harv_light = tsd_harv * harv_light, tsd_harv_heavy = tsd_harv * harv_heavy)]
  base <- base_terms(lb); cols <- unique(c(base, HARV_TERMS))
  bad <- d[, .(bad = anyNA(.SD)), by = .(id_pe, sojourn), .SDcols = cols][bad == TRUE]
  if (nrow(bad)) d <- d[!bad, on = .(id_pe, sojourn)]
  d[, soj := .GRP, by = .(id_pe, sojourn)]
  d[, gap := truncation - data.table::shift(Tstop), by = soj]
  n_gap <- d[!is.na(gap) & abs(gap) > 1e-5, .N]
  if (n_gap > 0) stop(sprintf("design %s: %d segments not contiguous within their sojourn", lb, n_gap))
  D <- d[, c("id_pe", "sojourn", "seg_id", "soj", "truncation", "Tstop", "harv", cols), with = FALSE]
  rr <- d[, .(events = sum(status), decades = sum(Tstop - truncation)), keyby = harv][, rate := events / decades]
  real_loghr <- c(light = log(rr[harv == 1, rate] / rr[harv == 0, rate]), heavy = log(rr[harv == 2, rate] / rr[harv == 0, rate]))
  list(D = D, label = lb, base = base, n_soj = max(D$soj), n_rows = nrow(D), n_soj_na_dropped = nrow(bad), real_events = sum(d$status),
       exposure = rr, real_loghr = real_loghr)
}
truth_for <- function(lb, eff, des) {
  if (identical(eff, "real")) { eff <- des$real_loghr; if (!all(is.finite(eff))) stop(sprintf("%s: real log rate ratio not finite (a level without events)", lb)) }
  p <- pooled[label == lb]; a <- p[parameter == "alpha", pooled_est]; base <- base_terms(lb)
  b <- setNames(p[match(base, parameter), pooled_est], base); stopifnot(!anyNA(b))
  h <- c(harv_light = eff[["light"]], harv_heavy = eff[["heavy"]], gdd_x_harv_light = TRUE_GDD_X_HARV, gdd_x_harv_heavy = TRUE_GDD_X_HARV,
         tsd_harv_light = TRUE_TSD_HARV, tsd_harv_heavy = TRUE_TSD_HARV) / a       # log HR -> beta
  list(alpha = a, beta = c(b, h))
}
# per-row cumulative hazard of a segment WITHOUT the intercept: exp(alpha * (log 10 + x'beta)) * (Tstop^a - trunc^a)
seg_parts <- function(D, tru) {
  X <- as.matrix(D[, names(tru$beta), with = FALSE]); eta0 <- drop(X %*% tru$beta) + LOG_TIME_SCALE
  list(eta0 = eta0, b0 = exp(tru$alpha * eta0) * (D$Tstop^tru$alpha - D$truncation^tru$alpha))
}
expected_counts <- function(c_yr, D, sp, alpha) {
  dH <- exp(alpha * c_yr) * sp$b0; cs <- ave(dH, D$soj, FUN = cumsum); p <- exp(-(cs - dH)) * (-expm1(-dH))
  c(d0 = sum(p[D$harv == 0]), d_light = sum(p[D$harv == 1]), d_heavy = sum(p[D$harv == 2]))
}
calibrate <- function(target, D, sp, alpha) {
  f <- function(c) log(expected_counts(c, D, sp, alpha)[["d0"]]) - log(target)
  r <- tryCatch(uniroot(f, c(-25, 0), extendInt = "upX", tol = 1e-10)$root, error = function(e) NA_real_)
  r
}
dataset_seed <- function(ti, ei, di, rep) MASTER_SEED + 100000L * ti + 10000L * ei + 1000L * di + rep
simulate_outcome <- function(D, sp, alpha, c_yr, seed) {
  set.seed(seed); n <- max(D$soj); u <- rexp(n)
  dH <- exp(alpha * c_yr) * sp$b0; cs <- ave(dH, D$soj, FUN = cumsum); prev <- cs - dH
  hit <- which(cs >= u[D$soj]); hit <- hit[!duplicated(D$soj[hit])]          # first segment where the drawn exit falls
  if (!length(hit)) return(data.table(soj = integer(0), row = integer(0), t = numeric(0)))
  rem <- u[D$soj[hit]] - prev[hit]; lam_a <- exp(alpha * (c_yr + sp$eta0[hit]))
  t <- (D$truncation[hit]^alpha + rem / lam_a)^(1 / alpha)
  t <- pmin(pmax(t, D$truncation[hit] + EPSILON), D$Tstop[hit])
  data.table(soj = D$soj[hit], row = hit, t = t)
}
regenerate_dataset <- function(D, outcome) {
  cut <- rep(.Machine$integer.max, max(D$soj)); cut[outcome$soj] <- outcome$row
  keep <- seq_len(nrow(D)) <= cut[D$soj]
  d <- D[keep]; d[, status := 0L]; ri <- match(outcome$row, which(keep))
  d[ri, `:=`(status = 1L, Tstop = outcome$t)]
  d[]
}
H_at <- function(d, tru, c_yr) { X <- as.matrix(d[, names(tru$beta), with = FALSE]); la <- exp(tru$alpha * (c_yr + LOG_TIME_SCALE + drop(X %*% tru$beta)))
  list(tr = la * d$truncation^tru$alpha, T = la * d$Tstop^tru$alpha) }

run_cell <- function(cell) {
  f <- file.path(OUT, "cells", sprintf("%s.rds", cell$cell_id)); if (file.exists(f)) return(invisible(NULL))
  des <- readRDS(file.path(OUT, "design", sprintf("design_%s.rds", gsub("->", "_to_", cell$template)))); D <- des$D
  tru <- truth_for(cell$template, EFFECTS[[cell$effect]], des); sp <- seg_parts(D, tru); c_yr <- cell$c_yr
  truth_tab <- data.table(parameter = c("alpha", "(Intercept)", names(tru$beta)), truth = c(tru$alpha, c_yr, tru$beta))
  res <- vector("list", N_REP); idx <- vector("list", N_REP)
  for (r in seq_len(N_REP)) {
    seed <- dataset_seed(cell$ti, cell$ei, cell$di, r)
    oc <- simulate_outcome(D, sp, tru$alpha, c_yr, seed)
    saveRDS(oc, file.path(OUT, "outcomes", sprintf("%s_rep%02d.rds", cell$cell_id, r)))
    d <- regenerate_dataset(D, oc); H <- H_at(d, tru, c_yr)
    ev <- d[status == 1L]
    info <- data.table(cell_id = cell$cell_id, rep = r, seed = seed, n_rows = nrow(d), n_events = nrow(ev), d0 = ev[harv == 0, .N],
                       d_light = ev[harv == 1, .N], d_heavy = ev[harv == 2, .N], maxH_tr = max(H$tr), maxH_T = max(H$T),
                       n_rows_above_cap = sum(H$tr > H_CAP | H$T > H_CAP), outcome_md5 = digest_md5(oc))
    idx[[r]] <- info
    out <- list()
    for (j in seq_len(nrow(ARMS))) {
      terms <- if (ARMS$harvest[j]) c(des$base, HARV_TERMS) else des$base
      t0 <- proc.time()[["elapsed"]]; m <- fit_one(d, terms, ARMS$beta_prec[j]); sec <- proc.time()[["elapsed"]] - t0
      tb <- tidy_one(m, terms); rm(m)
      out[[j]] <- cbind(info[, .(cell_id, rep, seed, n_events, d0, d_light, d_heavy, maxH_tr, n_rows_above_cap)], arm = ARMS$arm[j], sec = sec, tb)
    }
    res[[r]] <- merge(rbindlist(out), truth_tab, by = "parameter", all.x = TRUE, sort = FALSE)
    prog(FALSE, "  %s rep %d/%d: events %d (d0 %d, light %d, heavy %d), fits %s s, errors %d", cell$cell_id, r, N_REP, info$n_events, info$d0,
         info$d_light, info$d_heavy, paste(round(sapply(out, function(o) o$sec[1]), 1), collapse = "/"), sum(sapply(out, function(o) !is.na(o$error[1]))))
  }
  R <- rbindlist(res); I <- rbindlist(idx)
  R <- cbind(cell[rep(1, nrow(R)), .(template, effect, d0_target)], R)
  saveRDS(list(est = R, index = I), f)
  n_done <- sum(file.exists(file.path(OUT, "cells", paste0(CELL_IDS, ".rds"))))
  prog(TRUE, "[sim] cell %d/%d done: %s | fit errors %d of %d | %s", n_done, length(CELL_IDS), cell$cell_id, sum(sapply(res, function(x) sum(!is.na(x$error[x$parameter == "alpha"])))),
       N_REP * nrow(ARMS), eta_text(PROG_T0, n_done - N_DONE0, length(CELL_IDS) - n_done))
  invisible(NULL)
}
digest_md5 <- function(x) { tf <- tempfile(); on.exit(unlink(tf)); saveRDS(x, tf, compress = FALSE); unname(md5sum(tf)) }

if ("sim" %in% RUN_PARTS) {
  logl("---- TEST 1: simulation")
  rows_d <- rebuild(DESIGN_M)
  for (lb in TEMPLATES) {
    f <- file.path(OUT, "design", sprintf("design_%s.rds", gsub("->", "_to_", lb)))
    if (!file.exists(f)) saveRDS(make_design(rows_d, lb), f)
    des <- readRDS(f)
    logl("design %s: %d sojourns, %d rows, %d sojourns dropped for NA covariates, md5 %s; real events / decades at risk by harvest level: %s; real crude log HR light %.2f heavy %.2f",
         lb, des$n_soj, des$n_rows, des$n_soj_na_dropped, md5sum(f),
         paste(sprintf("%s: %d / %.0f", des$exposure$harv, des$exposure$events, des$exposure$decades), collapse = ", "), des$real_loghr[["light"]], des$real_loghr[["heavy"]])
  }
  rm(rows_d); invisible(gc())
  # calibration: intercept (years) so that E[events without harvest] = d0 target
  cal_f <- file.path(OUT, "calibration.csv")
  if (!file.exists(cal_f)) {
    cal <- rbindlist(lapply(seq_along(TEMPLATES), function(ti) { lb <- TEMPLATES[ti]; des <- readRDS(file.path(OUT, "design", sprintf("design_%s.rds", gsub("->", "_to_", lb))))
      rbindlist(lapply(seq_along(EFFECTS), function(ei) { tru <- truth_for(lb, EFFECTS[[ei]], des); sp <- seg_parts(des$D, tru)
        tg <- c(D0_TARGETS, D0_EXTRA[[lb]])
        rbindlist(lapply(seq_along(tg), function(di) { c_yr <- calibrate(tg[di], des$D, sp, tru$alpha)
          ec <- if (is.na(c_yr)) c(d0 = NA, d_light = NA, d_heavy = NA) else expected_counts(c_yr, des$D, sp, tru$alpha)
          data.table(template = lb, effect = names(EFFECTS)[ei], d0_target = tg[di], ti = ti, ei = ei, di = di, c_yr = c_yr,
                     alpha = tru$alpha, loghr_light = tru$alpha * tru$beta[["harv_light"]], loghr_heavy = tru$alpha * tru$beta[["harv_heavy"]],
                     E_d0 = ec[["d0"]], E_light = ec[["d_light"]], E_heavy = ec[["d_heavy"]]) })) })) }))
    cal[, cell_id := sprintf("%s__%s__d0_%s", gsub("->", "_to_", template), effect, d0_target)]
    fwrite(cal, cal_f)
  }
  cal <- fread(cal_f)
  if (cal[is.na(c_yr), .N]) logl("WARNING: calibration failed for %d cells (skipped): %s", cal[is.na(c_yr), .N], paste(cal[is.na(c_yr), cell_id], collapse = ", "))
  print(cal[, .(template, effect, d0_target, c_yr = round(c_yr, 3), E_d0 = round(E_d0, 2), E_light = round(E_light, 1), E_heavy = round(E_heavy, 1))])
  # reproducibility check: the same seed gives the identical outcome
  c1 <- cal[!is.na(c_yr)][1]; des1 <- readRDS(file.path(OUT, "design", sprintf("design_%s.rds", gsub("->", "_to_", c1$template))))
  tru1 <- truth_for(c1$template, EFFECTS[[c1$effect]], des1); sp1 <- seg_parts(des1$D, tru1); s1 <- dataset_seed(c1$ti, c1$ei, c1$di, 1L)
  stopifnot(identical(simulate_outcome(des1$D, sp1, tru1$alpha, c1$c_yr, s1), simulate_outcome(des1$D, sp1, tru1$alpha, c1$c_yr, s1)))
  logl("reproducibility check passed (same seed -> identical outcome)")
  todo <- cal[!is.na(c_yr)][!file.exists(file.path(OUT, "cells", paste0(cell_id, ".rds")))]
  logl("cells: %d total, %d to run (%d datasets x %d fits each)", cal[!is.na(c_yr), .N], nrow(todo), N_REP, nrow(ARMS))
  CELL_IDS <- cal[!is.na(c_yr), cell_id]; N_DONE0 <- length(CELL_IDS) - nrow(todo); PROG_T0 <- Sys.time()   # read by the workers
  prog(TRUE, "[sim] start: %d of %d cells to run on %d workers (progress per dataset in %s)", nrow(todo), length(CELL_IDS), N_WORKERS, progf)
  t0 <- proc.time()[["elapsed"]]
  err <- mclapply(split(todo, seq_len(nrow(todo))), function(cell) tryCatch({ run_cell(cell); NULL }, error = function(e) paste(cell$cell_id, conditionMessage(e))),
                  mc.cores = N_WORKERS, mc.preschedule = FALSE)
  err <- unlist(err); if (length(err)) { logl("ERRORS in cells:\n%s", paste(err, collapse = "\n")); stop("cell errors, see run_log.txt") }
  logl("simulation fits done in %.1f min", (proc.time()[["elapsed"]] - t0) / 60)

  # ---- collect and summarise -------------------------------------------------------------------------------------
  cells <- lapply(file.path(OUT, "cells", paste0(cal[!is.na(c_yr), cell_id], ".rds")), readRDS)
  E <- rbindlist(lapply(cells, `[[`, "est")); I <- rbindlist(lapply(cells, `[[`, "index"))
  fwrite(I, file.path(OUT, "datasets_index.csv"))
  E[, `:=`(err = mean - truth, covered = !is.na(q025) & q025 <= truth & truth <= q975, ci_width = q975 - q025, failed = !is.na(error) | is.na(mean))]
  al <- E[parameter == "alpha", .(cell_id, rep, arm, alpha_est = mean, alpha_sd = sd)]
  E <- merge(E, al, by = c("cell_id", "rep", "arm"), all.x = TRUE)
  E[, alpha_sd_ratio := alpha_sd / (0.78 * alpha_est / sqrt(pmax(n_events, 1)))]            # the pipeline's frozen screen quantity
  E[!parameter %in% c("alpha", "(Intercept)"), `:=`(loghr_est = alpha_est * mean, loghr_truth = cal$alpha[match(cell_id, cal$cell_id)] * truth)]
  fwrite(E, file.path(OUT, "sim_estimates.csv"))
  smry <- function(x) x[, .(n_fits = .N, n_failed = sum(failed), n_screen_fail = sum(!failed & alpha_sd_ratio < 0.1, na.rm = TRUE),
                            d0_mean = mean(d0), d_light_mean = mean(d_light), d_heavy_mean = mean(d_heavy), truth = truth[1],
                            mean_bias = mean(err, na.rm = TRUE), median_bias = median(err, na.rm = TRUE), sd_est = sd(mean, na.rm = TRUE),
                            rmse = sqrt(mean(err^2, na.rm = TRUE)), err_q05 = quantile(err, 0.05, na.rm = TRUE), err_q95 = quantile(err, 0.95, na.rm = TRUE),
                            coverage = mean(covered[!failed]), mean_ci_width = mean(ci_width, na.rm = TRUE),
                            mean_loghr_bias = mean(loghr_est - loghr_truth, na.rm = TRUE), n_rows_above_cap = sum(n_rows_above_cap > 0))]
  S <- E[, smry(.SD), keyby = .(template, effect, d0_target, arm, parameter)]
  fwrite(S, file.path(OUT, "sim_summary.csv"))
  # same summary WITHOUT datasets that contain a row above INLA's cap (H >= 22.85 at the truth): there the fit is
  # affected by the known INLA clamp, not by the question tested here
  S_cap <- E[n_rows_above_cap == 0, smry(.SD), keyby = .(template, effect, d0_target, arm, parameter)]
  fwrite(S_cap, file.path(OUT, "sim_summary_below_cap.csv"))
  logl("datasets with rows above the INLA cap (H >= %.2f): %d of %d", H_CAP, I[n_rows_above_cap > 0, .N], nrow(I))
  E[, d0_bin := cut(d0, c(-1, 0, 1, 4, 9, 19, Inf), labels = c("0", "1", "2-4", "5-9", "10-19", "20+"))]
  S2 <- E[, smry(.SD), keyby = .(template, effect, d0_bin, arm, parameter)]
  fwrite(S2, file.path(OUT, "sim_summary_by_d0.csv"))

  # ---- trend figure: bias and coverage vs d0 target, per template (rows) x effect (columns), lines = fits A/B/C ----
  pdf(file.path(OUT, "sim_trends.pdf"), width = 11, height = 8.5)
  cols <- c(A = "#d95f02", B = "#1b9e77", C = "#7570b3")
  for (par in c("harv_heavy", "harv_light", "gdd_x_harv_heavy", "tsd_harv_heavy", "z_gdd_3", "z_cmi_gs", "z_soil_pc1", "alpha")) {
    for (what in c("mean_bias", "coverage")) {
      par_s <- S[parameter == par]; if (!nrow(par_s)) next
      op <- par(mfrow = c(length(TEMPLATES), length(EFFECTS)), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2.5, 0))
      for (lb in TEMPLATES) for (ef in names(EFFECTS)) {
        z <- par_s[template == lb & effect == ef]
        yl <- if (what == "coverage") c(0, 1) else { v <- z[[what]]; v <- v[is.finite(v)]; if (length(v)) range(c(0, v)) else c(-1, 1) }
        xt <- sort(unique(cal$d0_target))
        plot(NA, xlim = range(log(xt)), ylim = yl, xaxt = "n", xlab = "expected events without harvest (d0 target)", ylab = what,
             main = sprintf("%s | %s", lb, ef), cex.main = 0.9)
        axis(1, at = log(xt), labels = xt); abline(h = if (what == "coverage") 0.95 else 0, lty = 2, col = "grey50")
        for (a in names(cols)) { w <- z[arm == a][order(d0_target)]; if (nrow(w)) lines(log(w$d0_target), w[[what]], type = "b", pch = 19, col = cols[[a]]) }
      }
      mtext(sprintf("%s: %s   (A harvest in, default prior | B harvest in, N(0,%g^2) | C harvest out)", par, what, SIGMA_B), outer = TRUE, cex = 1)
      legend("bottomright", legend = names(cols), col = cols, lty = 1, pch = 19, bty = "n")
      par(op)
    }
  }
  invisible(dev.off())
  logl("TEST 1 outputs written: sim_estimates.csv (%d rows), sim_summary.csv, sim_summary_by_d0.csv, sim_trends.pdf", nrow(E))
}

# ==================================================================================================================
# TEST 2: real data, default coefficient prior vs N(0, SIGMA_B^2)
# ==================================================================================================================
if ("real" %in% RUN_PARTS) {
  logl("---- TEST 2: real data, prior check (imputations %s)", paste(REAL_M, collapse = ", "))
  mi_raw <- fread(file.path(FIT_DIR, sprintf("mi_raw_%s.csv", RUN_TAG)))
  m_todo <- REAL_M[!file.exists(file.path(OUT, "real", sprintf("real_m%02d.rds", REAL_M)))]
  PROG_T0 <- Sys.time(); n_m_done <- 0L
  prog(TRUE, "[real] start: %d of %d imputations to run (%d transitions x 2 fits each, %d workers)", length(m_todo), length(REAL_M), nrow(inp$td), N_WORKERS)
  for (m in m_todo) {
    f <- file.path(OUT, "real", sprintf("real_m%02d.rds", m))
    rows <- rebuild(m); d_all <- rows[kept == TRUE]; d_all[, label := paste(from, to, sep = "->")]
    cols <- c("status", "Tstop", "truncation", inp$terms)
    res <- mclapply(seq_len(nrow(inp$td)), function(k) {
      lb <- inp$td$label[k]; tms <- colnames(inp$mask)[inp$mask[k, ]]; d <- d_all[label == lb, ..cols]
      out <- rbindlist(lapply(c(default = PREC_DEFAULT, sigma = 1 / SIGMA_B^2), function(pr) {
        t0 <- proc.time()[["elapsed"]]; mm <- fit_one(d, tms, pr); tb <- tidy_one(mm, tms)
        tb[, `:=`(label = lb, m = m, prior = if (pr == PREC_DEFAULT) "default" else sprintf("sigma%g", SIGMA_B), n_events = sum(d$status), sec = proc.time()[["elapsed"]] - t0)] }))
      prog(FALSE, "  real m%02d transition %d/%d %s: fits %s s, errors %d", m, k, nrow(inp$td), lb,
           paste(round(out[parameter == "alpha", sec], 1), collapse = "/"), out[parameter == "alpha" & !is.na(error), .N])
      out
    }, mc.cores = N_WORKERS, mc.preschedule = FALSE)
    bad <- vapply(res, inherits, logical(1), "try-error"); if (any(bad)) stop("TEST 2 worker error: ", paste(res[bad], collapse = " | "))
    saveRDS(rbindlist(res), f); rm(rows, d_all); invisible(gc())
    n_m_done <- n_m_done + 1L
    prog(TRUE, "[real] imputation m%02d done (%d/%d) | %s", m, n_m_done, length(m_todo), eta_text(PROG_T0, n_m_done, length(m_todo) - n_m_done))
  }
  Rr <- rbindlist(lapply(file.path(OUT, "real", sprintf("real_m%02d.rds", REAL_M)), readRDS))
  W <- data.table::dcast(Rr, label + m + parameter + n_events ~ prior, value.var = c("mean", "sd", "q025", "q975", "error"))
  setnames(W, function(x) sub(sprintf("sigma%g", SIGMA_B), "sigma", x))
  W <- merge(W, mi_raw[, .(label, m, parameter, est_main = est, sd_main = sd)], by = c("label", "m", "parameter"), all.x = TRUE)
  fv <- fread(file.path(RESULTS_DIR, sprintf("fit_validity_%s.csv", RUN_TAG)))
  W <- merge(W, fv[, .(label, m, main_fit_valid = valid)], by = c("label", "m"), all.x = TRUE)
  W[, `:=`(shift_in_sd = (mean_sigma - mean_default) / sd_default, repro_diff_in_sd = (mean_default - est_main) / sd_main)]
  fwrite(W, file.path(OUT, "real_prior_shift.csv"))
  # the prior comparison is only meaningful where the refit itself is a healthy fit: same frozen alpha-SD screen as 07/08
  W[parameter == "alpha", ok_fit := is.finite(mean_default) & sd_default / (0.78 * mean_default / sqrt(n_events)) >= 0.1]
  W[, ok_fit := any(ok_fit, na.rm = TRUE), by = .(label, m)]
  logl("TEST 2: %d of %d refits (default prior) pass the alpha-SD screen; shifts below use those only", W[parameter == "alpha" & ok_fit == TRUE, .N], W[parameter == "alpha", .N])
  Wc <- W[ok_fit == TRUE & !parameter %in% c("alpha", "(Intercept)") & is.finite(shift_in_sd)]
  sm <- Wc[, .(n = .N, median_abs_shift_sd = median(abs(shift_in_sd)), q95_abs_shift_sd = quantile(abs(shift_in_sd), 0.95),
               max_abs_shift_sd = max(abs(shift_in_sd)), n_shift_gt_0.1 = sum(abs(shift_in_sd) > 0.1), n_shift_gt_0.25 = sum(abs(shift_in_sd) > 0.25))]
  sm2 <- W[ok_fit == TRUE & parameter %in% c("alpha", "(Intercept)") & is.finite(shift_in_sd), .(median_abs_shift_sd = median(abs(shift_in_sd)), max_abs_shift_sd = max(abs(shift_in_sd))), by = parameter]
  rp <- W[main_fit_valid == TRUE & is.finite(repro_diff_in_sd), .(n = .N, median_abs = median(abs(repro_diff_in_sd)), max_abs = max(abs(repro_diff_in_sd)), n_gt_0.01 = sum(abs(repro_diff_in_sd) > 0.01))]
  fwrite(cbind(scope = "coefficients", sm), file.path(OUT, "real_prior_shift_summary.csv"))
  logl("TEST 2 coefficients: %s", paste(names(sm), round(unlist(sm), 3), sep = "=", collapse = " "))
  for (i in seq_len(nrow(sm2))) logl("TEST 2 %s: median |shift| %.3f SD, max %.3f SD", sm2$parameter[i], sm2$median_abs_shift_sd[i], sm2$max_abs_shift_sd[i])
  logl("TEST 2 reproduction of the main run (default prior vs stored estimates, in SDs): %s", paste(names(rp), round(unlist(rp), 4), sep = "=", collapse = " "))
}
logl("==== done")
