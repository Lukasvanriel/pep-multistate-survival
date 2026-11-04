# ==============================================================================
# R/simulation/sim03_variant0.R   (2026-11-01)   RUN ON THE MACHINE WITH INLA, AFTER sim01
#
# QUESTION: should the pipeline switch from INLA weibullsurv variant 1 (S = exp(-(lambda t)^alpha), coefficients
#  multiplied by alpha) to variant 0 (S = exp(-lambda t^alpha), coefficients = log hazard ratios)? sim02 showed that the
#  variant-1 POSTERIOR of alpha is biased low with few events per coefficient (INLA = exact MCMC), and variant 0 is not.
#  Variant 1 was adopted on D. Rustand's advice (r-inla group, 15 Nov 2025: "more stable numerically", from an
#  intercept-only, interval-censored simulation). What is left to know for OUR pipeline:
#  PART A (real data): is variant 0 as stable as variant 1 on the real transitions, and do hazard ratios / alpha differ?
#    All kept transitions x stored main-run imputations REAL_M x variant {0, 1} x coefficient prior {INLA default,
#    N(0, 5^2)}. Same rows, same pc.alphaw(ALPHA_PRIOR), one thread. Measures: INLA errors, frozen alpha-SD screen,
#    alpha, log hazard ratios (variant 1: alpha * beta; variant 0: beta), max cumulative hazard at truncation H(tr) at the
#    posterior means (INLA's cap 22.85 applies to both variants), run time. Variant 1 + default prior = the main run
#    (reproduction check against mi_raw).
#  PART B (simulation): sim01's 1,440 seeded datasets refitted under variant 0 (fits A default prior, B N(0, 5^2),
#    C harvest out N(0, 5^2)), compared with sim01's variant-1 fits ON THE LOG-HAZARD-RATIO SCALE. Redoes Q1 (harvest
#    effects with few events without harvest), Q2 (bias from leaving harvest out) and alpha recovery.
# Truth on the variant-0 scale (exact conversion of sim01's variant-1 truth, time in decades):
#   alpha unchanged; coefficient = alpha * beta_v1 (= log HR); intercept = alpha * (c_yr + log 10).
# READ-ONLY inputs: sim01 output folder, main-run tables/extras/results. Nothing else is changed.
# OUTPUT data-output/simulation/sim03_variant0/: realA/ realA_fits.csv realA_summary.csv realA_transitions.csv
#   cellsB/ simB_estimates.csv simB_summary.csv simB_vs_sim01.csv run_log.txt progress.log
# Run from the project root:  Rscript R/simulation/sim03_variant0.R      (re-running skips finished parts)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(tools) })
source(here("R", "config.R"))
options(width = 200)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# ---- settings -------------------------------------------------------------------------------------------------
FIT_ENGINE <- Sys.getenv("SIM03_ENGINE", "inla")              # "mle" = local plumbing check only
RUN_PARTS  <- strsplit(Sys.getenv("SIM03_PARTS", "A,B"), ",")[[1]]
N_WORKERS  <- as.integer(Sys.getenv("SIM03_WORKERS", "6"))
SIM01_DIR  <- here("data-output", "simulation", Sys.getenv("SIM03_SIM01", "sim01_harvest_reference"))
REAL_M     <- 1:3
SIGMA_B    <- 5; PREC_DEFAULT <- 0.001
HARV_TERMS <- c("harv_light", "harv_heavy", "gdd_x_harv_light", "gdd_x_harv_heavy", "tsd_harv_light", "tsd_harv_heavy")
EFFECTS    <- list(real = "real", small = c(light = 0.2, heavy = 0.3), negative = c(light = -0.5, heavy = -1.0))   # = sim01
ARMS       <- data.table(arm = c("A", "B", "C"), harvest = c(TRUE, TRUE, FALSE), beta_prec = c(PREC_DEFAULT, 1 / SIGMA_B^2, 1 / SIGMA_B^2))
H_CAP      <- -log(0.001 * 2^-23)

OUT <- here("data-output", "simulation", paste0("sim03_variant0", if (FIT_ENGINE != "inla") paste0("_", FIT_ENGINE, "_check")))
for (s in c("realA", "cellsB")) dir.create(file.path(OUT, s), recursive = TRUE, showWarnings = FALSE)
SIM_M <- 200L; M_MI <- 30L; N_COV_SAMPLES <- 0L
source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
if (FIT_ENGINE == "inla") suppressPackageStartupMessages(library(INLA))
logf <- file.path(OUT, "run_log.txt"); progf <- file.path(OUT, "progress.log")
logl <- function(...) { s <- sprintf(...); cat(s, "\n"); cat(format(Sys.time()), s, "\n", file = logf, append = TRUE) }
prog <- function(console, ...) { s <- sprintf(...); if (console) cat(s, "\n"); cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), s, "\n", file = progf, append = TRUE) }
eta_text <- function(t_start, n_done_run, n_left) {
  if (n_done_run < 1) return("ETA: after the first one finishes")
  el <- as.numeric(difftime(Sys.time(), t_start, units = "mins")); rem <- el / n_done_run * n_left
  sprintf("elapsed %.1f min, ~%.0f min left (finish ~%s)", el, rem, format(Sys.time() + rem * 60, "%H:%M"))
}
self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); if (!length(self)) self <- here("R", "simulation", "sim03_variant0.R")
logl("==== start | engine %s | R %s | data.table %s | INLA %s | parts %s", FIT_ENGINE, getRversion(), packageVersion("data.table"),
     if (FIT_ENGINE == "inla") as.character(packageVersion("INLA")) else "not used", paste(RUN_PARTS, collapse = ","))
logl("script %s md5 %s", self, md5sum(self))
for (f in c(file.path(TABLE_DIR, c("visit_intervals.rds", "transition_set.rds")), file.path(FIT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, REAL_M)),
            file.path(FIT_DIR, sprintf("mi_raw_%s.csv", RUN_TAG)), file.path(RESULTS_DIR, sprintf("pooled_%s_final.csv", RUN_TAG)),
            file.path(SIM01_DIR, c("calibration.csv", "datasets_index.csv", "sim_estimates.csv"))))
  if (file.exists(f)) logl("input %s md5 %s", sub(here(), "", f, fixed = TRUE), md5sum(f))
capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))

# ---- fitting: one function for both variants -------------------------------------------------------------------------
# returns alpha, intercept (decades scale, as fitted) and coefficients with mean/sd/quantiles; plus mode status / error
fit_inla_v <- function(d, terms, variant, beta_prec) {
  f <- as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", paste(c("1", terms), collapse = " + ")))
  m <- tryCatch(INLA::inla(f, family = "weibullsurv", data = d,
                           control.family = list(variant = variant, hyper = list(alpha = list(prior = "pc.alphaw", param = ALPHA_PRIOR))),
                           control.fixed = list(prec = beta_prec, prec.intercept = 0), control.compute = list(config = TRUE),
                           control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)),
                           num.threads = "1:1"), error = function(e) list(error = conditionMessage(e)))
  pars <- c("alpha", "(Intercept)", terms)
  if (!is.null(m$error)) return(data.table(parameter = pars, mean = NA_real_, sd = NA_real_, q025 = NA_real_, q500 = NA_real_, q975 = NA_real_,
                                           mode_status = NA_integer_, error = substr(m$error, 1, 150)))
  cols <- c("mean", "sd", "0.025quant", "0.5quant", "0.975quant"); fx <- m$summary.fixed; hy <- m$summary.hyperpar
  g <- function(nm) if (nm == "alpha") unlist(hy[ALPHA_NAME, cols]) else if (nm %in% rownames(fx)) unlist(fx[nm, cols]) else rep(NA_real_, 5)
  v <- t(vapply(pars, g, numeric(5)))
  data.table(parameter = pars, mean = v[, 1], sd = v[, 2], q025 = v[, 3], q500 = v[, 4], q975 = v[, 5],
             mode_status = tryCatch(as.integer(m$mode$mode.status)[1], error = function(e) NA_integer_), error = NA_character_)
}
# plumbing-only stand-in: MAP with the coefficient penalty, Wald intervals, either variant
fit_map_v <- function(d, terms, variant, beta_prec) {
  d <- d[complete.cases(d[, c("status", "Tstop", "truncation", terms), with = FALSE])]
  X <- if (length(terms)) cbind(1, as.matrix(d[, ..terms])) else matrix(1, nrow(d), 1); T <- d$Tstop; tr <- d$truncation; s <- d$status; lT <- log(T)
  nll <- function(th) { la <- th[1]; a <- exp(la); eta <- drop(X %*% th[-1])
    ll <- if (variant == 1) sum(s * (la + a * eta + (a - 1) * lT)) - sum(exp(a * eta) * (T^a - tr^a)) else sum(s * (la + eta + (a - 1) * lT)) - sum(exp(eta) * (T^a - tr^a))
    -(ll - if (length(terms)) beta_prec / 2 * sum(th[-(1:2)]^2) else 0) }
  o <- tryCatch(optim(c(0, log(max(sum(s), 0.5) / sum(T - tr)), rep(0, length(terms))), nll, method = "BFGS", control = list(maxit = 10000, reltol = 1e-12)), error = function(e) NULL)
  pars <- c("alpha", "(Intercept)", terms)
  if (is.null(o)) return(data.table(parameter = pars, mean = NA_real_, sd = NA_real_, q025 = NA_real_, q500 = NA_real_, q975 = NA_real_, mode_status = NA_integer_, error = "optim failed"))
  V <- tryCatch(solve(optimHess(o$par, nll)), error = function(e) NULL); se <- if (is.null(V)) rep(NA_real_, length(o$par)) else sqrt(pmax(diag(V), 0))
  est <- c(exp(o$par[1]), o$par[-1]); sdv <- c(exp(o$par[1]) * se[1], se[-1])
  lo <- c(exp(o$par[1] - 1.96 * se[1]), o$par[-1] - 1.96 * se[-1]); hi <- c(exp(o$par[1] + 1.96 * se[1]), o$par[-1] + 1.96 * se[-1])
  data.table(parameter = pars, mean = est, sd = sdv, q025 = lo, q500 = est, q975 = hi, mode_status = o$convergence, error = NA_character_)
}
fit_v <- function(d, terms, variant, beta_prec) if (FIT_ENGINE == "inla") fit_inla_v(d, terms, variant, beta_prec) else fit_map_v(d, terms, variant, beta_prec)
# max cumulative hazard at truncation and at Tstop, at the posterior means (decades; same row data as the fit)
maxH <- function(d, terms, tb, variant) {
  a <- tb[parameter == "alpha", mean]; b <- tb[parameter != "alpha", mean]; if (anyNA(c(a, b))) return(c(NA_real_, NA_real_))
  X <- cbind(1, as.matrix(d[, ..terms])); X[is.na(X)] <- 0; eta <- drop(X %*% b); la <- if (variant == 1) exp(a * eta) else exp(eta)
  c(max(la * d$truncation^a), max(la * d$Tstop^a))
}

# ==================================================================================================================
# PART A: real data
# ==================================================================================================================
if ("A" %in% RUN_PARTS) {
  logl("---- PART A: real data, variant 0 vs 1 (imputations %s)", paste(REAL_M, collapse = ", "))
  inp <- prepare_fit_inputs(TABLE_DIR)
  key_iv <- paste(inp$iv$id_pe, inp$iv$interval); b_abs <- inp$iv$t_end
  rebuild <- function(m) {                                 # identical to R/tests/08f_failed_vs_valid.R and sim01
    ex <- readRDS(file.path(FIT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, m))); p1 <- as.data.table(ex$p1); p2 <- as.data.table(ex$p2); rm(ex)
    ex_row <- rep(NA_real_, nrow(inp$iv)); ii <- match(paste(p2$id_pe, p2$interval), key_iv); stopifnot(!anyNA(ii)); ex_row[ii] <- p2$exit_abs
    p2f <- function(a_dur, b_dur, from, to, r) { E <- b_abs[r] - b_dur; min(max(ex_row[r] - E, a_dur), b_dur) }
    build_rows(inp$iv, setNames(p1$entry, as.character(p1$id_pe)), inp$states, inp$tr, p2_fun = p2f, keep_cols = inp$keep_cols)$rows[kept == TRUE]
  }
  SPECS <- CJ(variant = c(1L, 0L), prior = c("default", "sigma5"))
  m_todo <- REAL_M[!file.exists(file.path(OUT, "realA", sprintf("m%02d.rds", REAL_M)))]
  t0A <- Sys.time(); nd <- 0L
  prog(TRUE, "[A] start: %d of %d imputations to run (%d transitions x %d fits)", length(m_todo), length(REAL_M), nrow(inp$td), nrow(SPECS))
  for (m in m_todo) {
    rows <- rebuild(m); d_all <- rows[kept == TRUE]; d_all[, label := paste(from, to, sep = "->")]
    cols <- c("status", "Tstop", "truncation", inp$terms)
    res <- mclapply(seq_len(nrow(inp$td)), function(k) {
      lb <- inp$td$label[k]; tms <- colnames(inp$mask)[inp$mask[k, ]]; d <- d_all[label == lb, ..cols]
      out <- rbindlist(lapply(seq_len(nrow(SPECS)), function(j) {
        pr <- if (SPECS$prior[j] == "default") PREC_DEFAULT else 1 / SIGMA_B^2
        t1 <- proc.time()[["elapsed"]]; tb <- fit_v(d, tms, SPECS$variant[j], pr); sec <- proc.time()[["elapsed"]] - t1
        h <- maxH(d, tms, tb, SPECS$variant[j])
        tb[, `:=`(label = lb, m = m, variant = SPECS$variant[j], prior = SPECS$prior[j], n_events = sum(d$status), n_coef = length(tms), sec = sec,
                  maxH_tr = h[1], maxH_T = h[2])] }))
      prog(FALSE, "  A m%02d %d/%d %s: errors %d", m, k, nrow(inp$td), lb, out[parameter == "alpha" & !is.na(error), .N])
      out
    }, mc.cores = N_WORKERS, mc.preschedule = FALSE)
    bad <- vapply(res, inherits, logical(1), "try-error"); if (any(bad)) stop("PART A worker error: ", paste(res[bad], collapse = " | "))
    saveRDS(rbindlist(res), file.path(OUT, "realA", sprintf("m%02d.rds", m))); rm(rows, d_all); invisible(gc())
    nd <- nd + 1L; prog(TRUE, "[A] imputation m%02d done (%d/%d) | %s", m, nd, length(m_todo), eta_text(t0A, nd, length(m_todo) - nd))
  }
  A <- rbindlist(lapply(file.path(OUT, "realA", sprintf("m%02d.rds", REAL_M)), readRDS))
  al <- A[parameter == "alpha", .(label, m, variant, prior, a = mean, a_sd = sd)]
  A <- merge(A, al, by = c("label", "m", "variant", "prior"))
  A[, `:=`(failed = !is.na(error) | is.na(mean), screen_fail = !is.na(a_sd) & a_sd / (0.78 * a / sqrt(n_events)) < 0.1)]
  A[!parameter %in% c("alpha", "(Intercept)"), `:=`(loghr = fifelse(variant == 1L, a * mean, mean))]
  fwrite(A, file.path(OUT, "realA_fits.csv"))
  fits <- A[, .(n_events = n_events[1], n_coef = n_coef[1], failed = any(failed), screen_fail = any(screen_fail), maxH_tr = maxH_tr[1], maxH_T = maxH_T[1], sec = sec[1]),
            by = .(label, m, variant, prior)]   # per fit (fixed 2026-11-01: was computed over the whole table)
  smA <- fits[, .(fits = .N, inla_errors = sum(failed), screen_fail = sum(screen_fail & !failed), above_cap = sum(maxH_tr > H_CAP | maxH_T > H_CAP, na.rm = TRUE),
                  median_sec = round(median(sec), 1)), keyby = .(variant, prior)]
  # paired comparison variant 0 vs 1 (same data, same prior), where BOTH fits are healthy
  W <- dcast(A[!parameter %in% c("(Intercept)")], label + m + prior + parameter + n_events + n_coef ~ variant, value.var = c("mean", "sd", "loghr", "failed", "screen_fail"))
  W[, ok := !failed_0 & !failed_1 & !screen_fail_0 & !screen_fail_1]
  W[, `:=`(alpha_ratio = fifelse(parameter == "alpha", mean_0 / mean_1, NA_real_),
           loghr_diff_sd = fifelse(parameter != "alpha", (loghr_0 - loghr_1) / sd_0, NA_real_))]
  cmp <- W[ok == TRUE, .(n = .N, median_alpha_ratio_v0_over_v1 = round(median(alpha_ratio, na.rm = TRUE), 3),
                         median_abs_loghr_diff_sd = round(median(abs(loghr_diff_sd), na.rm = TRUE), 3),
                         q95_abs_loghr_diff_sd = round(quantile(abs(loghr_diff_sd), 0.95, na.rm = TRUE), 3)), keyby = prior]
  epc <- W[ok == TRUE & parameter == "alpha"][, epc := cut(n_events / n_coef, c(0, 4, 8, 16, Inf))][, .(n = .N, median_alpha_ratio_v0_over_v1 = round(median(alpha_ratio), 3)), keyby = .(prior, epc)]
  smA <- merge(smA, cmp, by = "prior", all.x = TRUE)
  fwrite(smA, file.path(OUT, "realA_summary.csv")); fwrite(epc, file.path(OUT, "realA_alpha_by_events_per_coef.csv"))
  tr_view <- fits[, .(errors = sum(failed), screen_fail = sum(screen_fail & !failed), above_cap = sum(maxH_tr > H_CAP | maxH_T > H_CAP, na.rm = TRUE)), keyby = .(label, variant, prior)]
  fwrite(dcast(tr_view, label ~ variant + prior, value.var = c("errors", "screen_fail", "above_cap"))[order(label)], file.path(OUT, "realA_transitions.csv"))
  mi_raw <- fread(file.path(FIT_DIR, sprintf("mi_raw_%s.csv", RUN_TAG)))
  rp <- merge(A[variant == 1L & prior == "default", .(label, m, parameter, est = fifelse(parameter == "(Intercept)", mean - LOG_TIME_SCALE, mean))],
              mi_raw[, .(label, m, parameter, est_main = est, sd_main = sd)], by = c("label", "m", "parameter"))
  logl("PART A reproduction (variant 1, default prior vs main run): median |diff| %.2g SD, %d of %d estimates > 0.01 SD",
       rp[, median(abs(est - est_main) / sd_main, na.rm = TRUE)], rp[abs(est - est_main) / sd_main > 0.01, .N], nrow(rp))
  cat("\nPART A summary\n"); print(smA); cat("\nalpha variant 0 / variant 1 by events per coefficient\n"); print(epc)
}

# ==================================================================================================================
# PART B: sim01 datasets under variant 0
# ==================================================================================================================
if ("B" %in% RUN_PARTS) {
  logl("---- PART B: sim01 datasets refitted under variant 0")
  pooled <- fread(file.path(RESULTS_DIR, sprintf("pooled_%s_final.csv", RUN_TAG)))
  cal <- fread(file.path(SIM01_DIR, "calibration.csv"))[!is.na(c_yr)]
  designs <- lapply(setNames(unique(cal$template), unique(cal$template)), function(lb) readRDS(file.path(SIM01_DIR, "design", sprintf("design_%s.rds", gsub("->", "_to_", lb)))))
  truth_for <- function(lb, eff, des) {                    # sim01's truth_for(), variant-1 scale
    if (identical(eff, "real")) eff <- des$real_loghr
    p <- pooled[label == lb]; a <- p[parameter == "alpha", pooled_est]
    b <- setNames(p[match(des$base, parameter), pooled_est], des$base)
    h <- c(harv_light = eff[["light"]], harv_heavy = eff[["heavy"]], gdd_x_harv_light = 0, gdd_x_harv_heavy = 0, tsd_harv_light = 0, tsd_harv_heavy = 0) / a
    list(alpha = a, beta = c(b, h))
  }
  regenerate_dataset <- function(D, outcome) {             # identical to sim01
    cut <- rep(.Machine$integer.max, max(D$soj)); cut[outcome$soj] <- outcome$row
    keep <- seq_len(nrow(D)) <= cut[D$soj]; d <- D[keep]; d[, status := 0L]; ri <- match(outcome$row, which(keep))
    d[ri, `:=`(status = 1L, Tstop = outcome$t)]; d[]
  }
  idx <- fread(file.path(SIM01_DIR, "datasets_index.csv"))
  run_cellB <- function(cell) {
    f <- file.path(OUT, "cellsB", paste0(cell$cell_id, ".rds")); if (file.exists(f)) return(invisible(NULL))
    des <- designs[[cell$template]]; tru <- truth_for(cell$template, EFFECTS[[cell$effect]], des)
    # truth on the variant-0 scale (decades): alpha, alpha * (c_yr + log 10), alpha * beta
    truth0 <- data.table(parameter = c("alpha", "(Intercept)", names(tru$beta)), truth = c(tru$alpha, tru$alpha * (cell$c_yr + LOG_TIME_SCALE), tru$alpha * tru$beta))
    reps <- idx[cell_id == cell$cell_id][order(rep)]
    R <- rbindlist(lapply(seq_len(nrow(reps)), function(r) {
      oc <- readRDS(file.path(SIM01_DIR, "outcomes", sprintf("%s_rep%02d.rds", cell$cell_id, reps$rep[r])))
      d <- regenerate_dataset(des$D, oc); stopifnot(sum(d$status) == reps$n_events[r])
      rbindlist(lapply(seq_len(nrow(ARMS)), function(j) {
        terms <- if (ARMS$harvest[j]) c(des$base, HARV_TERMS) else des$base
        t1 <- proc.time()[["elapsed"]]; tb <- fit_v(d, terms, 0L, ARMS$beta_prec[j])
        tb[, `:=`(arm = ARMS$arm[j], rep = reps$rep[r], sec = proc.time()[["elapsed"]] - t1)] }))
    }))
    R <- merge(R, truth0, by = "parameter", all.x = TRUE, sort = FALSE)
    R <- merge(R, reps[, .(rep, n_events, d0, d_light, d_heavy, n_rows_above_cap)], by = "rep")
    R[, `:=`(cell_id = cell$cell_id, template = cell$template, effect = cell$effect, d0_target = cell$d0_target)]
    saveRDS(R, f)
    n_done <- sum(file.exists(file.path(OUT, "cellsB", paste0(cal$cell_id, ".rds"))))
    prog(TRUE, "[B] cell %d/%d done: %s | fit errors %d of %d | %s", n_done, nrow(cal), cell$cell_id, R[parameter == "alpha" & !is.na(error), .N], nrow(reps) * nrow(ARMS),
         eta_text(T0B, n_done - N0B, nrow(cal) - n_done))
    invisible(NULL)
  }
  todo <- cal[!file.exists(file.path(OUT, "cellsB", paste0(cell_id, ".rds")))]
  N0B <- nrow(cal) - nrow(todo); T0B <- Sys.time()
  prog(TRUE, "[B] start: %d of %d cells to run (%d workers)", nrow(todo), nrow(cal), N_WORKERS)
  err <- unlist(mclapply(split(todo, seq_len(nrow(todo))), function(cell) tryCatch({ run_cellB(cell); NULL }, error = function(e) paste(cell$cell_id, conditionMessage(e))),
                         mc.cores = N_WORKERS, mc.preschedule = FALSE))
  if (length(err)) { logl("ERRORS:\n%s", paste(err, collapse = "\n")); stop("cell errors, see run_log.txt") }
  B <- rbindlist(lapply(file.path(OUT, "cellsB", paste0(cal$cell_id, ".rds")), readRDS))
  al <- B[parameter == "alpha", .(cell_id, rep, arm, a = mean, a_sd = sd)]
  B <- merge(B, al, by = c("cell_id", "rep", "arm"))
  B[, `:=`(err = mean - truth, covered = !is.na(q025) & q025 <= truth & truth <= q975, failed = !is.na(error) | is.na(mean),
           screen_ok = !is.na(a_sd) & a_sd / (0.78 * a / sqrt(pmax(n_events, 1))) >= 0.1)]
  fwrite(B, file.path(OUT, "simB_estimates.csv"))
  # comparison with sim01 (variant 1) on the log-HR scale for coefficients, relative error for alpha;
  # same filter as the sim01 analysis: >= 15 events, no row above INLA's cap, fit ok and passing the alpha-SD screen
  E1 <- fread(file.path(SIM01_DIR, "sim_estimates.csv"))
  E1[, `:=`(screen_ok = !failed & alpha_sd_ratio >= 0.1, v1_est = fifelse(parameter %in% c("alpha", "(Intercept)"), mean, loghr_est),
            v1_truth = fifelse(parameter %in% c("alpha", "(Intercept)"), truth, loghr_truth))]
  J <- merge(B[parameter != "(Intercept)", .(cell_id, rep, arm, parameter, template, effect, d0_target, n_events, d0, n_rows_above_cap,
                                            v0_est = mean, v0_truth = truth, v0_cov = covered, v0_ok = !failed & screen_ok)],
             E1[parameter != "(Intercept)", .(cell_id, rep, arm, parameter, v1_est, v1_truth, v1_cov_beta = covered, v1_ok = screen_ok)],
             by = c("cell_id", "rep", "arm", "parameter"))
  J <- J[n_events >= 15 & n_rows_above_cap == 0]
  fwrite(J, file.path(OUT, "simB_vs_sim01.csv"))
  smB <- J[v0_ok & v1_ok, .(n = .N,
    v0_median_err = round(median(v0_est - v0_truth), 3), v1_median_err = round(median(v1_est - v1_truth), 3),
    v0_median_rel_err_pct = round(100 * median((v0_est - v0_truth) / abs(v0_truth)), 1), v1_median_rel_err_pct = round(100 * median((v1_est - v1_truth) / abs(v1_truth)), 1),
    v0_coverage = round(mean(v0_cov), 2), v1_coverage_beta_scale = round(mean(v1_cov_beta), 2)), keyby = .(template, effect, arm, parameter)]
  fwrite(smB, file.path(OUT, "simB_summary.csv"))
  smB_alpha <- J[v0_ok & v1_ok & parameter == "alpha"][, tot := cut(n_events, c(14, 29, 59, 119, Inf), labels = c("15-29", "30-59", "60-119", "120+"))][
    , .(n = .N, v0_rel_bias_pct = round(100 * median((v0_est - v0_truth) / v0_truth), 1), v1_rel_bias_pct = round(100 * median((v1_est - v1_truth) / v1_truth), 1),
        v0_coverage = round(mean(v0_cov), 2), v1_coverage = round(mean(v1_cov_beta), 2)), keyby = .(arm, tot)]
  fwrite(smB_alpha, file.path(OUT, "simB_alpha_summary.csv"))
  logl("PART B: %d cells, %d fits, %d INLA errors", nrow(cal), B[parameter == "alpha", .N], B[parameter == "alpha" & failed == TRUE, .N])
  cat("\nPART B alpha, variant 0 vs sim01 variant 1 (same datasets)\n"); print(smB_alpha)
}
logl("==== done")
