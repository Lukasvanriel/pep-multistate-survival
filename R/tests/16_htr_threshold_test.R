# ==============================================================================
# 16_htr_threshold_test.R   RUN ON THE MACHINE WITH INLA (part A ~20-40 min, part B ~10-20 min with 6 workers)  -- DIAGNOSTIC ONLY
#
# QUESTION: is INLA's collapse caused by ONE row with an extreme cumulative hazard at truncation, H(tr) = exp(alpha*eta) * tr^alpha, above ~23?
#
# PART A (real data, regen->fir, the stored "int" imputations used in 15_failure_ablation.R):
#   The plot with the highest H(tr) (MLE) is identified per dataset (expected: 7409800101). One change per variant:
#     as_is          data as fitted (control)
#     drop_top       all rows of the top-H(tr) plot removed            -> if this cures the failures, one plot is the cause
#     drop_2nd       the plot with the 2nd-highest H(tr) removed      -> control: removing a plot that is not extreme
#     drop_rand1/2   a random plot with an event removed              -> control: does ANY small data change flip the failure?
#     censor_top     top plot kept, its event set to censored        -> tests "its event status matters"
#     nopest_top     top plot kept, its pest terms set to 0           -> tests "pest pushes its hazard up"
#     H_xx           top plot kept, its time axis shifted (truncation and Tstop of all its rows by the same delta) so that its maximum
#                    H(tr) at the MLE equals xx                       -> dose-response on real data, in failing AND working datasets
#   All fits with INLA_THREADS "1:1" (deterministic, no oversubscription).
# PART B (simulation, no real data): single-row Weibull with left truncation, covariates x1 ~ N(0,1), x2 ~ Bernoulli(0.2), N_SIM rows,
#   ~25% events, true alpha in {1, 1.8, 3}, time unit decades or years (exact reparametrisation). ONE extra row with x1 = x2 = 0 is added
#   whose H(tr) = H_TARGET (0 = no extra row), as an event (Tstop = tr + 2 yr) or censored (Tstop = tr + 10 yr). 3 replicate datasets.
#   If INLA fails above the same H(tr) for every alpha, unit and status, the trigger is H(tr) itself (survival at truncation exp(-H) < ~1e-10),
#   i.e. a numerical limit of INLA, independent of our data and pipeline.
# FAILED (both parts, same rule as 15): INLA error/NA, alpha-SD ratio < 0.1, or |alpha_INLA - alpha_MLE| > 3 MLE-SE (MLE on the same data).
# MUST be run BEFORE 04_prepare_inla.R is re-run (part A needs the int tables; checked). Part B runs on its own (set RUN_REAL <- FALSE).
# OUTPUT: data-output/inttest/fits/htr_real_fits.csv, data-output/diagnostics/htr_sim_fits.csv + printed summaries.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(INLA) }); options(width = 220)
RUN_REAL <- TRUE; RUN_SIM <- TRUE; N_WORKERS <- 6L
INLA_THREADS <- "1:1"; INLA::inla.setOption(num.threads = INLA_THREADS)
LABEL <- "regen->fir"; H_REAL <- c(15, 20, 22, 23, 24, 26, 30)
H_SIM <- c(0, 10, 15, 18, 20, 21, 22, 23, 24, 25, 26, 28, 30, 35, 40); ALPHA_SIM <- c(1, 1.8, 3); N_SIM <- 3000L; REPS <- 1:3
cat("INLA version:", as.character(packageVersion("INLA")), " threads:", INLA_THREADS, "\n")

fit_inla <- function(d, rhs) tryCatch(inla(as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", rhs)), family = "weibullsurv", data = d,
  control.family = list(variant = 1), control.compute = list(config = TRUE),
  control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE))), error = function(e) list(error = conditionMessage(e)))
fit_mle <- function(d, tms) {                                   # S(t) = exp(-(lambda t)^alpha), lambda = exp(eta); theta = (log alpha, intercept, betas)
  X <- cbind(1, as.matrix(d[, tms, drop = FALSE])); T <- d$Tstop; tr <- d$truncation; s <- d$status; lT <- log(T); ltr <- ifelse(tr > 0, log(pmax(tr, 1e-300)), 0)
  nll <- function(th) { a <- exp(th[1]); eta <- drop(X %*% th[-1]); -(sum(s * (th[1] + a * eta + (a - 1) * lT)) - sum(exp(a * eta) * (T^a - tr^a))) }
  gr <- function(th) { a <- exp(th[1]); eta <- drop(X %*% th[-1]); e <- exp(a * eta); Ta <- T^a; tra <- tr^a
    gb <- -(colSums(X * (s * a)) - colSums(X * (a * e * (Ta - tra)))); ga <- -(sum(s * (1 / a + eta + lT)) - sum(e * (eta * (Ta - tra) + Ta * lT - tra * ltr))); c(ga * a, gb) }
  o <- optim(c(0, log(sum(s) / sum(T - tr)), rep(0, length(tms))), nll, gr, method = "BFGS", control = list(maxit = 5000, reltol = 1e-12))
  V <- tryCatch(solve(optimHess(o$par, nll, gr)), error = function(e) NULL); a <- exp(o$par[1])
  list(par = o$par, alpha = a, alpha_se = if (is.null(V)) NA_real_ else a * sqrt(V[1, 1]), conv = o$convergence,
       Htr = exp(a * drop(X %*% o$par[-1])) * tr^a) }
fit_one <- function(d, tms) {                                   # INLA + MLE on the same data; terms that became constant are dropped
  tms <- tms[vapply(tms, function(tm) length(unique(d[[tm]])) > 1, logical(1))]
  dd <- as.data.frame(d[, c("status", "Tstop", "truncation", tms), drop = FALSE]); rhs <- if (length(tms)) paste(tms, collapse = " + ") else "1"
  mo <- fit_inla(dd, rhs); ml <- tryCatch(fit_mle(dd, tms), error = function(e) list(alpha = NA_real_, alpha_se = NA_real_, conv = NA_integer_, Htr = NA_real_))
  a <- if (is.null(mo$error)) { h <- mo$summary.hyperpar; h[grep("alpha", rownames(h))[1], ] } else NULL
  data.table(rows = nrow(dd), events = sum(dd$status == 1L), inla_error = !is.null(mo$error),
             alpha = if (is.null(a)) NA_real_ else a$mean, alpha_sd = if (is.null(a)) NA_real_ else a$sd,
             mle_alpha = ml$alpha, mle_alpha_se = ml$alpha_se, mle_conv = ml$conv, Htr_max = max(ml$Htr)) }
flag <- function(res) { res[, `:=`(ratio = alpha_sd / (0.78 * alpha / sqrt(events)), z_mle = (alpha - mle_alpha) / mle_alpha_se)]
                        res[, failed := inla_error | is.na(alpha) | ratio < 0.1 | abs(z_mle) > 3][] }

# ==== PART A: real data ======================================================================================================================================
if (RUN_REAL) {
  RUN_TAG <- "int"; IN_DIR <- here("data-output", "inttest", "fits"); SIM_M <- 200L; M_MI <- 30L; SD_MIN <- 1e-3; N_COV_SAMPLES <- 0L
  source(here("R", "covariate_spec.R")); COV_CONT <- c("gdd_3", "cmi_sum", "soil_pc1")      # the int run used annual CMI
  source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
  if (!"cmi_sum" %in% names(readRDS(here("data", "tables", "inttest", "visit_intervals.rds")))) stop("visit_intervals.rds is already the gs version: part A needs the int tables.")
  INLA::inla.setOption(num.threads = INLA_THREADS)
  inp <- prepare_fit_inputs(here("data", "tables", "inttest")); k <- match(LABEL, inp$td$label); stopifnot(!is.na(k)); tms0 <- colnames(inp$mask)[inp$mask[k, ]]
  key_iv <- paste(inp$iv$id_pe, inp$iv$interval); b_abs <- inp$iv$t_end
  rebuild_rows <- function(m) {
    ex <- readRDS(file.path(IN_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, m))); p1 <- as.data.table(ex$p1); p2 <- as.data.table(ex$p2); rm(ex)
    ex_row <- rep(NA_real_, nrow(inp$iv)); ii <- match(paste(p2$id_pe, p2$interval), key_iv); stopifnot(!anyNA(ii)); ex_row[ii] <- p2$exit_abs
    p2f <- function(a_dur, b_dur, from, to, r) { E <- b_abs[r] - b_dur; min(max(ex_row[r] - E, a_dur), b_dur) }
    r <- build_rows(inp$iv, setNames(p1$entry, as.character(p1$id_pe)), inp$states, inp$tr, p2_fun = p2f, keep_cols = inp$keep_cols)$rows[kept == TRUE]
    as.data.frame(r[paste(from, to, sep = "->") == LABEL]) }
  jobs <- unique(fread(file.path(IN_DIR, "ablation_fits.csv"))[variant == "as_is" & label == LABEL, .(m, int_failed)])
  cat(sprintf("Part A: %s, %d imputations (%d failed in the int run)\n", LABEL, nrow(jobs), sum(jobs$int_failed)))
  pest_cols <- intersect(c("pest_f", "gdd_x_pest", "tsd_pest"), tms0)

  resA <- rbindlist(mclapply(jobs$m, function(mi) tryCatch({
    d0 <- rebuild_rows(mi); ml0 <- fit_mle(d0, tms0); a <- ml0$alpha; e <- ml0$Htr / pmax(d0$truncation, 1e-300)^a   # exp(alpha*eta) per row
    ph <- data.table(id_pe = d0$id_pe, H = ml0$Htr)[, .(H = max(H)), by = id_pe][order(-H)]; p1 <- ph$id_pe[1]; p2 <- ph$id_pe[2]
    set.seed(1000L + mi); rnd <- sample(setdiff(unique(d0$id_pe[d0$status == 1L]), c(p1, p2)), 2)
    i1 <- which(d0$id_pe == p1); j <- i1[which.max(ml0$Htr[i1])]
    shift <- function(Ht) { dl <- d0$truncation[j] - (Ht / e[j])^(1 / a); d <- d0
      d$truncation[i1] <- pmax(d0$truncation[i1] - dl, 0); d$Tstop[i1] <- d0$Tstop[i1] - dl; d[!(seq_len(nrow(d)) %in% i1) | d$Tstop > d$truncation, ] }
    V <- list(as_is = d0, drop_top = d0[d0$id_pe != p1, ], drop_2nd = d0[d0$id_pe != p2, ], drop_rand1 = d0[d0$id_pe != rnd[1], ], drop_rand2 = d0[d0$id_pe != rnd[2], ],
              censor_top = within(d0, status[id_pe == p1] <- 0L),
              nopest_top = { d <- d0; for (cc in pest_cols) d[[cc]][i1] <- 0; d })
    for (h in H_REAL) V[[sprintf("H_%02d", h)]] <- shift(h)
    rbindlist(lapply(names(V), function(vn) cbind(data.table(m = mi, variant = vn, top_plot = p1, top_H = ph$H[1], second_H = ph$H[2]), fit_one(V[[vn]], tms0)))) },
    error = function(e) data.table(m = mi, variant = paste("run_error:", conditionMessage(e)))), mc.cores = N_WORKERS), fill = TRUE)
  resA <- flag(merge(resA, jobs, by = "m")); fwrite(resA, file.path(IN_DIR, "htr_real_fits.csv"))

  cat("\n== A1. Top-H(tr) plot per dataset (MLE) ==\n")
  print(resA[variant == "as_is", .(m, int_failed, failed_now = failed, top_plot, top_H = round(top_H, 1), second_H = round(second_H, 1))][order(-top_H)])
  cat("\n== A2. Failures per variant: datasets that failed in the int run | that worked ==\n")
  print(resA[!grepl("^H_", variant), .(fail_in_intfailed = sprintf("%d/%d", sum(failed[int_failed]), sum(int_failed)),
            fail_in_intok = sprintf("%d/%d", sum(failed[!int_failed]), sum(!int_failed)), mle_alpha = round(median(mle_alpha), 3), Htr_max = round(median(Htr_max), 1)), by = variant])
  cat("\n== A3. Dose-response: top plot's H(tr) set to H_xx (all datasets) ==\n")
  print(resA[grepl("^H_", variant), .(failed = sprintf("%d/%d", sum(failed), .N), Htr_max_achieved = round(median(Htr_max), 1)), by = variant][order(variant)])
}

# ==== PART B: simulation =====================================================================================================================================
if (RUN_SIM) {
  INLA::inla.setOption(num.threads = INLA_THREADS)
  sim <- function(alpha, unit, rep, Ht, status_x) {           # unit = 1 (decades) or 10 (years); lambda chosen for ~25% events
    set.seed(rep * 100L + round(alpha * 10)); lam <- 0.25 / unit; n <- N_SIM
    x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.2); eta <- log(lam) + 0.3 * x1 + 0.5 * x2; e <- exp(alpha * eta)
    tr <- runif(n, 0, 5 * unit); Tend <- tr + runif(n, 0.5, 1.5) * unit
    tev <- ((e * tr^alpha + rexp(n)) / e)^(1 / alpha); st <- as.integer(tev <= Tend)
    d <- data.frame(status = st, Tstop = pmin(tev, Tend), truncation = tr, x1 = x1, x2 = x2)
    if (Ht > 0) { trx <- (Ht / lam^alpha)^(1 / alpha)                # x1 = x2 = 0 -> exp(alpha*eta) = lam^alpha
      d <- rbind(d, data.frame(status = status_x, Tstop = trx + if (status_x == 1L) 0.2 * unit else 1 * unit, truncation = trx, x1 = 0, x2 = 0)) }
    d }
  grid <- CJ(a_true = ALPHA_SIM, unit = c(1, 10), rep = REPS, H = H_SIM, status_x = c(1L, 0L))[!(H == 0 & status_x == 0L)]
  resB <- rbindlist(mclapply(seq_len(nrow(grid)), function(i) { g <- grid[i]
    tryCatch(cbind(g, fit_one(sim(g$a_true, g$unit, g$rep, g$H, g$status_x), c("x1", "x2"))), error = function(e) cbind(g, data.table(inla_error = TRUE))) },
    mc.cores = N_WORKERS), fill = TRUE)
  resB <- flag(resB); resB[, z_true := (alpha - a_true) / alpha_sd]    # INLA alpha vs the true alpha, in posterior SDs
  dir.create(here("data-output", "diagnostics"), showWarnings = FALSE); fwrite(resB, here("data-output", "diagnostics", "htr_sim_fits.csv"))

  cat("\n== B1. Simulation: failed fits / replicates, by H(tr) of the extra row (rows) and alpha x unit x extra-row status (columns) ==\n")
  resB[, cell := sprintf("a%s_%s_%s", a_true, fifelse(unit == 1, "dec", "yr"), fifelse(status_x == 1L, "ev", "cens"))]
  print(dcast(resB, H ~ cell, value.var = "failed", fun.aggregate = function(x) sprintf("%d/%d", sum(x), length(x))))
  cat("\n== B2. Lowest H(tr) with failures in >= 2 of 3 replicates, per cell (NA = never) ==\n")
  print(resB[, .(f = sum(failed)), by = .(cell, H)][, .(threshold = if (any(f >= 2)) min(H[f >= 2]) else NA_real_, highest_ok = max(c(-Inf, H[f == 0]))), by = cell])
  cat("\n== B3. Healthy fits only: median (INLA alpha - MLE alpha) / MLE SE, by H ==\n")
  print(resB[failed == FALSE, .(n = .N, median_z_mle = round(median(z_mle), 2), median_z_true = round(median(z_true), 2)), by = H][order(H)])
}
