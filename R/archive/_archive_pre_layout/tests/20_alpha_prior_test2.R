# ==============================================================================
# 20_alpha_prior_test2.R   RUN ON THE MACHINE WITH INLA, FROM A TERMINAL (~10-20 min with 6 workers)   -- DIAGNOSTIC ONLY
#   cd ~/Desktop/RProjects/pep-multistate-survival && Rscript R/tests/20_alpha_prior_test2.R
#
# Follow-up of 19_alpha_prior_test.R (weaker PC prior on the Weibull shape alpha). Three parts, run in this order:
#   PART 1  red-flag checks
#     1a  Is the explicit prior specification identical to INLA's default? Prints INLA's own default for alpha (inla.models()), then fits
#         default vs explicit pc.alphaw(5) on the 30 HEALTHY datasets of part 3. STOPS if they differ (then parts 2-3 would be meaningless).
#     1b  other_temperate->white_birch m78 / m111 crashed in script 18 but not in 19 (default prior). Only difference between the scripts:
#         18 kept INLA's working files (working.directory + keep = TRUE). Each dataset is fitted twice with and twice without.
#   PART 2  failing datasets, priors default / pc1 / pc0.1: 10 int-failed imputations of black_spruce->regen and 10 of fir->sugar_maple.
#   PART 3  three healthy, quick transitions, 10 healthy imputations each, priors default / pc5 / pc1 / pc0.1:
#         jack_pine->black_spruce (~106 events, alpha ~1.6), yellow_birch->fir (~105 events, alpha ~1.0), sugar_maple->other_temperate (~118, ~1.0);
#         added (sparser, 22-57 events): other_pioneer->white_birch 22, cedar->yellow_birch 26, other_pioneer->other_temperate 31,
#         red_maple->yellow_birch 33, yellow_birch->red_maple 38, cedar->fir 57. Fits already saved are skipped on re-running.
#   Every dataset also gets a maximum-likelihood fit (no prior) as a reference.
# Priors: control.family = list(variant = 1, hyper = list(alpha = list(prior = "pc.alphaw", param = <strength>))); lower = weaker pull toward 1.
# FAILED: INLA error/NA, or alpha-SD ratio < 0.1. One thread per fit. SAVES AS IT GOES (one file per dataset); re-running skips finished ones.
# MUST be run BEFORE 04_prepare_inla.R is re-run (needs the int tables; checked).
# OUTPUT (data/processed/fit_int/prior_test2/):
#   all_fits.csv     every fit x parameter (alpha, intercept, betas): est, sd, log hazard ratio alpha*beta and its approximate SD (alpha*sd(beta))
#   pooled.csv       per transition x prior x parameter: Rubin-pooled over the imputations (alpha; log HR -> HR with 95% interval)
#   + printed summaries per transition (alpha, failures, hazard ratios side by side per prior)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(INLA) }); options(width = 230)
RUN_TAG <- "int"; IN_DIR <- here("data", "processed", "fit_int"); N_WORKERS <- 6L; INLA::inla.setOption(num.threads = "1:1")
N_M <- 10L
HEALTHY <- c("jack_pine->black_spruce", "yellow_birch->fir", "sugar_maple->other_temperate",          # ~105-118 events (first run)
             "other_pioneer->white_birch", "cedar->yellow_birch", "other_pioneer->other_temperate",     # 20-60 events (added 2026-10-23)
             "red_maple->yellow_birch", "yellow_birch->red_maple", "cedar->fir")
FAILING <- c("black_spruce->regen", "fir->sugar_maple")
OTWB <- "other_temperate->white_birch"; OTWB_M <- c(78L, 111L)
SIM_M <- 200L; M_MI <- 30L; SD_MIN <- 1e-3; N_COV_SAMPLES <- 0L
source(here("R", "covariate_spec.R")); COV_CONT <- c("gdd_3", "cmi_sum", "soil_pc1")      # the int run used annual CMI
source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
if (!"cmi_sum" %in% names(readRDS(here("data", "processed", "tables_int", "visit_intervals.rds")))) stop("visit_intervals.rds is already the gs version: run this before re-running 04.")
inp <- prepare_fit_inputs(here("data", "processed", "tables_int")); stopifnot(all(c(HEALTHY, FAILING, OTWB) %in% inp$td$label))
OUT <- file.path(IN_DIR, "prior_test2"); dir.create(OUT, showWarnings = FALSE)
cat("INLA version:", as.character(packageVersion("INLA")), "\n")

key_iv <- paste(inp$iv$id_pe, inp$iv$interval); b_abs <- inp$iv$t_end
rebuild_rows <- function(m, lb) {
  ex <- readRDS(file.path(IN_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, m))); p1 <- as.data.table(ex$p1); p2 <- as.data.table(ex$p2); rm(ex)
  ex_row <- rep(NA_real_, nrow(inp$iv)); ii <- match(paste(p2$id_pe, p2$interval), key_iv); stopifnot(!anyNA(ii)); ex_row[ii] <- p2$exit_abs
  p2f <- function(a_dur, b_dur, from, to, r) { E <- b_abs[r] - b_dur; min(max(ex_row[r] - E, a_dur), b_dur) }
  r <- build_rows(inp$iv, setNames(p1$entry, as.character(p1$id_pe)), inp$states, inp$tr, p2_fun = p2f, keep_cols = inp$keep_cols)$rows[kept == TRUE]
  d <- r[paste(from, to, sep = "->") == lb]; k <- match(lb, inp$td$label); tms <- colnames(inp$mask)[inp$mask[k, ]]
  list(dd = as.data.frame(d[, c("status", "Tstop", "truncation", tms), with = FALSE]), tms = tms, rhs = inp$rhs[[k]]) }
fit_inla <- function(dd, rhs, strength, wdir = NULL) {
  cf <- list(variant = 1); if (!is.null(strength)) cf$hyper <- list(alpha = list(prior = "pc.alphaw", param = strength))
  args <- list(formula = as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", rhs)), family = "weibullsurv", data = dd,
               control.family = cf, control.compute = list(config = TRUE),
               control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)))
  if (!is.null(wdir)) { unlink(wdir, recursive = TRUE); args$working.directory <- wdir; args$keep <- TRUE }
  mo <- tryCatch(do.call(inla, args), error = function(e) list(error = conditionMessage(e)))
  if (!is.null(wdir)) unlink(wdir, recursive = TRUE)
  if (!is.null(mo$error)) return(data.table(parameter = "alpha", est = NA_real_, sd = NA_real_, inla_error = TRUE))
  h <- mo$summary.hyperpar; f <- mo$summary.fixed
  data.table(parameter = c("alpha", rownames(f)), est = c(h[ALPHA_NAME, "mean"], f$mean), sd = c(h[ALPHA_NAME, "sd"], f$sd), inla_error = FALSE) }
fit_mle <- function(d, tms) {                                   # S(t) = exp(-(lambda t)^alpha), lambda = exp(eta)
  X <- cbind(1, as.matrix(d[, tms, drop = FALSE])); T <- d$Tstop; tr <- d$truncation; s <- d$status; lT <- log(T); ltr <- ifelse(tr > 0, log(pmax(tr, 1e-300)), 0)
  nll <- function(th) { a <- exp(th[1]); eta <- drop(X %*% th[-1]); -(sum(s * (th[1] + a * eta + (a - 1) * lT)) - sum(exp(a * eta) * (T^a - tr^a))) }
  gr <- function(th) { a <- exp(th[1]); eta <- drop(X %*% th[-1]); e <- exp(a * eta); Ta <- T^a; tra <- tr^a
    gb <- -(colSums(X * (s * a)) - colSums(X * (a * e * (Ta - tra)))); ga <- -(sum(s * (1 / a + eta + lT)) - sum(e * (eta * (Ta - tra) + Ta * lT - tra * ltr))); c(ga * a, gb) }
  o <- optim(c(0, log(sum(s) / sum(T - tr)), rep(0, length(tms))), nll, gr, method = "BFGS", control = list(maxit = 5000, reltol = 1e-12))
  V <- tryCatch(solve(optimHess(o$par, nll, gr)), error = function(e) NULL); se <- if (is.null(V)) rep(NA_real_, length(o$par)) else sqrt(diag(V))
  data.table(parameter = c("alpha", "(Intercept)", tms), est = c(exp(o$par[1]), o$par[-1]), sd = c(exp(o$par[1]) * se[1], se[-1]), inla_error = FALSE) }

fv <- fread(file.path(IN_DIR, sprintf("fit_validity_%s.csv", RUN_TAG)))
avail <- as.integer(sub(".*_m(\\d+)\\.rds$", "\\1", list.files(IN_DIR, pattern = sprintf("^extras_%s_m\\d+\\.rds$", RUN_TAG))))
pick <- function(lb, valid_, n = N_M) head(fv[label == lb & valid == valid_ & m %in% avail, sort(m)], n)
run_jobs <- function(J) invisible(mclapply(seq_len(nrow(J)), function(i) { j <- J[i]
  res <- tryCatch({ x <- rebuild_rows(j$m, j$label); nev <- sum(x$dd$status == 1L)
    fits <- rbindlist(lapply(j$priors[[1]], function(p) cbind(data.table(prior = p), fit_inla(x$dd, x$rhs, if (p == "default") NULL else as.numeric(sub("pc", "", p))))), fill = TRUE)
    cbind(j[, .(part, label, m)], events = nev, rbind(fits, cbind(data.table(prior = "MLE"), fit_mle(x$dd, x$tms)), fill = TRUE)) },
    error = function(e) cbind(j[, .(part, label, m)], data.table(prior = "run_error", parameter = substr(conditionMessage(e), 1, 150))))
  saveRDS(res, j$file); NULL }, mc.cores = N_WORKERS, mc.preschedule = FALSE))
flag <- function(r) { a <- r[parameter == "alpha"]; a[, failed := prior != "MLE" & (inla_error %in% TRUE | is.na(est) | sd / (0.78 * est / sqrt(events)) < 0.1)]; a }

# ==== PART 1a: INLA's default prior, and default vs explicit pc5 on the healthy datasets ====================================================================
cat("\n== 1a. INLA's own default for the Weibull shape (weibullsurv) ==\n")
hy <- inla.models()$likelihood$weibullsurv$hyper; str(hy[[1]][c("name", "short.name", "initial", "prior", "param", "to.theta")])
jobs3 <- rbindlist(lapply(HEALTHY, function(lb) data.table(part = "3_healthy", label = lb, m = pick(lb, TRUE))))
jobs3[, priors := list(list(c("default", "pc5", "pc1", "pc0.1")))]
jobs3[, file := file.path(OUT, sprintf("p3_%s_m%03d.rds", gsub("->", "_to_", label), m))]
t0 <- proc.time()["elapsed"]; run_jobs(jobs3[!file.exists(file)])
r3 <- rbindlist(lapply(jobs3[file.exists(file), file], readRDS), fill = TRUE)
chk <- merge(r3[prior == "default", .(label, m, parameter, d = est, dsd = sd)], r3[prior == "pc5", .(label, m, parameter, p5 = est)], by = c("label", "m", "parameter"))
chk[, z := abs(p5 - d) / dsd]                                     # difference in posterior SDs (tolerance 1%; rounding-level differences do occur)
cat(sprintf("Default vs explicit pc5: %d parameter estimates on %d datasets; identical: %d; largest |difference| = %.3g = %.3g posterior SD\n",
            nrow(chk), uniqueN(chk[, .(label, m)]), sum(chk$p5 == chk$d, na.rm = TRUE), max(abs(chk$d - chk$p5), na.rm = TRUE), max(chk$z, na.rm = TRUE)))
if (nrow(chk) == 0 || max(chk$z, na.rm = TRUE) > 0.01) {
  print(chk[order(-z)][1:20]); stop("Explicit pc.alphaw(5) does NOT reproduce the default: the prior specification is wrong. Parts 2-3 not run/summarised.") }
cat("OK: the explicit specification reproduces INLA's default (all differences < 1% of a posterior SD).\n")

# ==== PART 1b: m78 / m111 with and without INLA's working files kept ====================================================================================
cat("\n== 1b. other_temperate->white_birch m78 / m111: default prior, with vs without keep = TRUE (2 fits each) ==\n")
f1b <- file.path(OUT, "p1b.rds")
if (!file.exists(f1b)) {
  j1 <- CJ(m = OTWB_M, keep = c(FALSE, TRUE), rep = 1:2)
  r1 <- rbindlist(mclapply(seq_len(nrow(j1)), function(i) { g <- j1[i]; x <- rebuild_rows(g$m, OTWB)
    f <- fit_inla(x$dd, x$rhs, NULL, if (g$keep) file.path(OUT, sprintf("wd_m%d_%d", g$m, g$rep)) else NULL)
    cbind(g, f[parameter == "alpha", .(alpha = est, alpha_sd = sd, inla_error)]) }, mc.cores = N_WORKERS))
  saveRDS(r1, f1b) }
print(readRDS(f1b))

# ==== PART 2: failing datasets ======================================================================================================================
jobs2 <- rbindlist(lapply(FAILING, function(lb) data.table(part = "2_failing", label = lb, m = pick(lb, FALSE))))
jobs2[, priors := list(list(c("default", "pc1", "pc0.1")))]
jobs2[, file := file.path(OUT, sprintf("p2_%s_m%03d.rds", gsub("->", "_to_", label), m))]
run_jobs(jobs2[!file.exists(file)])
cat(sprintf("\nRun time parts 1-3: %.1f min\n", (proc.time()["elapsed"] - t0) / 60))

# ==== Summaries ============================================================================================================================================
res <- rbindlist(lapply(c(jobs3[file.exists(file), file], jobs2[file.exists(file), file]), readRDS), fill = TRUE)
a <- flag(res)
res <- merge(res, a[, .(label, m, prior, alpha = est, failed)], by = c("label", "m", "prior"), all.x = TRUE)
res[!parameter %in% c("alpha", "(Intercept)"), `:=`(logHR = alpha * est, logHR_sd = alpha * sd)]
fwrite(res[order(part, label, m, prior, parameter)], file.path(OUT, "all_fits.csv"))

cat("\n== 2. Failing datasets: failed fits per prior (int-failed imputations) ==\n")
print(dcast(a[part == "2_failing" & prior != "MLE", .(failed = sprintf("%d/%d", sum(failed), .N)), by = .(label, prior)], label ~ prior, value.var = "failed"))

pool <- function(x, y) {                                        # Rubin: x = estimates, y = within-imputation SDs
  M <- length(x); W <- mean(y^2); B <- if (M > 1) var(x) else 0; list(M = M, est = mean(x), se = sqrt(W + (1 + 1 / M) * B)) }
PR <- c("MLE", "default", "pc5", "pc1", "pc0.1")
pa <- a[failed == FALSE, pool(est, sd), by = .(part, label, prior)][, parameter := "alpha"]
pb <- res[failed == FALSE & !is.na(logHR), pool(logHR, logHR_sd), by = .(part, label, prior, parameter)]
pooled <- rbind(pa, pb, fill = TRUE)[, prior := factor(prior, PR)][order(part, label, parameter, prior)]
pooled[parameter != "alpha", `:=`(HR = exp(est), HR_lo = exp(est - 1.96 * se), HR_hi = exp(est + 1.96 * se))]
fwrite(pooled, file.path(OUT, "pooled.csv"))

cat("\n== 3. Per transition: alpha and hazard ratios per prior (pooled over the imputations; failed fits excluded; M = fits used) ==\n")
for (lb in unique(pooled$label)) {
  p <- pooled[label == lb]; cat(sprintf("\n----- %s (%s; events %d) -----\n", lb, p$part[1], res[label == lb, events[1]]))
  al <- p[parameter == "alpha", .(prior, M, alpha = round(est, 3), se = round(se, 3))]; print(al)
  hr <- p[parameter != "alpha", .(parameter, prior, cell = sprintf("%.2f (%.2f-%.2f)", HR, HR_lo, HR_hi))]
  print(dcast(hr, parameter ~ prior, value.var = "cell")) }
cat("\nHazard ratios use exp(alpha*beta); their interval ignores the alpha-beta covariance (diagnostic approximation, not the pipeline's 09 method).\n")
