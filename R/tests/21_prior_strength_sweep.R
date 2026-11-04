# ==============================================================================
# 21_prior_strength_sweep.R   RUN ON THE MACHINE WITH INLA, FROM A TERMINAL (~10-15 min with 6 workers)   -- DIAGNOSTIC ONLY
#   cd ~/Desktop/RProjects/pep-multistate-survival && Rscript R/tests/21_prior_strength_sweep.R
#
# QUESTION: how do alpha, its uncertainty and the hazard ratios change as the PC prior on the Weibull shape is weakened step by step,
#   and where does further weakening stop mattering (or start causing trouble: crashes, implausible alpha tails)?
# PRIORS: default (= pc.alphaw(5)) | 1 | 0.1 | 0.05 | 0.01   (+ maximum likelihood = data only, as reference)
# TRANSITIONS (healthy in the int run; 10 healthy imputations each, rows rebuilt from extras):
#   higher N: fir->white_birch (~276 events), jack_pine->black_spruce (~106)
#   lower N : other_pioneer->white_birch (~22; alpha snaps to ~1 under the default), cedar->yellow_birch (~26; alpha > 1)
# FAILED: INLA error/NA, or alpha-SD ratio < 0.1. One thread per fit. SAVES AS IT GOES (one file per dataset); re-running skips finished ones.
# MUST be run BEFORE 04_prepare_inla.R is re-run (needs the int tables; checked).
# OUTPUT (data-output/inttest/fits/prior_sweep/): all_fits.csv (every fit x parameter, incl. alpha's 2.5% / 97.5% posterior quantiles),
#   pooled.csv (Rubin-pooled alpha and hazard ratios per transition x prior) + printed summary.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(INLA) }); options(width = 230)
RUN_TAG <- "int"; IN_DIR <- here("data-output", "inttest", "fits"); N_WORKERS <- 6L; INLA::inla.setOption(num.threads = "1:1")
N_M <- 10L
LABELS <- c("fir->white_birch", "jack_pine->black_spruce", "other_pioneer->white_birch", "cedar->yellow_birch")
PRIORS <- c(default = NA, pc1 = 1, pc0.1 = 0.1, pc0.05 = 0.05, pc0.01 = 0.01)
SIM_M <- 200L; M_MI <- 30L; SD_MIN <- 1e-3; N_COV_SAMPLES <- 0L
source(here("R", "covariate_spec.R")); COV_CONT <- c("gdd_3", "cmi_sum", "soil_pc1")      # the int run used annual CMI
source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
if (!"cmi_sum" %in% names(readRDS(here("data", "tables", "inttest", "visit_intervals.rds")))) stop("visit_intervals.rds is already the gs version: run this before re-running 04.")
inp <- prepare_fit_inputs(here("data", "tables", "inttest")); stopifnot(all(LABELS %in% inp$td$label))
OUT <- file.path(IN_DIR, "prior_sweep"); dir.create(OUT, showWarnings = FALSE)
cat("INLA version:", as.character(packageVersion("INLA")), "\n")

key_iv <- paste(inp$iv$id_pe, inp$iv$interval); b_abs <- inp$iv$t_end
rebuild_rows <- function(m, lb) {
  ex <- readRDS(file.path(IN_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, m))); p1 <- as.data.table(ex$p1); p2 <- as.data.table(ex$p2); rm(ex)
  ex_row <- rep(NA_real_, nrow(inp$iv)); ii <- match(paste(p2$id_pe, p2$interval), key_iv); stopifnot(!anyNA(ii)); ex_row[ii] <- p2$exit_abs
  p2f <- function(a_dur, b_dur, from, to, r) { E <- b_abs[r] - b_dur; min(max(ex_row[r] - E, a_dur), b_dur) }
  r <- build_rows(inp$iv, setNames(p1$entry, as.character(p1$id_pe)), inp$states, inp$tr, p2_fun = p2f, keep_cols = inp$keep_cols)$rows[kept == TRUE]
  d <- r[paste(from, to, sep = "->") == lb]; k <- match(lb, inp$td$label); tms <- colnames(inp$mask)[inp$mask[k, ]]
  list(dd = as.data.frame(d[, c("status", "Tstop", "truncation", tms), with = FALSE]), tms = tms, rhs = inp$rhs[[k]]) }
fit_inla <- function(dd, rhs, strength) {
  cf <- list(variant = 1); if (!is.na(strength)) cf$hyper <- list(alpha = list(prior = "pc.alphaw", param = strength))
  mo <- tryCatch(inla(as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", rhs)), family = "weibullsurv", data = dd,
    control.family = cf, control.compute = list(config = TRUE),
    control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE))), error = function(e) list(error = conditionMessage(e)))
  if (!is.null(mo$error)) return(data.table(parameter = "alpha", est = NA_real_, sd = NA_real_, q025 = NA_real_, q975 = NA_real_, inla_error = TRUE))
  h <- mo$summary.hyperpar; f <- mo$summary.fixed
  data.table(parameter = c("alpha", rownames(f)), est = c(h[ALPHA_NAME, "mean"], f$mean), sd = c(h[ALPHA_NAME, "sd"], f$sd),
             q025 = c(h[ALPHA_NAME, "0.025quant"], f[, "0.025quant"]), q975 = c(h[ALPHA_NAME, "0.975quant"], f[, "0.975quant"]), inla_error = FALSE) }
fit_mle <- function(d, tms) {                                   # S(t) = exp(-(lambda t)^alpha), lambda = exp(eta)
  X <- cbind(1, as.matrix(d[, tms, drop = FALSE])); T <- d$Tstop; tr <- d$truncation; s <- d$status; lT <- log(T); ltr <- ifelse(tr > 0, log(pmax(tr, 1e-300)), 0)
  nll <- function(th) { a <- exp(th[1]); eta <- drop(X %*% th[-1]); -(sum(s * (th[1] + a * eta + (a - 1) * lT)) - sum(exp(a * eta) * (T^a - tr^a))) }
  gr <- function(th) { a <- exp(th[1]); eta <- drop(X %*% th[-1]); e <- exp(a * eta); Ta <- T^a; tra <- tr^a
    gb <- -(colSums(X * (s * a)) - colSums(X * (a * e * (Ta - tra)))); ga <- -(sum(s * (1 / a + eta + lT)) - sum(e * (eta * (Ta - tra) + Ta * lT - tra * ltr))); c(ga * a, gb) }
  o <- optim(c(0, log(sum(s) / sum(T - tr)), rep(0, length(tms))), nll, gr, method = "BFGS", control = list(maxit = 5000, reltol = 1e-12))
  V <- tryCatch(solve(optimHess(o$par, nll, gr)), error = function(e) NULL); se <- if (is.null(V)) rep(NA_real_, length(o$par)) else sqrt(diag(V))
  est <- c(exp(o$par[1]), o$par[-1]); sdv <- c(exp(o$par[1]) * se[1], se[-1])
  data.table(parameter = c("alpha", "(Intercept)", tms), est = est, sd = sdv, q025 = est - 1.96 * sdv, q975 = est + 1.96 * sdv, inla_error = FALSE) }

fv <- fread(file.path(IN_DIR, sprintf("fit_validity_%s.csv", RUN_TAG)))
avail <- as.integer(sub(".*_m(\\d+)\\.rds$", "\\1", list.files(IN_DIR, pattern = sprintf("^extras_%s_m\\d+\\.rds$", RUN_TAG))))
jobs <- rbindlist(lapply(LABELS, function(lb) data.table(label = lb, m = head(fv[label == lb & valid == TRUE & m %in% avail, sort(m)], N_M))))
jobs[, file := file.path(OUT, sprintf("ds_%s_m%03d.rds", gsub("->", "_to_", label), m))]
todo <- jobs[!file.exists(file)]; cat(sprintf("Datasets: %d, already done: %d\n", nrow(jobs), nrow(jobs) - nrow(todo)))
t0 <- proc.time()["elapsed"]
invisible(mclapply(seq_len(nrow(todo)), function(i) { j <- todo[i]
  res <- tryCatch({ x <- rebuild_rows(j$m, j$label); nev <- sum(x$dd$status == 1L)
    fits <- rbindlist(lapply(names(PRIORS), function(p) cbind(data.table(prior = p), fit_inla(x$dd, x$rhs, PRIORS[[p]]))), fill = TRUE)
    cbind(j[, .(label, m)], events = nev, rbind(fits, cbind(data.table(prior = "MLE"), fit_mle(x$dd, x$tms)), fill = TRUE)) },
    error = function(e) cbind(j[, .(label, m)], data.table(prior = "run_error", parameter = substr(conditionMessage(e), 1, 150))))
  saveRDS(res, j$file); NULL }, mc.cores = N_WORKERS, mc.preschedule = FALSE))
cat(sprintf("Run time: %.1f min\n", (proc.time()["elapsed"] - t0) / 60))

# ---- summary --------------------------------------------------------------------------------------------------------------------------------------------
res <- rbindlist(lapply(jobs[file.exists(file), file], readRDS), fill = TRUE)
a <- res[parameter == "alpha"][, failed := prior != "MLE" & (inla_error %in% TRUE | is.na(est) | sd / (0.78 * est / sqrt(events)) < 0.1)]
res <- merge(res, a[, .(label, m, prior, alpha = est, failed)], by = c("label", "m", "prior"), all.x = TRUE)
res[!parameter %in% c("alpha", "(Intercept)"), `:=`(logHR = alpha * est, logHR_sd = alpha * sd)]
fwrite(res, file.path(OUT, "all_fits.csv"))
PR <- c("MLE", names(PRIORS)); ev <- a[, .(events = events[1]), by = label]
pool <- function(x, y) { M <- length(x); W <- mean(y^2); B <- if (M > 1) var(x) else 0; list(M = M, est = mean(x), se = sqrt(W + (1 + 1 / M) * B)) }
pa <- a[failed == FALSE, c(pool(est, sd), list(q025 = median(q025), q975 = median(q975))), by = .(label, prior)][, parameter := "alpha"]
pb <- res[failed == FALSE & !is.na(logHR), pool(logHR, logHR_sd), by = .(label, prior, parameter)]
pooled <- merge(rbind(pa, pb, fill = TRUE), ev, by = "label")[, prior := factor(prior, PR)][order(-events, label, parameter, prior)]
pooled[parameter != "alpha", `:=`(HR = exp(est), HR_lo = exp(est - 1.96 * se), HR_hi = exp(est + 1.96 * se))]
fwrite(pooled, file.path(OUT, "pooled.csv"))

cat("\n== Failed fits per prior (out of", N_M, ") ==\n")
print(dcast(a[prior != "MLE", .(failed = sum(failed)), by = .(label, prior)][, prior := factor(prior, PR)], label ~ prior, value.var = "failed"))
cat("\n== alpha: pooled estimate (SE) and median 95% interval per prior ==\n")
al <- pooled[parameter == "alpha", .(label, events, prior, cell = sprintf("%.3f (%.3f) [%.2f-%.2f]", est, se, q025, q975))]
print(dcast(al, label + events ~ prior, value.var = "cell")[order(-events)])
cat("\n== alpha SE relative to the data-only (MLE) SE ==\n")
r <- pooled[parameter == "alpha"]; r <- merge(r, r[prior == "MLE", .(label, se_mle = se)], by = "label")
print(dcast(r[prior != "MLE", .(label, events, prior, ratio = round(se / se_mle, 2))], label + events ~ prior, value.var = "ratio")[order(-events)])
cat("\n== Hazard ratios: shift of the log HR vs the default prior, in default posterior SDs (median / max over terms and imputations) ==\n")
b <- res[!is.na(logHR) & failed == FALSE & prior != "MLE"]
bd <- merge(b[prior != "default"], b[prior == "default", .(label, m, parameter, l0 = logHR, s0 = logHR_sd)], by = c("label", "m", "parameter"))
print(dcast(bd[, .(cell = sprintf("%.2f / %.2f", median(abs((logHR - l0) / s0)), max(abs((logHR - l0) / s0)))), by = .(label, prior)][, prior := factor(prior, PR)],
            label ~ prior, value.var = "cell"))
cat("\n== Hazard ratios per transition and prior (pooled; interval ignores alpha-beta covariance) ==\n")
for (lb in unique(pooled$label)) { cat("\n---", lb, "---\n")
  print(dcast(pooled[label == lb & parameter != "alpha", .(parameter, prior, cell = sprintf("%.2f (%.2f-%.2f)", HR, HR_lo, HR_hi))], parameter ~ prior, value.var = "cell")) }
