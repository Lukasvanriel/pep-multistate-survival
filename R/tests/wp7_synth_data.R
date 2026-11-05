# ==============================================================================
# R/tests/wp7_synth_data.R   (WP7, cloud tests only; NEVER run in the real project root)
# Writes SYNTHETIC 05/07 outputs (mi_raw_<tag>.csv, extras_<tag>_mXX.rds with stored covariances, screen_ref_<tag>.csv) for the 90 kept transitions x 32 imputations into a
# scratch root, so that 08 -> 09 -> 12 -> 10 can be run without INLA. The scratch root needs data/tables/<tag>/transition_set.rds and visit_intervals.rds (copies of the real ones)
# and the R/ scripts. Usage:  WP7_SYNTH_OK=1 Rscript R/tests/wp7_synth_data.R <scratch_root>
# Truth: random coefficients per transition; fits = truth + noise; the stored covariance is a random positive-definite matrix over (intercept, terms, alpha).
# ==============================================================================
suppressPackageStartupMessages(library(data.table))
if (Sys.getenv("WP7_SYNTH_OK") != "1") stop("synthetic data only: set WP7_SYNTH_OK=1 and pass a scratch root")
root <- commandArgs(TRUE)[1]; stopifnot(!is.na(root), dir.exists(root))
setwd(root); stopifnot("refusing: real fits found in this root" = !file.exists("data-output/final/fits/mi_raw_final.csv"))
source("R/config.R"); source("R/covariate_spec.R"); set.seed(7)
ts <- readRDS(file.path("data", "tables", RUN_TAG, "transition_set.rds")); setDT(ts); mk <- term_mask(ts)
fd <- file.path("data-output", RUN_TAG, "fits"); dir.create(fd, recursive = TRUE, showWarnings = FALSE)
k_ts <- ts[kept == TRUE]; labs <- rownames(mk); stopifnot(identical(labs, paste(k_ts$from, k_ts$to, sep = "->")))
HALFNAMES <- c("alpha parameter for weibullsurv")
raw <- list(); ex <- vector("list", 32)   # m31, m32 = replacements for the two deliberately invalid fits (k = 5, m = 3; k = 9, m = 7)
truth <- lapply(seq_along(labs), function(k) { tm <- colnames(mk)[mk[k, ]]; list(tm = tm, b0 = rnorm(1, -3, .5), b = setNames(rnorm(length(tm), 0, .4), tm), a = runif(1, .8, 1.5)) })
A <- lapply(seq_along(labs), function(k) { p <- length(truth[[k]]$tm) + 2L; M <- matrix(rnorm(p * p, 0, .06), p); M %*% t(M) + diag(.01, p) })
for (m in 1:32) { ex[[m]] <- list(fits = list())
  for (k in seq_along(labs)) { tr <- truth[[k]]; nm <- c("(Intercept)", tr$tm, "alpha"); V <- A[[k]]; dimnames(V) <- list(nm, nm)
    est <- c(tr$b0 + rnorm(1, 0, .05) + if ((k == 5 && m == 3) || (k == 9 && m == 7)) 3 else 0, tr$b + rnorm(length(tr$tm), 0, .05), tr$a + rnorm(1, 0, .03)); names(est) <- nm
    sdv <- sqrt(diag(V)) * runif(length(nm), .95, 1.05)
    ne <- k_ts$n_events[k]
    cap <- c(maxH_inla = runif(1, 5, 25), n_above_inla = rpois(1, .5), maxH_map = runif(1, 5, 25), n_above_map = rpois(1, .5), z_map_max = abs(rnorm(1)), z_map_alpha = abs(rnorm(1)), z_map_intercept = abs(rnorm(1)))
    raw[[length(raw) + 1L]] <- data.table(k = k, parameter = c("alpha", setdiff(nm, "alpha")), est = est[c("alpha", setdiff(nm, "alpha"))], sd = sdv[match(c("alpha", setdiff(nm, "alpha")), nm)],
                                          fit_error = FALSE, error_msg = NA_character_, mode_status = 0L, mlik = -100, log_warn = "", label = labs[k], m = m, n_events = ne)[, names(cap) := as.list(cap)]
    sf <- cbind(mean = est[setdiff(nm, "alpha")], sd = sdv[match(setdiff(nm, "alpha"), nm)]); rownames(sf) <- setdiff(nm, "alpha")
    sh <- matrix(c(est["alpha"], sdv[nm == "alpha"]), 1, dimnames = list(HALFNAMES, c("mean", "sd")))
    ex[[m]]$fits[[labs[k]]] <- list(summary_fixed = sf, summary_hyper = sh, cov = list(mean = est, cov = V, n = 1000L)) }
  saveRDS(ex[[m]], file.path(fd, sprintf("extras_%s_m%02d.rds", RUN_TAG, m))) }
raw <- rbindlist(raw); fwrite(raw, file.path(fd, sprintf("mi_raw_%s.csv", RUN_TAG)))
ref <- raw[parameter == "(Intercept)", .(icpt_ref = median(est), icpt_sd_ref = median(sd), n_ref = .N), by = k]; fwrite(ref, file.path(fd, sprintf("screen_ref_%s.csv", RUN_TAG)))
saveRDS(list(truth = truth, labs = labs), file.path(fd, "synthetic_truth.rds"))
cat("synthetic fits written:", length(labs), "transitions x 32 imputations ->", fd, "\n")
