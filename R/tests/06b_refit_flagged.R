# ==============================================================================
# 06b_refit_flagged.R   RUN ON THE MACHINE WITH INLA (a few minutes)
#
# Are the flagged fits optimiser failures that a better INLA start value repairs?
# For every flagged fit (flagged_fits_diag.csv) the saved fitting data
# (baddata_diag_mXX.rds) are refitted with the pipeline's exact settings:
#   A  as in the pipeline (INLA defaults)             -> does the failure repeat?
#   B  start value for alpha_intern = median of that transition's converged fits
#      in the other imputations, control.mode(restart = TRUE)
#   C  single-threaded (num.threads = "1:1")
#   D  verbose = TRUE (only for the INLA-error fits, output kept in the console)
# For each: alpha mean/sd, log marginal likelihood, function evaluations, internal theta.
# The settings replicate fit_inla_models() in the (read-only) simulation folder.
# Output: data/processed/fit/refit_flagged_<tag>.csv
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(INLA) })
RUN_TAG <- "diag"; OUT_DIR <- here("data", "processed", "fit")
ALPHA_NAME <- "alpha parameter for weibullsurv"
flagged <- fread(file.path(OUT_DIR, sprintf("flagged_fits_%s.csv", RUN_TAG)))

# reference internal theta per transition: median over NON-flagged imputations (from the stored extras)
ref <- rbindlist(lapply(1:30, function(m) { f <- file.path(OUT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, m)); if (!file.exists(f)) return(NULL)
  ex <- readRDS(f)$fits; rbindlist(lapply(names(ex), function(lb) if (is.null(ex[[lb]]$error)) data.table(m = m, label = lb, theta = as.numeric(ex[[lb]]$mode_theta[1]), nfunc = as.numeric(ex[[lb]]$nfunc)))) }))
ref <- merge(ref, flagged[, .(m, label, bad = TRUE)], by = c("m", "label"), all.x = TRUE); ref <- ref[is.na(bad)]
refstart <- ref[, .(theta_ref = median(theta), nfunc_ref = median(nfunc)), by = label]

fit1 <- function(d, rhs, start = NULL, threads = NULL, verbose = FALSE) {
  a <- list(formula = as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", rhs)),
            family = "weibullsurv", data = d, control.family = list(variant = 1), control.compute = list(config = TRUE),
            control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)), verbose = verbose)
  if (!is.null(start))   a$control.mode <- list(theta = start, restart = TRUE)
  if (!is.null(threads)) a$num.threads <- threads
  tryCatch(do.call(inla, a), error = function(e) list(error = conditionMessage(e)))
}
summ <- function(m, variant, mm, lb) {
  if (!is.null(m$error)) return(data.table(m = mm, label = lb, variant = variant, alpha = NA_real_, sd = NA_real_, mlik = NA_real_, nfunc = NA_real_, theta = NA_real_, mode_status = NA_integer_, error = substr(m$error, 1, 200)))
  data.table(m = mm, label = lb, variant = variant, alpha = m$summary.hyperpar[ALPHA_NAME, "mean"], sd = m$summary.hyperpar[ALPHA_NAME, "sd"],
             mlik = as.numeric(m$mlik[1, 1]), nfunc = as.numeric(m$misc$nfunc), theta = as.numeric(m$mode$theta[1]), mode_status = as.integer(m$mode$mode.status), error = NA_character_)
}

out <- list()
for (i in seq_len(nrow(flagged))) {
  mm <- flagged$m[i]; lb <- flagged$label[i]
  b <- readRDS(file.path(OUT_DIR, sprintf("baddata_%s_m%02d.rds", RUN_TAG, mm)))[[lb]]
  th_ref <- refstart[label == lb, theta_ref]
  cat(sprintf("\n[%d/%d] m=%d %s (rows %d, events %d, reference theta %.3f)\n", i, nrow(flagged), mm, lb, nrow(b$data), sum(b$data$status), th_ref))
  ex <- readRDS(file.path(OUT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, mm)))$fits[[lb]]
  out[[length(out) + 1]] <- if (is.null(ex$error)) data.table(m = mm, label = lb, variant = "orig", alpha = ex$summary_hyper[1, "mean"], sd = ex$summary_hyper[1, "sd"], mlik = as.numeric(ex$mlik[1, 1]),
                                                              nfunc = as.numeric(ex$nfunc), theta = as.numeric(ex$mode_theta[1]), mode_status = as.integer(ex$mode_status), error = NA_character_) else
                                                    data.table(m = mm, label = lb, variant = "orig", alpha = NA_real_, sd = NA_real_, mlik = NA_real_, nfunc = NA_real_, theta = NA_real_, mode_status = NA_integer_, error = substr(ex$error, 1, 200))
  out[[length(out) + 1]] <- summ(fit1(b$data, b$rhs), "A_default", mm, lb)
  out[[length(out) + 1]] <- summ(fit1(b$data, b$rhs, start = th_ref), "B_start_ref", mm, lb)
  out[[length(out) + 1]] <- summ(fit1(b$data, b$rhs, threads = "1:1"), "C_1thread", mm, lb)
  if (is.na(out[[length(out) - 3]]$alpha)) { cat("  -- verbose refit (INLA error in the pipeline) --\n"); r <- fit1(b$data, b$rhs, verbose = TRUE); out[[length(out) + 1]] <- summ(r, "D_verbose", mm, lb) }
  print(rbindlist(tail(out, 4))[, .(variant, alpha = round(alpha, 3), sd = signif(sd, 3), mlik = round(mlik, 1), nfunc, theta = round(theta, 2), error)])
}
res <- rbindlist(out); fwrite(res, file.path(OUT_DIR, sprintf("refit_flagged_%s.csv", RUN_TAG)))
cat("\n==== SUMMARY: variant A vs B (A should repeat the failure; B = repaired?) ====\n")
w <- dcast(res[variant %in% c("orig", "A_default", "B_start_ref")], m + label ~ variant, value.var = c("alpha", "sd", "mlik"))
print(w[, .(m, label, alpha_orig = round(alpha_orig, 3), alpha_A = round(alpha_A_default, 3), alpha_B = round(alpha_B_start_ref, 3),
            sd_B = signif(sd_B_start_ref, 3), mlik_A = round(mlik_A_default, 1), mlik_B = round(mlik_B_start_ref, 1), dmlik_B_minus_A = round(mlik_B_start_ref - mlik_A_default, 1))])
cat("Repaired = B has alpha SD ~0.05 (not ~1e-4) and a higher mlik than A. Not repaired = B repeats the collapsed fit.\n")
