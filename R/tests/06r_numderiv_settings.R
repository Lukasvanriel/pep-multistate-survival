# ==============================================================================
# 06r_numderiv_settings.R   RUN ON THE MACHINE WITH INLA (~ 30 min)
# INLA's control.inla documentation: the Newton step for the latent field uses NUMERICAL derivatives of the log-likelihood
# (step.len = step length, stencil = 3/5/7/9 points; default stencil 5) and a damping factor step.factor; tolerance.step controls
# the inner convergence.  The truncated Weibull log-likelihood has huge terms (lam*Tstop)^alpha - (lam*trunc)^alpha at the cold start,
# so inaccurate finite-difference derivatives could (a) leave the intercept at its start value and (b) make mlik jagged.
# These are native, global INLA settings (identical for every transition) - no change of model or method.
# TEST on the 8 saved collapsed black_spruce->regen datasets (pipeline formula, free alpha) + cold-start fixed-alpha 1.7 on m02 and m08:
#   base | stencil9 | steplen_div10 | steplen_x10 | stepfactor0.5 | tolstep_div100 | combo (stencil9 + steplen_div10 + tolstep_div100)
# Healthy = SD ratio ~ 1 (collapse < 0.1); fixed cold start healthy = mlik near -2170 (stuck ~ -25,000 or lower).
# Output: data/processed/fit/numderiv_settings_<tag>.csv
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(INLA) })
RUN_TAG <- "diag"; OUT_DIR <- here("data", "processed", "fit"); options(width = 200)
LABEL <- "black_spruce->regen"; ALPHA_NAME <- "alpha parameter for weibullsurv"; M_SET <- c(2, 8, 14, 15, 19, 23, 25, 29)
D <- INLA::inla.set.control.inla.default()
cat("INLA defaults: step.len =", format(D$step.len), "| stencil =", D$stencil, "| step.factor =", D$step.factor, "| tolerance.step =", format(D$tolerance.step), "| cmin =", D$cmin, "| h =", D$h, "\n")
sl <- if (is.null(D$step.len) || is.na(D$step.len)) 1e-3 else D$step.len
ts <- if (is.null(D$tolerance.step) || is.na(D$tolerance.step)) 1e-4 else D$tolerance.step
V <- list(base = list(), stencil9 = list(stencil = 9L), steplen_div10 = list(step.len = sl / 10), steplen_x10 = list(step.len = sl * 10),
          stepfactor0.5 = list(step.factor = 0.5), tolstep_div100 = list(tolerance.step = ts / 100),
          combo = list(stencil = 9L, step.len = sl / 10, tolerance.step = ts / 100))
run <- function(d, rhs, extra, fixed_alpha = NULL) {
  ci <- modifyList(list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)), extra)
  fam <- list(variant = 1); if (!is.null(fixed_alpha)) fam$hyper <- list(alpha = list(initial = 10 * log(fixed_alpha), fixed = TRUE))
  t0 <- proc.time()["elapsed"]
  r <- tryCatch(inla(as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", rhs)), family = "weibullsurv", data = d,
                     control.family = fam, control.compute = list(config = TRUE), control.inla = ci), error = function(e) conditionMessage(e))
  list(r = r, sec = proc.time()["elapsed"] - t0) }
out <- list()
for (mm in M_SET) {
  b <- readRDS(file.path(OUT_DIR, sprintf("baddata_%s_m%02d.rds", RUN_TAG, mm)))[[LABEL]]; ev <- sum(b$data$status == 1L)
  for (v in names(V)) {
    f <- run(b$data, b$rhs, V[[v]]); r <- f$r
    row <- if (is.character(r)) data.table(m = mm, setting = v, alpha = NA_real_, sd = NA_real_, ratio = NA_real_, mlik = NA_real_, fixed17_mlik = NA_real_, sec = f$sec, error = substr(r, 1, 60)) else {
      a <- r$summary.hyperpar[ALPHA_NAME, "mean"]; s <- r$summary.hyperpar[ALPHA_NAME, "sd"]
      data.table(m = mm, setting = v, alpha = a, sd = s, ratio = s / (0.78 * a / sqrt(ev)), mlik = as.numeric(r$mlik[1, 1]), fixed17_mlik = NA_real_, sec = f$sec, error = NA_character_) }
    if (mm %in% c(2, 8)) { g <- run(b$data, b$rhs, V[[v]], fixed_alpha = 1.7)$r; row$fixed17_mlik <- if (is.character(g)) NA_real_ else as.numeric(g$mlik[1, 1]) }
    out[[length(out) + 1]] <- row
    print(row[, .(m, setting, alpha = round(alpha, 3), sd = signif(sd, 3), ratio = round(ratio, 3), mlik = round(mlik, 1), fixed17_mlik = round(fixed17_mlik, 1), sec = round(sec), error)])
    fwrite(rbindlist(out), file.path(OUT_DIR, sprintf("numderiv_settings_%s.csv", RUN_TAG)))
  }
}
R <- rbindlist(out)
cat("\n==== SD ratio by setting (healthy ~1.2; collapse < 0.1) ====\n"); print(dcast(R, m ~ setting, value.var = "ratio")[, lapply(.SD, function(v) round(v, 3))])
cat("\nalpha by setting (MLE of these data: 1.81-1.96):\n"); print(dcast(R, m ~ setting, value.var = "alpha")[, lapply(.SD, function(v) round(v, 3))])
cat("\nfixed-alpha 1.7 cold start, mlik (healthy ~ -2170; stuck <= -25,000):\n"); print(dcast(R[m %in% c(2, 8)], m ~ setting, value.var = "fixed17_mlik")[, lapply(.SD, function(v) round(v, 1))])
