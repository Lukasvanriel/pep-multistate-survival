# ==============================================================================
# 06s_laplace_settings.R   RUN ON THE MACHINE WITH INLA (~ 10 min)
# 06r: no step.len / stencil / step.factor / tolerance.step setting repairs the free-fit collapse (stencil = 9 repairs the fixed-alpha
# cold start of m02 only) -> the stuck intercept and the free-fit collapse are different problems.  The collapse comes with a jagged
# log marginal likelihood (+-0.5 over alpha steps of 1e-4, 06m).  control.inla documentation: `fast` (default TRUE) replaces conditional
# MODES in the Laplace approximation by conditional EXPECTATIONS; `correct` (default FALSE) adds a correction to the Laplace approximation;
# `cmin` clips the negative likelihood Hessian (default -Inf; diagnostic only, can bias).
# TEST on the 8 saved collapsed black_spruce->regen datasets (pipeline formula, free alpha) + cold-start fixed-alpha 1.7 on m02, m08:
#   base | fast_FALSE | fast_FALSE_stencil9 | correct_TRUE | fast_FALSE_correct_TRUE | cmin_0 (diagnostic)
# Healthy = SD ratio ~ 1 (collapse < 0.1).  Output: data/processed/fit/laplace_settings_<tag>.csv
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(INLA) })
RUN_TAG <- "diag"; OUT_DIR <- here("data", "processed", "fit"); options(width = 200)
LABEL <- "black_spruce->regen"; ALPHA_NAME <- "alpha parameter for weibullsurv"; M_SET <- c(2, 8, 14, 15, 19, 23, 25, 29)
V <- list(base = list(), fast_FALSE = list(fast = FALSE), fast_FALSE_stencil9 = list(fast = FALSE, stencil = 9L), correct_TRUE = list(correct = TRUE),
          fast_FALSE_correct_TRUE = list(fast = FALSE, correct = TRUE), cmin_0 = list(cmin = 0))
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
    fwrite(rbindlist(out), file.path(OUT_DIR, sprintf("laplace_settings_%s.csv", RUN_TAG)))
  }
}
R <- rbindlist(out)
cat("\n==== SD ratio by setting (healthy ~1.2; collapse < 0.1) ====\n"); print(dcast(R, m ~ setting, value.var = "ratio")[, lapply(.SD, function(v) round(v, 3))])
cat("\nalpha by setting (MLE of these data: 1.81-1.96):\n"); print(dcast(R, m ~ setting, value.var = "alpha")[, lapply(.SD, function(v) round(v, 3))])
cat("\nfixed-alpha 1.7 cold start, mlik (healthy ~ -2170; stuck <= -25,000):\n"); print(dcast(R[m %in% c(2, 8)], m ~ setting, value.var = "fixed17_mlik")[, lapply(.SD, function(v) round(v, 1))])
