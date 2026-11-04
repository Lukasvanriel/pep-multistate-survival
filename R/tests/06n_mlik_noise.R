# ==============================================================================
# 06n_mlik_noise.R   RUN ON THE MACHINE WITH INLA (~ 15 min)
# 06m: INLA's log posterior looks jagged (+-0.2-0.5) over alpha steps of 1e-4, where a smooth function would move ~1e-4.
# TEST: INLA fixed-alpha log marginal likelihood on a fine grid (alpha 1.690-1.710, step 0.001; the smooth second difference
# there is ~4e-4) for the saved collapsed data of m02 and m08 (pipeline settings, full formula), in two versions:
#   with_trunc  = data as they are
#   no_trunc    = truncation set to 0 (Tstop unchanged)
# noise = second differences d2 = mlik[i+1] - 2*mlik[i] + mlik[i-1]; reported: SD(d2)/sqrt(6) = per-evaluation noise, and the
# intercept INLA found (stuck ones sit near 0 instead of ~ -2.3).
# Output: data/processed/fit/mlik_noise_<tag>.csv
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(INLA) })
RUN_TAG <- "diag"; OUT_DIR <- here("data", "processed", "fit"); options(width = 200)
LABEL <- "black_spruce->regen"; GRID <- seq(1.690, 1.710, by = 0.001)
fit_fixed <- function(d, rhs, a) {
  args <- list(formula = as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", rhs)), family = "weibullsurv", data = d,
               control.family = list(variant = 1, hyper = list(alpha = list(initial = 10 * log(a), fixed = TRUE))), control.compute = list(config = TRUE),
               control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)))
  r <- tryCatch(do.call(inla, args), error = function(e) NULL)
  if (is.null(r)) return(c(NA_real_, NA_real_))
  c(as.numeric(r$mlik[1, 1]), r$summary.fixed["(Intercept)", "mean"])
}
out <- list()
for (mm in c(2, 8)) {
  b <- readRDS(file.path(OUT_DIR, sprintf("baddata_%s_m%02d.rds", RUN_TAG, mm)))[[LABEL]]
  for (ver in c("with_trunc", "no_trunc")) {
    d <- b$data; if (ver == "no_trunc") d$truncation <- 0
    for (a in GRID) { r <- fit_fixed(d, b$rhs, a); out[[length(out) + 1]] <- data.table(m = mm, version = ver, alpha = a, mlik = r[1], intercept = r[2]) }
    cat(sprintf("m%02d %s done\n", mm, ver)); fwrite(rbindlist(out), file.path(OUT_DIR, sprintf("mlik_noise_%s.csv", RUN_TAG)))
  }
}
R <- rbindlist(out); setorder(R, m, version, alpha)
R[, d2 := shift(mlik, -1) - 2 * mlik + shift(mlik, 1), by = .(m, version)]
cat("\n==== per-evaluation noise of INLA's mlik (smooth function: d2 ~ 4e-4) ====\n")
print(R[, .(n_ok = sum(!is.na(mlik)), mlik_range = round(diff(range(mlik, na.rm = TRUE)), 2), sd_d2 = round(sd(d2, na.rm = TRUE), 3), noise_sd = round(sd(d2, na.rm = TRUE) / sqrt(6), 3),
            n_intercept_near0 = sum(abs(intercept) < 0.5, na.rm = TRUE), n_fail = sum(is.na(mlik))), by = .(m, version)])
cat("\nmlik by alpha:\n"); print(dcast(R, alpha ~ m + version, value.var = "mlik")[, lapply(.SD, round, 2)])
