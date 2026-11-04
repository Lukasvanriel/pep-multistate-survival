# ==============================================================================
# 06l_inla_vs_own_loglik.R   RUN ON THE MACHINE WITH INLA (~ 5 min)
# Where does INLA go wrong at fixed alpha past the cliff?  Intercept-only model on the saved collapsed datasets (m02, m08).
# For alpha in ALPHAS, INLA fit with alpha FIXED (dic = TRUE), then, against our own Weibull loglik at the SAME alpha:
#   (1) is INLA's intercept the true optimum?   (intercept_inla vs intercept_opt; our loglik at both -> loss)
#   (2) does INLA's per-row posterior-mean deviance (local.dic - local.p.eff) agree with our -2*loglik_i at INLA's intercept?
#       constants may differ, so rows are compared as CHANGES from the reference alpha (1.5): delta_i = dD_INLA_i - dD_own_i
#       (sum over rows should be ~0 if INLA's likelihood is evaluated correctly; a few huge rows or a pattern in truncation/Tstop = the bug)
# Output: data/processed/fit/inla_vs_own_summary_<tag>.csv, inla_vs_own_rows_<tag>.csv
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(INLA) })
RUN_TAG <- "diag"; OUT_DIR <- here("data", "processed", "fit"); options(width = 200)
LABEL <- "black_spruce->regen"; REF <- 1.5; ALPHAS <- c(1.3, 1.5, 1.6, 1.7, 1.8, 1.9)
dev_own <- function(b, a, x) { eta <- b; lam <- exp(eta); y <- x$Tstop; tr <- x$truncation; st <- x$status
  -2 * (st * (log(a) + a * eta + (a - 1) * log(y)) - (lam * y)^a + (lam * tr)^a) }
opt_b <- function(a, x) optimize(function(b) sum(dev_own(b, a, x)), c(-8, 3), tol = 1e-10)$minimum
fit_fixed <- function(x, a) {
  args <- list(formula = inla.surv(time = Tstop, event = status, truncation = truncation) ~ 1, family = "weibullsurv", data = x,
               control.family = list(variant = 1, hyper = list(alpha = list(initial = 10 * log(a), fixed = TRUE))),
               control.compute = list(dic = TRUE, config = TRUE),
               control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)))
  tryCatch(do.call(inla, args), error = function(e) { args$num.threads <- "1:1"; tryCatch(do.call(inla, args), error = function(e2) NULL) })
}
summ <- list(); rows <- list()
for (mm in c(2, 8)) {
  x <- as.data.frame(readRDS(file.path(OUT_DIR, sprintf("baddata_%s_m%02d.rds", RUN_TAG, mm)))[[LABEL]]$data)[, c("status", "Tstop", "truncation")]
  ref <- NULL
  for (a in ALPHAS) {
    r <- fit_fixed(x, a); if (is.null(r)) { cat("m", mm, "alpha", a, ": fit failed\n"); next }
    b_inla <- r$summary.fixed["(Intercept)", "mean"]; b_opt <- opt_b(a, x)
    Dbar_i <- r$dic$local.dic - r$dic$local.p.eff; own_i <- dev_own(b_inla, a, x)
    if (a == REF) ref <- list(Dbar = Dbar_i, own = own_i)
    summ[[length(summ) + 1]] <- data.table(m = mm, alpha = a, intercept_inla = b_inla, intercept_opt = b_opt, inla_sd = r$summary.fixed["(Intercept)", "sd"],
      loglik_own_at_inla = -sum(own_i) / 2, loglik_own_at_opt = -sum(dev_own(b_opt, a, x)) / 2, mlik = as.numeric(r$mlik[1, 1]),
      inla_mean_deviance = r$dic$mean.deviance, own_deviance_at_inla = sum(own_i), mode_status = r$mode$mode.status)
    print(tail(rbindlist(summ), 1)[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)])
    if (!is.null(ref) && a != REF) {
      d <- data.table(m = mm, alpha = a, row = seq_len(nrow(x)), status = x$status, Tstop = x$Tstop, truncation = x$truncation,
                      delta = (Dbar_i - ref$Dbar) - (own_i - ref$own))
      rows[[length(rows) + 1]] <- d
    }
  }
}
S <- rbindlist(summ); fwrite(S, file.path(OUT_DIR, sprintf("inla_vs_own_summary_%s.csv", RUN_TAG)))
R <- rbindlist(rows); fwrite(R, file.path(OUT_DIR, sprintf("inla_vs_own_rows_%s.csv", RUN_TAG)))
cat("\n==== (1) INLA intercept vs true optimum; loglik lost by INLA's intercept ====\n")
print(S[, .(m, alpha, intercept_inla = round(intercept_inla, 4), intercept_opt = round(intercept_opt, 4), loglik_lost = round(loglik_own_at_opt - loglik_own_at_inla, 2), mlik = round(mlik, 1), mode_status)])
cat("\n==== (2) per-row disagreement vs reference alpha (delta_i; sum should be ~0) ====\n")
print(R[, .(sum_delta = round(sum(delta), 1), max_abs = round(max(abs(delta)), 2), n_abs_gt_1 = sum(abs(delta) > 1), n_rows = .N), by = .(m, alpha)])
cat("\nby status / truncation>0 (alpha 1.7):\n"); print(R[alpha == 1.7, .(sum_delta = round(sum(delta), 1), mean_delta = round(mean(delta), 4), n = .N), by = .(m, status, trunc_pos = truncation > 0)])
cat("\nby Tstop decile (alpha 1.7, m02):\n"); print(R[alpha == 1.7 & m == 2, .(sum_delta = round(sum(delta), 1), mean_delta = round(mean(delta), 4), Tstop_max = round(max(Tstop), 2)), by = .(decile = cut(Tstop, quantile(Tstop, 0:10 / 10), include.lowest = TRUE, labels = FALSE))][order(decile)])
