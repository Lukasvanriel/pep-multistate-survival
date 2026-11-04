# ==============================================================================
# 06m_integration_points.R   RUN ON THE MACHINE WITH INLA (~ 3 min)
# 06l: with alpha FIXED, INLA's intercept search sometimes never leaves its start (0) -> log-lik lost ~ 20,000.
# Hypothesis: in the FREE fit the integration points next to the mode hit such failed states, so the log posterior of alpha
# has huge fake curvature -> alpha SD ~ 0.0002 and the mode sits at the edge. Test on the saved collapsed fitting data
# (pipeline settings, original formula incl. covariates), m02 and m08: for every alpha integration point INLA used print
#   theta, alpha, log.posterior, intercept (posterior mean in the config), and our intercept-only optimum at that alpha
# as a rough yardstick (covariates are z-scaled, so intercepts should be close). Failed points = intercept near 0 and a
# log.posterior drop of thousands.
# Output: data/processed/fit/integration_points_<tag>.csv
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(INLA) })
RUN_TAG <- "diag"; OUT_DIR <- here("data", "processed", "fit"); options(width = 200)
LABEL <- "black_spruce->regen"; ALPHA_NAME <- "alpha parameter for weibullsurv"
dev_own <- function(b, a, x) { lam <- exp(b); y <- x$Tstop; tr <- x$truncation; st <- x$status
  -2 * (st * (log(a) + a * b + (a - 1) * log(y)) - (lam * y)^a + (lam * tr)^a) }
opt_b <- function(a, x) optimize(function(b) sum(dev_own(b, a, x)), c(-8, 3), tol = 1e-10)$minimum
out <- list()
for (mm in c(2, 8)) {
  b <- readRDS(file.path(OUT_DIR, sprintf("baddata_%s_m%02d.rds", RUN_TAG, mm)))[[LABEL]]
  r <- inla(as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", b$rhs)), family = "weibullsurv", data = b$data,
            control.family = list(variant = 1), control.compute = list(config = TRUE),
            control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)))
  cat(sprintf("\nm%02d free fit: alpha %.4f sd %.2e | mode theta %.3f | mode_status %s | nfunc %s | mlik %.1f\n", mm,
              r$summary.hyperpar[ALPHA_NAME, "mean"], r$summary.hyperpar[ALPHA_NAME, "sd"], r$mode$theta[1], r$mode$mode.status, r$misc$nfunc, r$mlik[1, 1]))
  cf <- r$misc$configs; ct <- cf$contents; ii <- which(ct$tag == "(Intercept)"); i0 <- ct$start[ii]
  x <- b$data
  d <- rbindlist(lapply(seq_len(cf$nconfig), function(i) { cc <- cf$config[[i]]; th <- cc$theta[1]; a <- exp(th / 10)
    data.table(m = mm, config = i, theta = th, alpha = a, log_post = cc$log.posterior, intercept_inla = cc$mean[i0], intercept_own = opt_b(a, x)) }))
  d[, `:=`(drop_vs_best = log_post - max(log_post), intercept_gap = intercept_inla - intercept_own)]
  setorder(d, theta); print(d[, lapply(.SD, function(v) round(v, 3)), .SDcols = c("config", "theta", "alpha", "log_post", "drop_vs_best", "intercept_inla", "intercept_own", "intercept_gap")])
  out[[length(out) + 1]] <- d
}
fwrite(rbindlist(out), file.path(OUT_DIR, sprintf("integration_points_%s.csv", RUN_TAG)))
