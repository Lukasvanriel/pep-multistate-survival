# ==============================================================================
# 06q_truncation_dose.R   RUN ON THE MACHINE WITH INLA (~ 20 min)
# 06n: truncation = 0 (inla.surv's default, same code path) -> INLA fine; actual truncation -> intercept stuck at 0.
# Is it the PRESENCE of a truncation term, its SIZE, or the NUMBER of truncated rows?  Saved collapsed data of m02 and m08, variants of the
# truncation column (Tstop, status, covariates untouched; the likelihood changes, this is a diagnostic only):
#   zero          0 for all rows (= no truncation)
#   const_1e-7    1e-7 for all rows
#   const_0.01 / const_0.5      constants
#   scale_0.01 / 0.1 / 0.25 / 0.5 / 0.75 / 1    actual truncation x factor   (1 = as in the data)
#   cap_1 / cap_2 / cap_4       pmin(actual, cap)   (cap in decades; limits the SIZE, keeps all truncated rows)
#   frac_0.05 / 0.2 / 0.5       actual truncation kept on a random fraction of rows, 0 elsewhere (limits the NUMBER)
# For each variant: (a) fixed-alpha (1.7) cold-start fit: does INLA's intercept reach our optimum? stuck = gap > 0.5
#                   (b) free-alpha fit: alpha, SD ratio (collapse = ratio < 0.1)
# Output: data/processed/fit/truncation_dose_<tag>.csv
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(INLA) })
RUN_TAG <- "diag"; OUT_DIR <- here("data", "processed", "fit"); options(width = 200); set.seed(2026)
LABEL <- "black_spruce->regen"; ALPHA_NAME <- "alpha parameter for weibullsurv"; A_FIX <- 1.7
dev_own <- function(b0, a, x) { lam <- exp(b0); y <- x$Tstop; tr <- x$truncation; st <- x$status
  -2 * (st * (log(a) + a * b0 + (a - 1) * log(y)) - (lam * y)^a + (lam * tr)^a) }
opt_b <- function(a, x) optimize(function(b0) sum(dev_own(b0, a, x)), c(-8, 3), tol = 1e-10)$minimum
mk <- function(d, rhs, fam_hyper) {
  r <- tryCatch(inla(as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", rhs)), family = "weibullsurv", data = d,
                     control.family = c(list(variant = 1), fam_hyper), control.compute = list(config = TRUE),
                     control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE))), error = function(e) NULL); r }
variants <- function(tr) {
  n <- length(tr); v <- list(zero = rep(0, n), `const_1e-7` = rep(1e-7, n), const_0.01 = rep(0.01, n), const_0.5 = rep(0.5, n))
  for (s in c(0.01, 0.1, 0.25, 0.5, 0.75, 1)) v[[sprintf("scale_%g", s)]] <- tr * s
  for (cp in c(1, 2, 4)) v[[sprintf("cap_%g", cp)]] <- pmin(tr, cp)
  for (f in c(0.05, 0.2, 0.5)) v[[sprintf("frac_%g", f)]] <- ifelse(runif(n) < f, tr, 0)
  v }
out <- list()
for (mm in c(2, 8)) {
  b <- readRDS(file.path(OUT_DIR, sprintf("baddata_%s_m%02d.rds", RUN_TAG, mm)))[[LABEL]]; V <- variants(b$data$truncation)
  for (nm in names(V)) {
    d <- b$data; d$truncation <- V[[nm]]; ev <- sum(d$status == 1L)
    rf <- mk(d, b$rhs, list(hyper = list(alpha = list(initial = 10 * log(A_FIX), fixed = TRUE))))
    gap <- if (is.null(rf)) NA_real_ else rf$summary.fixed["(Intercept)", "mean"] - opt_b(A_FIX, d)
    ml <- if (is.null(rf)) NA_real_ else as.numeric(rf$mlik[1, 1])
    rr <- mk(d, b$rhs, list())
    a <- if (is.null(rr)) NA_real_ else rr$summary.hyperpar[ALPHA_NAME, "mean"]; s <- if (is.null(rr)) NA_real_ else rr$summary.hyperpar[ALPHA_NAME, "sd"]
    out[[length(out) + 1]] <- data.table(m = mm, variant = nm, share_trunc_pos = mean(d$truncation > 0), trunc_max = max(d$truncation), fixed_gap = gap, fixed_stuck = !is.na(gap) & abs(gap) > 0.5,
                                         fixed_mlik = ml, free_alpha = a, free_sd = s, free_ratio = s / (0.78 * a / sqrt(ev)))
    print(tail(rbindlist(out), 1)[, .(m, variant, share_trunc_pos = round(share_trunc_pos, 2), trunc_max = round(trunc_max, 2), fixed_gap = round(fixed_gap, 3), fixed_stuck, free_alpha = round(free_alpha, 3), free_ratio = round(free_ratio, 3))])
    fwrite(rbindlist(out), file.path(OUT_DIR, sprintf("truncation_dose_%s.csv", RUN_TAG)))
  }
}
R <- rbindlist(out)
cat("\n==== fixed-alpha stuck (1) / ok (0) and free-fit SD ratio by variant ====\n")
print(dcast(R, variant ~ m, value.var = c("fixed_stuck", "free_ratio"))[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)])
