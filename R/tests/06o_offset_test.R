# ==============================================================================
# 06o_offset_test.R   RUN ON THE MACHINE WITH INLA (~ 20 min)
# 06n: with left truncation INLA's intercept search stays at its start value 0 (true intercept ~ -3.8) in fixed-alpha fits.
# TEST of an exact reparametrisation: linear predictor = off + (intercept' + covariates). The default intercept prior is flat
# (prec.intercept = 0), so posterior(intercept) = posterior(intercept') + off; likelihood, alpha, betas unchanged.
# PART 1  free-alpha fits of the 8 saved collapsed black_spruce->regen datasets with off in {0 (as pipeline), -2, -3, rule}
#         rule = argmax of the intercept-only Weibull loglik (alpha 1.5, truncation included), computed from the data.
#         Expect if the fix works: alpha near the MLE (1.81-1.96) and SD ratio ~ 1 (collapse = ratio < 0.1).
# PART 2  invariance: on 3 random 20% row subsets of m02, fit with off = 0 and off = rule; where off = 0 does NOT fail, the
#         posteriors must agree (alpha, sd, intercept_total = intercept' + off).
# Output: data/processed/fit/offset_test_<tag>.csv
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(INLA) })
RUN_TAG <- "diag"; OUT_DIR <- here("data", "processed", "fit"); options(width = 200)
LABEL <- "black_spruce->regen"; ALPHA_NAME <- "alpha parameter for weibullsurv"; M_SET <- c(2, 8, 14, 15, 19, 23, 25, 29)
dev_own <- function(b, a, x) { lam <- exp(b); y <- x$Tstop; tr <- x$truncation; st <- x$status
  -2 * (st * (log(a) + a * b + (a - 1) * log(y)) - (lam * y)^a + (lam * tr)^a) }
off_rule <- function(x, a = 1.5) optimize(function(b) sum(dev_own(b, a, x)), c(-8, 3), tol = 1e-6)$minimum
fit_off <- function(d, rhs, off) {
  d$off <- off
  f <- as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s + offset(off)", rhs))
  r <- tryCatch(inla(f, family = "weibullsurv", data = d, control.family = list(variant = 1), control.compute = list(config = TRUE),
                     control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE))), error = function(e) conditionMessage(e))
  ev <- sum(d$status == 1L)
  if (is.character(r)) return(data.table(off = off, alpha = NA_real_, sd = NA_real_, ratio = NA_real_, intercept_total = NA_real_, mlik = NA_real_, mode_status = NA_integer_, error = substr(r, 1, 60)))
  a <- r$summary.hyperpar[ALPHA_NAME, "mean"]; s <- r$summary.hyperpar[ALPHA_NAME, "sd"]
  data.table(off = off, alpha = a, sd = s, ratio = s / (0.78 * a / sqrt(ev)), intercept_total = r$summary.fixed["(Intercept)", "mean"] + off,
             mlik = as.numeric(r$mlik[1, 1]), mode_status = as.integer(r$mode$mode.status), error = NA_character_)
}
tu <- fread(file.path(OUT_DIR, sprintf("time_unit_test_%s.csv", RUN_TAG)))[, .(mle_alpha = mle_alpha[1]), by = .(m = as.integer(sub("^m(\\d+) .*", "\\1", set)), lab = sub("^m\\d+ ", "", set))]
res <- list()
for (mm in M_SET) {
  b <- readRDS(file.path(OUT_DIR, sprintf("baddata_%s_m%02d.rds", RUN_TAG, mm)))[[LABEL]]; r0 <- off_rule(b$data)
  for (o in c(0, -2, -3, r0)) {
    x <- fit_off(b$data, b$rhs, o); res[[length(res) + 1]] <- cbind(data.table(part = 1L, m = mm, variant = if (o == r0) "rule" else if (o == 0) "none" else as.character(o)), x)
    print(tail(rbindlist(res), 1)[, .(m, variant, off = round(off, 2), alpha = round(alpha, 3), sd = signif(sd, 3), ratio = round(ratio, 3), intercept_total = round(intercept_total, 3), mlik = round(mlik, 1), error)])
    fwrite(rbindlist(res, fill = TRUE), file.path(OUT_DIR, sprintf("offset_test_%s.csv", RUN_TAG)))
  }
}
set.seed(2026); b <- readRDS(file.path(OUT_DIR, sprintf("baddata_%s_m02.rds", RUN_TAG)))[[LABEL]]
for (i in 1:3) {
  d <- b$data[sample.int(nrow(b$data), round(0.2 * nrow(b$data))), ]; r0 <- off_rule(d)
  for (o in c(0, r0)) { x <- fit_off(d, b$rhs, o); res[[length(res) + 1]] <- cbind(data.table(part = 2L, m = i, variant = if (o == 0) "none" else "rule"), x)
    print(tail(rbindlist(res), 1)[, .(subset = m, variant, off = round(off, 2), alpha = round(alpha, 4), sd = signif(sd, 4), intercept_total = round(intercept_total, 4), mlik = round(mlik, 2), error)]) }
}
R <- rbindlist(res, fill = TRUE); fwrite(R, file.path(OUT_DIR, sprintf("offset_test_%s.csv", RUN_TAG)))
cat("\n==== PART 1: SD ratio by offset (collapse = ratio < 0.1; expected ~1) ====\n")
print(dcast(R[part == 1L], m ~ variant, value.var = "ratio")[, lapply(.SD, function(v) round(v, 3))])
cat("\nalpha by offset (MLE of these data: 1.81-1.96):\n"); print(dcast(R[part == 1L], m ~ variant, value.var = "alpha")[, lapply(.SD, function(v) round(v, 3))])
cat("\n==== PART 2: invariance (off = 0 vs rule, same data) ====\n")
print(dcast(R[part == 2L], m ~ variant, value.var = c("alpha", "sd", "intercept_total"))[, lapply(.SD, function(v) round(v, 5))])
