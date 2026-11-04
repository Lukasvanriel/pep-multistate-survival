# ==============================================================================
# 06g_truncation_hazard_test.R   RUN ON THE MACHINE WITH INLA (a few minutes)
# Hypothesis: INLA's truncated-Weibull likelihood collapses when the cumulative hazard at the truncation time,
#   H_tr = (lambda_i * truncation_i)^alpha  (lambda from the maximum-likelihood fit), is large (black_spruce->regen:
#   largest of 91 transitions, up to ~19, i.e. truncation probability ~ exp(-19)).
# Test: black_spruce->regen (imputations 2, 8). Rows with H_tr > cutoff are DROPPED (cutoffs 10, 5, 3, 1), the same
# fixed-alpha grid as 06f is run on the reduced data, and INLA's curve is compared with the MLE profile of the SAME reduced data.
#   cliff metric = [INLA mlik(a) - INLA mlik(1.5)] - [MLE profile(a) - MLE profile(1.5)]; < -50 = INLA collapsed.
# Negative control: the same NUMBER of randomly chosen rows is dropped instead (dropping rows alone should not cure it).
# Reading: collapse disappears below some cutoff, and not in the control -> large truncation hazard is the trigger.
# Output: data/processed/fit/trunc_hazard_test_<tag>.csv. Nothing else is modified. Dropping rows is a diagnostic only.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(INLA) })
RUN_TAG <- "diag"; OUT_DIR <- here("data", "processed", "fit"); options(width = 200); set.seed(1)
ALPHAS <- c(1.3, 1.5, 1.6, 1.7, 1.9, 2.1); REF <- 1.5

prof_fun <- function(x, tms) {
  X <- if (length(tms)) cbind(1, as.matrix(x[, tms, drop = FALSE])) else matrix(1, nrow(x), 1)
  y <- x$Tstop; tr <- x$truncation; st <- x$status
  nll <- function(b, a) { eta <- drop(X %*% b); lam <- exp(eta); -sum(st * (log(a) + a * eta + (a - 1) * log(y)) - (lam * y)^a + (lam * tr)^a) }
  function(a) { best <- Inf; for (i0 in c(-2, -4)) { o <- tryCatch(optim(c(i0, rep(0, ncol(X) - 1)), nll, a = a, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12)), error = function(e) NULL)
    if (!is.null(o) && o$value < best) best <- o$value }; -best }
}
mlik_fixed <- function(x, tms, a) {
  rhs <- if (length(tms)) paste(tms, collapse = " + ") else "1"
  args <- list(formula = as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", rhs)), family = "weibullsurv", data = x,
               control.family = list(variant = 1, hyper = list(alpha = list(initial = 10 * log(a), fixed = TRUE))), control.compute = list(config = TRUE),
               control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)))
  r <- tryCatch(do.call(inla, args), error = function(e) { args$num.threads <- "1:1"; tryCatch(do.call(inla, args), error = function(e2) NULL) })
  if (is.null(r) || !isTRUE(r$all.hyper$family[[1]]$hyper[[1]]$fixed)) NA_real_ else as.numeric(r$mlik[1, 1])
}


mle_H <- function(x, tms) {   # full MLE, returns cumulative hazard at truncation per row
  X <- cbind(1, as.matrix(x[, tms, drop = FALSE])); y <- x$Tstop; tr <- x$truncation; st <- x$status
  nll <- function(p) { a <- exp(p[1]); eta <- drop(X %*% p[-1]); lam <- exp(eta); -sum(st * (log(a) + a * eta + (a - 1) * log(y)) - (lam * y)^a + (lam * tr)^a) }
  best <- NULL; for (s in c(0.5, 1.5)) { o <- tryCatch(optim(c(log(s + 0.5), -2, rep(0, ncol(X) - 1)), nll, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12)), error = function(e) NULL)
    if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o }
  a <- exp(best$par[1]); list(alpha = a, H = (exp(drop(X %*% best$par[-1])) * tr)^a) }

out <- list()
for (m in c(2, 8)) {
  b <- readRDS(file.path(OUT_DIR, sprintf("baddata_%s_m%02d.rds", RUN_TAG, m)))[["black_spruce->regen"]]; tms <- strsplit(b$rhs, " \\+ ")[[1]]
  fh <- mle_H(b$data, tms); d0 <- b$data; d0$H <- fh$H
  cat(sprintf("\n===== m%02d black_spruce->regen: MLE alpha %.3f, rows %d, events %d; H_tr quantiles: %s =====\n", m, fh$alpha, nrow(d0), sum(d0$status),
              paste(round(quantile(fh$H, c(.5, .9, .99, 1)), 2), collapse = " / ")))
  variants <- list(all_rows = d0)
  for (cut in c(10, 5, 3, 1)) {
    keep <- d0$H <= cut; variants[[sprintf("drop_H>%g", cut)]] <- d0[keep, ]
    variants[[sprintf("control_%g", cut)]] <- d0[sample(nrow(d0), sum(keep)), ]    # same number of rows kept, chosen at random
  }
  for (vn in names(variants)) { v <- variants[[vn]][, setdiff(names(variants[[vn]]), "H")]; if (sum(v$status) < 20) next; pf <- prof_fun(v, tms)
    ml <- sapply(ALPHAS, function(a) mlik_fixed(v, tms, a)); pr <- sapply(ALPHAS, pf)
    metric <- (ml - ml[ALPHAS == REF]) - (pr - pr[ALPHAS == REF]); cliff <- ALPHAS[which(metric < -50)[1]]
    cat(sprintf("%-14s rows %5d events %4d | metric by alpha [%s] = [%s]  CLIFF_AT = %s\n", vn, nrow(v), sum(v$status), paste(ALPHAS, collapse = " "), paste(round(metric), collapse = " "), ifelse(is.na(cliff), "none", cliff)))
    out[[length(out) + 1]] <- data.table(m = m, variant = vn, rows = nrow(v), events = sum(v$status), alpha = ALPHAS, mlik = ml, prof = pr, metric = metric, cliff_at = cliff)
    fwrite(rbindlist(out), file.path(OUT_DIR, sprintf("trunc_hazard_test_%s.csv", RUN_TAG))) }
}
