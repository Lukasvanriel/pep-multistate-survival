# ==============================================================================
# 06h_random_subsets.R   RUN ON THE MACHINE WITH INLA (several minutes)
# Does the INLA collapse on black_spruce->regen depend on sample size / on particular rows?
# Data: black_spruce->regen of imputation 2 (8651 rows, 809 events). Rows are independent likelihood terms
# (split-row likelihood), so a random subset of rows is a valid smaller data set (no plot id is stored in the fit data).
# 10 random subsets at 50% of the rows and 10 at 20% (fixed seed). For each subset and the full data:
#   (1) free INLA fit (pipeline settings): alpha, SD            vs   MLE alpha, SD (same data)
#   (2) fixed-alpha grid 1.3...2.1: cliff metric = [INLA mlik(a)-mlik(1.5)] - [MLE profile(a)-profile(1.5)]; < -50 = collapse
# Reading: collapse in all / most subsets = a property of the whole transition; only at some subsets = particular rows;
# collapse only in large subsets = a size effect.
# Output: data/processed/fit/random_subsets_<tag>.csv (+ _grid). Nothing else is modified.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(INLA) })
RUN_TAG <- "diag"; OUT_DIR <- here("data", "processed", "fit"); options(width = 200); set.seed(2026)
ALPHAS <- c(1.3, 1.5, 1.6, 1.7, 1.9, 2.1); REF <- 1.5; ALPHA_NAME <- "alpha parameter for weibullsurv"

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


mle_full <- function(x, tms) {
  X <- cbind(1, as.matrix(x[, tms, drop = FALSE])); y <- x$Tstop; tr <- x$truncation; st <- x$status
  nll <- function(p) { a <- exp(p[1]); eta <- drop(X %*% p[-1]); lam <- exp(eta); -sum(st * (log(a) + a * eta + (a - 1) * log(y)) - (lam * y)^a + (lam * tr)^a) }
  best <- NULL; for (s in c(0.5, 1.5)) { o <- tryCatch(optim(c(log(s + 0.5), -2, rep(0, ncol(X) - 1)), nll, method = "BFGS", hessian = TRUE, control = list(maxit = 2000, reltol = 1e-12)), error = function(e) NULL)
    if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o }
  V <- tryCatch(solve(best$hessian), error = function(e) matrix(NA, nrow(best$hessian), ncol(best$hessian))); a <- exp(best$par[1]); c(alpha = a, sd = a * sqrt(V[1, 1]))
}
free_fit <- function(x, tms) {
  args <- list(formula = as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", paste(tms, collapse = " + "))), family = "weibullsurv", data = x,
               control.family = list(variant = 1), control.compute = list(config = TRUE), control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)))
  r <- tryCatch(do.call(inla, args), error = function(e) { args$num.threads <- "1:1"; tryCatch(do.call(inla, args), error = function(e2) NULL) })
  if (is.null(r)) c(alpha = NA_real_, sd = NA_real_) else c(alpha = r$summary.hyperpar[ALPHA_NAME, "mean"], sd = r$summary.hyperpar[ALPHA_NAME, "sd"])
}

b <- readRDS(file.path(OUT_DIR, sprintf("baddata_%s_m02.rds", RUN_TAG)))[["black_spruce->regen"]]; tms <- strsplit(b$rhs, " \\+ ")[[1]]; full <- b$data
subsets <- c(list(full = full),
             setNames(lapply(1:10, function(i) full[sample(nrow(full), round(0.5 * nrow(full))), ]), sprintf("half_%02d", 1:10)),
             setNames(lapply(1:10, function(i) full[sample(nrow(full), round(0.2 * nrow(full))), ]), sprintf("fifth_%02d", 1:10)))
res <- list(); grid <- list()
for (sn in names(subsets)) { v <- subsets[[sn]]; pf <- prof_fun(v, tms)
  ml <- sapply(ALPHAS, function(a) mlik_fixed(v, tms, a)); pr <- sapply(ALPHAS, pf); metric <- (ml - ml[ALPHAS == REF]) - (pr - pr[ALPHAS == REF]); cliff <- ALPHAS[which(metric < -50)[1]]
  fi <- free_fit(v, tms); mo <- mle_full(v, tms)
  cat(sprintf("%-9s rows %5d ev %4d | free INLA alpha %.3f sd %.4g | MLE alpha %.3f sd %.3f | z %.1f sd-ratio %.2f | cliff metric [%s] CLIFF_AT %s\n", sn, nrow(v), sum(v$status),
              fi["alpha"], fi["sd"], mo["alpha"], mo["sd"], (fi["alpha"] - mo["alpha"]) / mo["sd"], fi["sd"] / mo["sd"], paste(round(metric), collapse = " "), ifelse(is.na(cliff), "none", cliff)))
  res[[sn]] <- data.table(subset = sn, rows = nrow(v), events = sum(v$status), inla_alpha = fi["alpha"], inla_sd = fi["sd"], mle_alpha = mo["alpha"], mle_sd = mo["sd"], cliff_at = cliff)
  grid[[sn]] <- data.table(subset = sn, alpha = ALPHAS, mlik = ml, prof = pr, metric = metric)
  fwrite(rbindlist(res), file.path(OUT_DIR, sprintf("random_subsets_%s.csv", RUN_TAG))); fwrite(rbindlist(grid), file.path(OUT_DIR, sprintf("random_subsets_grid_%s.csv", RUN_TAG))) }
r <- rbindlist(res); r[, kind := sub("_.*", "", subset)]
cat("\n==== SUMMARY: share of subsets with a collapse ====\n"); print(r[, .(n = .N, n_collapse = sum(!is.na(cliff_at)), median_events = median(events)), by = kind])
