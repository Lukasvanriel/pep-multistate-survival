# ==============================================================================
# 06f_find_cliff_trigger.R   RUN ON THE MACHINE WITH INLA (a few minutes)
# What in the DATA of black_spruce->regen / other_temperate->white_birch triggers INLA's numerical
# collapse (marginal likelihood falling by thousands once alpha is fixed above ~1.5)?
# Variants of the same data: full model; each covariate dropped in turn; intercept only; left truncation
# removed (truncation = 0). For every variant, alpha is fixed (key "alpha", theta = 10*log(alpha)) at 1.3...2.1 and
#   cliff metric = [INLA mlik(a) - INLA mlik(1.5)] - [MLE profile(a) - MLE profile(1.5)]
# The model is the same Weibull likelihood, so the two differences must agree (up to a smooth offset of a few units).
# Large negative metric (< -50) = INLA numerically collapsed at alpha = a. CLIFF_AT = smallest such alpha.
# Output: data/processed/fit/cliff_trigger_<tag>.csv + printed table. Nothing else is modified.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(INLA) })
RUN_TAG <- "diag"; OUT_DIR <- here("data", "processed", "fit"); options(width = 200)
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

targets <- list(list(m = 2, label = "black_spruce->regen"), list(m = 8, label = "black_spruce->regen"), list(m = 24, label = "other_temperate->white_birch"))
out <- list()
for (tg in targets) {
  b <- readRDS(file.path(OUT_DIR, sprintf("baddata_%s_m%02d.rds", RUN_TAG, tg$m)))[[tg$label]]; full <- strsplit(b$rhs, " \\+ ")[[1]]
  variants <- c(list(full = list(d = b$data, tms = full)),
                setNames(lapply(full, function(v) list(d = b$data, tms = setdiff(full, v))), paste0("drop_", full)),
                list(intercept_only = list(d = b$data, tms = character(0)),
                     no_truncation  = list(d = transform(b$data, truncation = 0), tms = full)))
  cat(sprintf("\n===== m%02d %s (rows %d, events %d) =====\n", tg$m, tg$label, nrow(b$data), sum(b$data$status)))
  for (vn in names(variants)) { v <- variants[[vn]]; pf <- prof_fun(v$d, v$tms)
    ml <- sapply(ALPHAS, function(a) mlik_fixed(v$d, v$tms, a)); pr <- sapply(ALPHAS, pf)
    metric <- (ml - ml[ALPHAS == REF]) - (pr - pr[ALPHAS == REF]); cliff <- ALPHAS[which(metric < -50)[1]]
    cat(sprintf("%-18s metric by alpha [%s] = [%s]  CLIFF_AT = %s\n", vn, paste(ALPHAS, collapse = " "), paste(round(metric), collapse = " "), ifelse(is.na(cliff), "none", cliff)))
    out[[length(out) + 1]] <- data.table(m = tg$m, label = tg$label, variant = vn, alpha = ALPHAS, mlik = ml, prof = pr, metric = metric, cliff_at = cliff)
    fwrite(rbindlist(out), file.path(OUT_DIR, sprintf("cliff_trigger_%s.csv", RUN_TAG))) }
}
