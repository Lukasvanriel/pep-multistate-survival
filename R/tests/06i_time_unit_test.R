# ==============================================================================
# 06i_time_unit_test.R   RUN ON THE MACHINE WITH INLA (several minutes)
# Finding from 06h: the fixed-alpha collapse appears in ALL random subsets, its size grows with the number of rows, and free fits
# collapse more often in larger subsets (fifths 1/10, halves 7/10, full data most). Free fits that work are warm-started; fixed-alpha
# fits start from default latent values (intercept 0 = hazard 1 per time unit) while the true intercept is ~ -3.1 (hazard 0.04).
# With intercept 0 the truncation term (lambda*truncation)^alpha is huge for long truncation times: (14.7)^1.6 ~ 75, (14.7)^2.4 ~ 650.
# HYPOTHESIS: INLA's inner optimisation diverges from that start. TEST: change the TIME UNIT so the fitted hazard is ~1 (intercept ~ 0):
#   time' = time / s, s = median event time of the data set (data-driven; alpha is unchanged; intercept' = intercept + log(s)).
# For every flagged data set (+ black_spruce->regen of imputation 1) and s in {1 (as now), median event time}:
#   free INLA fit alpha, SD vs MLE alpha, SD (the MLE of the original data; alpha is invariant to the time unit) and the fixed-alpha cliff metric.
# Output: data/processed/fit/time_unit_test_<tag>.csv. Nothing else is modified.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(INLA) })
RUN_TAG <- "diag"; OUT_DIR <- here("data", "processed", "fit"); options(width = 200)
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


flagged <- fread(file.path(OUT_DIR, sprintf("flagged_fits_%s.csv", RUN_TAG))); sets <- list()
for (i in seq_len(nrow(flagged))) { b <- readRDS(file.path(OUT_DIR, sprintf("baddata_%s_m%02d.rds", RUN_TAG, flagged$m[i])))[[flagged$label[i]]]
  sets[[sprintf("m%02d %s", flagged$m[i], flagged$label[i])]] <- list(data = b$data, rhs = b$rhs) }
meta <- readRDS(file.path(OUT_DIR, sprintf("run_meta_%s.rds", RUN_TAG))); lab <- meta$td$label
M1 <- readRDS(file.path(OUT_DIR, sprintf("models_%s_m01.rds", RUN_TAG))); k <- match("black_spruce->regen", lab)
sets[["m01 black_spruce->regen"]] <- list(data = as.data.frame(M1[[k]]$.args$data), rhs = meta$rhs[[k]]); rm(M1); invisible(gc())

out <- list()
for (nm in names(sets)) { st <- sets[[nm]]; tms <- strsplit(st$rhs, " \\+ ")[[1]]; mo <- mle_full(st$data, tms); smed <- median(st$data$Tstop[st$data$status == 1])
  cat(sprintf("\n== %s | MLE alpha %.3f sd %.3f | median event time (decades) %.2f ==\n", nm, mo["alpha"], mo["sd"], smed))
  pt_ev <- sum(st$data$Tstop - st$data$truncation) / sum(st$data$status)   # person-time per event: baseline hazard ~ 1 in this unit
  for (sc in c(1, smed, pt_ev, 2 * pt_ev)) { v <- st$data; v$Tstop <- v$Tstop / sc; v$truncation <- v$truncation / sc; pf <- prof_fun(v, tms)
    ml <- sapply(ALPHAS, function(a) mlik_fixed(v, tms, a)); pr <- sapply(ALPHAS, pf); metric <- (ml - ml[ALPHAS == REF]) - (pr - pr[ALPHAS == REF]); cliff <- ALPHAS[which(metric < -50)[1]]
    fi <- free_fit(v, tms)
    cat(sprintf("  unit %-6.2f free INLA alpha %.3f sd %.4g (z %.1f, sd-ratio %.2f) | cliff metric [%s] CLIFF_AT %s\n", sc, fi["alpha"], fi["sd"], (fi["alpha"] - mo["alpha"]) / mo["sd"], fi["sd"] / mo["sd"],
                paste(round(metric), collapse = " "), ifelse(is.na(cliff), "none", cliff)))
    out[[length(out) + 1]] <- data.table(set = nm, time_unit = sc, inla_alpha = fi["alpha"], inla_sd = fi["sd"], mle_alpha = mo["alpha"], mle_sd = mo["sd"], cliff_at = cliff)
    fwrite(rbindlist(out), file.path(OUT_DIR, sprintf("time_unit_test_%s.csv", RUN_TAG))) } }
r <- rbindlist(out); r[, ok := abs(inla_alpha - mle_alpha) / mle_sd < 2 & inla_sd / mle_sd > 0.5 & inla_sd / mle_sd < 2]
cat("\n==== SUMMARY: fits agreeing with the MLE (|z|<2, SD ratio 0.5-2) and collapse found, by time unit ====\n")
r[, unit := c("1 (as now)", "median event time", "person-time/event", "2 x person-time/event")[ave(seq_len(.N), set, FUN = seq_along)]]
print(r[, .(n = .N, n_ok = sum(ok, na.rm = TRUE), n_cliff = sum(!is.na(cliff_at))), by = unit])
