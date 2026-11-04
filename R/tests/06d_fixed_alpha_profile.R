# ==============================================================================
# 06d_fixed_alpha_profile.R   RUN ON THE MACHINE WITH INLA
#
# Is INLA's collapsed alpha a failure to find the maximum of ITS OWN likelihood,
# or does INLA's likelihood differ from the Weibull likelihood used in the MLE?
# For each dataset (as in 06c) alpha is FIXED on a grid and:
#   INLA: log marginal likelihood mlik(alpha) (alpha fixed via control.family hyper)
#   MLE : profile log-likelihood prof(alpha) (beta maximised at each alpha)
# Reading: curves should have the same shape (differences between alphas agree within
# ~1-2 log units, since mlik additionally integrates beta under its vague prior).
#   * INLA curve peaks near the MLE alpha  -> INLA's likelihood is fine; the free fit's
#     optimiser/integration lands on the wrong alpha (optimiser/software problem).
#   * INLA curve peaks near the collapsed alpha -> INLA evaluates a different likelihood
#     (specification / data handling / prior problem).
# Sanity check: fixing alpha at the free fit's own mode must reproduce the free fit's mlik.
# Output: data/processed/fit/fixed_alpha_profile_<tag>.csv + printed tables. Nothing else modified.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(INLA) })
RUN_TAG <- "diag"; OUT_DIR <- here("data", "processed", "fit"); ALPHA_NAME <- "alpha parameter for weibullsurv"
options(width = 200)

design <- function(x, rhs) { tms <- strsplit(rhs, " \\+ ")[[1]]; cbind(1, as.matrix(x[, tms])) }
mle_prof <- function(x, rhs) {
  X <- design(x, rhs); y <- x$Tstop; tr <- x$truncation; st <- x$status
  nll <- function(b, a) { eta <- drop(X %*% b); lam <- exp(eta); -sum(st * (log(a) + a * eta + (a - 1) * log(y)) - (lam * y)^a + (lam * tr)^a) }
  prof <- function(a) { best <- Inf
    for (i0 in c(-2, -4)) { o <- tryCatch(optim(c(i0, rep(0, ncol(X) - 1)), nll, a = a, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12)), error = function(e) NULL)
      if (!is.null(o) && o$value < best) best <- o$value }
    -best }
  o <- optimize(function(a) -prof(a), c(0.3, 4), tol = 1e-4); list(prof = prof, alpha = o$minimum, ll = -o$objective)
}

# ---- datasets (same as 06c) --------------------------------------------------
flagged <- fread(file.path(OUT_DIR, sprintf("flagged_fits_%s.csv", RUN_TAG))); sets <- list()
for (i in seq_len(nrow(flagged))) {
  b <- readRDS(file.path(OUT_DIR, sprintf("baddata_%s_m%02d.rds", RUN_TAG, flagged$m[i])))[[flagged$label[i]]]
  sets[[sprintf("m%02d %s", flagged$m[i], flagged$label[i])]] <- list(m = flagged$m[i], label = flagged$label[i], data = b$data, rhs = b$rhs, kind = "flagged") }
meta <- readRDS(file.path(OUT_DIR, sprintf("run_meta_%s.rds", RUN_TAG))); lab <- meta$td$label
M1 <- readRDS(file.path(OUT_DIR, sprintf("models_%s_m01.rds", RUN_TAG))); k <- match("black_spruce->regen", lab)
sets[["m01 black_spruce->regen"]] <- list(m = 1L, label = "black_spruce->regen", data = as.data.frame(M1[[k]]$.args$data), rhs = meta$rhs[[k]], kind = "unflagged_off")
rm(M1); invisible(gc())

inla_args <- function(d, rhs, a = NULL) {
  cf <- list(variant = 1); if (!is.null(a)) cf$hyper <- list(alpha = list(initial = 10 * log(a), fixed = TRUE))   # NB key "alpha" (or "theta1"); key "theta" is silently ignored; theta = 10*log(alpha)
  list(formula = as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", rhs)), family = "weibullsurv", data = d,
       control.family = cf, control.compute = list(config = TRUE),
       control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)))
}
run <- function(args) { r <- tryCatch(do.call(inla, args), error = function(e) NULL)
  if (is.null(r)) { args$num.threads <- "1:1"; r <- tryCatch(do.call(inla, args), error = function(e) NULL) }; r }
mlik_fixed <- function(d, rhs, a) { r <- run(inla_args(d, rhs, a)); if (is.null(r)) return(NA_real_)
  if (!isTRUE(r$all.hyper$family[[1]]$hyper[[1]]$fixed) || NROW(r$summary.hyperpar) > 0) stop("alpha fixing was NOT honoured")
  as.numeric(r$mlik[1, 1]) }

out <- list(); summ <- list()
for (nm in names(sets)) {
  s <- sets[[nm]]; cat(sprintf("\n== %s (%s) ==\n", nm, s$kind))
  mp <- mle_prof(s$data, s$rhs)
  free <- run(inla_args(s$data, s$rhs))
  a_free <- if (is.null(free)) NA_real_ else free$summary.hyperpar[ALPHA_NAME, "mean"]; mlik_free <- if (is.null(free)) NA_real_ else as.numeric(free$mlik[1, 1])
  th_free <- if (is.null(free)) NA_real_ else 0.1 * as.numeric(free$mode$theta[1])   # = log(alpha) at the mode (alpha = exp(0.1*theta))
  chk <- if (is.na(th_free)) NA_real_ else mlik_fixed(s$data, s$rhs, exp(th_free))
  cat(sprintf("  MLE alpha %.3f | free INLA alpha %.3f (mode %.3f), mlik %.2f | fixed at INLA mode: mlik %.2f (sanity: should be ~ free mlik)\n", mp$alpha, a_free, exp(th_free), mlik_free, chk))
  grid <- sort(unique(round(c(seq(0.8, 2.6, by = 0.1), mp$alpha, if (!is.na(th_free)) exp(th_free)), 4)))
  tab <- rbindlist(lapply(grid, function(a) data.table(alpha = a, inla_mlik = mlik_fixed(s$data, s$rhs, a), mle_prof = mp$prof(a))))
  tab[, `:=`(inla_rel = inla_mlik - max(inla_mlik, na.rm = TRUE), mle_rel = mle_prof - max(mle_prof))]
  print(tab[, .(alpha = round(alpha, 3), inla_rel = round(inla_rel, 2), mle_rel = round(mle_rel, 2), diff = round(inla_rel - mle_rel, 2))])
  ia <- tab$alpha[which.max(tab$inla_mlik)]
  cat(sprintf("  INLA fixed-alpha curve peaks at grid alpha %.3f; MLE profile peaks at %.3f; free INLA alpha %.3f\n", ia, mp$alpha, a_free))
  out[[nm]] <- cbind(data.table(m = s$m, label = s$label, kind = s$kind), tab)
  summ[[nm]] <- data.table(m = s$m, label = s$label, kind = s$kind, alpha_mle = mp$alpha, alpha_inla_free = a_free, alpha_inla_fixedcurve_peak = ia,
                           mlik_free = mlik_free, mlik_fixed_at_mode = chk)
  fwrite(rbindlist(out), file.path(OUT_DIR, sprintf("fixed_alpha_profile_%s.csv", RUN_TAG)))
}
cat("\n==== SUMMARY ====\n"); print(rbindlist(summ)[, .(m, label, alpha_mle = round(alpha_mle, 3), alpha_inla_free = round(alpha_inla_free, 3), alpha_fixedcurve_peak = alpha_inla_fixedcurve_peak,
                                                       mlik_free = round(mlik_free, 2), mlik_fixed_at_mode = round(mlik_fixed_at_mode, 2))])
fwrite(rbindlist(summ), file.path(OUT_DIR, sprintf("fixed_alpha_summary_%s.csv", RUN_TAG)))
