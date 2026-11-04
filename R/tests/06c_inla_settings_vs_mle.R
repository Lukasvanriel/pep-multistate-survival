# ==============================================================================
# 06c_inla_settings_vs_mle.R   RUN ON THE MACHINE WITH INLA
#
# Do other INLA control settings reproduce the maximum-likelihood alpha/SD on the
# datasets where the pipeline's INLA fit looked wrong?
# Datasets: all saved flagged fits (baddata_diag_mXX.rds) + the unflagged-but-off
#           black_spruce->regen fit of imputation 1 (from models_diag_m01.rds).
# Reference: independent Weibull MLE (same likelihood: split rows, left truncation,
#           variant 1 S = exp(-(lambda t)^alpha)); SD from the inverse Hessian.
# Settings (everything else as in fit_inla_models()): see SETTINGS below.
# Output: data/processed/fit/settings_vs_mle_<tag>.csv  (+ printed summary)
# No pipeline file is modified.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(INLA) })
RUN_TAG <- "diag"; OUT_DIR <- here("data", "processed", "fit"); ALPHA_NAME <- "alpha parameter for weibullsurv"
options(width = 200)

mle <- function(x, rhs) {
  tms <- strsplit(rhs, " \\+ ")[[1]]; X <- cbind(1, as.matrix(x[, tms])); y <- x$Tstop; tr <- x$truncation; st <- x$status
  nll <- function(p) { a <- exp(p[1]); eta <- drop(X %*% p[-1]); lam <- exp(eta)
    -sum(st * (log(a) + a * eta + (a - 1) * log(y)) - (lam * y)^a + (lam * tr)^a) }
  best <- NULL
  for (s in c(0.5, 1.5)) { st0 <- c(log(s + 0.5), -2, rep(0, ncol(X) - 1))
    o <- tryCatch(optim(st0, nll, method = "BFGS", hessian = TRUE, control = list(maxit = 2000, reltol = 1e-12)), error = function(e) NULL)
    if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o }
  V <- tryCatch(solve(best$hessian), error = function(e) matrix(NA, nrow(best$hessian), ncol(best$hessian)))
  a <- exp(best$par[1]); list(alpha = a, sd = a * sqrt(V[1, 1]), ll = -best$value)
}

# ---- datasets ---------------------------------------------------------------
flagged <- fread(file.path(OUT_DIR, sprintf("flagged_fits_%s.csv", RUN_TAG)))
sets <- list()
for (i in seq_len(nrow(flagged))) {
  b <- readRDS(file.path(OUT_DIR, sprintf("baddata_%s_m%02d.rds", RUN_TAG, flagged$m[i])))[[flagged$label[i]]]
  sets[[sprintf("m%02d %s", flagged$m[i], flagged$label[i])]] <- list(m = flagged$m[i], label = flagged$label[i], data = b$data, rhs = b$rhs, kind = "flagged")
}
meta <- readRDS(file.path(OUT_DIR, sprintf("run_meta_%s.rds", RUN_TAG))); lab <- meta$td$label
M1 <- readRDS(file.path(OUT_DIR, sprintf("models_%s_m01.rds", RUN_TAG))); k <- match("black_spruce->regen", lab)
sets[["m01 black_spruce->regen"]] <- list(m = 1L, label = "black_spruce->regen", data = as.data.frame(M1[[k]]$.args$data), rhs = meta$rhs[[k]], kind = "unflagged_off")
rm(M1); invisible(gc())

# ---- settings ---------------------------------------------------------------
SETTINGS <- list(
  S0_default      = list(ci = list()),
  S1_tol1e6       = list(ci = list(tolerance = 1e-6)),
  S2_tol_h        = list(ci = list(tolerance = 1e-6, h = 0.001)),
  S3_ccd          = list(ci = list(int.strategy = "ccd")),
  S4_eb           = list(ci = list(int.strategy = "eb")),
  S5_laplace      = list(ci = list(strategy = "laplace")),
  S6_tight_ccd    = list(ci = list(tolerance = 1e-6, h = 0.001, int.strategy = "ccd")),
  S7_default_1thr = list(ci = list(), threads = "1:1"),
  S8_tight_1thr   = list(ci = list(tolerance = 1e-6, h = 0.001), threads = "1:1"))

fit1 <- function(d, rhs, s) {
  ci <- modifyList(list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)), s$ci)
  a <- list(formula = as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", rhs)),
            family = "weibullsurv", data = d, control.family = list(variant = 1), control.compute = list(config = TRUE), control.inla = ci)
  if (!is.null(s$threads)) a$num.threads <- s$threads
  t0 <- Sys.time(); r <- tryCatch(do.call(inla, a), error = function(e) list(error = conditionMessage(e)))
  sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (!is.null(r$error)) return(data.table(alpha = NA_real_, sd = NA_real_, mlik = NA_real_, nfunc = NA_real_, sec = sec, error = substr(r$error, 1, 120)))
  data.table(alpha = r$summary.hyperpar[ALPHA_NAME, "mean"], sd = r$summary.hyperpar[ALPHA_NAME, "sd"], mlik = as.numeric(r$mlik[1, 1]),
             nfunc = as.numeric(r$misc$nfunc), sec = sec, error = NA_character_)
}

out <- list()
for (nm in names(sets)) {
  s <- sets[[nm]]; ref <- mle(s$data, s$rhs)
  cat(sprintf("\n== %s (%s)  rows %d events %d | MLE alpha %.3f sd %.4f ==\n", nm, s$kind, nrow(s$data), sum(s$data$status), ref$alpha, ref$sd))
  for (sn in names(SETTINGS)) {
    r <- fit1(s$data, s$rhs, SETTINGS[[sn]])
    out[[length(out) + 1]] <- cbind(data.table(m = s$m, label = s$label, kind = s$kind, setting = sn, mle_alpha = ref$alpha, mle_sd = ref$sd), r)
    cat(sprintf("  %-16s alpha %s sd %s  z=%s  nfunc %s  %.0fs %s\n", sn, format(round(r$alpha, 3)), format(signif(r$sd, 3)),
                format(round((r$alpha - ref$alpha) / ref$sd, 1)), format(r$nfunc), r$sec, ifelse(is.na(r$error), "", paste("ERR:", r$error))))
  }
  fwrite(rbindlist(out), file.path(OUT_DIR, sprintf("settings_vs_mle_%s.csv", RUN_TAG)))   # written after each dataset
}
res <- rbindlist(out); res[, z := (alpha - mle_alpha) / mle_sd]; res[, sd_ratio := sd / mle_sd]
cat("\n==== SUMMARY per setting: datasets where INLA agrees with MLE (|z|<1 & SD ratio in 0.7-1.4) ====\n")
print(res[, .(n = .N, n_error = sum(!is.na(error)), n_agree = sum(abs(z) < 1 & sd_ratio > 0.7 & sd_ratio < 1.4, na.rm = TRUE), median_sec = round(median(sec), 1)), by = setting])
fwrite(res, file.path(OUT_DIR, sprintf("settings_vs_mle_%s.csv", RUN_TAG)))
