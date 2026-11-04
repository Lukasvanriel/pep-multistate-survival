# ==============================================================================
# 17_htr_sim_grid.R   RUN ON THE MACHINE WITH INLA (~3-5 min with 6 workers)   -- DIAGNOSTIC ONLY, uses no project data
# Run from a terminal, not RStudio (parallel workers inside RStudio can hang):  Rscript R/tests/17_htr_sim_grid.R
#
# QUESTION: rough failure curve: how often does ONE extreme row make INLA fail, as a function of its cumulative hazard at truncation H(tr)?
# DESIGN (small, for a reproducible example): simple Weibull data, N_SIM rows, one row per subject, left truncation, x1 ~ N(0,1) capped at
#   +-3, x2 ~ Bernoulli(0.2), true alpha 1.8, time in years, ~30% events; base rows have H(tr) below ~10. Plus ONE extra event row with
#   x1 = x2 = 0 and H(tr) = H (Tstop = tr + 2 years). H = 0 (no extra row: control) and 15, 16, ..., 30; N_REP replicate datasets.
#   Within a replicate the base data are identical across H, so only the extra row changes.
#   INLA call identical to the pipeline (weibullsurv variant 1, adaptive, grid, VB off), one thread per call.
# FAILED: INLA error/NA, alpha-SD ratio < 0.1, or |alpha_INLA - alpha_MLE| > 3 MLE-SE (MLE = exact optimum of the same likelihood, reference only).
# SAVES AS IT GOES: every fit writes its own file to data-output/diagnostics/htr_grid/ when it finishes, and a line to progress.log
#   (start and end of each fit). Re-running skips fits already saved, so an interrupted run loses nothing.
# OUTPUT: data-output/diagnostics/htr_sim_grid.csv (all fits) + printed failures per H.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(INLA) }); options(width = 200)
N_WORKERS <- 6L; INLA::inla.setOption(num.threads = "1:1")
H_GRID <- c(0, 15:30); A_TRUE <- 1.8; UNIT <- 10; N_REP <- 5L; N_SIM <- 3000L
OUT <- here("data-output", "diagnostics", "htr_grid"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(OUT, "progress.log")
cat("INLA version:", as.character(packageVersion("INLA")), "\n")

sim <- function(alpha, unit, rep, H) {                      # unit = 10: time in years
  set.seed(rep * 100L + round(alpha * 10)); lam <- 0.25 / unit; n <- N_SIM
  x1 <- pmin(pmax(rnorm(n), -3), 3); x2 <- rbinom(n, 1, 0.2); eta <- log(lam) + 0.3 * x1 + 0.5 * x2; e <- exp(alpha * eta)
  tr <- runif(n, 0, 3.5 * unit)                            # base rows stay below the grid (H(tr) <= ~10)
  Tend <- tr + runif(n, 0.5, 1.5) * unit
  tev <- ((e * tr^alpha + rexp(n)) / e)^(1 / alpha); st <- as.integer(tev <= Tend)
  d <- data.frame(status = st, Tstop = pmin(tev, Tend), truncation = tr, x1 = x1, x2 = x2)
  base_maxH <- max(e * tr^alpha)
  if (H > 0) { trx <- (H / lam^alpha)^(1 / alpha)               # x1 = x2 = 0 -> exp(alpha*eta) = lam^alpha
    d <- rbind(d, data.frame(status = 1L, Tstop = trx + 0.2 * unit, truncation = trx, x1 = 0, x2 = 0)) }
  list(d = d, base_maxH = base_maxH) }
fit_mle <- function(d) {
  X <- cbind(1, d$x1, d$x2); T <- d$Tstop; tr <- d$truncation; s <- d$status; lT <- log(T); ltr <- ifelse(tr > 0, log(pmax(tr, 1e-300)), 0)
  nll <- function(th) { a <- exp(th[1]); eta <- drop(X %*% th[-1]); -(sum(s * (th[1] + a * eta + (a - 1) * lT)) - sum(exp(a * eta) * (T^a - tr^a))) }
  gr <- function(th) { a <- exp(th[1]); eta <- drop(X %*% th[-1]); e <- exp(a * eta); Ta <- T^a; tra <- tr^a
    gb <- -(colSums(X * (s * a)) - colSums(X * (a * e * (Ta - tra)))); ga <- -(sum(s * (1 / a + eta + lT)) - sum(e * (eta * (Ta - tra) + Ta * lT - tra * ltr))); c(ga * a, gb) }
  o <- optim(c(0, log(sum(s) / sum(T - tr)), 0, 0), nll, gr, method = "BFGS", control = list(maxit = 5000, reltol = 1e-12))
  V <- tryCatch(solve(optimHess(o$par, nll, gr)), error = function(e) NULL); a <- exp(o$par[1])
  c(mle_alpha = a, mle_alpha_se = if (is.null(V)) NA_real_ else a * sqrt(V[1, 1])) }
fit_inla <- function(d) tryCatch(inla(inla.surv(time = Tstop, event = status, truncation = truncation) ~ x1 + x2, family = "weibullsurv", data = d,
  control.family = list(variant = 1), control.compute = list(config = TRUE),
  control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE))), error = function(e) list(error = conditionMessage(e)))

grid <- CJ(rep = seq_len(N_REP), H = H_GRID)
grid[, file := file.path(OUT, sprintf("fit_H%02d_rep%02d.csv", H, rep))]
todo <- grid[!file.exists(file)]
cat(sprintf("Fits: %d in total, %d already saved, %d to run\n", nrow(grid), nrow(grid) - nrow(todo), nrow(todo)))
t0 <- proc.time()["elapsed"]
invisible(mclapply(seq_len(nrow(todo)), function(i) { g <- todo[i]; ts <- proc.time()["elapsed"]
  cat(sprintf("%s start H=%d rep=%d\n", format(Sys.time(), "%H:%M:%S"), g$H, g$rep), file = LOG, append = TRUE)
  r <- tryCatch({ s <- sim(A_TRUE, UNIT, g$rep, g$H); mo <- fit_inla(s$d); ml <- fit_mle(s$d)
    a <- if (is.null(mo$error)) { h <- mo$summary.hyperpar; h[grep("alpha", rownames(h))[1], ] } else NULL
    data.table(H = g$H, rep = g$rep, base_maxH = s$base_maxH, events = sum(s$d$status), inla_error = !is.null(mo$error),
               alpha = if (is.null(a)) NA_real_ else a$mean, alpha_sd = if (is.null(a)) NA_real_ else a$sd,
               mle_alpha = ml[["mle_alpha"]], mle_alpha_se = ml[["mle_alpha_se"]]) },
    error = function(e) data.table(H = g$H, rep = g$rep, inla_error = TRUE))
  r[, seconds := round(proc.time()["elapsed"] - ts)]
  fwrite(r, g$file)
  cat(sprintf("%s done  H=%d rep=%d (%d s)\n", format(Sys.time(), "%H:%M:%S"), g$H, g$rep, r$seconds), file = LOG, append = TRUE)
  NULL }, mc.cores = N_WORKERS, mc.preschedule = FALSE))
cat(sprintf("Run time: %.1f min\n", (proc.time()["elapsed"] - t0) / 60))

res <- rbindlist(lapply(grid[file.exists(file), file], fread), fill = TRUE)
res[, `:=`(ratio = alpha_sd / (0.78 * alpha / sqrt(events)), z_mle = (alpha - mle_alpha) / mle_alpha_se)]
res[, failed := inla_error | is.na(alpha) | ratio < 0.1 | abs(z_mle) > 3]
fwrite(res, here("data-output", "diagnostics", "htr_sim_grid.csv"))
cat(sprintf("Largest H(tr) among the base rows: %.1f (should be below 15)\n", max(res$base_maxH, na.rm = TRUE)))
cat("\n== Failed fits per H(tr) of the extra event row (H = 0: no extra row) ==\n")
print(res[, .(failed = sprintf("%d/%d", sum(failed), .N), slowest_fit_s = max(seconds)), by = H][order(H)])
