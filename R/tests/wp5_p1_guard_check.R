# ==============================================================================
# R/tests/wp5_p1_guard_check.R   (WP5 P1 guards; Mac; needs INLA only for the posterior sampler)
#   Rscript R/tests/wp5_p1_guard_check.R [fit_dir] [table_dir]
#   defaults: fit_dir = data-output/final/fits_subsample (stage2.rds of the smoke test), table_dir = TABLE_DIR of R/config.R
# Reproduces the hang diagnosed on the Mac (smoke test, SUBSAMPLE_N = 500): the theta DRAWN for m = 1 has non-finite hazards, so a P1 path exits at
# time 0 forever. Same inputs as the MI loop: RNGkind L'Ecuyer-CMRG (as 05_fit_inla.R), prepare_fit_inputs(TABLE_DIR, subsample_n = 500, seed = 2026),
# set.seed(2026 * 1000 + 1), draw_theta_m(stage 2 models, stage 2 theta), then p1_entries(th, inp, 20, "draw", 1, m = 1).
#   PASS 1  p1_entries on the Stage 2 theta runs (mean mode, M = 200) in < 60 s (diagnosed: 0.5 s)
#   PASS 2  p1_entries on the drawn theta stops within 60 s with the guard-A message; printed as "GUARD OK: ..."
# A time-out (hang), a different error or no stop at all is a FAIL and is reported, nothing is changed.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("R", "config.R")); source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
if (requireNamespace("INLA", quietly = TRUE)) suppressPackageStartupMessages(library(INLA))
args <- commandArgs(TRUE)
fit_dir <- if (length(args) >= 1) args[1] else paste0(FIT_DIR, "_subsample"); table_dir <- if (length(args) >= 2) args[2] else TABLE_DIR
sampler <- if (exists("TEST_SAMPLER", inherits = FALSE)) TEST_SAMPLER else default_sampler        # TEST_SAMPLER: cloud test of this script only
SEED <- 2026L
RNGkind("L'Ecuyer-CMRG")
inp <- prepare_fit_inputs(table_dir, subsample_n = 500L, seed = SEED)
s2 <- readRDS(file.path(fit_dir, "stage2.rds"))
cat(sprintf("Plots %d | intervals %d | kept transitions %d | P1_MAX_JUMPS %d\n", nrow(inp$p1), nrow(inp$iv), nrow(inp$td), P1_MAX_JUMPS))

# ---- PASS 1: Stage 2 theta ---------------------------------------------------------------------------------------------------------------------
t0 <- proc.time()["elapsed"]
pe <- p1_entries(s2$theta, inp, 200L, "mean", 1L)
el <- proc.time()["elapsed"] - t0
cat(sprintf("1. p1_entries on the Stage 2 theta (mean, M = 200): %.2f s, %d plots, %d fallbacks\n", el, nrow(pe), sum(pe$fallback)))
if (nrow(pe) != nrow(inp$p1) || el >= 60) stop("FAIL 1: p1_entries on the Stage 2 theta (", round(el, 1), " s)")

# ---- PASS 2: drawn theta, m = 1 --------------------------------------------------------------------------------------------------------------------
set.seed(SEED * 1000L + 1L)
th <- draw_theta_m(s2$models, inp, s2$theta, sampler)
L0 <- lambda_matrix(th, inp$X0); Liv <- lambda_matrix(th, inp$Xiv)
cat(sprintf("2. drawn theta m = 1: non-finite hazards %d in X0 + %d in Xiv (diagnosed on the Mac: 12,291 infinite) | drawn |beta| > 5: %d | largest |beta| %.1f | alpha %.2f-%.2f\n",
            sum(!is.finite(L0)), sum(!is.finite(Liv)), sum(abs(th$beta_v0) > 5), max(abs(th$beta_v0)), min(th$alpha, na.rm = TRUE), max(th$alpha, na.rm = TRUE)))
t0 <- proc.time()["elapsed"]
res <- tryCatch({ setTimeLimit(elapsed = 60, transient = TRUE); p1_entries(th, inp, 20L, "draw", 1L, m = 1L); "NO_STOP" }, error = function(e) conditionMessage(e))
setTimeLimit(elapsed = Inf); el <- proc.time()["elapsed"] - t0
if (grepl("^P1 guard A", res)) cat(sprintf("GUARD OK: (%.2f s) %s\n", el, res)) else
  if (identical(res, "NO_STOP")) stop("FAIL 2: p1_entries on the drawn theta did not stop (no non-finite hazard in this draw? compare the counts printed above with the diagnosed 12,291)") else
  stop("FAIL 2: p1_entries on the drawn theta ended with a different error / time-out after ", round(el, 1), " s: ", res)
cat("\nALL CHECKS PASSED\n")
