# ==============================================================================
# 06p_transition_comparison.R   RUN ON THE MACHINE WITH INLA (~ 10 min)
# What is special about black_spruce->regen?  Using the Stage-1 rows (entry W_j/2, midpoint exits; no draws involved):
#  (1) DATA FEATURES for all kept transitions: rows, events, share of rows with truncation > 0, truncation median / p95 / max,
#      Tstop max, and the "cold-start burden" = sum over rows of (Tstop^1.5 - truncation^1.5) (what INLA faces at intercept 0, alpha 1.5).
#  (2) COLD-START TEST for every transition with >= MIN_EV events: intercept-only INLA fit with alpha FIXED (1.5 and 1.8),
#      pipeline settings. Does INLA's intercept leave 0 (compare with our own optimum)?  loglik_lost ~ 0: fine; thousands: stuck.
#  If other large transitions also get stuck, the INLA problem is general and black_spruce->regen is just where the free fit fails.
#  If only black_spruce->regen (or the transitions with the largest cold-start burden) get stuck, the driver is that data property.
# Output: data/processed/fit/transition_comparison_<tag>.csv  (+ _features.csv)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(INLA) })
SIM_DIR <- "/Users/lukas/Desktop/RProjects/Quebec-forest-succession/R-scripts/Sim4stateExact/p1p2_v2_test"
RUN_TAG <- "diag"; OUT_DIR <- here("data", "processed", "fit"); SEED <- 2026L; options(width = 200)
MIN_EV <- 100L; ALPHAS <- c(1.5, 1.8)
source(file.path(SIM_DIR, "model_fitting.R"))
source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
inp <- prepare_fit_inputs(here("data", "processed"), subsample_n = NULL, seed = SEED)
b <- build_rows(inp$iv, setNames(inp$p1$W_j / 2, inp$p1$id_pe), inp$states, inp$tr, keep_cols = inp$keep_cols)
rows <- b$rows[kept == TRUE]
dev_own <- function(b0, a, x) { lam <- exp(b0); y <- x$Tstop; tr <- x$truncation; st <- x$status
  -2 * (st * (log(a) + a * b0 + (a - 1) * log(y)) - (lam * y)^a + (lam * tr)^a) }
opt_b <- function(a, x) optimize(function(b0) sum(dev_own(b0, a, x)), c(-8, 3), tol = 1e-10)$minimum
feat <- rows[, .(rows = .N, events = sum(status == 1L), share_trunc_pos = mean(truncation > 0), trunc_p50 = median(truncation),
                 trunc_p95 = quantile(truncation, .95), trunc_max = max(truncation), Tstop_max = max(Tstop),
                 cold_burden = sum(Tstop^1.5 - truncation^1.5), cold_burden_per_row = mean(Tstop^1.5 - truncation^1.5)), by = .(label = paste(from, to, sep = "->"))]
setorder(feat, -events); fwrite(feat, file.path(OUT_DIR, sprintf("transition_comparison_features_%s.csv", RUN_TAG)))
cat("\n==== (1) data features, top 15 by events ====\n"); print(feat[1:15, lapply(.SD, function(v) if (is.numeric(v)) signif(v, 3) else v)])
fit_fixed <- function(x, a) {
  args <- list(formula = inla.surv(time = Tstop, event = status, truncation = truncation) ~ 1, family = "weibullsurv", data = x,
               control.family = list(variant = 1, hyper = list(alpha = list(initial = 10 * log(a), fixed = TRUE))), control.compute = list(config = TRUE),
               control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)))
  tryCatch(do.call(inla, args), error = function(e) NULL)
}
out <- list(); todo <- feat[events >= MIN_EV, label]
for (lb in todo) {
  x <- as.data.frame(rows[paste(from, to, sep = "->") == lb, .(status, Tstop, truncation)])
  for (a in ALPHAS) {
    r <- fit_fixed(x, a); bo <- opt_b(a, x)
    bi <- if (is.null(r)) NA_real_ else r$summary.fixed["(Intercept)", "mean"]
    out[[length(out) + 1]] <- data.table(label = lb, events = sum(x$status), alpha = a, intercept_inla = bi, intercept_opt = bo,
      loglik_lost = if (is.na(bi)) NA_real_ else (-sum(dev_own(bi, a, x)) / 2) - (-sum(dev_own(bo, a, x)) / 2) |> abs())
    cat(sprintf("%-40s events %4d alpha %.1f  inla %.3f  opt %.3f\n", lb, sum(x$status), a, bi, bo))
  }
  fwrite(rbindlist(out), file.path(OUT_DIR, sprintf("transition_comparison_%s.csv", RUN_TAG)))
}
R <- rbindlist(out); R[, stuck := is.na(intercept_inla) | abs(intercept_inla - intercept_opt) > 0.3]
R <- merge(R, feat[, .(label, rows, share_trunc_pos, trunc_p95, trunc_max, Tstop_max, cold_burden_per_row)], by = "label")
cat("\n==== (2) cold-start test: stuck = INLA intercept > 0.3 from the optimum ====\n")
print(R[order(-events), .(label, events, alpha, intercept_inla = round(intercept_inla, 3), intercept_opt = round(intercept_opt, 3), stuck, share_trunc_pos = round(share_trunc_pos, 2),
                          trunc_p95 = round(trunc_p95, 1), Tstop_max = round(Tstop_max, 1), cold_burden_per_row = round(cold_burden_per_row, 1))])
cat("\nstuck by alpha:\n"); print(R[, .(n = .N, n_stuck = sum(stuck)), by = alpha])
