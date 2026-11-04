# ==============================================================================
# 06t_smoothness_cmin.R   RUN ON THE MACHINE WITH INLA (~ 10 min)
# 06s: cmin = 0 (clip the negative likelihood Hessian; the truncated Weibull loglik is concave in the linear predictor, so this only removes
# numerical artefacts) and stencil = 9 cure the stuck cold start, but the FREE fit still collapses.
# Q1  Is INLA's log marginal likelihood SMOOTH in alpha once the cold start is cured?  Fixed-alpha fits on a fine grid (1.690-1.710, step 0.001)
#     for m02, variants base | cmin0 | stencil9 | cmin0_stencil9; noise = SD(second difference)/sqrt(6)  (smooth function: ~0).
# Q2  Does the hyperparameter step h (default 0.01) matter?  Free fits of the 8 collapsed datasets with cmin = 0 and h in {0.001, 0.01, 0.05, 0.1}.
# Output: data/processed/fit/smoothness_cmin_<tag>.csv (grid), smoothness_cmin_h_<tag>.csv (free fits)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(INLA) })
RUN_TAG <- "diag"; OUT_DIR <- here("data", "processed", "fit"); options(width = 200)
LABEL <- "black_spruce->regen"; ALPHA_NAME <- "alpha parameter for weibullsurv"; M_SET <- c(2, 8, 14, 15, 19, 23, 25, 29); GRID <- seq(1.690, 1.710, by = 0.001)
run <- function(d, rhs, extra, fixed_alpha = NULL) {
  ci <- modifyList(list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)), extra)
  fam <- list(variant = 1); if (!is.null(fixed_alpha)) fam$hyper <- list(alpha = list(initial = 10 * log(fixed_alpha), fixed = TRUE))
  tryCatch(inla(as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", rhs)), family = "weibullsurv", data = d,
                control.family = fam, control.compute = list(config = TRUE), control.inla = ci), error = function(e) conditionMessage(e)) }
# Q1
V1 <- list(base = list(), cmin0 = list(cmin = 0), stencil9 = list(stencil = 9L), cmin0_stencil9 = list(cmin = 0, stencil = 9L))
b <- readRDS(file.path(OUT_DIR, sprintf("baddata_%s_m02.rds", RUN_TAG)))[[LABEL]]; g1 <- list()
for (v in names(V1)) { for (a in GRID) { r <- run(b$data, b$rhs, V1[[v]], a)
    g1[[length(g1) + 1]] <- data.table(variant = v, alpha = a, mlik = if (is.character(r)) NA_real_ else as.numeric(r$mlik[1, 1])) }
  cat("grid", v, "done\n"); fwrite(rbindlist(g1), file.path(OUT_DIR, sprintf("smoothness_cmin_%s.csv", RUN_TAG))) }
G <- rbindlist(g1); setorder(G, variant, alpha); G[, d2 := shift(mlik, -1) - 2 * mlik + shift(mlik, 1), by = variant]
cat("\n==== Q1: fixed-alpha mlik on the fine grid, m02 (noise = SD(d2)/sqrt(6); smooth ~ 1e-3) ====\n")
print(G[, .(n_ok = sum(!is.na(mlik)), mlik_min = round(min(mlik, na.rm = TRUE), 2), mlik_max = round(max(mlik, na.rm = TRUE), 2), noise = signif(sd(d2, na.rm = TRUE) / sqrt(6), 3)), by = variant])
# Q2
H <- c(0.001, 0.01, 0.05, 0.1); out <- list()
for (mm in M_SET) { bb <- readRDS(file.path(OUT_DIR, sprintf("baddata_%s_m%02d.rds", RUN_TAG, mm)))[[LABEL]]; ev <- sum(bb$data$status == 1L)
  for (h in H) { r <- run(bb$data, bb$rhs, list(cmin = 0, h = h))
    out[[length(out) + 1]] <- if (is.character(r)) data.table(m = mm, h = h, alpha = NA_real_, sd = NA_real_, ratio = NA_real_, mlik = NA_real_, error = substr(r, 1, 50)) else {
      a <- r$summary.hyperpar[ALPHA_NAME, "mean"]; s <- r$summary.hyperpar[ALPHA_NAME, "sd"]
      data.table(m = mm, h = h, alpha = a, sd = s, ratio = s / (0.78 * a / sqrt(ev)), mlik = as.numeric(r$mlik[1, 1]), error = NA_character_) }
    print(tail(rbindlist(out), 1)[, .(m, h, alpha = round(alpha, 3), sd = signif(sd, 3), ratio = round(ratio, 3), mlik = round(mlik, 1), error)]) }
  fwrite(rbindlist(out), file.path(OUT_DIR, sprintf("smoothness_cmin_h_%s.csv", RUN_TAG))) }
R <- rbindlist(out)
cat("\n==== Q2: cmin = 0, SD ratio by h (healthy ~1.2) ====\n"); print(dcast(R, m ~ h, value.var = "ratio")[, lapply(.SD, function(v) round(v, 3))])
cat("\nalpha by h (MLE of these data: 1.81-1.96):\n"); print(dcast(R, m ~ h, value.var = "alpha")[, lapply(.SD, function(v) round(v, 3))])
