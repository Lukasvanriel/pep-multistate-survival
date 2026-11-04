# ==============================================================================
# 06y_warm_start.R   RUN ON THE MACHINE WITH INLA (~ 10 min), default (compact) mode
# 06u trace: in the free m08 fit the SAME theta (the mode) gives -2166.3 during optimisation but -23417 during integration; 06x: classic mode does not help.
# HYPOTHESIS: every evaluation of the latent-field solve starts cold (x = 0 = hazard 1/decade, truncation term ~1e4-1e5) and sometimes ends in a wrong point.
# TEST (INLA-native, documented: control.mode / inla.rerun; no data, model or prior change):
#  A  m02 fixed-alpha grid 1.690-1.710 (cold start fails everywhere in default settings): start every fit from the converged latent field of ONE good fit
#     (alpha 1.699, cmin = 0 + stencil 9)  -> control.mode = list(result = good, restart = TRUE).  Smooth, ~ -2165 = cured.
#  B  free fits of the 8 collapsed black_spruce->regen data sets: fit 1 (as pipeline) then inla.rerun() 1x and 2x  (mode of the previous fit as start).
# Output: data/processed/fit/warm_start_grid_<tag>.csv, warm_start_free_<tag>.csv
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(INLA) })
RUN_TAG <- "diag"; OUT_DIR <- here("data", "processed", "fit"); options(width = 200)
LABEL <- "black_spruce->regen"; ALPHA_NAME <- "alpha parameter for weibullsurv"; M_SET <- c(2, 8, 14, 15, 19, 23, 25, 29)
cat("inla.mode before:", tryCatch(inla.getOption("inla.mode"), error = function(e) "n/a"), "\n")
inla.setOption(inla.mode = "compact")      # 06x set "classic" for the whole R session; this script must run in the default (compact) mode
cat("inla.mode now   :", inla.getOption("inla.mode"), "\n"); stopifnot(identical(inla.getOption("inla.mode"), "compact"))
FORM <- function(rhs) as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", rhs))
fit <- function(d, rhs, extra = list(), fixed_alpha = NULL, cm = NULL) {
  ci <- modifyList(list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)), extra)
  fam <- list(variant = 1); if (!is.null(fixed_alpha)) fam$hyper <- list(alpha = list(initial = 10 * log(fixed_alpha), fixed = TRUE))
  a <- list(formula = FORM(rhs), family = "weibullsurv", data = d, control.family = fam, control.inla = ci); if (!is.null(cm)) a$control.mode <- cm
  tryCatch(do.call(inla, a), error = function(e) conditionMessage(e)) }
mk <- function(r) if (is.character(r)) NA_real_ else as.numeric(r$mlik[1, 1])
# Part A
b <- readRDS(file.path(OUT_DIR, sprintf("baddata_%s_m02.rds", RUN_TAG)))[[LABEL]]
good <- fit(b$data, b$rhs, list(cmin = 0, stencil = 9L), 1.699); cat("good start mlik:", mk(good), "(expect -2165.16)\n"); stopifnot(!is.character(good), abs(mk(good) + 2165.16) < 1)   # must reproduce the 06t value, else stop
g <- list()
for (a in seq(1.690, 1.710, by = 0.001)) {
  g[[length(g) + 1]] <- data.table(alpha = a, cold = mk(fit(b$data, b$rhs, list(), a)), warm = mk(fit(b$data, b$rhs, list(), a, list(result = good, restart = TRUE))),
                                   warm_cmin0 = mk(fit(b$data, b$rhs, list(cmin = 0), a, list(result = good, restart = TRUE))))
  fwrite(rbindlist(g), file.path(OUT_DIR, sprintf("warm_start_grid_%s.csv", RUN_TAG))) }
G <- rbindlist(g); cat("\n==== A: m02 fixed-alpha grid (cured = about -2165 and smooth; stuck = about -25,000) ====\n"); print(G[, lapply(.SD, round, 2), .SDcols = names(G)])
cat("stuck (< -10000):", sapply(G[, -1], function(v) sum(v < -10000, na.rm = TRUE)), "of", nrow(G), "\n")
# Part B
out <- list()
for (mm in M_SET) { bb <- readRDS(file.path(OUT_DIR, sprintf("baddata_%s_m%02d.rds", RUN_TAG, mm)))[[LABEL]]; ev <- sum(bb$data$status == 1L)
  r <- fit(bb$data, bb$rhs); stage <- "fit1"
  for (pass in 1:3) {
    if (!is.character(r)) { a <- r$summary.hyperpar[ALPHA_NAME, "mean"]; s <- r$summary.hyperpar[ALPHA_NAME, "sd"]
      out[[length(out) + 1]] <- data.table(m = mm, pass = stage, alpha = a, sd = s, ratio = s / (0.78 * a / sqrt(ev)), mlik = mk(r), error = NA_character_) }
    else out[[length(out) + 1]] <- data.table(m = mm, pass = stage, alpha = NA_real_, sd = NA_real_, ratio = NA_real_, mlik = NA_real_, error = substr(r, 1, 60))
    print(tail(rbindlist(out), 1)[, .(m, pass, alpha = round(alpha, 3), sd = signif(sd, 3), ratio = round(ratio, 3), mlik = round(mlik, 1), error)])
    if (pass == 3 || is.character(r)) break
    r <- tryCatch(inla.rerun(r), error = function(e) conditionMessage(e)); stage <- paste0("rerun", pass) }
  fwrite(rbindlist(out), file.path(OUT_DIR, sprintf("warm_start_free_%s.csv", RUN_TAG))) }
R <- rbindlist(out); cat("\n==== B: SD ratio by pass (collapse < 0.1, healthy ~1.0-1.2) ====\n"); W <- dcast(R, m ~ pass, value.var = "ratio"); num <- setdiff(names(W), "m"); W[, (num) := lapply(.SD, round, 3), .SDcols = num]; print(W)
cat("\nalpha by pass:\n"); W <- dcast(R, m ~ pass, value.var = "alpha"); W[, (num) := lapply(.SD, round, 3), .SDcols = num]; print(W)
