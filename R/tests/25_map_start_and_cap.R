# ==============================================================================
# R/tests/25_map_start_and_cap.R  (2026-11-05) -- RUN ON THE MAC (INLA). Reads the WP6 run, edits nothing.
# Saved collapsed black_spruce->regen datasets (baddata_final_mXX.rds), default m = 1 37 7 2 8 5 (user choice 2026-11-05); other m as command-line arguments.
# Per dataset, four INLA/MAP results:
#   MAP   fit_map_v0 (05a): same log-posterior, no cap
#   cold  INLA exactly as in the run (fit_inla_models_pc, exact = TRUE); must reproduce the stored mi_raw alpha / SD (printed)
#   Q2    same INLA call (replicated, checked once against fit_inla_models_pc with an empty control.mode) with
#         control.mode = list(theta = <alpha MAP on INLA's internal scale>, x = <latent mode with eta and the fixed effects set to the MAP>, restart = TRUE)
#   Q3    fit_inla_models_pc on the rows with MAP H(tr) > CAP (22.85) removed (one pass, MAP of the full data); no rows above the cap = identical to cold, not refitted
# healthy (as in Q1) = ratio >= 0.1 and z_map_alpha <= 3 and z_map_intercept <= 3, z against the MAP of the data that fit used; ratio uses the events of those data.
# Stops (does not guess) if INLA's internal alpha scale or latent layout does not match what the script expects.
# Run from the project root:  Rscript R/tests/25_map_start_and_cap.R 2>&1 | tee data-output/diagnostics/25_map_start_and_cap.log     (resumable: one rds per m)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(INLA) })
source(here::here("R", "config.R")); options(width = 220)
INLA::inla.setOption(num.threads = "1:1")
SIM_M <- 200L; M_MI <- 30L; SD_MIN <- 1e-3; N_COV_SAMPLES <- 0L
source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
LABEL <- "black_spruce->regen"; RATIO_MIN <- 0.1; Z_HEALTHY <- 3
args <- suppressWarnings(as.integer(commandArgs(TRUE))); M_SET <- if (length(args) && !anyNA(args)) args else c(1L, 37L, 7L, 2L, 8L, 5L)
OUT <- here("data-output", "diagnostics", "25_map_start_and_cap"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
mr <- rbindlist(lapply(file.path(FIT_DIR, sprintf(c("mi_raw_%s.csv", "mi_raw_redraw_%s.csv"), RUN_TAG)), fread), fill = TRUE)[label == LABEL]

fit_cm <- function(dd, rhs, cm) tryCatch(   # = the INLA call of fit_inla_models_pc (exact = TRUE) plus control.mode
  INLA::inla(as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", rhs)), family = "weibullsurv", data = dd,
             control.family = list(variant = WEIBULL_VARIANT, hyper = list(alpha = list(prior = "pc.alphaw", param = ALPHA_PRIOR))),
             control.fixed = list(mean = 0, prec = 1 / COEF_PRIOR_SD^2, mean.intercept = 0, prec.intercept = 0),
             control.compute = list(config = TRUE), control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)),
             control.mode = cm, num.threads = "1:1"), error = function(e) list(error = conditionMessage(e)))
tab_of <- function(r, tms) {   # as tidy_fit
  pars <- c("alpha", "(Intercept)", tms)
  if (!is.null(r$error)) return(data.table(parameter = pars, est = NA_real_, sd = NA_real_, fit_error = TRUE))
  fx <- r$summary.fixed; hy <- r$summary.hyperpar; g <- function(nm, col) if (nm %in% rownames(fx)) fx[nm, col] else NA_real_
  data.table(parameter = pars, est = c(hy[ALPHA_NAME, "mean"], g("(Intercept)", "mean"), vapply(tms, g, numeric(1), col = "mean")),
             sd = c(hy[ALPHA_NAME, "sd"], g("(Intercept)", "sd"), vapply(tms, g, numeric(1), col = "sd")), fit_error = FALSE) }
summ <- function(variant, r, dd, tms, secs) {
  tab <- tab_of(r, tms); cap <- cap_diagnostics(dd, tms, tab); dg <- fit_diagnostics(r); al <- tab[parameter == "alpha"]; ic <- tab[parameter == "(Intercept)"]
  ev <- sum(dd$status == 1L); ratio <- al$sd / (0.78 * al$est / sqrt(ev))
  list(summary = data.table(variant = variant, rows = nrow(dd), events = ev, alpha = al$est, alpha_sd = al$sd, ratio = ratio, icpt = ic$est, icpt_sd = ic$sd,
         z_alpha = cap[["z_map_alpha"]], z_icpt = cap[["z_map_intercept"]], z_terms_max = cap[["z_map_max"]], maxH_inla = cap[["maxH_inla"]], n_above_inla = cap[["n_above_inla"]],
         maxH_map = cap[["maxH_map"]], n_above_map = cap[["n_above_map"]], mlik = dg$mlik, mode_status = dg$mode_status, secs = secs, err = if (is.null(r$error)) NA_character_ else substr(r$error, 1, 80),
         healthy = isTRUE(ratio >= RATIO_MIN && cap[["z_map_alpha"]] <= Z_HEALTHY && cap[["z_map_intercept"]] <= Z_HEALTHY)),
       est = tab[, .(variant = variant, parameter, est, sd)]) }

run_m <- function(mm, first) {
  b <- readRDS(file.path(FIT_DIR, sprintf("baddata_%s_m%02d.rds", RUN_TAG, mm)))[[LABEL]]; stopifnot("no baddata for this label" = !is.null(b))
  d <- as.data.frame(b$data); rhs <- b$rhs; st <- mr[m == mm]; tms <- setdiff(st$parameter, c("alpha", "(Intercept)"))
  stopifnot(setequal(all.vars(as.formula(paste("~", rhs))), tms), all(c("status", "Tstop", "truncation", tms) %in% names(d)))
  cat(sprintf("\n######## m = %d: %d rows, %d events, %d terms\n", mm, nrow(d), sum(d$status == 1L), length(tms)))
  # MAP, cold fit
  fm <- fit_map_v0(as.data.table(d), tms); mapv <- setNames(fm$tab$est, fm$tab$parameter)
  t0 <- proc.time()[["elapsed"]]; r0 <- fit_inla_models_pc(list(d), rhs, exact = TRUE)[[1]]; s0 <- proc.time()[["elapsed"]] - t0
  stopifnot("cold fit returned an error" = is.null(r0$error))
  sa <- st[parameter == "alpha"]; sc <- summ("cold", r0, d, tms, s0)
  cat(sprintf("Cold fit vs stored mi_raw: alpha %.6f vs %.6f, alpha SD %.3e vs %.3e (%.0f s)\n", sc$summary$alpha, sa$est, sc$summary$alpha_sd, sa$sd, s0))
  # INLA internal alpha scale and latent layout: checked, not assumed
  hy1 <- INLA::inla.models()$likelihood$weibullsurv$hyper[[1]]; stopifnot("to.theta / from.theta not available" = is.function(hy1$to.theta) && is.function(hy1$from.theta))
  a_mode <- hy1$from.theta(unname(r0$mode$theta[1]))
  cat(sprintf("alpha from mode$theta via from.theta: %.6f (posterior mean %.6f); to.theta(MAP alpha %.4f) = %.4f; 10*log(alpha) = %.4f\n", a_mode, sc$summary$alpha, mapv[["alpha"]], hy1$to.theta(mapv[["alpha"]]), 10 * log(mapv[["alpha"]])))
  stopifnot("theta scale check failed: from.theta(mode$theta) is not the fitted alpha (diagnose, do not guess)" = abs(a_mode / sc$summary$alpha - 1) < 0.01)
  cn <- r0$misc$configs$contents; x0 <- r0$mode$x
  cat("latent layout (tag, start, length):\n"); print(data.table(tag = cn$tag, start = cn$start, length = cn$length), row.names = FALSE)
  stopifnot("latent layout: sum of lengths != length(mode$x)" = length(x0) == sum(cn$length), "no 'Predictor' block of nrow(d)" = any(cn$tag == "Predictor" & cn$length == nrow(d)))
  eta <- unname(mapv[["(Intercept)"]] + drop(as.matrix(d[, tms, drop = FALSE]) %*% mapv[tms])); x1 <- x0
  for (i in seq_along(cn$tag)) { idx <- cn$start[i] - 1L + seq_len(cn$length[i]); tg <- cn$tag[i]
    if (tg == "Predictor") x1[idx] <- eta else if (tg %in% names(mapv) && length(idx) == 1L) x1[idx] <- mapv[[tg]] else stop("latent tag not understood: ", tg) }
  cat(sprintf("start x: max |x1 - x0| over the fixed effects %.3f; over eta %.3f\n", max(abs(x1 - x0)[cn$start[cn$tag != "Predictor"]]), max(abs(x1 - x0)[seq_len(nrow(d))])))
  if (first) { tw <- fit_cm(d, rhs, list()); cat(sprintf("Replicated call with empty control.mode vs fit_inla_models_pc: alpha diff %.2e, mlik diff %.2e\n", abs(tw$summary.hyperpar[ALPHA_NAME, "mean"] - sc$summary$alpha), abs(tw$mlik[1, 1] - sc$summary$mlik)))
    stopifnot("replicated INLA call differs from fit_inla_models_pc" = abs(tw$summary.hyperpar[ALPHA_NAME, "mean"] / sc$summary$alpha - 1) < 1e-8) }
  # Q2
  t0 <- proc.time()[["elapsed"]]; r2 <- fit_cm(d, rhs, list(theta = hy1$to.theta(mapv[["alpha"]]), x = x1, restart = TRUE)); s2 <- proc.time()[["elapsed"]] - t0
  if (!is.null(r2$error)) cat("Q2 error:", r2$error, "\n"); s2s <- summ("Q2_warm", r2, d, tms, s2)
  # Q3
  Hrow <- exp(eta) * d$truncation^mapv[["alpha"]]; above <- Hrow > CAP
  cat(sprintf("Q3: rows with MAP H(tr) > %.2f: %d (events among them %d); H of those: %s\n", CAP, sum(above), sum(d$status[above] == 1L), paste(round(Hrow[above], 1), collapse = ", ")))
  if (any(above)) { d3 <- d[!above, , drop = FALSE]; t0 <- proc.time()[["elapsed"]]; r3 <- fit_inla_models_pc(list(d3), rhs, exact = TRUE)[[1]]; s3 <- proc.time()[["elapsed"]] - t0
    if (!is.null(r3$error)) cat("Q3 error:", r3$error, "\n"); s3s <- summ("Q3_cap", r3, d3, tms, s3); f3 <- fit_map_v0(as.data.table(d3), tms)
    m3 <- data.table(variant = "MAP_Q3data", f3$tab[, .(parameter, est, sd)]) } else { s3s <- summ("Q3_cap", r0, d, tms, s0); s3s$summary[, `:=`(secs = 0, err = "no row above the cap: identical to cold, not refitted")]; s3s$est[, variant := "Q3_cap"]; m3 <- NULL }
  mp <- data.table(variant = "MAP", fm$tab[, .(parameter, est, sd)])
  Ssum <- rbindlist(list(data.table(variant = "MAP", rows = nrow(d), events = sum(d$status == 1L), alpha = mapv[["alpha"]], icpt = mapv[["(Intercept)"]], maxH_map = max(Hrow), n_above_map = sum(above)),
                         sc$summary, s2s$summary, s3s$summary), fill = TRUE)
  Eest <- rbindlist(list(mp, sc$est, s2s$est, s3s$est, m3), fill = TRUE)
  res <- list(m = mm, S = cbind(m = mm, Ssum), E = cbind(m = mm, Eest), removed = data.table(m = mm, status = d$status[above], Tstop = d$Tstop[above], truncation = d$truncation[above], H_map = Hrow[above]))
  saveRDS(res, file.path(OUT, sprintf("m%03d.rds", mm))); res }

first <- TRUE; RES <- list()
for (mm in M_SET) { f <- file.path(OUT, sprintf("m%03d.rds", mm))
  RES[[as.character(mm)]] <- if (file.exists(f)) { cat("m =", mm, ": loaded", f, "\n"); readRDS(f) } else tryCatch({ r <- run_m(mm, first); first <- FALSE; r }, error = function(e) { cat("!! m =", mm, "stopped:", conditionMessage(e), "\n"); NULL }) }
RES <- Filter(Negate(is.null), RES)
S <- rbindlist(lapply(RES, `[[`, "S"), fill = TRUE); E <- rbindlist(lapply(RES, `[[`, "E"), fill = TRUE); fwrite(S, file.path(OUT, "summary.csv")); fwrite(E, file.path(OUT, "estimates.csv"))
cat("\n=========== One row per fit (healthy = ratio >= 0.1 and z_alpha, z_icpt <= 3; z against the MAP of the data that fit used) ===========\n")
print(S[, .(m, variant, rows, events, alpha = round(alpha, 4), alpha_sd = signif(alpha_sd, 3), ratio = signif(ratio, 3), icpt = round(icpt, 3), z_alpha = round(z_alpha, 2), z_icpt = round(z_icpt, 2),
            z_terms_max = round(z_terms_max, 2), maxH_inla = round(maxH_inla, 2), n_above_inla, maxH_map = round(maxH_map, 2), n_above_map, mode_status, secs = round(secs), healthy, err)], row.names = FALSE)
cat("\n=========== Estimates per dataset (INLA mean, MAP estimate) ===========\n")
for (mm in names(RES)) { cat("\n--- m =", mm, "---\n"); W <- dcast(E[m == as.integer(mm)], parameter ~ variant, value.var = "est"); W <- W[match(unique(E[m == as.integer(mm)]$parameter), parameter)]
  print(W, digits = 4, row.names = FALSE); rmv <- RES[[mm]]$removed; if (nrow(rmv)) { cat("rows removed in Q3:\n"); print(rmv, row.names = FALSE) } }
cat("\nSummary: healthy fits per variant:\n"); print(S[variant %in% c("cold", "Q2_warm", "Q3_cap"), .(n = .N, healthy = sum(healthy)), by = variant])
