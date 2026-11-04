# ==============================================================================
# R/simulation/sim05_tsd_form.R   (2026-11-01)   RUN ON THE MACHINE WITH INLA
#
# QUESTION: is the time-since-disturbance (TSD) effect better described as LINEAR in TSD (current pipeline) or as
#  log(TSD + 1)? Both forms have the same number of parameters (main disturbance term + one TSD coefficient).
#  Linear:  log HR(TSD) = b + c * TSD[decades]          (constant rate of change)
#  Log:     log HR(TSD) = b + c * log(TSD[years] + 1)   (fast change in the first years, then levelling off)
#  log(0 + 1) = 0, so undisturbed rows keep their structural 0 and nothing else in the model changes.
# REAL DATA, 8 transitions chosen to cover the expected ecological shapes (decided with L, 2026-11-01):
#   pulse at disturbance   : fir->regen, jack_pine->regen
#   recruitment window     : regen->trembling_aspen, regen->white_birch
#   delayed onset          : white_birch->fir, trembling_aspen->fir
#   clocks differ (disturbance during the sojourn): fir->white_birch, fir->yellow_birch
#  Imputations 1-5 of the main run (rows rebuilt exactly as in sim03 Part A / R/tests/08f).
# MODEL (decided 2026-11-01): weibullsurv VARIANT 0 (coefficients = log hazard ratios), coefficient prior N(0, 5^2),
#  flat intercept, pc.alphaw(ALPHA_PRIOR from config.R), one thread. Current term set of the main build (harvest TSD is
#  still ONE term for any harvest; separate light/heavy TSD is not yet in the pipeline).
# FITS per transition x imputation: "linear" (all TSD terms linear), "log_all" (all TSD terms log), and one fit per
#  disturbance with ONLY that disturbance's TSD term in log form ("log_fire", "log_pest", "log_harv"), so each
#  disturbance gets its own linear-vs-log comparison with everything else held linear.
# MEASURES: WAIC and DIC per fit; paired WAIC difference (log minus linear; NEGATIVE = log fits better) with its standard
#  error from the row-wise WAIC contributions (same rows in both fits); disturbance effect curves (main term + TSD
#  coefficient * form(TSD)) at mean GDD (z_gdd_3 = 0, so the GDD x disturbance term is 0), from 250 joint posterior draws
#  per fit, pooled over imputations by stacking draws; observed TSD range and events per TSD band.
# OUTPUT data-output/simulation/sim05_tsd_form/: fits/ (one .rds per transition x imputation), sim05_waic.csv,
#  sim05_waic_summary.csv, sim05_coef.csv, sim05_support.csv, sim05_curves.csv, sim05_curves.png, run_log.txt, progress.log
# READ-ONLY inputs: data/tables/main/, data-output/main/fits/extras_main_m01-05.rds. Nothing else is changed.
# Run from the project root, with nothing else running:
#   Rscript R/simulation/sim05_tsd_form.R 2>&1 | tee data-output/simulation/sim05_run.log   (re-running skips finished fits)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(tools) })
source(here("R", "config.R"))
options(width = 200)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

FIT_ENGINE <- Sys.getenv("SIM05_ENGINE", "inla")              # "mle" = local plumbing check only
N_WORKERS  <- as.integer(Sys.getenv("SIM05_WORKERS", "6"))
REAL_M     <- eval(parse(text = Sys.getenv("SIM05_M", "1:5")))
LABELS     <- c("fir->regen", "jack_pine->regen", "regen->trembling_aspen", "regen->white_birch",
                "white_birch->fir", "trembling_aspen->fir", "fir->white_birch", "fir->yellow_birch")
SIGMA_B    <- 5; VARIANT <- 0L; N_DRAW <- 250L; MASTER_SEED <- 20261105L
DIST       <- list(fire = list(tsd = "tsd_fire", main = "fire_f"), pest = list(tsd = "tsd_pest", main = "pest_f"),
                   harv = list(tsd = "tsd_harv", main = c("harv_light", "harv_heavy", "harv_any")))
TSD_BANDS  <- c(-Inf, 10, 20, 40, Inf)                        # years

OUT <- here("data-output", "simulation", paste0("sim05_tsd_form", if (FIT_ENGINE != "inla") paste0("_", FIT_ENGINE, "_check")))
dir.create(file.path(OUT, "fits"), recursive = TRUE, showWarnings = FALSE)
SIM_M <- 200L; M_MI <- 30L; N_COV_SAMPLES <- 0L
source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
if (FIT_ENGINE == "inla") suppressPackageStartupMessages(library(INLA))
logf <- file.path(OUT, "run_log.txt"); progf <- file.path(OUT, "progress.log")
logl <- function(...) { s <- sprintf(...); cat(s, "\n"); cat(format(Sys.time()), s, "\n", file = logf, append = TRUE) }
prog <- function(console, ...) { s <- sprintf(...); if (console) cat(s, "\n"); cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), s, "\n", file = progf, append = TRUE) }
eta_text <- function(t_start, n_done_run, n_left) {
  if (n_done_run < 1) return("ETA: after the first one finishes")
  el <- as.numeric(difftime(Sys.time(), t_start, units = "mins")); rem <- el / n_done_run * n_left
  sprintf("elapsed %.1f min, ~%.0f min left (finish ~%s)", el, rem, format(Sys.time() + rem * 60, "%H:%M"))
}
self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); if (!length(self)) self <- here("R", "simulation", "sim05_tsd_form.R")
logl("==== start | engine %s | R %s | data.table %s | INLA %s | imputations %s | variant %d | sigma %g | alpha prior pc.alphaw(%g)",
     FIT_ENGINE, getRversion(), packageVersion("data.table"), if (FIT_ENGINE == "inla") as.character(packageVersion("INLA")) else "not used",
     paste(range(REAL_M), collapse = "-"), VARIANT, SIGMA_B, ALPHA_PRIOR)
logl("script %s md5 %s", self, md5sum(self))
for (f in c(file.path(TABLE_DIR, c("visit_intervals.rds", "transition_set.rds", "plot_info.rds")), file.path(FIT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, REAL_M)),
            here("R", c("covariate_spec.R", "04b_build_rows.R", "05a_fit_funs.R", "config.R"))))
  if (file.exists(f)) logl("input %s md5 %s", sub(here(), "", f, fixed = TRUE), md5sum(f)) else logl("input %s MISSING", f)
capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))

# ---- fitting ----------------------------------------------------------------------------------------------------------
# returns list(tab = summary table, waic, dic, local = row-wise WAIC, draws = N_DRAW x (disturbance terms) matrix, error)
fit_inla_f <- function(d, terms, draw_terms, seed) {
  f <- as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", paste(c("1", terms), collapse = " + ")))
  m <- tryCatch(INLA::inla(f, family = "weibullsurv", data = d,
                           control.family = list(variant = VARIANT, hyper = list(alpha = list(prior = "pc.alphaw", param = ALPHA_PRIOR))),
                           control.fixed = list(prec = 1 / SIGMA_B^2, prec.intercept = 0),
                           control.compute = list(config = TRUE, waic = TRUE, dic = TRUE),
                           control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)),
                           num.threads = "1:1"), error = function(e) list(error = conditionMessage(e)))
  pars <- c("alpha", "(Intercept)", terms)
  if (!is.null(m$error)) return(list(tab = data.table(parameter = pars, mean = NA_real_, sd = NA_real_, q025 = NA_real_, q975 = NA_real_),
                                     waic = NA_real_, dic = NA_real_, local = NULL, draws = NULL, error = substr(m$error, 1, 150)))
  cols <- c("mean", "sd", "0.025quant", "0.975quant"); fx <- m$summary.fixed; hy <- m$summary.hyperpar
  g <- function(nm) if (nm == "alpha") unlist(hy[ALPHA_NAME, cols]) else if (nm %in% rownames(fx)) unlist(fx[nm, cols]) else rep(NA_real_, 4)
  v <- t(vapply(pars, g, numeric(4)))
  sel <- setNames(as.list(rep(1L, length(draw_terms))), draw_terms)
  smp <- tryCatch(INLA::inla.posterior.sample(N_DRAW, m, selection = sel, seed = seed, num.threads = "1:1"), error = function(e) NULL)
  dr <- if (is.null(smp)) NULL else t(vapply(smp, function(s) s$latent[paste0(draw_terms, ":1"), 1], numeric(length(draw_terms))))
  if (!is.null(dr)) { dr <- matrix(dr, ncol = length(draw_terms)); colnames(dr) <- draw_terms }
  list(tab = data.table(parameter = pars, mean = v[, 1], sd = v[, 2], q025 = v[, 3], q975 = v[, 4]),
       waic = m$waic$waic, dic = m$dic$dic, local = m$waic$local.waic, draws = dr, error = NA_character_)
}
fit_map_f <- function(d, terms, draw_terms, seed) {            # plumbing-only stand-in: MAP, Wald, Gaussian draws, -2 loglik as "WAIC"
  X <- cbind(1, as.matrix(d[, ..terms])); X[is.na(X)] <- 0; T <- d$Tstop; tr <- d$truncation; s <- d$status; lT <- log(T)   # NA covariate = 0, as INLA does
  lli <- function(th) { la <- th[1]; a <- exp(la); eta <- drop(X %*% th[-1]); s * (la + eta + (a - 1) * lT) - exp(eta) * (T^a - tr^a) }
  nll <- function(th) -(sum(lli(th)) - 1 / SIGMA_B^2 / 2 * sum(th[-(1:2)]^2))
  o <- optim(c(0, log(max(sum(s), 0.5) / sum(T - tr)), rep(0, length(terms))), nll, method = "BFGS", control = list(maxit = 10000, reltol = 1e-12))
  V <- solve(optimHess(o$par, nll)); se <- sqrt(pmax(diag(V), 0)); pars <- c("alpha", "(Intercept)", terms)
  set.seed(seed); idx <- match(draw_terms, pars); L <- chol(V[idx, idx, drop = FALSE])
  dr <- sweep(matrix(rnorm(N_DRAW * length(idx)), N_DRAW) %*% L, 2, o$par[idx], "+"); colnames(dr) <- draw_terms
  est <- c(exp(o$par[1]), o$par[-1])
  list(tab = data.table(parameter = pars, mean = est, sd = c(est[1] * se[1], se[-1]), q025 = est - 1.96 * c(est[1] * se[1], se[-1]), q975 = est + 1.96 * c(est[1] * se[1], se[-1])),
       waic = -2 * sum(lli(o$par)), dic = NA_real_, local = -2 * lli(o$par), draws = dr, error = NA_character_)
}
fit_f <- function(...) if (FIT_ENGINE == "inla") fit_inla_f(...) else fit_map_f(...)

# ---- data --------------------------------------------------------------------------------------------------------------
inp <- prepare_fit_inputs(TABLE_DIR)
stopifnot(all(LABELS %in% inp$td$label))
key_iv <- paste(inp$iv$id_pe, inp$iv$interval); b_abs <- inp$iv$t_end
rebuild <- function(m) {                                       # identical to sim03 Part A / R/tests/08f_failed_vs_valid.R
  ex <- readRDS(file.path(FIT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, m))); p1 <- as.data.table(ex$p1); p2 <- as.data.table(ex$p2); rm(ex)
  ex_row <- rep(NA_real_, nrow(inp$iv)); ii <- match(paste(p2$id_pe, p2$interval), key_iv); stopifnot(!anyNA(ii)); ex_row[ii] <- p2$exit_abs
  p2f <- function(a_dur, b_dur, from, to, r) { E <- b_abs[r] - b_dur; min(max(ex_row[r] - E, a_dur), b_dur) }
  build_rows(inp$iv, setNames(p1$entry, as.character(p1$id_pe)), inp$states, inp$tr, p2_fun = p2f, keep_cols = inp$keep_cols)$rows[kept == TRUE]
}
dist_in <- function(tms) Filter(function(nm) DIST[[nm]]$tsd %in% tms, names(DIST))
mains_of <- function(nm, tms) intersect(DIST[[nm]]$main, tms)

# ---- run: one job = one transition x one imputation (all forms, so WAIC differences are paired on identical rows) -------
run_m <- function(m, d_all) {
  jobs <- LABELS[!file.exists(file.path(OUT, "fits", sprintf("m%02d_%s.rds", m, gsub("->", "_to_", LABELS))))]
  res <- mclapply(jobs, function(lb) {
    k <- match(lb, inp$td$label); tms <- colnames(inp$mask)[inp$mask[k, ]]; ds <- dist_in(tms)
    d <- d_all[label == lb, c("status", "Tstop", "truncation", tms), with = FALSE]
    draw_terms <- unlist(lapply(ds, function(nm) c(mains_of(nm, tms), DIST[[nm]]$tsd)))
    sup <- rbindlist(lapply(ds, function(nm) { fl <- d[[DIST[[nm]]$tsd]]; on <- if (nm == "harv") rowSums(d[, mains_of(nm, tms), with = FALSE]) > 0 else d[[DIST[[nm]]$main]] == 1
      yrs <- 10 * fl[on]; ev <- d$status[on] == 1
      data.table(dist = nm, rows = sum(on), events = sum(ev), tsd_min = min(yrs), tsd_q05 = quantile(yrs, 0.05), tsd_q95 = quantile(yrs, 0.95), tsd_max = max(yrs),
                 ev_tsd_q05 = if (any(ev)) quantile(yrs[ev], 0.05) else NA_real_, ev_tsd_q95 = if (any(ev)) quantile(yrs[ev], 0.95) else NA_real_,
                 band = list(table(cut(yrs[ev], TSD_BANDS, right = FALSE)))) }))
    forms <- c("linear", "log_all", if (length(ds) > 1) paste0("log_", ds))
    fits <- lapply(setNames(forms, forms), function(fm) {
      dd <- copy(d); logd <- if (fm == "linear") character(0) else if (fm == "log_all") ds else sub("^log_", "", fm)
      for (nm in logd) { v <- DIST[[nm]]$tsd; set(dd, j = v, value = log1p(10 * dd[[v]])) }   # decades -> years -> log(years + 1); 0 stays 0
      seed <- MASTER_SEED + 1000L * m + 10L * k + match(fm, c("linear", "log_all", paste0("log_", names(DIST))))
      t1 <- proc.time()[["elapsed"]]; r <- fit_f(dd, tms, draw_terms, seed); r$sec <- proc.time()[["elapsed"]] - t1; r$log_dist <- logd; r$seed <- seed; r })
    out <- list(label = lb, m = m, terms = tms, dists = ds, draw_terms = draw_terms, support = sup, n_events = sum(d$status), n_rows = nrow(d), fits = fits)
    saveRDS(out, file.path(OUT, "fits", sprintf("m%02d_%s.rds", m, gsub("->", "_to_", lb))))
    n_done <- sum(file.exists(file.path(OUT, "fits", sprintf("m%02d_%s.rds", rep(REAL_M, each = length(LABELS)), gsub("->", "_to_", LABELS)))))
    prog(TRUE, "[sim05] %d/%d done: m%02d %s | %d fits, INLA errors %d | WAIC log_all - linear = %.1f | %s", n_done, length(REAL_M) * length(LABELS), m, lb,
         length(fits), sum(!is.na(vapply(fits, `[[`, "", "error"))), fits$log_all$waic - fits$linear$waic, eta_text(T0, n_done - N0, length(REAL_M) * length(LABELS) - n_done))
    NULL
  }, mc.cores = N_WORKERS, mc.preschedule = FALSE)
  bad <- vapply(res, inherits, logical(1), "try-error"); if (any(bad)) stop("worker error: ", paste(res[bad], collapse = " | "))
}
all_files <- file.path(OUT, "fits", sprintf("m%02d_%s.rds", rep(REAL_M, each = length(LABELS)), gsub("->", "_to_", rep(LABELS, length(REAL_M)))))
N0 <- sum(file.exists(all_files)); T0 <- Sys.time()
prog(TRUE, "[sim05] start: %d of %d transition x imputation jobs to run (%d workers)", length(all_files) - N0, length(all_files), N_WORKERS)
for (m in REAL_M) {
  if (all(file.exists(file.path(OUT, "fits", sprintf("m%02d_%s.rds", m, gsub("->", "_to_", LABELS)))))) next
  rows <- rebuild(m); rows[, label := paste(from, to, sep = "->")]; d_all <- rows[label %in% LABELS]; rm(rows); invisible(gc())
  run_m(m, d_all); rm(d_all); invisible(gc())
}

# ---- summaries ---------------------------------------------------------------------------------------------------------
R <- lapply(all_files, readRDS)
W <- rbindlist(lapply(R, function(r) {
  lin <- r$fits$linear
  rbindlist(lapply(setdiff(names(r$fits), "linear"), function(fm) { f <- r$fits[[fm]]
    ok <- is.na(lin$error) && is.na(f$error) && !is.null(lin$local) && !is.null(f$local)
    dl <- if (ok) f$local - lin$local else NA_real_
    data.table(label = r$label, m = r$m, comparison = fm, n_events = r$n_events, n_rows = r$n_rows,
               waic_linear = lin$waic, waic_log = f$waic, dwaic = f$waic - lin$waic, dwaic_se = if (ok) sqrt(length(dl)) * sd(dl) else NA_real_,
               ddic = f$dic - lin$dic, error_linear = lin$error, error_log = f$error) }))
}))
fwrite(W, file.path(OUT, "sim05_waic.csv"))
WS <- W[, .(imputations = .N, errors = sum(!is.na(error_linear) | !is.na(error_log)), dwaic_median = round(median(dwaic, na.rm = TRUE), 1),
            dwaic_min = round(min(dwaic, na.rm = TRUE), 1), dwaic_max = round(max(dwaic, na.rm = TRUE), 1), se_median = round(median(dwaic_se, na.rm = TRUE), 1),
            ddic_median = round(median(ddic, na.rm = TRUE), 1)), by = .(label, comparison)]
WS[, verdict := fifelse(dwaic_max < -2 * se_median & dwaic_max < -4, "log better (all imputations)",
                fifelse(dwaic_min > 2 * se_median & dwaic_min > 4, "linear better (all imputations)", "no clear difference"))]
fwrite(WS, file.path(OUT, "sim05_waic_summary.csv"))
C <- rbindlist(lapply(R, function(r) rbindlist(lapply(names(r$fits), function(fm) r$fits[[fm]]$tab[, `:=`(label = r$label, m = r$m, form = fm, error = r$fits[[fm]]$error, sec = r$fits[[fm]]$sec)]))))
fwrite(C, file.path(OUT, "sim05_coef.csv"))
S <- rbindlist(lapply(R, function(r) r$support[, `:=`(label = r$label, m = r$m, ev_tsd_0_10 = sapply(band, `[`, 1), ev_tsd_10_20 = sapply(band, `[`, 2),
                                                      ev_tsd_20_40 = sapply(band, `[`, 3), ev_tsd_40plus = sapply(band, `[`, 4))][, band := NULL]))
SS <- S[, lapply(.SD, function(x) round(median(x), 1)), by = .(label, dist), .SDcols = setdiff(names(S), c("label", "dist", "m"))]
fwrite(SS, file.path(OUT, "sim05_support.csv"))
# curves: effect(TSD) = main + c * form(TSD), GDD at its mean; grid from 0 to the 95th percentile of TSD among disturbed rows
CV <- rbindlist(lapply(LABELS, function(lb) {
  rr <- Filter(function(r) r$label == lb, R); tms <- rr[[1]]$terms
  rbindlist(lapply(rr[[1]]$dists, function(nm) {
    tmax <- SS[label == lb & dist == nm, tsd_q95]; grid <- seq(0, max(tmax, 1), length.out = 61); tv <- DIST[[nm]]$tsd
    rbindlist(lapply(mains_of(nm, tms), function(mn) rbindlist(lapply(c("linear", "log_all"), function(fm) {
      D <- do.call(rbind, lapply(rr, function(r) r$fits[[fm]]$draws)); if (is.null(D)) return(NULL)
      x <- if (fm == "linear") grid / 10 else log1p(grid)
      E <- outer(D[, mn], rep(1, length(grid))) + outer(D[, tv], x)
      data.table(label = lb, dist = nm, main = mn, form = fifelse(fm == "linear", "linear", "log"), tsd_years = grid,
                 median = apply(E, 2, median), q025 = apply(E, 2, quantile, 0.025), q975 = apply(E, 2, quantile, 0.975), draws = nrow(D)) }))))
  }))
}))
fwrite(CV, file.path(OUT, "sim05_curves.csv"))
ok_plot <- tryCatch({
  suppressPackageStartupMessages(library(ggplot2))
  CV[, panel := sprintf("%s | %s", label, fifelse(main %in% c("harv_light", "harv_heavy", "harv_any"), main, dist))]
  p <- ggplot(CV, aes(tsd_years, median, colour = form, fill = form, linetype = form)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.15, colour = NA) + geom_line(linewidth = 0.8) +
    facet_wrap(~ panel, scales = "free", ncol = 4) +
    scale_colour_manual(values = c(linear = "#5a5a57", log = "#2a78d6")) + scale_fill_manual(values = c(linear = "#5a5a57", log = "#2a78d6")) +
    scale_linetype_manual(values = c(linear = "22", log = "solid")) +
    labs(x = "Time since disturbance (years; 0 to the 95th percentile of disturbed rows)", y = "Disturbance effect, log hazard ratio (at mean GDD)",
         title = "sim05: linear vs log(TSD + 1), real data, variant 0, N(0, 5^2), imputations pooled",
         subtitle = "Line = posterior median, band = 95% interval. Panel = transition | disturbance (harvest: per harvest level).") +
    theme_minimal(base_size = 9) + theme(legend.position = "top", strip.text = element_text(size = 7.5, hjust = 0), plot.background = element_rect(fill = "white", colour = NA))
  ggsave(file.path(OUT, "sim05_curves.png"), p, width = 14, height = 2.6 * ceiling(uniqueN(CV$panel) / 4) + 1, dpi = 130, limitsize = FALSE); TRUE
}, error = function(e) { logl("plot failed: %s", conditionMessage(e)); FALSE })
cat("\nWAIC log minus linear (negative = log better), per transition x comparison\n"); print(WS)
cat("\nTSD support (medians over imputations)\n"); print(SS)
logl("outputs: sim05_waic.csv, sim05_waic_summary.csv, sim05_coef.csv, sim05_support.csv, sim05_curves.csv%s", if (ok_plot) ", sim05_curves.png" else "")
logl("==== done")
