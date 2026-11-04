# ==============================================================================
# R/simulation/sim05b_tsd_check.R   (2026-11-01)   RUN ON THE MACHINE WITH INLA, AFTER sim05
#
# QUESTION: what do the DATA say about the shape of the time-since-disturbance (TSD) effect, without assuming a form?
#  Observed / expected check (functional-form diagnostic): for each transition x disturbance of sim05, refit the model
#  WITHOUT that disturbance's TSD term (everything else as in sim05: variant 0, N(0, 5^2), pc.alphaw(ALPHA_PRIOR), other
#  TSD terms linear). The disturbance main term then carries the AVERAGE disturbance effect. For every disturbed row the
#  fit gives the expected number of events, E = exp(eta) * (Tstop^alpha - truncation^alpha) at the posterior means.
#  Grouping disturbed rows by TSD band, log(observed / expected) shows how the hazard departs from the average
#  disturbance effect at that TSD: above 0 = more transitions than the average effect predicts.
#  Overlay: the sim05 linear and log(TSD + 1) TSD curves (posterior median of the TSD coefficient, imputations pooled),
#  centred so that their expected-weighted average over the disturbed rows is 0 (same reference as the points).
#  Observed and expected are averaged over imputations 1-5; bands start at 0-5-10-20-30-45-60-inf years and are merged
#  with a neighbour until each has >= MIN_EV events (average per imputation). 95% interval: log(O/E) +- 1.96 / sqrt(O).
#  CAVEAT: where the state starts at the disturbance (regen->..., ...->regen), TSD ~ time in state, so the Weibull shape
#  alpha absorbs part of the TSD pattern; the points then show what is left after alpha.
# OUTPUT data-output/simulation/sim05b_tsd_check/: fits/ (one .rds per transition x imputation), sim05b_bands.csv,
#  sim05b_curves.csv, sim05b_check.png, run_log.txt, progress.log
# READ-ONLY inputs: as sim05, plus data-output/simulation/sim05_tsd_form/fits/.
# Run from the project root, with nothing else running:
#   Rscript R/simulation/sim05b_tsd_check.R 2>&1 | tee data-output/simulation/sim05b_run.log   (re-running skips finished fits)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(tools) })
source(here("R", "config.R"))
options(width = 200)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

FIT_ENGINE <- Sys.getenv("SIM05B_ENGINE", "inla")             # "mle" = local plumbing check only
N_WORKERS  <- as.integer(Sys.getenv("SIM05B_WORKERS", "6"))
REAL_M     <- eval(parse(text = Sys.getenv("SIM05B_M", "1:5")))
LABELS     <- c("fir->regen", "jack_pine->regen", "regen->trembling_aspen", "regen->white_birch",
                "white_birch->fir", "trembling_aspen->fir", "fir->white_birch", "fir->yellow_birch")      # = sim05
SIGMA_B    <- 5; VARIANT <- 0L; MIN_EV <- 10
BREAKS     <- c(0, 5, 10, 20, 30, 45, 60, Inf)               # years
DIST       <- list(fire = list(tsd = "tsd_fire", main = "fire_f"), pest = list(tsd = "tsd_pest", main = "pest_f"),
                   harv = list(tsd = "tsd_harv", main = c("harv_light", "harv_heavy", "harv_any")))
SIM05_DIR  <- here("data-output", "simulation", paste0("sim05_tsd_form", if (FIT_ENGINE != "inla") paste0("_", FIT_ENGINE, "_check")))

OUT <- here("data-output", "simulation", paste0("sim05b_tsd_check", if (FIT_ENGINE != "inla") paste0("_", FIT_ENGINE, "_check")))
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
self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); if (!length(self)) self <- here("R", "simulation", "sim05b_tsd_check.R")
logl("==== start | engine %s | R %s | INLA %s | imputations %s | variant %d | sigma %g | alpha prior pc.alphaw(%g) | min events per band %g",
     FIT_ENGINE, getRversion(), if (FIT_ENGINE == "inla") as.character(packageVersion("INLA")) else "not used", paste(range(REAL_M), collapse = "-"),
     VARIANT, SIGMA_B, ALPHA_PRIOR, MIN_EV)
logl("script %s md5 %s", self, md5sum(self))
for (f in c(file.path(TABLE_DIR, c("visit_intervals.rds", "transition_set.rds", "plot_info.rds")), file.path(FIT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, REAL_M)),
            here("R", c("covariate_spec.R", "04b_build_rows.R", "05a_fit_funs.R", "config.R")), file.path(SIM05_DIR, "sim05_waic_summary.csv")))
  if (file.exists(f)) logl("input %s md5 %s", sub(here(), "", f, fixed = TRUE), md5sum(f)) else logl("input %s MISSING", f)
capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))

# ---- fitting: posterior means of alpha and the fixed effects ---------------------------------------------------------
fit_inla_b <- function(d, terms) {
  f <- as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", paste(c("1", terms), collapse = " + ")))
  m <- tryCatch(INLA::inla(f, family = "weibullsurv", data = d,
                           control.family = list(variant = VARIANT, hyper = list(alpha = list(prior = "pc.alphaw", param = ALPHA_PRIOR))),
                           control.fixed = list(prec = 1 / SIGMA_B^2, prec.intercept = 0),
                           control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)),
                           num.threads = "1:1"), error = function(e) list(error = conditionMessage(e)))
  if (!is.null(m$error)) return(list(alpha = NA_real_, beta = NULL, error = substr(m$error, 1, 150)))
  fx <- m$summary.fixed
  list(alpha = m$summary.hyperpar[ALPHA_NAME, "mean"], beta = setNames(fx[c("(Intercept)", terms), "mean"], c("(Intercept)", terms)), error = NA_character_)
}
fit_map_b <- function(d, terms) {                              # plumbing-only stand-in (MAP)
  X <- cbind(1, as.matrix(d[, ..terms])); X[is.na(X)] <- 0; T <- d$Tstop; tr <- d$truncation; s <- d$status; lT <- log(T)
  nll <- function(th) { la <- th[1]; a <- exp(la); eta <- drop(X %*% th[-1])
    -(sum(s * (la + eta + (a - 1) * lT) - exp(eta) * (T^a - tr^a)) - 1 / SIGMA_B^2 / 2 * sum(th[-(1:2)]^2)) }
  o <- optim(c(0, log(max(sum(s), 0.5) / sum(T - tr)), rep(0, length(terms))), nll, method = "BFGS", control = list(maxit = 10000, reltol = 1e-12))
  list(alpha = exp(o$par[1]), beta = setNames(o$par[-1], c("(Intercept)", terms)), error = NA_character_)
}
fit_b <- function(...) if (FIT_ENGINE == "inla") fit_inla_b(...) else fit_map_b(...)

# ---- data (identical to sim05) -----------------------------------------------------------------------------------------
inp <- prepare_fit_inputs(TABLE_DIR)
stopifnot(all(LABELS %in% inp$td$label))
key_iv <- paste(inp$iv$id_pe, inp$iv$interval); b_abs <- inp$iv$t_end
rebuild <- function(m) {
  ex <- readRDS(file.path(FIT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, m))); p1 <- as.data.table(ex$p1); p2 <- as.data.table(ex$p2); rm(ex)
  ex_row <- rep(NA_real_, nrow(inp$iv)); ii <- match(paste(p2$id_pe, p2$interval), key_iv); stopifnot(!anyNA(ii)); ex_row[ii] <- p2$exit_abs
  p2f <- function(a_dur, b_dur, from, to, r) { E <- b_abs[r] - b_dur; min(max(ex_row[r] - E, a_dur), b_dur) }
  build_rows(inp$iv, setNames(p1$entry, as.character(p1$id_pe)), inp$states, inp$tr, p2_fun = p2f, keep_cols = inp$keep_cols)$rows[kept == TRUE]
}
dist_in <- function(tms) Filter(function(nm) DIST[[nm]]$tsd %in% tms, names(DIST))
mains_of <- function(nm, tms) intersect(DIST[[nm]]$main, tms)

# ---- run: one job = one transition x one imputation; one refit per disturbance (its TSD term removed) ------------------
fname <- function(m, lb) file.path(OUT, "fits", sprintf("m%02d_%s.rds", m, gsub("->", "_to_", lb)))
all_files <- unlist(lapply(REAL_M, function(m) fname(m, LABELS)))
N0 <- sum(file.exists(all_files)); T0 <- Sys.time()
prog(TRUE, "[sim05b] start: %d of %d transition x imputation jobs to run (%d workers)", length(all_files) - N0, length(all_files), N_WORKERS)
for (m in REAL_M) {
  jobs <- LABELS[!file.exists(fname(m, LABELS))]; if (!length(jobs)) next
  rows <- rebuild(m); rows[, label := paste(from, to, sep = "->")]; d_all <- rows[label %in% LABELS]; rm(rows); invisible(gc())
  res <- mclapply(jobs, function(lb) {
    k <- match(lb, inp$td$label); tms <- colnames(inp$mask)[inp$mask[k, ]]
    d <- d_all[label == lb, c("status", "Tstop", "truncation", tms), with = FALSE]
    out <- rbindlist(lapply(dist_in(tms), function(nm) {
      tv <- DIST[[nm]]$tsd; tms_b <- setdiff(tms, tv); t1 <- proc.time()[["elapsed"]]; fb <- fit_b(d, tms_b); sec <- proc.time()[["elapsed"]] - t1
      if (!is.na(fb$error)) return(data.table(dist = nm, main = NA_character_, tsd_years = NA_real_, obs = NA_real_, exp = NA_real_, error = fb$error, sec = sec))
      X <- cbind(1, as.matrix(d[, ..tms_b])); X[is.na(X)] <- 0
      E <- exp(drop(X %*% fb$beta)) * (d$Tstop^fb$alpha - d$truncation^fb$alpha)
      rbindlist(lapply(mains_of(nm, tms), function(mn) { on <- d[[mn]] == 1
        data.table(dist = nm, main = mn, tsd_years = 10 * d[[tv]][on], obs = d$status[on], exp = E[on], error = NA_character_, sec = sec) }))
    }))
    out[, `:=`(label = lb, m = m)]
    saveRDS(out, fname(m, lb))
    n_done <- sum(file.exists(all_files))
    prog(TRUE, "[sim05b] %d/%d done: m%02d %s | refits %d, INLA errors %d | %s", n_done, length(all_files), m, lb, uniqueN(out$dist),
         uniqueN(out[!is.na(error), dist]), eta_text(T0, n_done - N0, length(all_files) - n_done))
    NULL
  }, mc.cores = N_WORKERS, mc.preschedule = FALSE)
  bad <- vapply(res, inherits, logical(1), "try-error"); if (any(bad)) stop("worker error: ", paste(res[bad], collapse = " | "))
  rm(d_all); invisible(gc())
}

# ---- bands -------------------------------------------------------------------------------------------------------------
P <- rbindlist(lapply(all_files, readRDS))
nerr <- P[!is.na(error), uniqueN(paste(label, m, dist))]; logl("refits with INLA errors: %d (their imputation is left out of that panel)", nerr)
P <- P[is.na(error)]
nimp <- P[, .(n_imp = uniqueN(m)), by = .(label, dist, main)]
P[, band0 := findInterval(tsd_years, BREAKS, rightmost.closed = TRUE)]
merge_bands <- function(ev) {                                    # merge adjacent bands until each has >= MIN_EV events (left to right, last merged back)
  grp <- integer(length(ev)); g <- 1L; acc <- 0
  for (i in seq_along(ev)) { grp[i] <- g; acc <- acc + ev[i]; if (acc >= MIN_EV && i < length(ev)) { g <- g + 1L; acc <- 0 } }
  if (sum(ev[grp == max(grp)]) < MIN_EV && max(grp) > 1L) grp[grp == max(grp)] <- max(grp) - 1L
  grp
}
B0 <- P[, .(obs = sum(obs), exp = sum(exp), w_tsd = sum(exp * tsd_years), lo = min(tsd_years), hi = max(tsd_years)), by = .(label, dist, main, band0)]
B0 <- merge(B0, nimp, by = c("label", "dist", "main"))
B0[, `:=`(obs = obs / n_imp, exp = exp / n_imp, w_tsd = w_tsd / n_imp)]                      # averages per imputation
setorder(B0, label, dist, main, band0)
B0[, band := merge_bands(obs), by = .(label, dist, main)]
B <- B0[, .(obs = sum(obs), exp = sum(exp), tsd_mid = sum(w_tsd) / sum(exp), tsd_from = BREAKS[min(band0)], tsd_to = BREAKS[max(band0) + 1L],
            n_imp = n_imp[1]), by = .(label, dist, main, band)]
B[, `:=`(log_oe = log(pmax(obs, 0.5) / exp), se = 1 / sqrt(pmax(obs, 0.5)))][, `:=`(lo95 = log_oe - 1.96 * se, hi95 = log_oe + 1.96 * se)]
fwrite(B, file.path(OUT, "sim05b_bands.csv"))

# ---- curves from sim05 (TSD part only), centred on the disturbed rows with the no-TSD expected counts as weights -----------
S5 <- lapply(list.files(file.path(SIM05_DIR, "fits"), full.names = TRUE, pattern = "\\.rds$"), readRDS)
coef5 <- function(lb, tv, fm) { D <- do.call(rbind, lapply(Filter(function(r) r$label == lb, S5), function(r) r$fits[[fm]]$draws)); if (is.null(D)) NA_real_ else median(D[, tv]) }
CV <- rbindlist(lapply(split(P, by = c("label", "dist", "main")), function(x) {
  lb <- x$label[1]; nm <- x$dist[1]; tv <- DIST[[nm]]$tsd; tmax <- max(B[label == lb & dist == nm & main == x$main[1], tsd_mid], quantile(x$tsd_years, 0.95))
  grid <- seq(0, tmax, length.out = 81)
  rbindlist(lapply(c("linear", "log"), function(fm) {
    cc <- coef5(lb, tv, if (fm == "linear") "linear" else "log_all"); g <- function(t) if (fm == "linear") t / 10 else log1p(t)
    ctr <- log(sum(x$exp * exp(cc * g(x$tsd_years))) / sum(x$exp))
    data.table(label = lb, dist = nm, main = x$main[1], form = fm, tsd_years = grid, effect = cc * g(grid) - ctr, coef = cc)
  }))
}))
fwrite(CV, file.path(OUT, "sim05b_curves.csv"))

ok_plot <- tryCatch({
  suppressPackageStartupMessages(library(ggplot2))
  pn <- function(dt) dt[, panel := sprintf("%s | %s", label, fifelse(main %in% c("harv_light", "harv_heavy", "harv_any"), main, dist))]
  pn(B); pn(CV)
  B[, lab := sprintf("%.0f", obs)]
  p <- ggplot() + geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_line(data = CV, aes(tsd_years, effect, colour = form, linetype = form), linewidth = 0.8) +
    geom_errorbar(data = B, aes(x = tsd_mid, ymin = lo95, ymax = hi95), width = 0, colour = "#1a1a19", linewidth = 0.4) +
    geom_point(data = B, aes(tsd_mid, log_oe), colour = "#1a1a19", size = 1.8) +
    geom_text(data = B, aes(tsd_mid, hi95, label = lab), vjust = -0.5, size = 2.3, colour = "#5a5a57") +
    facet_wrap(~ panel, scales = "free", ncol = 4) +
    scale_colour_manual(values = c(linear = "#5a5a57", log = "#2a78d6")) + scale_linetype_manual(values = c(linear = "22", log = "solid")) +
    labs(x = "Time since disturbance (years); point at the expected-weighted mean TSD of its band",
         y = "log(observed / expected) vs the model WITHOUT this TSD term",
         title = "sim05b: what the data say about the TSD shape (points) vs the fitted linear and log curves (lines)",
         subtitle = sprintf("Points: disturbed rows grouped by TSD band (>= %g events each, averaged over %d imputations; number = events), 95%% Poisson interval. Lines: sim05 TSD part, centred on the same rows.", MIN_EV, length(REAL_M))) +
    theme_minimal(base_size = 9) + theme(legend.position = "top", strip.text = element_text(size = 7.5, hjust = 0), plot.subtitle = element_text(size = 7.5),
                                         plot.background = element_rect(fill = "white", colour = NA))
  ggsave(file.path(OUT, "sim05b_check.png"), p, width = 14, height = 2.6 * ceiling(uniqueN(B$panel) / 4) + 1, dpi = 130, limitsize = FALSE); TRUE
}, error = function(e) { logl("plot failed: %s", conditionMessage(e)); FALSE })
cat("\nbands\n"); print(B[, .(label, main, tsd_from, tsd_to, obs = round(obs, 1), exp = round(exp, 1), log_oe = round(log_oe, 2), lo95 = round(lo95, 2), hi95 = round(hi95, 2))])
logl("outputs: sim05b_bands.csv, sim05b_curves.csv%s", if (ok_plot) ", sim05b_check.png" else "")
logl("==== done")
