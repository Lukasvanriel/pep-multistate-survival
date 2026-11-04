# ==============================================================================
# R/tests/08f_failed_vs_valid.R  (2026-10-27) -- DIAGNOSTIC, no INLA needed.
# Question: are the imputed datasets whose INLA fit failed the screen (05/07_redraw) different from those whose fit was valid?
# If failing and valid datasets give the same maximum-likelihood (MLE) estimates and have the same data features, the failures are
# INLA numerical problems on ordinary datasets, and replacing them by new imputations does not select the data.
# For every transition with at least one failed fit, every tried imputation m: the fitting rows are rebuilt from the saved P1 entries
# and P2 exits (extras_<tag>_mXX.rds; same build_rows call as 05; checked against the saved failing data), then
#   (1) MLE fit (Weibull with left truncation, the same parameterisation as INLA: S(t) = exp(-(lambda t)^alpha), lambda = exp(eta));
#   (2) data features: rows, events, max truncation time, max exit time, max cumulative hazard at the truncation time H(tr) under the
#       dataset's own MLE (the failure mechanism found for the int run: one row with an extreme H(tr)).
# Output (data-output/<tag>/fits/failcheck_mle/): one rds per m (resumable), failcheck_<tag>.csv (all), failcheck_summary_<tag>.csv.
# NOTE: rows with NA covariates (climate NA for 9 intervals in 5 plots; INLA treats an NA covariate as 0) are left out of the MLE;
# the rebuild check therefore compares row counts including them.
# Run from the project root after 08:  Rscript R/tests/08f_failed_vs_valid.R   (re-running skips finished m)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel) })
source(here::here("R", "config.R")); options(width = 220)
N_CORES <- 2L
SIM_M <- 200L; M_MI <- 30L; SD_MIN <- 1e-3; N_COV_SAMPLES <- 0L
source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
inp <- prepare_fit_inputs(TABLE_DIR)
OUT <- file.path(FIT_DIR, "failcheck_mle2"); dir.create(OUT, showWarnings = FALSE)
fv <- fread(file.path(RESULTS_DIR, sprintf("fit_validity_%s.csv", RUN_TAG)))
LABS <- fv[valid == FALSE, unique(label)]
jobs <- fv[label %in% LABS, .(label, m, valid, reason)]
cat("Transitions with failed fits:", paste(LABS, collapse = ", "), "| datasets:", nrow(jobs), "\n")

key_iv <- paste(inp$iv$id_pe, inp$iv$interval); b_abs <- inp$iv$t_end
rebuild <- function(m) {
  ex <- readRDS(file.path(FIT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, m))); p1 <- as.data.table(ex$p1); p2 <- as.data.table(ex$p2); rm(ex)
  ex_row <- rep(NA_real_, nrow(inp$iv)); ii <- match(paste(p2$id_pe, p2$interval), key_iv); stopifnot(!anyNA(ii)); ex_row[ii] <- p2$exit_abs
  p2f <- function(a_dur, b_dur, from, to, r) { E <- b_abs[r] - b_dur; min(max(ex_row[r] - E, a_dur), b_dur) }
  build_rows(inp$iv, setNames(p1$entry, as.character(p1$id_pe)), inp$states, inp$tr, p2_fun = p2f, keep_cols = inp$keep_cols)$rows[kept == TRUE]
}
fit_mle <- function(d, tms) {                       # as R/tests/20_alpha_prior_test2.R
  X <- cbind(1, as.matrix(d[, tms, drop = FALSE])); T <- d$Tstop; tr <- d$truncation; s <- d$status; lT <- log(T); ltr <- ifelse(tr > 0, log(pmax(tr, 1e-300)), 0)
  nll <- function(th) { a <- exp(th[1]); eta <- drop(X %*% th[-1]); -(sum(s * (th[1] + a * eta + (a - 1) * lT)) - sum(exp(a * eta) * (T^a - tr^a))) }
  gr <- function(th) { a <- exp(th[1]); eta <- drop(X %*% th[-1]); e <- exp(a * eta); Ta <- T^a; tra <- tr^a
    gb <- -(colSums(X * (s * a)) - colSums(X * (a * e * (Ta - tra)))); ga <- -(sum(s * (1 / a + eta + lT)) - sum(e * (eta * (Ta - tra) + Ta * lT - tra * ltr))); c(ga * a, gb) }
  o <- optim(c(0, log(sum(s) / sum(T - tr)), rep(0, length(tms))), nll, gr, method = "BFGS", control = list(maxit = 5000, reltol = 1e-12))
  V <- tryCatch(solve(optimHess(o$par, nll, gr)), error = function(e) NULL); se <- if (is.null(V)) rep(NA_real_, length(o$par)) else sqrt(diag(V))
  a <- exp(o$par[1]); eta <- drop(X %*% o$par[-1])
  list(est = data.table(parameter = c("alpha", "(Intercept)", tms), est = c(a, o$par[2] - log(TIME_SCALE), o$par[-(1:2)]), se = c(a * se[1], se[-1])),
       conv = o$convergence, maxH_tr = max(exp(a * eta) * d$truncation^a)) }

ms <- sort(unique(jobs$m)); todo <- ms[!file.exists(file.path(OUT, sprintf("m%03d.rds", ms)))]
cat("Imputations to rebuild:", length(todo), "of", length(ms), "\n")
invisible(mclapply(todo, function(m) {
  res <- tryCatch({ rows <- rebuild(m)
    mm <- m; rbindlist(lapply(jobs[m == mm, label], function(lb) {
      d <- as.data.frame(rows[paste(from, to, sep = "->") == lb]); k <- match(lb, inp$td$label); tms <- colnames(inp$mask)[inp$mask[k, ]]
      ok <- complete.cases(d[, c("status", "Tstop", "truncation", tms)])   # rows with NA climate (9 intervals in 5 plots) are left out of the MLE
      n_na <- sum(!ok); d <- d[ok, , drop = FALSE]
      f <- fit_mle(d, tms)
      cbind(data.table(label = lb, m = m, rows = nrow(d), rows_na_dropped = n_na, events = sum(d$status == 1L), max_trunc = max(d$truncation), max_Tstop = max(d$Tstop),
                       maxH_tr = f$maxH_tr, conv = f$conv), f$est) })) },
    error = function(e) data.table(label = NA_character_, m = m, error = conditionMessage(e)))
  saveRDS(res, file.path(OUT, sprintf("m%03d.rds", m))); NULL }, mc.cores = N_CORES, mc.preschedule = FALSE))

R <- rbindlist(lapply(file.path(OUT, sprintf("m%03d.rds", ms)), readRDS), fill = TRUE)
if ("error" %in% names(R) && R[!is.na(error), .N]) { print(R[!is.na(error)]); stop("rebuild errors") }
R <- merge(R, jobs, by = c("label", "m"))
# check: the rebuilt data of failed fits equal the data saved by 05/07 (baddata) for the same m
chk <- rbindlist(lapply(R[valid == FALSE & reason != "error/NA", unique(m)], function(m) { f <- file.path(FIT_DIR, sprintf("baddata_%s_m%02d.rds", RUN_TAG, m))
  if (!file.exists(f)) return(NULL); b <- readRDS(f); rbindlist(lapply(names(b), function(lb) { if (!lb %in% LABS) return(NULL)
    data.table(label = lb, m = m, rows_saved = nrow(b[[lb]]$data), events_saved = sum(b[[lb]]$data$status == 1L), max_Tstop_saved = max(b[[lb]]$data$Tstop)) })) }))
if (nrow(chk)) { cc <- merge(chk, unique(R[, .(label, m, rows = rows + rows_na_dropped, events, max_Tstop)]), by = c("label", "m"))
  cat(sprintf("Rebuild check vs saved failing data: %d datasets, identical rows/events/max exit: %d\n", nrow(cc), cc[rows == rows_saved & events == events_saved & abs(max_Tstop - max_Tstop_saved) < 1e-9, .N])) }
fwrite(R, file.path(FIT_DIR, sprintf("failcheck_%s.csv", RUN_TAG)))

# ---- comparison failed vs valid -----------------------------------------------------------------------------------------------------------------
R[, group := fifelse(valid, "valid", "failed")]
z <- R[, { v <- est[group == "valid"]; f <- est[group == "failed"]
  .(n_valid = length(v), n_failed = length(f), mean_valid = mean(v), sd_valid = sd(v), mean_failed = mean(f),
    shift_in_sd = (mean(f) - mean(v)) / sd(v), median_abs_z_failed = median(abs((f - mean(v)) / sd(v))),
    pct_failed_outside_2sd = 100 * mean(abs((f - mean(v)) / sd(v)) > 2),
    p_wilcox = if (length(f) >= 2) suppressWarnings(wilcox.test(f, v)$p.value) else NA_real_,
    mean_mle_se = mean(se, na.rm = TRUE)) }, by = .(label, parameter)]
feat <- unique(R[, .(label, m, group, rows, rows_na_dropped, events, max_trunc, max_Tstop, maxH_tr, conv)])
fz <- melt(feat, id.vars = c("label", "m", "group"), measure.vars = c("max_trunc", "max_Tstop", "maxH_tr"))[
  , .(median_valid = median(value[group == "valid"]), median_failed = median(value[group == "failed"]),
      p_wilcox = suppressWarnings(wilcox.test(value[group == "failed"], value[group == "valid"])$p.value)), by = .(label, variable)]
fwrite(z, file.path(FIT_DIR, sprintf("failcheck_summary_%s.csv", RUN_TAG)))
cat("\nMLE convergence problems:", feat[conv != 0, .N], "of", nrow(feat), "datasets\n")
cat("\nData features, failed vs valid (median; Wilcoxon p):\n"); print(fz, digits = 3)
cat("\nMLE estimates, failed vs valid datasets (shift of the failed mean in SDs of the valid MLEs; % of failed outside +/- 2 SD):\n")
print(z[, .(label, parameter, n_valid, n_failed, shift_in_sd = round(shift_in_sd, 2), median_abs_z_failed = round(median_abs_z_failed, 2),
            pct_failed_outside_2sd = round(pct_failed_outside_2sd), p_wilcox = signif(p_wilcox, 2))], row.names = FALSE)
cat(sprintf("\nOverall: %d parameter comparisons; |shift| > 0.5 SD: %d; Wilcoxon p < 0.05: %d (expected by chance ~%.0f)\n",
            nrow(z), z[abs(shift_in_sd) > 0.5, .N], z[p_wilcox < 0.05, .N], 0.05 * z[!is.na(p_wilcox), .N]))
