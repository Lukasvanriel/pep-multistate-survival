# ==============================================================================
# R/tests/08f_failed_vs_valid.R  (2026-10-27) -- DIAGNOSTIC, no INLA needed.
# Question: are the imputed datasets whose INLA fit failed the screen (05/07_redraw) different from those whose fit was valid?
# If failing and valid datasets give the same penalised maximum-likelihood (MAP) estimates and have the same data features, the failures are
# INLA numerical problems on ordinary datasets, and replacing them by new imputations does not select the data.
# For every transition with at least one failed fit, every tried imputation m: the fitting rows are rebuilt from the saved P1 entries
# and P2 exits (extras_<tag>_mXX.rds; same build_rows call as 05; checked against the saved failing data; extras$p2 holds the final exits, rule or draw), then
#   (1) penalised maximum-likelihood fit with the FINAL model: fit_map_v0() from 05a (variant 0, N(0, COEF_PRIOR_SD^2) on the coefficients, flat intercept, no cap), terms from the term mask;
#   (2) data features: rows, events, max truncation time, max exit time, maxH_map = max cumulative hazard at the truncation time H(tr) under the dataset's own MAP fit
#       (H_at() from 05a; the failure mechanism found for the int run: one row with an extreme H(tr)), n_above_map (rows above the cap);
#   (3) (WP7, D6) the B6 per-fit columns of the INLA fit of that same dataset (mi_raw / mi_raw_redraw: maxH_inla, n_above_inla, maxH_map, n_above_map, z_map_max, z_map_alpha,
#       z_map_intercept; prefix mi_) are compared failed vs valid as well, and the 08f maxH_map is checked against the stored one.
# Estimates are on the INLA scale (decade scale intercept; no years conversion). There are no rows with NA covariates after A1: the script stops if there are.
# Output (data-output/<tag>/fits/failcheck_map/): one rds per m (resumable), failcheck_<tag>.csv (all), failcheck_summary_<tag>.csv.
# Run from the project root after 08:  Rscript R/tests/08f_failed_vs_valid.R   (re-running skips finished m)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel) })
source(here::here("R", "config.R")); options(width = 220)
N_CORES <- 2L
SIM_M <- 200L; M_MI <- 30L; SD_MIN <- 1e-3; N_COV_SAMPLES <- 0L
source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
inp <- prepare_fit_inputs(TABLE_DIR)
OUT <- file.path(FIT_DIR, "failcheck_map"); dir.create(OUT, showWarnings = FALSE)
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
fit_map <- function(d, tms) {                       # d: data.table of the rows of one transition; MAP fit and H(tr) from 05a (B6)
  stopifnot("NA in the fitting columns (none expected after A1)" = !anyNA(d[, c("status", "Tstop", "truncation", tms), with = FALSE]))
  f <- fit_map_v0(d, tms); H <- H_at(d, tms, f$tab)
  list(est = f$tab[, .(parameter, est, se = sd)], conv = f$conv, maxH_map = H[["maxH"]], n_above_map = H[["n_above"]]) }

ms <- sort(unique(jobs$m)); todo <- ms[!file.exists(file.path(OUT, sprintf("m%03d.rds", ms)))]
cat("Imputations to rebuild:", length(todo), "of", length(ms), "\n")
invisible(mclapply(todo, function(m) {
  res <- tryCatch({ rows <- rebuild(m)
    mm <- m; rbindlist(lapply(jobs[m == mm, label], function(lb) {
      d <- rows[paste(from, to, sep = "->") == lb]; k <- match(lb, inp$td$label); tms <- colnames(inp$mask)[inp$mask[k, ]]
      f <- fit_map(d, tms)
      cbind(data.table(label = lb, m = m, rows = nrow(d), events = sum(d$status == 1L), max_trunc = max(d$truncation), max_Tstop = max(d$Tstop),
                       maxH_map = f$maxH_map, n_above_map = f$n_above_map, conv = f$conv), f$est) })) },
    error = function(e) data.table(label = NA_character_, m = m, error = conditionMessage(e)))
  saveRDS(res, file.path(OUT, sprintf("m%03d.rds", m))); NULL }, mc.cores = N_CORES, mc.preschedule = FALSE))

R <- rbindlist(lapply(file.path(OUT, sprintf("m%03d.rds", ms)), readRDS), fill = TRUE)
if ("error" %in% names(R) && R[!is.na(error), .N]) { print(R[!is.na(error)]); stop("rebuild errors") }
R <- merge(R, jobs, by = c("label", "m"))
# check: the rebuilt data of failed fits equal the data saved by 05/07 (baddata) for the same m
chk <- rbindlist(lapply(R[valid == FALSE & reason != "error/NA", unique(m)], function(m) { f <- file.path(FIT_DIR, sprintf("baddata_%s_m%02d.rds", RUN_TAG, m))
  if (!file.exists(f)) return(NULL); b <- readRDS(f); rbindlist(lapply(names(b), function(lb) { if (!lb %in% LABS) return(NULL)
    data.table(label = lb, m = m, rows_saved = nrow(b[[lb]]$data), events_saved = sum(b[[lb]]$data$status == 1L), max_Tstop_saved = max(b[[lb]]$data$Tstop)) })) }))
if (nrow(chk)) { cc <- merge(chk, unique(R[, .(label, m, rows, events, max_Tstop)]), by = c("label", "m"))
  cat(sprintf("Rebuild check vs saved failing data: %d datasets, identical rows/events/max exit: %d\n", nrow(cc), cc[rows == rows_saved & events == events_saved & abs(max_Tstop - max_Tstop_saved) < 1e-9, .N])) }
fwrite(R, file.path(FIT_DIR, sprintf("failcheck_%s.csv", RUN_TAG)))

# ---- comparison failed vs valid -----------------------------------------------------------------------------------------------------------------
R[, group := fifelse(valid, "valid", "failed")]
z <- R[, { v <- est[group == "valid"]; f <- est[group == "failed"]
  .(n_valid = length(v), n_failed = length(f), mean_valid = mean(v), sd_valid = sd(v), mean_failed = mean(f),
    shift_in_sd = (mean(f) - mean(v)) / sd(v), median_abs_z_failed = median(abs((f - mean(v)) / sd(v))),
    pct_failed_outside_2sd = 100 * mean(abs((f - mean(v)) / sd(v)) > 2),
    p_wilcox = if (length(f) >= 2) suppressWarnings(wilcox.test(f, v)$p.value) else NA_real_,
    mean_map_se = mean(se, na.rm = TRUE)) }, by = .(label, parameter)]
# B6 columns of the INLA fit of the same dataset (one value per fit)
rf <- file.path(FIT_DIR, sprintf(c("mi_raw_%s.csv", "mi_raw_redraw_%s.csv"), RUN_TAG)); mr <- rbindlist(lapply(rf[file.exists(rf)], fread), fill = TRUE)
stopifnot("mi_raw lacks the B6 columns" = all(CAP_COLS %in% names(mr)))
b6 <- unique(mr[parameter == "alpha", c("label", "m", CAP_COLS), with = FALSE]); setnames(b6, CAP_COLS, paste0("mi_", CAP_COLS))
feat <- unique(R[, .(label, m, group, rows, events, max_trunc, max_Tstop, maxH_map, n_above_map, conv)]); feat <- merge(feat, b6, by = c("label", "m"), all.x = TRUE)
mcols <- c("max_trunc", "max_Tstop", "maxH_map", "n_above_map", paste0("mi_", CAP_COLS)); feat[, (mcols) := lapply(.SD, as.numeric), .SDcols = mcols]
cmp <- feat[!is.na(mi_maxH_map)]; cat(sprintf("\nB6 consistency: 08f maxH_map vs mi_raw maxH_map for %d datasets: max relative difference %.2e; n_above_map equal in %d\n", nrow(cmp), max(abs(cmp$maxH_map - cmp$mi_maxH_map) / cmp$mi_maxH_map), cmp[n_above_map == mi_n_above_map, .N]))
fz <- melt(feat, id.vars = c("label", "m", "group"), measure.vars = c("max_trunc", "max_Tstop", "maxH_map", "n_above_map", paste0("mi_", CAP_COLS)))[
  , .(median_valid = median(value[group == "valid"], na.rm = TRUE), median_failed = median(value[group == "failed"], na.rm = TRUE),
      p_wilcox = tryCatch(suppressWarnings(wilcox.test(value[group == "failed"], value[group == "valid"])$p.value), error = function(e) NA_real_)), by = .(label, variable)]
fwrite(z, file.path(FIT_DIR, sprintf("failcheck_summary_%s.csv", RUN_TAG)))
cat("\nMAP convergence problems:", feat[conv != 0, .N], "of", nrow(feat), "datasets\n")
cat("\nData features, failed vs valid (median; Wilcoxon p):\n"); print(fz, digits = 3)
cat("\nMAP estimates, failed vs valid datasets (shift of the failed mean in SDs of the valid MAP estimates; % of failed outside +/- 2 SD):\n")
print(z[, .(label, parameter, n_valid, n_failed, shift_in_sd = round(shift_in_sd, 2), median_abs_z_failed = round(median_abs_z_failed, 2),
            pct_failed_outside_2sd = round(pct_failed_outside_2sd), p_wilcox = signif(p_wilcox, 2))], row.names = FALSE)
cat(sprintf("\nOverall: %d parameter comparisons; |shift| > 0.5 SD: %d; Wilcoxon p < 0.05: %d (expected by chance ~%.0f)\n",
            nrow(z), z[abs(shift_in_sd) > 0.5, .N], z[p_wilcox < 0.05, .N], 0.05 * z[!is.na(p_wilcox), .N]))
