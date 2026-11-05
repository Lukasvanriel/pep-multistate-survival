# ==============================================================================
# R/tests/24_map_selection_regen.R  (2026-11-05) -- DIAGNOSTIC, no INLA needed. Reads the WP6 run, edits nothing.
# Question: is the *selection* of healthy INLA fits what drives the X->regen estimates? Same likelihood and prior, no INLA, no cap:
# for every imputation m with an INLA fit of black_spruce->regen, fir->regen, jack_pine->regen, trembling_aspen->regen (mi_raw + mi_raw_redraw,
# m1 to the latest candidate at the time of the run) the fitting data are rebuilt as the run did (08f rebuild: extras_<tag>_mXX.rds -> p1, p2 -> build_rows)
# and fit_map_v0() is fitted. Groups by the INLA fit of the same imputation:
#   healthy   = ratio >= 0.1 and z_map_alpha <= 3 and z_map_intercept <= 3     (ratio = sd_alpha / (0.78 alpha / sqrt(events)), as in 07)
#   collapsed = INLA fit without error, not healthy
#   error     = INLA error or NA in ratio / z (in 'all' only)
# Identity check of the rebuild (stops before pooling if it fails): events == mi_raw n_events, n_above_map == stored, data identical to baddata
# where a baddata file exists; maxH_map vs stored and z recomputed from the stored INLA estimates are printed (information).
# Pooling: Rubin (pool_rubin of 05a) over all / healthy / collapsed; alpha on the natural scale (delta-method SD) and on the log scale.
# Run from the project root:  PEP_DIAG_CORES=2 Rscript R/tests/24_map_selection_regen.R 2>&1 | tee data-output/diagnostics/24_map_selection.log
# Output: data-output/diagnostics/24_map_selection/ (snapshots of mi_raw, one rds per m, map_by_imputation.csv, map_pooled.csv). Resumable.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel) })
source(here::here("R", "config.R")); options(width = 220)
N_CORES <- as.integer(Sys.getenv("PEP_DIAG_CORES", "2"))
SIM_M <- 200L; M_MI <- 30L; SD_MIN <- 1e-3; N_COV_SAMPLES <- 0L
source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
LABS <- c("black_spruce->regen", "fir->regen", "jack_pine->regen", "trembling_aspen->regen"); RATIO_MIN <- 0.1; Z_HEALTHY <- 3
OUT <- here("data-output", "diagnostics", "24_map_selection"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
inp <- prepare_fit_inputs(TABLE_DIR)

# ---- snapshot of the INLA results (07 appends to the redraw file while it runs) -------------------------------------------------------------------
src <- file.path(FIT_DIR, sprintf(c("mi_raw_%s.csv", "mi_raw_redraw_%s.csv"), RUN_TAG)); snap <- file.path(OUT, paste0("snap_", basename(src)))
stopifnot(all(file.copy(src, snap, overwrite = TRUE)))
mr <- rbindlist(lapply(snap, fread), fill = TRUE)[label %in% LABS]
n_par <- vapply(LABS, function(lb) 2L + sum(inp$mask[match(lb, inp$td$label), ]), integer(1))
cnt <- mr[, .N, by = .(label, m)][, ok := N == unname(n_par[label])]
if (any(!cnt$ok)) { cat("Dropping incomplete blocks (file copied while 07 was writing?):\n"); print(cnt[ok == FALSE]); mr <- merge(mr, cnt[ok == TRUE, .(label, m)], by = c("label", "m")) }
a <- mr[parameter == "alpha", .(label, m, n_events, alpha_inla = est, sd_inla = sd, fit_error = fit_error %in% TRUE, z_a = z_map_alpha, z_i = z_map_intercept,
                                maxH_map_st = maxH_map, nab_map_st = n_above_map)]
i0 <- mr[parameter == "(Intercept)", .(label, m, icpt_inla = est)]; a <- merge(a, i0, by = c("label", "m"))
a[, ratio := sd_inla / (0.78 * alpha_inla / sqrt(n_events))]
a[, group := fifelse(fit_error | is.na(ratio) | is.na(z_a) | is.na(z_i), "error",
             fifelse(ratio >= RATIO_MIN & abs(z_a) <= Z_HEALTHY & abs(z_i) <= Z_HEALTHY, "healthy", "collapsed"))]
a[, reason := fifelse(group != "collapsed", group, fifelse(ratio < RATIO_MIN, "ratio<0.1", "z>3 only"))]
cat("Imputations per transition and INLA group (snapshot of", format(Sys.time()), "):\n"); print(dcast(a[, .N, by = .(label, reason)], label ~ reason, value.var = "N", fill = 0))
cat("m range per transition:\n"); print(a[, .(n = .N, m_min = min(m), m_max = max(m)), by = label])
jobs <- a[, .(label, m)]; stopifnot(all(file.exists(file.path(FIT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, unique(jobs$m))))))

# ---- rebuild (as 08f) and MAP -------------------------------------------------------------------------------------------------------------------
key_iv <- paste(inp$iv$id_pe, inp$iv$interval); b_abs <- inp$iv$t_end
rebuild <- function(m) {
  ex <- readRDS(file.path(FIT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, m))); p1 <- as.data.table(ex$p1); p2 <- as.data.table(ex$p2); rm(ex)
  ex_row <- rep(NA_real_, nrow(inp$iv)); ii <- match(paste(p2$id_pe, p2$interval), key_iv); stopifnot(!anyNA(ii)); ex_row[ii] <- p2$exit_abs
  p2f <- function(a_dur, b_dur, from, to, r) { E <- b_abs[r] - b_dur; min(max(ex_row[r] - E, a_dur), b_dur) }
  build_rows(inp$iv, setNames(p1$entry, as.character(p1$id_pe)), inp$states, inp$tr, p2_fun = p2f, keep_cols = inp$keep_cols)$rows[kept == TRUE]
}
done_keys <- function() { f <- list.files(OUT, "^m[0-9]+_.*\\.rds$", full.names = TRUE)
  if (!length(f)) return(character(0)); r <- rbindlist(lapply(f, readRDS), fill = TRUE); if ("label" %in% names(r)) unique(paste(r$label[!is.na(r$label)], r$m[!is.na(r$label)])) else character(0) }
dk <- done_keys(); todo <- jobs[!paste(label, m) %in% dk]; ms <- sort(unique(todo$m)); cat("Fits to do:", nrow(todo), "in", length(ms), "imputations (", length(dk), "done before )\n")
work <- function(mm) {
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch({
    rows <- rebuild(mm); labs_m <- todo$label[todo$m == mm]; bf <- file.path(FIT_DIR, sprintf("baddata_%s_m%02d.rds", RUN_TAG, mm)); bd <- if (file.exists(bf)) readRDS(bf) else NULL
    rbindlist(lapply(labs_m, function(lb) {
      d <- rows[paste(from, to, sep = "->") == lb]; k <- match(lb, inp$td$label); tms <- colnames(inp$mask)[inp$mask[k, ]]
      stopifnot("NA in the fitting columns" = !anyNA(d[, c("status", "Tstop", "truncation", tms), with = FALSE]))
      f <- fit_map_v0(d, tms); H <- H_at(d, tms, f$tab); cc <- c("status", "Tstop", "truncation", tms)
      ident <- if (!is.null(bd) && lb %in% names(bd)) isTRUE(all.equal(as.matrix(as.data.frame(bd[[lb]]$data)[, cc]), as.matrix(as.data.frame(d)[, cc]), check.attributes = FALSE, tolerance = 1e-12)) else NA
      cbind(data.table(label = lb, m = mm, rows = nrow(d), events = sum(d$status == 1L), conv = f$conv, maxH_map = H[["maxH"]], n_above_map = H[["n_above"]],
                       baddata_identical = ident, secs_m = NA_real_), f$tab) })) },
    error = function(e) data.table(label = NA_character_, m = mm, error = conditionMessage(e)))
  res[, secs_m := proc.time()[["elapsed"]] - t0]
  saveRDS(res, file.path(OUT, sprintf("m%03d_%s.rds", mm, format(Sys.time(), "%H%M%S")))); cat(sprintf("m=%d done (%.0f s)\n", mm, proc.time()[["elapsed"]] - t0)); NULL }
invisible(mclapply(ms, work, mc.cores = N_CORES, mc.preschedule = FALSE))

R <- rbindlist(lapply(list.files(OUT, "^m[0-9]+_.*\\.rds$", full.names = TRUE), readRDS), fill = TRUE)
if ("error" %in% names(R) && R[!is.na(error), .N]) { print(R[!is.na(error), .(m, error)]); stop("rebuild errors (see above)") }
R <- unique(R, by = c("label", "m", "parameter")); stopifnot(nrow(unique(R[, .(label, m)])) == nrow(jobs))
# ---- identity check -----------------------------------------------------------------------------------------------------------------------------
f1 <- unique(R[, .(label, m, rows, events, conv, maxH_map, n_above_map, baddata_identical)]); chk <- merge(f1, a, by = c("label", "m"))
za <- merge(R[parameter == "alpha", .(label, m, a_map = est, a_se = sd)], R[parameter == "(Intercept)", .(label, m, i_map = est, i_se = sd)], by = c("label", "m"))
chk <- merge(chk, za, by = c("label", "m")); chk[, `:=`(z_a_new = abs(alpha_inla - a_map) / a_se, z_i_new = abs(icpt_inla - i_map) / i_se)]
cat("\nIdentity check of the rebuild (all imputations):\n")
print(chk[, .(n = .N, events_equal = sum(events == n_events), n_above_equal = sum(n_above_map == nab_map_st, na.rm = TRUE),
              maxH_rel_diff_max = max(abs(maxH_map - maxH_map_st) / maxH_map_st, na.rm = TRUE), n_maxH_diff_gt_1e6 = sum(abs(maxH_map - maxH_map_st) / maxH_map_st > 1e-6, na.rm = TRUE),
              baddata_checked = sum(!is.na(baddata_identical)), baddata_identical = sum(baddata_identical, na.rm = TRUE),
              z_alpha_max_abs_diff = max(abs(z_a_new - z_a), na.rm = TRUE), z_icpt_max_abs_diff = max(abs(z_i_new - z_i), na.rm = TRUE), map_conv_problems = sum(conv != 0)), by = label])
bad <- chk[events != n_events | n_above_map != nab_map_st | baddata_identical %in% FALSE]
if (nrow(bad)) { print(bad[, .(label, m, events, n_events, n_above_map, nab_map_st, baddata_identical)]); stop("Identity check failed for the imputations above: pooling not done. Diagnose before using any MAP number.") }

# ---- pooling ------------------------------------------------------------------------------------------------------------------------------------
E <- merge(R[, .(label, m, parameter, est, sd)], a[, .(label, m, group, reason)], by = c("label", "m"))
E <- rbind(E, E[parameter == "alpha", .(label, m, parameter = "log_alpha", est = log(est), sd = sd / est, group, reason)])   # sd of alpha = alpha * se(log alpha)
porder <- unique(E$parameter); porder <- c("alpha", "log_alpha", setdiff(porder, c("alpha", "log_alpha")))
setsof <- list(all = c("healthy", "collapsed", "error"), healthy = "healthy", collapsed = "collapsed")
P <- rbindlist(lapply(names(setsof), function(s) E[group %in% setsof[[s]], { p <- pool_rubin(est, sd^2)
  .(n = p$m_eff, pooled = p$pooled_est, se = sqrt(p$total_var), sd_within = sqrt(p$W_bar), sd_between = sqrt(p$B), median = median(est), min = min(est), max = max(est)) }, by = .(label, parameter)][, set := s]))
P[, parameter := factor(parameter, porder)]; setorder(P, label, parameter, set)
fwrite(merge(R, a, by = c("label", "m")), file.path(OUT, "map_by_imputation.csv")); fwrite(P, file.path(OUT, "map_pooled.csv"))
cat("\n=========== MAP (fit_map_v0, no cap), Rubin-pooled; est = pooled mean, se = sqrt(total variance); sd_btw = SD across imputations ===========\n")
for (lb in LABS) { cat("\n#####", lb, "#####\n"); print(a[label == lb, .N, by = group])
  W <- dcast(P[label == lb], parameter ~ set, value.var = c("n", "pooled", "se", "sd_between"))
  setcolorder(W, c("parameter", paste0(rep(c("n_", "pooled_", "se_", "sd_between_"), each = 3), c("all", "healthy", "collapsed"))))
  W[, `:=`(diff_all_minus_healthy = pooled_all - pooled_healthy, in_se_healthy = (pooled_all - pooled_healthy) / se_healthy)]
  print(W, digits = 3, row.names = FALSE) }
cat("\nNote: 'all' includes error fits (MAP exists for them). Rubin SEs with n_healthy small are rough; sd_between is the SD of the MAP estimates across imputations.\n")
