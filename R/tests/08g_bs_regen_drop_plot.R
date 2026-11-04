# ==============================================================================
# R/tests/08g_bs_regen_drop_plot.R  (2026-10-27) -- QUICK SANITY CHECK, not part of the pipeline.
# Refit black_spruce->regen WITHOUT plot 7401600201 on the ORIGINAL imputations m = 1..30 (no redraw), same INLA call as 05
# (fit_inla_models_pc, pc.alphaw prior), and apply the frozen 08 screen. Compare with
#   (a) the main pooled result (redraw-selected valid fits) and (b) the same m = 1..30 fits WITH the plot (mi_raw_main.csv).
# Rows are rebuilt from the saved P1/P2 of each m (as 08f). Output: FIT_DIR/dropplot_bs_regen/ (one rds per m, resumable).
# Run from the project root:  Rscript R/tests/08g_bs_regen_drop_plot.R
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(INLA) })
source(here::here("R", "config.R")); options(width = 200)
INLA::inla.setOption(num.threads = "1:1")
N_CORES <- 2L; SIM_M <- 200L; M_MI <- 30L; SD_MIN <- 1e-3; N_COV_SAMPLES <- 0L
source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
DROP <- "7401600201"; LB <- "black_spruce->regen"; MS <- 1:30
inp <- prepare_fit_inputs(TABLE_DIR); K <- match(LB, inp$td$label)
OUT <- file.path(FIT_DIR, "dropplot_bs_regen"); dir.create(OUT, showWarnings = FALSE)

key_iv <- paste(inp$iv$id_pe, inp$iv$interval); b_abs <- inp$iv$t_end
rebuild <- function(m) {                                  # identical to 08f
  ex <- readRDS(file.path(FIT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, m))); p1 <- as.data.table(ex$p1); p2 <- as.data.table(ex$p2); rm(ex)
  ex_row <- rep(NA_real_, nrow(inp$iv)); ii <- match(paste(p2$id_pe, p2$interval), key_iv); stopifnot(!anyNA(ii)); ex_row[ii] <- p2$exit_abs
  p2f <- function(a_dur, b_dur, from, to, r) { E <- b_abs[r] - b_dur; min(max(ex_row[r] - E, a_dur), b_dur) }
  build_rows(inp$iv, setNames(p1$entry, as.character(p1$id_pe)), inp$states, inp$tr, p2_fun = p2f, keep_cols = inp$keep_cols)$rows
}

todo <- MS[!file.exists(file.path(OUT, sprintf("m%03d.rds", MS)))]
cat("Imputations to fit:", length(todo), "of", length(MS), "\n")
invisible(mclapply(todo, function(m) {
  res <- tryCatch({ rows <- rebuild(m)
    n_drop <- rows[kept == TRUE & paste(from, to, sep = "->") == LB & as.character(id_pe) == DROP, .N]
    rows <- rows[as.character(id_pe) != DROP]             # the plot leaves only this transition's data (only K is fitted)
    f <- fit_rows(rows, inp, fit_inla_models_pc, only_k = K)
    tidy_fit(f$models, inp, ks = f$ks)[, `:=`(m = m, n_events = f$n_events[1], rows_dropped = n_drop)] },
    error = function(e) data.table(k = K, m = m, fit_error = TRUE, error_msg = conditionMessage(e)))
  saveRDS(res, file.path(OUT, sprintf("m%03d.rds", m))); NULL }, mc.cores = N_CORES, mc.preschedule = FALSE))

raw <- rbindlist(lapply(file.path(OUT, sprintf("m%03d.rds", MS)), readRDS), fill = TRUE)
cat("Rows of plot", DROP, "removed from", LB, "per m:", paste(unique(raw$rows_dropped), collapse = ", "), "\n")

# ---- frozen 08 screen ---------------------------------------------------------------------------------------------------------------------------
screen <- function(r) { ref <- fread(file.path(FIT_DIR, sprintf("screen_ref_%s.csv", RUN_TAG)))
  a <- r[parameter == "alpha", .(k, m, alpha = est, a_sd = sd, n_events, fit_error)]; i <- r[parameter == "(Intercept)", .(k, m, icpt = est)]
  x <- merge(merge(a, i, by = c("k", "m"), all.x = TRUE), ref, by = "k", all.x = TRUE)
  x[, `:=`(ratio = a_sd / (0.78 * alpha / sqrt(n_events)), z_icpt = (icpt - icpt_ref) / icpt_sd_ref)]
  x[, reason := fifelse(fit_error %in% TRUE | is.na(alpha) | is.na(a_sd) | is.na(icpt), "error/NA",
                 fifelse(ratio < 0.1, "alpha_sd_ratio", fifelse(abs(z_icpt) > 4, "intercept", "ok")))][, .(m, reason)] }
sc_drop <- screen(raw)
with <- rbindlist(lapply(file.path(FIT_DIR, sprintf(c("mi_raw_%s.csv", "mi_raw_redraw_%s.csv"), RUN_TAG)), function(f) if (file.exists(f)) fread(f)), fill = TRUE)[k == K & m %in% MS]
sc_with <- screen(with)
cat("\nScreen on m = 1..30 --  WITH plot:"); print(table(sc_with$reason)); cat("WITHOUT plot:"); print(table(sc_drop$reason))

# ---- pooled estimates ---------------------------------------------------------------------------------------------------------------------------
pool <- function(r, ok_m) r[m %in% ok_m & fit_error == FALSE, { p <- pool_rubin(est, sd^2); .(est = p$pooled_est, se = sqrt(p$total_var), m_eff = p$m_eff,
                                                                  fmi_approx = (1 + 1 / p$m_eff) * p$B / p$total_var) }, by = parameter]
main <- fread(file.path(RESULTS_DIR, sprintf("pooled_%s_final.csv", RUN_TAG)))[label == LB, .(parameter, est_main = pooled_est, se_main = se, fmi_main = fmi)]
a <- pool(raw, sc_drop[reason == "ok", m])                # without plot, valid fits only
b <- pool(raw, MS)                                        # without plot, all 30 that ran
cmp <- Reduce(function(x, y) merge(x, y, by = "parameter", sort = FALSE), list(main,
  a[, .(parameter, est_drop_valid = est, se_drop_valid = se, m_valid = m_eff, fmi_drop = fmi_approx)],
  b[, .(parameter, est_drop_all30 = est)]))
cmp[, z_main_vs_drop := (est_drop_valid - est_main) / sqrt(se_main^2 + se_drop_valid^2)]
cat("\n", LB, ": main (with plot, redraw-selected) vs without plot", DROP, "(m = 1..30)\n", sep = "")
print(cmp[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)], row.names = FALSE)
fwrite(cmp, file.path(FIT_DIR, "dropplot_bs_regen_compare.csv"))
