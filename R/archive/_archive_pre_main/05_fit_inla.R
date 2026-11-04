# ==============================================================================
# 05_fit_inla.R   RUN ON THE MACHINE WITH INLA
#
# Part B of the pipeline (real PEP data), in the order of fit_p1p2_combined.R:
#   Stage 1  naive rows (entry W_j/2, midpoint exits)        -> fit -> theta_stage1
#   Stage 2  P1 weighted-mean entry + moment-matched P2 exit
#            (both driven by theta_stage1)                   -> fit -> seed models, theta_seed
#   MI loop  M_MI times: draw theta_m from the seed fits, P1 draw, chained P2 -> fit
#   Pool     Rubin's rules per transition and term
# No bias correction, no truth, no coverage: the pooled table is the result.
#
# Needs: 04_prepare_inla.R run (data/processed/*.rds); covariate_spec.R, 04b_build_rows.R,
#        05a_fit_funs.R in this folder. The simulation folder is only SOURCED (read-only) for
#        model_fitting.R (fit_inla_models).
# Each stage is checkpointed in data/processed/fit/; delete a file to redo that stage.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(INLA) })

# ---- config ------------------------------------------------------------------
SIM_DIR       <- "/Users/lukas/Desktop/RProjects/Quebec-forest-succession/R-scripts/Sim4stateExact/p1p2_v2_test"
SUBSAMPLE_N   <- NULL          # e.g. 500 for a smoke test (fits will be sparse); NULL = all plots
SEED          <- 2026L
N_WORKERS     <- 6L            # MI loop workers; Stage 2 P1 uses the same number
M_MI          <- 30L           # simulation value
SIM_M         <- 200L          # P1 candidate paths per plot (simulation value)
RUN_TAG       <- "gs"          # 2026-10-22: growing-season CMI run ("int" = annual-CMI run, archived in fit_int/)
SD_MIN <- 1e-3                 # alpha SD below this = near-degenerate fit: screened, its fitting data are saved
N_COV_SAMPLES <- 1000L         # posterior samples per fit for the joint covariance (0 = skip)
SAVE_MODELS_FOR <- 1L          # full INLA objects saved for these imputations only (size check; ~200 MB each)
OUT_DIR       <- here("data", "processed", if (is.null(SUBSAMPLE_N)) paste0("fit_", RUN_TAG) else paste0("fit_", RUN_TAG, "_subsample"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(SIM_DIR, "model_fitting.R"))      # fit_inla_models(), read-only source
source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED)
fit_fun <- fit_inla_models

inp <- prepare_fit_inputs(here("data", "processed"), subsample_n = SUBSAMPLE_N, seed = SEED)
saveRDS(inp$scaling, file.path(OUT_DIR, "scaling.rds"))     # needed to convert z-scaled betas back to raw units
cat(sprintf("Plots %d | intervals %d | kept transitions %d | continuous: %s\n", nrow(inp$p1), nrow(inp$iv), nrow(inp$td), paste(COV_CONT, collapse = ", ")))
report_fits <- function(tidy, theta, n_events, tag) {
  cat(sprintf("%s: fits failed %d/%d | theta not ok %d | alpha range %.2f-%.2f | events per fit min %d median %d\n", tag,
              tidy[parameter == "alpha" & fit_error == TRUE, .N], nrow(inp$td), sum(!theta$ok),
              min(theta$alpha, na.rm = TRUE), max(theta$alpha, na.rm = TRUE), min(n_events), as.integer(median(n_events))))
}

# ---- Stage 1 -----------------------------------------------------------------
f1 <- file.path(OUT_DIR, "stage1.rds")
if (file.exists(f1)) { s1 <- readRDS(f1); cat("Stage 1: loaded checkpoint\n") } else {
  t0 <- proc.time()["elapsed"]; s1 <- run_stage1(inp, fit_fun); saveRDS(s1, f1)
  cat(sprintf("Stage 1 done (%.0f s), floored segments %d\n", proc.time()["elapsed"] - t0, s1$n_floored))
}
report_fits(s1$tidy, s1$theta, s1$n_events, "Stage 1")
if (any(!s1$theta$ok)) { cat("  transitions not ok (hazard 0 in P1, midpoint in P2):\n"); print(inp$td$label[!s1$theta$ok]) }

# ---- Stage 2 -----------------------------------------------------------------
f2 <- file.path(OUT_DIR, "stage2.rds")
if (file.exists(f2)) { s2 <- readRDS(f2); cat("Stage 2: loaded checkpoint\n") } else {
  t0 <- proc.time()["elapsed"]; s2 <- run_stage2(inp, s1$theta, fit_fun, M = SIM_M, n_cores = N_WORKERS); saveRDS(s2, f2)
  cat(sprintf("Stage 2 done (%.0f s)\n", proc.time()["elapsed"] - t0))
}
report_fits(s2$tidy, s2$theta, s2$n_events, "Stage 2")
cat(sprintf("  P1: %d of %d plots without a usable candidate -> entry W_j/2 (marked)\n", sum(s2$p1$fallback), nrow(s2$p1)))
cat(sprintf("  P2: %d exits, %d by midpoint (dropped / not-ok transition)\n", s2$p2_event, s2$p2_midpoint))
fb <- cbind(inp$p1[, .(id_pe, W_j, first_state)], s2$p1[, .(n_cand, fallback)])[fallback == TRUE]
fwrite(fb, file.path(OUT_DIR, "p1_fallback_plots_stage2.csv")); print(fb[, .N, by = first_state][order(-N)])

# first time only: look at one posterior sample's structure before trusting the name matching
k_ok <- which(vapply(s2$models, function(m) is.null(m$error), logical(1)))[1]
smp <- INLA::inla.posterior.sample(1, s2$models[[k_ok]])[[1]]
cat("posterior sample structure (", inp$td$label[k_ok], "):\n", sep = ""); print(rownames(smp$latent)[1:min(6, nrow(smp$latent))]); print(names(smp$hyperpar))

# run metadata once: software versions and the exact INLA settings (priors etc.) of a seed fit
meta_path <- file.path(OUT_DIR, sprintf("run_meta_%s.rds", RUN_TAG))
if (!file.exists(meta_path)) saveRDS(list(time = Sys.time(), SEED = SEED, M_MI = M_MI, SIM_M = SIM_M, N_COV_SAMPLES = N_COV_SAMPLES,
  COV_CONT = COV_CONT, inla_version = as.character(utils::packageVersion("INLA")), session = utils::sessionInfo(),
  inla_args = s2$models[[k_ok]]$.args[c("family", "control.fixed", "control.family", "control.inla", "control.compute", "control.mode")],
  fit_inla_models_src = deparse(fit_inla_models), rhs = inp$rhs, td = inp$td), meta_path)

# ---- MI loop (resumable: finished m are skipped) -------------------------------
raw_path  <- file.path(OUT_DIR, sprintf("mi_raw_%s.csv", RUN_TAG)); diag_path <- file.path(OUT_DIR, sprintf("mi_diag_%s.csv", RUN_TAG)); err_path <- file.path(OUT_DIR, sprintf("mi_errors_%s.log", RUN_TAG))
done <- if (file.exists(diag_path)) fread(diag_path)$m else integer(0)
todo <- setdiff(seq_len(M_MI), done)
cat(sprintf("MI loop: %d done, %d to do (%d workers)\n", length(done), length(todo), N_WORKERS))
t0 <- proc.time()["elapsed"]
for (batch in split(todo, ceiling(seq_along(todo) / N_WORKERS))) {
  out <- mclapply(batch, function(m) tryCatch({
    set.seed(SEED * 1000L + m)
    run_one_imputation(m, inp, s2$models, s2$theta, fit_fun, M = SIM_M, keep_models = m %in% SAVE_MODELS_FOR)
  }, error = function(e) list(error = conditionMessage(e), m = m)), mc.cores = N_WORKERS)
  for (r in out) {
    if (is.null(r$res)) { cat(sprintf("[%s] m=%s FAILED: %s\n", format(Sys.time()), r$m, if (is.null(r$error)) "worker returned nothing" else r$error), file = err_path, append = TRUE); next }
    if (!is.null(r$models)) saveRDS(r$models, file.path(OUT_DIR, sprintf("models_%s_m%02d.rds", RUN_TAG, r$diag$m)))
    if (!is.null(r$bad_data)) saveRDS(r$bad_data, file.path(OUT_DIR, sprintf("baddata_%s_m%02d.rds", RUN_TAG, r$diag$m)))
    saveRDS(r$extras, file.path(OUT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, r$diag$m)))
    fwrite(r$res,  raw_path,  append = file.exists(raw_path))
    fwrite(r$diag, diag_path, append = file.exists(diag_path))
  }
  cat(sprintf("  batch up to m=%d done (%.0f s)\n", max(batch), proc.time()["elapsed"] - t0))
}

# ---- Pool --------------------------------------------------------------------
raw <- fread(raw_path); dg <- fread(diag_path)
cat(sprintf("\nImputations in the pool: %d of %d\n", uniqueN(raw$m), M_MI))
cat("Per-imputation diagnostics (P1 fallbacks, P2 midpoints, failed fits):\n"); print(dg)
pooled <- pool_all(raw)                                   # regular Rubin pooling: every fit that ran
fwrite(pooled, file.path(OUT_DIR, sprintf("pooled_%s.csv", RUN_TAG)))
fl <- flag_fits(raw); fwrite(fl[flagged == TRUE], file.path(OUT_DIR, sprintf("flagged_fits_%s.csv", RUN_TAG)))
pooled_x <- pool_all(raw, exclude_flagged = TRUE)         # sensitivity: flagged fits (error or alpha SD < SD_MIN) removed
fwrite(pooled_x, file.path(OUT_DIR, sprintf("pooled_%s_excl.csv", RUN_TAG)))
cat(sprintf("\nFlagged fits: %d of %d (errors %d, near-degenerate %d) -> flagged_fits_%s.csv\n", fl[flagged == TRUE, .N], nrow(fl),
            fl[, sum(is.na(alpha_sd))], fl[!is.na(alpha_sd) & alpha_sd < SD_MIN, .N], RUN_TAG))
cat(sprintf("Pooled (all) %d rows and (excluding flagged) %d rows\n", nrow(pooled), nrow(pooled_x)))
cat("Parameters with m_eff < M_MI (some imputations failed or were NA):", pooled[m_eff < uniqueN(raw$m), .N], "\n")
