# ==============================================================================
# 06v_p1_tail_check.R   RUN ON THE MACHINE WITH INLA (no fitting, ~ 5-10 min)
# Is the black_spruce->regen collapse tied to the P1 entry draws (heavy tails / many long-truncated rows)?
# For ALL 30 imputations rebuild the transition's rows from the saved P1/P2 draws (variant A of 06j, no refit) and compute per imputation:
#   P1 entry draws (plots in this transition): quantiles, max, share beyond Stage-1 maximum, entry / (W_j/2)
#   rows: truncation quantiles/max, rows and plots beyond 5 / 8 / 10 decades, sum(truncation), sum(truncation^1.8) (burden), exposure
# Then compare the 8 collapsed imputations (m = 2, 8, 14, 15, 19, 23, 25, 29) with the 22 healthy ones: group medians and AUC per feature.
# Output: data/processed/fit/p1_tail_check_<tag>.csv  (+ _auc_<tag>.csv)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here) })
SIM_DIR <- "/Users/lukas/Desktop/RProjects/Quebec-forest-succession/R-scripts/Sim4stateExact/p1p2_v2_test"
RUN_TAG <- "diag"; OUT_DIR <- here("data", "processed", "fit"); SEED <- 2026L; options(width = 200)
FLAG <- c(2, 8, 14, 15, 19, 23, 25, 29); LABEL <- "black_spruce->regen"
suppressPackageStartupMessages(library(INLA))
source(file.path(SIM_DIR, "model_fitting.R"))
source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
inp <- prepare_fit_inputs(here("data", "processed"), subsample_n = NULL, seed = SEED)
k <- match(LABEL, inp$td$label); stopifnot(!is.na(k)); fr <- strsplit(LABEL, "->", fixed = TRUE)[[1]]
iv <- inp$iv; iv_key <- paste(iv$id_pe, iv$interval); W_half <- setNames(inp$p1$W_j / 2, inp$p1$id_pe)
make_p2_saved <- function(p2, first_entry) {
  ex <- setNames(p2$exit_abs, paste(p2$id_pe, p2$interval)); last <- NA; E <- NA
  function(a_dur, b_dur, from, to, r) {
    id <- iv$id_pe[r]; if (!identical(id, last)) { last <<- id; E <<- unname(first_entry[as.character(id)]) }
    Tc <- min(max(ex[[iv_key[r]]] - E, a_dur), b_dur); E <<- E + Tc; Tc } }
Q <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE)); out <- list()
for (mm in 1:30) {
  f <- file.path(OUT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, mm)); if (!file.exists(f)) next
  ex <- readRDS(f); fe <- setNames(ex$p1$entry, ex$p1$id_pe)
  b <- build_rows(iv, fe, inp$states, inp$tr, p2_fun = make_p2_saved(ex$p2, fe), keep_cols = inp$keep_cols)
  nm <- names(b$rows); idc <- intersect(c("id_pe", "id"), nm)[1]
  d <- b$rows[kept == TRUE & from == fr[1] & to == fr[2]]
  tr <- d$truncation; ts <- d$Tstop; ex_t <- ts - tr
  ids <- if (!is.na(idc)) unique(d[[idc]]) else NULL
  e1 <- if (!is.null(ids)) fe[as.character(ids)] else NA_real_; w1 <- if (!is.null(ids)) W_half[as.character(ids)] else NA_real_
  out[[length(out) + 1]] <- data.table(m = mm, flagged = mm %in% FLAG, n_rows = nrow(d), events = sum(d$status == 1L), n_plots = length(ids),
    entry_q50 = Q(e1, .5), entry_q90 = Q(e1, .9), entry_q99 = Q(e1, .99), entry_max = max(e1, na.rm = TRUE), entry_over_Wj2_q90 = Q(e1 / w1, .9), entry_over_Wj2_max = max(e1 / w1, na.rm = TRUE),
    plots_entry_gt8 = sum(e1 > 8, na.rm = TRUE), plots_entry_gt10 = sum(e1 > 10, na.rm = TRUE),
    trunc_q50 = Q(tr, .5), trunc_q90 = Q(tr, .9), trunc_q99 = Q(tr, .99), trunc_max = max(tr), rows_trunc_gt5 = sum(tr > 5), rows_trunc_gt8 = sum(tr > 8), rows_trunc_gt10 = sum(tr > 10),
    sum_trunc = sum(tr), burden_18 = sum(tr^1.8), exposure_sum = sum(ex_t), exposure_min = min(ex_t), rows_exposure_lt_0.01 = sum(ex_t < 0.01))
  cat("m", mm, "done\n") }
R <- rbindlist(out); fwrite(R, file.path(OUT_DIR, sprintf("p1_tail_check_%s.csv", RUN_TAG)))
feat <- setdiff(names(R), c("m", "flagged"))
auc <- rbindlist(lapply(feat, function(v) { x <- R[[v]]; g <- R$flagged; r <- rank(x); n1 <- sum(g); n0 <- sum(!g)
  data.table(feature = v, median_collapsed = round(median(x[g]), 3), median_healthy = round(median(x[!g]), 3), min_collapsed = round(min(x[g]), 3), max_healthy = round(max(x[!g]), 3),
             AUC = round((sum(r[g]) - n1 * (n1 + 1) / 2) / (n1 * n0), 3)) }))
auc[, sep := abs(AUC - 0.5)]; setorder(auc, -sep); fwrite(auc, file.path(OUT_DIR, sprintf("p1_tail_check_auc_%s.csv", RUN_TAG)))
cat("\n==== features separating collapsed (n=8) from healthy (n=22); AUC 1 or 0 = perfect separation, 0.5 = none ====\n"); print(auc[, !"sep"])
cat("\nper imputation:\n"); print(R[, .(m, flagged, entry_max, plots_entry_gt10, rows_trunc_gt10, trunc_max, burden_18 = round(burden_18))])
