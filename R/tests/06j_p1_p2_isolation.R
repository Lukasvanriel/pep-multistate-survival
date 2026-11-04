# ==============================================================================
# 06j_p1_p2_isolation.R   RUN ON THE MACHINE WITH INLA (~ 1 h; shorten M_SET for a quick look)
# Is the black_spruce->regen collapse caused by the P1 entries, the P2 exits, or neither?
# For imputations in M_SET (8 flagged + 4 unflagged) the black_spruce->regen data are REBUILT from the saved
# P1 / P2 draws (extras_diag_mXX.rds) and refitted (free alpha, pipeline settings):
#   A  P1 draw  + saved P2 draw   (control: must reproduce the original fit)
#   B  P1 draw  + midpoint exits
#   C  W_j/2    + saved P2 draw   (exit TIMES kept, clamped into the bracket)
#   D  W_j/2    + midpoint exits  (= Stage-1 construction)
# Output: data/processed/fit/p1p2_isolation_<tag>.csv   (ratio = alpha_sd / (0.78 * alpha / sqrt(events)); collapse = ratio < 0.1)
# Nothing in the pipeline or the simulation folder is modified.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(INLA) })
SIM_DIR <- "/Users/lukas/Desktop/RProjects/Quebec-forest-succession/R-scripts/Sim4stateExact/p1p2_v2_test"
RUN_TAG <- "diag"; OUT_DIR <- here("data", "processed", "fit"); SEED <- 2026L; options(width = 200)
M_SET <- c(2, 8, 14, 15, 19, 23, 25, 29,  1, 3, 4, 5)       # first 8 flagged, last 4 unflagged
LABEL <- "black_spruce->regen"
ALPHA_NAME <- "alpha parameter for weibullsurv"
source(file.path(SIM_DIR, "model_fitting.R"))
source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
fit_fun <- fit_inla_models
inp <- prepare_fit_inputs(here("data", "processed"), subsample_n = NULL, seed = SEED)
k <- match(LABEL, inp$td$label); stopifnot(!is.na(k))
fr <- strsplit(LABEL, "->", fixed = TRUE)[[1]]
cols <- c("from", "to", "status", "Tstop", "truncation", inp$terms)
iv <- inp$iv; iv_key <- paste(iv$id_pe, iv$interval)
W_half <- setNames(inp$p1$W_j / 2, inp$p1$id_pe)

# P2 from saved absolute exit times: Td = exit_abs - current entry, clamped into the bracket (counts clamps)
make_p2_saved <- function(p2, first_entry, cnt) {
  ex <- setNames(p2$exit_abs, paste(p2$id_pe, p2$interval)); last <- NA; E <- NA
  function(a_dur, b_dur, from, to, r) {
    id <- iv$id_pe[r]
    if (!identical(id, last)) { last <<- id; E <<- unname(first_entry[as.character(id)]) }
    Td <- ex[[iv_key[r]]] - E
    Tc <- min(max(Td, a_dur), b_dur)
    if (abs(Tc - Td) > 1e-6) cnt$clamped <- cnt$clamped + 1L
    E <<- E + Tc; Tc
  }
}
fit_variant <- function(fe, p2f) {
  b <- build_rows(iv, fe, inp$states, inp$tr, p2_fun = p2f, keep_cols = inp$keep_cols)
  d <- as.data.frame(b$rows[kept == TRUE & from == fr[1] & to == fr[2], ..cols])
  m <- fit_fun(list(d), rhs_formula = inp$rhs[k], exact = TRUE)[[1]]
  ev <- sum(d$status == 1L)
  if (!is.null(m$error)) return(data.table(alpha = NA_real_, sd = NA_real_, ratio = NA_real_, events = ev, mlik = NA_real_, error = substr(as.character(m$error)[1], 1, 80)))
  a <- m$summary.hyperpar[ALPHA_NAME, "mean"]; s <- m$summary.hyperpar[ALPHA_NAME, "sd"]
  data.table(alpha = a, sd = s, ratio = s / (0.78 * a / sqrt(ev)), events = ev, mlik = as.numeric(m$mlik[1, 1]), error = NA_character_)
}

res <- list(); outf <- file.path(OUT_DIR, sprintf("p1p2_isolation_%s.csv", RUN_TAG))
for (mm in M_SET) {
  f <- file.path(OUT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, mm)); if (!file.exists(f)) next
  ex <- readRDS(f); fe_draw <- setNames(ex$p1$entry, ex$p1$id_pe); cnt <- new.env()
  for (v in c("A_P1draw_P2draw", "B_P1draw_midpoint", "C_Wj2_P2draw", "D_Wj2_midpoint")) {
    cnt$clamped <- 0L
    fe <- if (v %in% c("A_P1draw_P2draw", "B_P1draw_midpoint")) fe_draw else W_half
    p2f <- if (v %in% c("A_P1draw_P2draw", "C_Wj2_P2draw")) make_p2_saved(ex$p2, fe, cnt) else p2_midpoint
    r <- tryCatch(fit_variant(fe, p2f), error = function(e) data.table(alpha = NA_real_, sd = NA_real_, ratio = NA_real_, events = NA_integer_, mlik = NA_real_, error = substr(conditionMessage(e), 1, 80)))
    res[[length(res) + 1]] <- cbind(data.table(m = mm, variant = v, clamped = cnt$clamped), r)
    print(tail(rbindlist(res), 1)[, .(m, variant, clamped, alpha = round(alpha, 3), sd = signif(sd, 3), ratio = round(ratio, 3), events)])
    fwrite(rbindlist(res), outf)
  }
}
R <- rbindlist(res)
cat("\n==== collapse (ratio < 0.1) by variant, flagged vs unflagged imputations ====\n")
R[, grp := ifelse(m %in% M_SET[1:8], "flagged", "unflagged")]
print(R[, .(n = .N, collapsed = sum(ratio < 0.1, na.rm = TRUE), errors = sum(!is.na(error)), median_alpha = round(median(alpha, na.rm = TRUE), 3)), by = .(grp, variant)])
cat("\nControl A should reproduce the pipeline fit for flagged imputations. Ratio by imputation:\n")
W <- dcast(R, m + grp ~ variant, value.var = "ratio"); num <- setdiff(names(W), c("m", "grp")); W[, (num) := lapply(.SD, round, 3), .SDcols = num]; print(W)
cat("Written:", outf, "\n")
