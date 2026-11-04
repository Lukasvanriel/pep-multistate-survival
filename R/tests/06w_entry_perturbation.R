# ==============================================================================
# 06w_entry_perturbation.R   RUN ON THE MACHINE WITH INLA (~ 1 h; progress is written after every fit)
# GOAL: move the 8 collapsed black_spruce->regen imputations OUT of the "twilight zone" by a SMALL, UNIFORM change of the P1 entry draws.
# Every rule below is applied to ALL plots alike (no plot- or imputation-specific tuning); P2 exits stay as saved (clamped into the bracket, as 06j-A).
#   scale_0.90 / scale_0.75 : entry' = s * entry
#   cap_100 / cap_80        : entry' = min(entry, C) years
#   blend_0.25 / blend_0.50 : entry' = (1 - L) * entry + L * W_j/2   (toward the Stage-1 value)
# Rows: the 8 collapsed imputations (m = 2, 8, 14, 15, 19, 23, 25, 29) + healthy control m = 4 with the strongest dose of each family.
# Fits use the pipeline settings (free alpha). Output: data/processed/fit/entry_perturbation_<tag>.csv   (collapse = ratio < 0.1; healthy ~ 1.0-1.2)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(INLA) })
SIM_DIR <- "/Users/lukas/Desktop/RProjects/Quebec-forest-succession/R-scripts/Sim4stateExact/p1p2_v2_test"
RUN_TAG <- "diag"; OUT_DIR <- here("data", "processed", "fit"); SEED <- 2026L; options(width = 200)

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

FLAG <- c(2, 8, 14, 15, 19, 23, 25, 29); CTRL <- 4
RULES <- list(scale_0.90 = function(e, w) 0.90 * e, scale_0.75 = function(e, w) 0.75 * e, cap_100 = function(e, w) pmin(e, 100), cap_80 = function(e, w) pmin(e, 80),
              blend_0.25 = function(e, w) 0.75 * e + 0.25 * w, blend_0.50 = function(e, w) 0.5 * e + 0.5 * w)
STRONG <- c("scale_0.75", "cap_80", "blend_0.50")
res <- list(); outf <- file.path(OUT_DIR, sprintf("entry_perturbation_%s.csv", RUN_TAG))
for (mm in c(FLAG, CTRL)) {
  f <- file.path(OUT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, mm)); if (!file.exists(f)) next
  ex <- readRDS(f); fe0 <- setNames(ex$p1$entry, ex$p1$id_pe); w0 <- W_half[names(fe0)]; cnt <- new.env()
  for (v in names(RULES)) {
    if (mm == CTRL && !(v %in% STRONG)) next
    fe <- setNames(RULES[[v]](unname(fe0), unname(w0)), names(fe0)); cnt$clamped <- 0L
    r <- tryCatch(fit_variant(fe, make_p2_saved(ex$p2, fe, cnt)), error = function(e) data.table(alpha = NA_real_, sd = NA_real_, ratio = NA_real_, events = NA_integer_, mlik = NA_real_, error = substr(conditionMessage(e), 1, 80)))
    res[[length(res) + 1]] <- cbind(data.table(m = mm, rule = v, mean_shift = mean(fe - fe0), clamped = cnt$clamped), r)
    print(tail(rbindlist(res), 1)[, .(m, rule, mean_shift = round(mean_shift, 1), clamped, alpha = round(alpha, 3), sd = signif(sd, 3), ratio = round(ratio, 3), error)])
    fwrite(rbindlist(res), outf)
  }
}
R <- rbindlist(res)
cat("\n==== SD ratio by rule (collapse < 0.1; healthy ~1.0-1.2) ====\n"); W <- dcast(R, m ~ rule, value.var = "ratio"); num <- setdiff(names(W), "m"); W[, (num) := lapply(.SD, round, 3), .SDcols = num]; print(W)
cat("\n==== alpha by rule (original pipeline alpha for these imputations: 1.68-1.83) ====\n"); W <- dcast(R, m ~ rule, value.var = "alpha"); W[, (num) := lapply(.SD, round, 3), .SDcols = num]; print(W)
cat("\ncollapsed count by rule (of the 8):\n"); print(R[m %in% FLAG, .(n = .N, collapsed = sum(ratio < 0.1, na.rm = TRUE), errors = sum(!is.na(error))), by = rule])
