# ==============================================================================
# 06k_exit_structure_test.R   RUN ON THE MACHINE WITH INLA (~ 20 min)
# 06j: the collapse needs P1 entries AND P2 exits together (A); midpoint exits (B) never collapse.
# Which property of the P2 exits triggers it?  For m02, m08 (flagged) and m04 (unflagged), black_spruce->regen:
#  1. FEATURES of the A and B data: gap = Tstop - truncation (events and all rows), floored rows, distinct Tstop values.
#  2. INTERPOLATION: exit = (1-w)*midpoint + w*moment-matched(saved), w = 0, .25, .5, .75, 1   (w=1 is A, w=0 is B)
#     -> gradual or abrupt onset of the collapse?
#  3. JITTER of A: exit + u*bracket_width*Uniform(-.5,.5) clamped to the bracket, u = 0.001, 0.01
#     -> does tiny noise remove the collapse (exact structure) or not (overall distribution)?
# Output: data/processed/fit/exit_structure_features_<tag>.csv, exit_structure_fits_<tag>.csv
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(INLA) })
SIM_DIR <- "/Users/lukas/Desktop/RProjects/Quebec-forest-succession/R-scripts/Sim4stateExact/p1p2_v2_test"
RUN_TAG <- "diag"; OUT_DIR <- here("data", "processed", "fit"); SEED <- 2026L; options(width = 200)
M_SET <- c(2, 8, 4); LABEL <- "black_spruce->regen"; ALPHA_NAME <- "alpha parameter for weibullsurv"
W_GRID <- c(0, 0.25, 0.5, 0.75, 1); JIT <- c(0.001, 0.01)
source(file.path(SIM_DIR, "model_fitting.R"))
source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
fit_fun <- fit_inla_models
inp <- prepare_fit_inputs(here("data", "processed"), subsample_n = NULL, seed = SEED)
k <- match(LABEL, inp$td$label); stopifnot(!is.na(k)); fr <- strsplit(LABEL, "->", fixed = TRUE)[[1]]
cols <- c("from", "to", "status", "Tstop", "truncation", inp$terms)
iv <- inp$iv; iv_key <- paste(iv$id_pe, iv$interval)

# exit = (1-w)*midpoint + w*saved (+ jitter), always clamped into the bracket; E is tracked with the exit actually used
make_p2_blend <- function(p2, first_entry, w, jit = 0) {
  ex <- setNames(p2$exit_abs, paste(p2$id_pe, p2$interval)); last <- NA; E <- NA
  function(a_dur, b_dur, from, to, r) {
    id <- iv$id_pe[r]
    if (!identical(id, last)) { last <<- id; E <<- unname(first_entry[as.character(id)]) }
    Ts <- min(max(ex[[iv_key[r]]] - E, a_dur), b_dur)
    Td <- (1 - w) * (a_dur + b_dur) / 2 + w * Ts
    if (jit > 0) Td <- Td + jit * (b_dur - a_dur) * (runif(1) - 0.5)
    Td <- min(max(Td, a_dur), b_dur)
    E <<- E + Td; Td
  }
}
build_d <- function(fe, p2f) { b <- build_rows(iv, fe, inp$states, inp$tr, p2_fun = p2f, keep_cols = inp$keep_cols)
  list(d = as.data.frame(b$rows[kept == TRUE & from == fr[1] & to == fr[2], ..cols]), n_floored = b$n_floored) }
fit_d <- function(d) {
  m <- fit_fun(list(d), rhs_formula = inp$rhs[k], exact = TRUE)[[1]]; ev <- sum(d$status == 1L)
  if (!is.null(m$error)) return(data.table(alpha = NA_real_, sd = NA_real_, ratio = NA_real_, events = ev, mlik = NA_real_, error = substr(as.character(m$error)[1], 1, 80)))
  a <- m$summary.hyperpar[ALPHA_NAME, "mean"]; s <- m$summary.hyperpar[ALPHA_NAME, "sd"]
  data.table(alpha = a, sd = s, ratio = s / (0.78 * a / sqrt(ev)), events = ev, mlik = as.numeric(m$mlik[1, 1]), error = NA_character_)
}
feat <- function(d, m, v) { g <- d$Tstop - d$truncation; e <- d$status == 1L
  q <- function(x, p) as.numeric(quantile(x, p, names = FALSE))
  data.table(m = m, variant = v, rows = nrow(d), events = sum(e), n_gap_eps_events = sum(g[e] < 2 * EPSILON), min_gap_events = min(g[e]),
             ev_gap_p05 = q(g[e], .05), ev_gap_p25 = q(g[e], .25), ev_gap_p50 = q(g[e], .5), ev_gap_p75 = q(g[e], .75),
             all_gap_p05 = q(g, .05), all_gap_p50 = q(g, .5), distinct_Tstop_events = uniqueN(round(d$Tstop[e], 10)),
             Tstop_ev_p50 = q(d$Tstop[e], .5), Tstop_ev_p95 = q(d$Tstop[e], .95), trunc_p50 = q(d$truncation, .5), trunc_p95 = q(d$truncation, .95)) }

fits <- list(); feats <- list(); set.seed(SEED)
ff <- file.path(OUT_DIR, sprintf("exit_structure_fits_%s.csv", RUN_TAG)); fe_f <- file.path(OUT_DIR, sprintf("exit_structure_features_%s.csv", RUN_TAG))
for (mm in M_SET) {
  ex <- readRDS(file.path(OUT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, mm))); fe <- setNames(ex$p1$entry, ex$p1$id_pe)
  specs <- c(lapply(W_GRID, function(w) list(v = sprintf("w%.2f", w), w = w, jit = 0)),
             lapply(JIT, function(u) list(v = sprintf("A_jit%g", u), w = 1, jit = u)))
  for (s in specs) {
    b <- build_d(fe, make_p2_blend(ex$p2, fe, s$w, s$jit))
    if (s$v %in% c("w0.00", "w1.00")) { feats[[length(feats) + 1]] <- cbind(feat(b$d, mm, s$v), n_floored_all_transitions = b$n_floored); fwrite(rbindlist(feats), fe_f) }
    r <- tryCatch(fit_d(b$d), error = function(e) data.table(alpha = NA_real_, sd = NA_real_, ratio = NA_real_, events = NA_integer_, mlik = NA_real_, error = substr(conditionMessage(e), 1, 80)))
    fits[[length(fits) + 1]] <- cbind(data.table(m = mm, variant = s$v), r); fwrite(rbindlist(fits), ff)
    print(tail(rbindlist(fits), 1)[, .(m, variant, alpha = round(alpha, 3), sd = signif(sd, 3), ratio = round(ratio, 3), mlik = round(mlik, 1))])
  }
}
cat("\n==== FEATURES: B (w0.00) vs A (w1.00) ====\n"); print(melt(rbindlist(feats), id.vars = c("m", "variant"))[, .(m, variant, variable, value = signif(value, 4))] |> dcast(m + variable ~ variant, value.var = "value"))
cat("\n==== FITS (collapse = ratio < 0.1) ====\n"); print(dcast(rbindlist(fits), m ~ variant, value.var = "ratio")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
cat("Written:", ff, "and", fe_f, "\n")
