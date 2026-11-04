# ==============================================================================
# 06u_inner_trace.R   RUN ON THE MACHINE WITH INLA (~ 3 min):  Rscript R/tests/06u_inner_trace.R > data/processed/fit/inner_trace_all.txt 2>&1
# 06t: with cmin=0 (+stencil 9) the mlik surface is smooth where INLA works (noise 0.003-0.04, vs ~2e-4 per step of curvature signal),
# but fixed-alpha fits fail intermittently (stuck at -25,000 or partly off by -1.4..-43). Free fits collapse for every h.
# TEST: record INLA's own verbose log (inner Newton iterations, gradients, step sizes) for m02, cmin=0, stencil 9:
#   alpha 1.699 (works), 1.700 (partial, -1.5), 1.701 (stuck), 1.704 (partial, -14); and one free fit of m08 (collapsed).
# Output: the redirected inner_trace_all.txt (INLA writes verbose output from its C process, capture.output cannot catch it; cases are marked ######## CASE)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(INLA) })
RUN_TAG <- "diag"; OUT_DIR <- here("data", "processed", "fit"); LABEL <- "black_spruce->regen"
run <- function(d, rhs, extra, fixed_alpha = NULL, tag = "") {
  ci <- modifyList(list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)), extra)
  fam <- list(variant = 1); if (!is.null(fixed_alpha)) fam$hyper <- list(alpha = list(initial = 10 * log(fixed_alpha), fixed = TRUE))
  cat("\n######## CASE", tag, "########\n")
  r <- tryCatch(inla(as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", rhs)), family = "weibullsurv", data = d,
                control.family = fam, control.inla = ci, verbose = TRUE), error = function(e) conditionMessage(e))
  cat("######## END", tag, "mlik =", if (is.character(r)) NA else as.numeric(r$mlik[1, 1]), "########\n")
  list(r = r, log = character(0)) }
b <- readRDS(file.path(OUT_DIR, sprintf("baddata_%s_m02.rds", RUN_TAG)))[[LABEL]]; S <- list()
for (a in c(1.699, 1.700, 1.701, 1.704)) { o <- run(b$data, b$rhs, list(cmin = 0, stencil = 9L), a, tag = sprintf("m02_fixed_%.3f", a))
  writeLines(o$log, file.path(OUT_DIR, sprintf("inner_trace_%.3f_%s.txt", a, RUN_TAG)))
  S[[length(S) + 1]] <- data.table(case = sprintf("m02_fixed_%.3f", a), n_lines = length(o$log), mlik = if (is.character(o$r)) NA_real_ else as.numeric(o$r$mlik[1, 1])) }
b8 <- readRDS(file.path(OUT_DIR, sprintf("baddata_%s_m08.rds", RUN_TAG)))[[LABEL]]; o <- run(b8$data, b8$rhs, list(cmin = 0, stencil = 9L), tag = "m08_free")
writeLines(o$log, file.path(OUT_DIR, sprintf("inner_trace_free_m08_%s.txt", RUN_TAG)))
S[[length(S) + 1]] <- data.table(case = "m08_free", n_lines = length(o$log), mlik = if (is.character(o$r)) NA_real_ else as.numeric(o$r$mlik[1, 1]))
R <- rbindlist(S); print(R); fwrite(R, file.path(OUT_DIR, sprintf("inner_trace_summary_%s.csv", RUN_TAG)))
