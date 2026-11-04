# 06e_check_fixing.R  RUN ON THE MACHINE WITH INLA (< 1 min). Diagnostic only.
# First run: hyper = list(theta = list(initial=, fixed=TRUE)) was IGNORED (all.hyper stayed initial -2, fixed FALSE).
# Docs: alpha = exp(0.1*theta). This run tries other key names / nesting for the fixing and reports which one is honoured.
suppressPackageStartupMessages({ library(data.table); library(here); library(INLA) }); options(width = 200)
OUT_DIR <- here("data", "processed", "fit")
b <- readRDS(file.path(OUT_DIR, "baddata_diag_m02.rds"))[["black_spruce->regen"]]
hy <- inla.models()$likelihood$weibullsurv$hyper
cat("names of weibullsurv hyper in this INLA version:", paste(names(hy), collapse = ", "), "\n")
for (nm in names(hy)[1:min(2, length(hy))]) cat(sprintf("  %s: short.name=%s name=%s initial=%s fixed=%s prior=%s\n", nm, hy[[nm]]$short.name, hy[[nm]]$name, hy[[nm]]$initial, hy[[nm]]$fixed, hy[[nm]]$prior))
f <- function(cf) inla(as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", b$rhs)), family = "weibullsurv", data = b$data,
                       control.family = cf, control.compute = list(config = TRUE), control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)))
# hyper key "theta" is the right name (names above). It was still ignored -> try other routes:
f2 <- function(cf, cm = NULL) { a <- list(formula = as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", b$rhs)), family = "weibullsurv", data = b$data,
                       control.family = cf, control.compute = list(config = TRUE), control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)))
                       if (!is.null(cm)) a$control.mode <- cm; do.call(inla, a) }
# Earlier runs printed NO line for some specs (theta1, alpha, mode_fixed): most likely sprintf() with a zero-length argument
# (character(0)) silently prints nothing - i.e. those fits may have WORKED. Printing is now zero-length safe.
S <- function(x, d = 3) if (length(x) == 0 || is.null(x)) "<none>" else paste(if (is.numeric(x)) round(x, d) else x, collapse = ",")
specs <- list(
  hyper_theta1  = function(a) list(cf = list(variant = 1, hyper = list(theta1 = list(initial = 10 * log(a), fixed = TRUE)))),
  hyper_alpha   = function(a) list(cf = list(variant = 1, hyper = list(alpha = list(initial = 10 * log(a), fixed = TRUE)))),
  mode_fixed    = function(a) list(cf = list(variant = 1), cm = list(theta = 10 * log(a), fixed = TRUE)),
  hyper_theta   = function(a) list(cf = list(variant = 1, hyper = list(theta = list(initial = 10 * log(a), fixed = TRUE)))))
for (sn in names(specs)) for (a in c(1.0, 2.2)) {
  sp <- specs[[sn]](a); r <- tryCatch(f2(sp$cf, sp$cm), error = function(e) e)
  if (inherits(r, "error")) { cat(sn, "alpha", a, ": ERROR", substr(conditionMessage(r), 1, 120), "\n"); next }
  h <- r$all.hyper$family[[1]]$hyper[[1]]
  cat(sn, "alpha", a, "| all.hyper initial", S(h$initial), "fixed", S(h$fixed), "| hyperpar rows", S(nrow(r$summary.hyperpar)),
      "| mode$theta", S(r$mode$theta), "(alpha", S(exp(0.1 * r$mode$theta)), ") | mlik", S(r$mlik[1, 1], 2), "| names(all.hyper$family[[1]]$hyper):", S(names(r$all.hyper$family[[1]]$hyper)), "\n")
}
cat("Honoured = all.hyper shows fixed TRUE and initial = 10*log(alpha) (0 for 1.0, 7.885 for 2.2), 0 hyperpar rows, mlik differs between the two alphas.\n")
