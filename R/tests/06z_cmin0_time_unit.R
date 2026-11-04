# ==============================================================================
# 06z_cmin0_time_unit.R   RUN ON THE MACHINE WITH INLA, fresh R session (~ 15-20 min), compact mode
# 06y: control.mode warm start has NO effect (bit-identical numbers); inla.rerun makes fits worse (singular-matrix crashes).
#      With cmin = 0 the fixed-alpha surface is smooth up to alpha 1.706 and fails from 1.707 (same in classic and compact mode): a WALL in alpha.
# 06i changed the time unit WITHOUT cmin = 0 (wall moved to ~1.9-2.1, free fits still erratic). cmin = 0 and the time unit have never been combined.
# TEST: for 5 collapsed data sets and time units s in {1, 4, 11, 22} (time' = time / s, an exact reparametrisation; alpha is unchanged), all with cmin = 0:
#   (a) the WALL: fixed-alpha grid 1.5-2.5 (step 0.05) for m02 and m08 -> where does mlik crash (drop > 100 below the best)?
#   (b) FREE fits for m02, m08, m14, m19, m23 -> alpha, SD, ratio  (healthy: ratio ~1.0-1.2, alpha 1.8-1.96)
# Output: data/processed/fit/cmin0_time_unit_wall_<tag>.csv, cmin0_time_unit_free_<tag>.csv
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(INLA) })
RUN_TAG <- "diag"; OUT_DIR <- here("data", "processed", "fit"); options(width = 200)
inla.setOption(inla.mode = "compact"); stopifnot(identical(inla.getOption("inla.mode"), "compact"))
LABEL <- "black_spruce->regen"; ALPHA_NAME <- "alpha parameter for weibullsurv"; S <- c(1, 4, 11, 22)
FORM <- function(rhs) as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", rhs))
fit <- function(d, rhs, fixed_alpha = NULL) {
  fam <- list(variant = 1); if (!is.null(fixed_alpha)) fam$hyper <- list(alpha = list(initial = 10 * log(fixed_alpha), fixed = TRUE))
  tryCatch(inla(FORM(rhs), family = "weibullsurv", data = d, control.family = fam,
                control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE), cmin = 0)), error = function(e) conditionMessage(e)) }
mk <- function(r) if (is.character(r)) NA_real_ else as.numeric(r$mlik[1, 1])
sc <- function(d, s) { d$Tstop <- d$Tstop / s; d$truncation <- d$truncation / s; d }
ld <- function(mm) readRDS(file.path(OUT_DIR, sprintf("baddata_%s_m%02d.rds", RUN_TAG, mm)))[[LABEL]]
# (b) free fits first (the decisive part)
fr <- list()
for (mm in c(2, 8, 14, 19, 23)) { b <- ld(mm); ev <- sum(b$data$status == 1L)
  for (s in S) { r <- fit(sc(b$data, s), b$rhs)
    fr[[length(fr) + 1]] <- if (is.character(r)) data.table(m = mm, s = s, alpha = NA_real_, sd = NA_real_, ratio = NA_real_, mlik = NA_real_, error = substr(r, 1, 50)) else {
      a <- r$summary.hyperpar[ALPHA_NAME, "mean"]; sdv <- r$summary.hyperpar[ALPHA_NAME, "sd"]
      data.table(m = mm, s = s, alpha = a, sd = sdv, ratio = sdv / (0.78 * a / sqrt(ev)), mlik = mk(r), error = NA_character_) }
    print(tail(rbindlist(fr), 1)[, .(m, s, alpha = round(alpha, 3), sd = signif(sd, 3), ratio = round(ratio, 3), mlik = round(mlik, 1), error)])
    fwrite(rbindlist(fr), file.path(OUT_DIR, sprintf("cmin0_time_unit_free_%s.csv", RUN_TAG))) } }
# (a) wall
wl <- list(); GR <- seq(1.5, 2.5, by = 0.05)
for (mm in c(2, 8)) { b <- ld(mm)
  for (s in S) { d <- sc(b$data, s); for (a in GR) wl[[length(wl) + 1]] <- data.table(m = mm, s = s, alpha = a, mlik = mk(fit(d, b$rhs, a)))
    fwrite(rbindlist(wl), file.path(OUT_DIR, sprintf("cmin0_time_unit_wall_%s.csv", RUN_TAG))); cat("wall grid done: m", mm, "s", s, "\n") } }
R <- rbindlist(fr); cat("\n==== (b) FREE fits with cmin = 0: SD ratio by time unit s (collapse < 0.1; healthy ~1.0-1.2) ====\n")
W <- dcast(R, m ~ s, value.var = "ratio"); num <- setdiff(names(W), "m"); W[, (num) := lapply(.SD, round, 3), .SDcols = num]; print(W)
cat("\nalpha (MLE of these data: 1.81-1.96):\n"); W <- dcast(R, m ~ s, value.var = "alpha"); W[, (num) := lapply(.SD, round, 3), .SDcols = num]; print(W)
G <- rbindlist(wl); G[, best := max(mlik, na.rm = TRUE), by = .(m, s)]; G[, fail := is.na(mlik) | mlik < best - 100]
cat("\n==== (a) WALL: failing alphas (mlik > 100 below best) per m and s ====\n")
print(G[, .(n_fail = sum(fail), first_fail = if (any(fail)) min(alpha[fail]) else NA_real_, fails_above_first = if (any(fail)) all(fail[alpha >= min(alpha[fail])]) else NA), by = .(m, s)])
